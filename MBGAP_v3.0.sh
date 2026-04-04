#!/usr/bin/env bash
set -euo pipefail

PROGRAM_NAME="$(basename "$0")"
VERSION="3.0"

# =========================
# DEFAULTS
# =========================
INPUT_DIR=""
SAMPLESHEET=""
OUTPUT_DIR=""
ASSEMBLER="spades"
TRIM="no"
ADAPTERS_FILE=""

GENUS=""
SPECIES=""
GRAM="?"
GENETIC_CODE="11"
USE_PROKKA_GENUS="no"
BAKTA_COMPLETE="no"

GAMBIT_DB=""
BAKTA_DB=""
PLASMIDFINDER_DB=""

SHOVILL_GSIZE=""
SHOVILL_GSIZE_DEFAULT=""

CONDA_ENV="bioinfo"
GTDBTK_ENV="gtdbtk"

THREADS_DEFAULT=""
THREADS_FASTQC="8"
THREADS_TRIMMO="8"
THREADS_SPADES="10"
THREADS_SHOVILL="10"
THREADS_UNICYCLER="10"
THREADS_QUAST="8"
THREADS_PROKKA="8"
THREADS_BAKTA="12"
THREADS_AMRFINDER="8"
THREADS_GTDBTK="10"

MEMORY_SPADES="50"
RAM_SHOVILL="50"

RUN_GTDBTK="yes"
RUN_GAMBIT="yes"
RUN_MULTIQC="yes"

LOG_FILE=""
CHECKPOINT_DIR=""

# =========================
# GLOBAL ARRAYS
# =========================
declare -a SAMPLE_ORDER=()

declare -A SAMPLE_R1
declare -A SAMPLE_R2
declare -A SAMPLE_GENUS
declare -A SAMPLE_SPECIES
declare -A SAMPLE_GRAMPHENO
declare -A SAMPLE_GSIZE
declare -A SAMPLE_SOURCE

declare -A SAMPLE_FINAL_GENUS
declare -A SAMPLE_FINAL_SPECIES
declare -A SAMPLE_GENUS_SOURCE
declare -A SAMPLE_SPECIES_SOURCE

# =========================
# HELP
# =========================
usage() {
    cat <<EOF
MBGAP 3.0 - Mykaella's Bacterial Genome Analysis Pipeline
Versão com metadata por amostra + inferência via Gambit + parser CSV robusto

Uso:
  $PROGRAM_NAME --output DIR --assembler {spades|shovill|unicycler} [opções]

Modos de entrada (escolha um):
  --input DIR                 Diretório contendo FASTQ(.gz) paired-end
  --samplesheet FILE          CSV com metadados por amostra

Obrigatórios:
  --output DIR                Diretório de saída
  --assembler NAME            spades | shovill | unicycler
  --gambit-db DIR             Caminho para database do Gambit
  --bakta-db DIR              Caminho para database do Bakta
  --plasmidfinder-db DIR      Caminho para database do PlasmidFinder

Condicionais:
  --adapters FILE             Arquivo de adaptadores do Trimmomatic (obrigatório se --trim yes)
  --gsize VALUE               Tamanho estimado global do genoma para Shovill
                              (obrigatório se --assembler shovill e não vier por amostra no CSV)

Anotação / organismo:
  --genus TEXT                Gênero bacteriano global
  --species TEXT              Espécie bacteriana global
  --gram {+|-|?}              Tipo de Gram global [default: ?]
  --genetic-code INT          Código genético para Prokka [default: 11]
  --use-prokka-genus yes|no   Usar --usegenus no Prokka [default: no]
  --bakta-complete yes|no     Usar --complete no Bakta [default: no]

Pipeline:
  --trim yes|no               Executar trimagem [default: no]
  --run-gtdbtk yes|no         Executar GTDB-Tk [default: yes]
  --run-gambit yes|no         Executar Gambit [default: yes]
  --run-multiqc yes|no        Executar MultiQC [default: yes]

Ambientes Conda:
  --conda-env NAME            Ambiente conda principal [default: bioinfo]
  --gtdbtk-env NAME           Ambiente conda do GTDB-Tk [default: gtdbtk]

Threads e recursos:
  --threads INT               Define todas as threads de uma vez
  --threads-fastqc INT
  --threads-trimmomatic INT
  --threads-spades INT
  --threads-shovill INT
  --threads-unicycler INT
  --threads-quast INT
  --threads-prokka INT
  --threads-bakta INT
  --threads-amrfinder INT
  --threads-gtdbtk INT
  --memory-spades INT         Memória em GB para SPAdes [default: 50]
  --ram-shovill INT           RAM em GB para Shovill [default: 50]

Outros:
  --help                      Mostra esta ajuda
  --version                   Mostra a versão

Formato do samplesheet CSV:
  sample,r1,r2,genus,species,gram,gsize

Exemplo:
  sample,r1,r2,genus,species,gram,gsize
  ISO001,"/dados/ISO001_R1.fastq.gz","/dados/ISO001_R2.fastq.gz",Escherichia,coli,-,5.2M
  ISO002,"/dados/ISO002_R1.fastq.gz","/dados/ISO002_R2.fastq.gz",Klebsiella,pneumoniae,-,5.5M
  ISO003,"/dados/ISO003_R1.fastq.gz","/dados/ISO003_R2.fastq.gz",,,,4.9M

Regras:
  - se --samplesheet for informado, ele substitui o modo --input
  - genus/species/gram/gsize são opcionais por amostra
  - gram é usado somente no Bakta
  - genus/species faltantes podem ser inferidos a partir de query + predicted.name do Gambit
  - se --assembler shovill, gsize deve existir por amostra ou via --gsize global

EOF
}

# =========================
# LOG / UTILS
# =========================
log() {
    local msg="$1"
    local entry="[ $(date '+%Y-%m-%d %H:%M:%S') ] $msg"
    echo "$entry" >&2
    [[ -n "${LOG_FILE:-}" ]] && echo "$entry" >> "$LOG_FILE"
}

die() {
    local msg="$1"
    echo "Erro: $msg" >&2
    [[ -n "${LOG_FILE:-}" ]] && echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Erro: $msg" >> "$LOG_FILE"
    exit 1
}

check_command() {
    command -v "$1" >/dev/null 2>&1 || die "Dependência não encontrada no PATH: $1"
}

check_file() {
    [[ -f "$1" ]] || die "Arquivo não encontrado: $1"
    [[ -s "$1" ]] || die "Arquivo vazio: $1"
}

check_dir() {
    [[ -d "$1" ]] || die "Diretório não encontrado: $1"
}

is_yes_no() {
    [[ "$1" == "yes" || "$1" == "no" ]]
}

is_assembler() {
    [[ "$1" == "spades" || "$1" == "shovill" || "$1" == "unicycler" ]]
}

is_gram() {
    [[ "$1" == "+" || "$1" == "-" || "$1" == "?" || -z "$1" ]]
}

trim_spaces() {
    local s="$1"
    s="${s#"${s%%[![:space:]]*}"}"
    s="${s%"${s##*[![:space:]]}"}"
    printf '%s' "$s"
}

strip_quotes() {
    local s="$1"
    s="$(trim_spaces "$s")"
    s="${s%\"}"
    s="${s#\"}"
    printf '%s' "$s"
}

normalize_na() {
    local s
    s="$(strip_quotes "$1")"
    case "${s,,}" in
        ""|"na"|"n/a"|"null"|"none"|".")
            printf ''
            ;;
        *)
            printf '%s' "$s"
            ;;
    esac
}

safe_sample_name() {
    local s="$1"
    [[ "$s" =~ ^[A-Za-z0-9._-]+$ ]] || die "Nome de amostra inválido: '$s'. Use apenas letras, números, ponto, underscore e hífen."
}


# =========================
# CHECKPOINT
# =========================
step_done() {
    local tag="$1"
    [[ -n "$CHECKPOINT_DIR" && -f "$CHECKPOINT_DIR/${tag}.done" ]]
}

mark_done() {
    local tag="$1"
    mkdir -p "$CHECKPOINT_DIR"
    touch "$CHECKPOINT_DIR/${tag}.done"
}

run_step() {
    local tag="$1"; shift
    local desc="$1"; shift
    if step_done "$tag"; then
        log "SKIP (já concluído): $desc"
        return 0
    fi
    log "$desc"
    "$@"
    mark_done "$tag"
}

trap 'die "Pipeline interrompido por erro."' ERR


# =========================
# ARG PARSING
# =========================
[[ $# -eq 0 ]] && { usage; exit 1; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input) INPUT_DIR="$2"; shift 2 ;;
        --samplesheet) SAMPLESHEET="$2"; shift 2 ;;
        --output) OUTPUT_DIR="$2"; shift 2 ;;
        --assembler) ASSEMBLER="$2"; shift 2 ;;
        --trim) TRIM="$2"; shift 2 ;;
        --adapters) ADAPTERS_FILE="$2"; shift 2 ;;
        --genus) GENUS="$2"; shift 2 ;;
        --species) SPECIES="$2"; shift 2 ;;
        --gram) GRAM="$2"; shift 2 ;;
        --genetic-code) GENETIC_CODE="$2"; shift 2 ;;
        --use-prokka-genus) USE_PROKKA_GENUS="$2"; shift 2 ;;
        --bakta-complete) BAKTA_COMPLETE="$2"; shift 2 ;;
        --gambit-db) GAMBIT_DB="$2"; shift 2 ;;
        --bakta-db) BAKTA_DB="$2"; shift 2 ;;
        --plasmidfinder-db) PLASMIDFINDER_DB="$2"; shift 2 ;;
        --gsize) SHOVILL_GSIZE_DEFAULT="$2"; shift 2 ;;
        --conda-env) CONDA_ENV="$2"; shift 2 ;;
        --gtdbtk-env) GTDBTK_ENV="$2"; shift 2 ;;
        --threads) THREADS_DEFAULT="$2"; shift 2 ;;
        --threads-fastqc) THREADS_FASTQC="$2"; shift 2 ;;
        --threads-trimmomatic) THREADS_TRIMMO="$2"; shift 2 ;;
        --threads-spades) THREADS_SPADES="$2"; shift 2 ;;
        --threads-shovill) THREADS_SHOVILL="$2"; shift 2 ;;
        --threads-unicycler) THREADS_UNICYCLER="$2"; shift 2 ;;
        --threads-quast) THREADS_QUAST="$2"; shift 2 ;;
        --threads-prokka) THREADS_PROKKA="$2"; shift 2 ;;
        --threads-bakta) THREADS_BAKTA="$2"; shift 2 ;;
        --threads-amrfinder) THREADS_AMRFINDER="$2"; shift 2 ;;
        --threads-gtdbtk) THREADS_GTDBTK="$2"; shift 2 ;;
        --memory-spades) MEMORY_SPADES="$2"; shift 2 ;;
        --ram-shovill) RAM_SHOVILL="$2"; shift 2 ;;
        --run-gtdbtk) RUN_GTDBTK="$2"; shift 2 ;;
        --run-gambit) RUN_GAMBIT="$2"; shift 2 ;;
        --run-multiqc) RUN_MULTIQC="$2"; shift 2 ;;
        --help|-h) usage; exit 0 ;;
        --version) echo "$VERSION"; exit 0 ;;
        *)
            die "Parâmetro desconhecido: $1. Use --help."
            ;;
    esac
done


# =========================
# GLOBAL THREAD OVERRIDE
# =========================
if [[ -n "$THREADS_DEFAULT" ]]; then
    THREADS_FASTQC="$THREADS_DEFAULT"
    THREADS_TRIMMO="$THREADS_DEFAULT"
    THREADS_SPADES="$THREADS_DEFAULT"
    THREADS_SHOVILL="$THREADS_DEFAULT"
    THREADS_UNICYCLER="$THREADS_DEFAULT"
    THREADS_QUAST="$THREADS_DEFAULT"
    THREADS_PROKKA="$THREADS_DEFAULT"
    THREADS_BAKTA="$THREADS_DEFAULT"
    THREADS_AMRFINDER="$THREADS_DEFAULT"
    THREADS_GTDBTK="$THREADS_DEFAULT"
fi


# =========================
# VALIDATION
# =========================
[[ -n "$OUTPUT_DIR" ]] || die "Informe --output"
[[ -n "$GAMBIT_DB" ]] || die "Informe --gambit-db"
[[ -n "$BAKTA_DB" ]] || die "Informe --bakta-db"
[[ -n "$PLASMIDFINDER_DB" ]] || die "Informe --plasmidfinder-db"

is_assembler "$ASSEMBLER" || die "--assembler deve ser: spades, shovill ou unicycler"
is_yes_no "$TRIM" || die "--trim deve ser yes ou no"
is_yes_no "$USE_PROKKA_GENUS" || die "--use-prokka-genus deve ser yes ou no"
is_yes_no "$BAKTA_COMPLETE" || die "--bakta-complete deve ser yes ou no"
is_yes_no "$RUN_GTDBTK" || die "--run-gtdbtk deve ser yes ou no"
is_yes_no "$RUN_GAMBIT" || die "--run-gambit deve ser yes ou no"
is_yes_no "$RUN_MULTIQC" || die "--run-multiqc deve ser yes ou no"
is_gram "$GRAM" || die "--gram deve ser +, - ou ?"

check_dir "$GAMBIT_DB"
check_dir "$BAKTA_DB"
check_dir "$PLASMIDFINDER_DB"

if [[ -n "$INPUT_DIR" && -n "$SAMPLESHEET" ]]; then
    die "Use apenas um modo de entrada: --input OU --samplesheet"
fi

if [[ -z "$INPUT_DIR" && -z "$SAMPLESHEET" ]]; then
    die "Informe --input ou --samplesheet"
fi

if [[ -n "$INPUT_DIR" ]]; then
    check_dir "$INPUT_DIR"
fi

if [[ -n "$SAMPLESHEET" ]]; then
    check_file "$SAMPLESHEET"
fi

if [[ "$TRIM" == "yes" ]]; then
    [[ -n "$ADAPTERS_FILE" ]] || die "--adapters é obrigatório quando --trim yes"
    check_file "$ADAPTERS_FILE"
fi

mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/Assembly_final"
LOG_FILE="$OUTPUT_DIR/pipeline_log.txt"
CHECKPOINT_DIR="$OUTPUT_DIR/.checkpoints"
mkdir -p "$CHECKPOINT_DIR"


# =========================
# CONDA
# =========================
check_command conda
eval "$(conda shell.bash hook)" || die "Falha ao inicializar o conda shell hook"
conda activate "$CONDA_ENV" || die "Não foi possível ativar o ambiente conda: $CONDA_ENV"

# =========================
# DEPENDENCIES (após ativação do env)
# =========================
check_command fastqc
check_command trimmomatic
check_command spades.py
check_command shovill
check_command unicycler
check_command quast.py
check_command prokka
check_command bakta
check_command amrfinder
check_command plasmidfinder.py
check_command assembly-scan
check_command multiqc
check_command conda
check_command awk
check_command sed
check_command grep
check_command python3

if [[ "$RUN_GAMBIT" == "yes" ]]; then
    check_command gambit
fi


# =========================
# SUMMARY
# =========================
log "MBGAP $VERSION iniciado"
log "Modo de entrada: $( [[ -n "$SAMPLESHEET" ]] && echo "samplesheet" || echo "input-dir" )"
[[ -n "$INPUT_DIR" ]] && log "Input: $INPUT_DIR"
[[ -n "$SAMPLESHEET" ]] && log "Samplesheet: $SAMPLESHEET"
log "Output: $OUTPUT_DIR"
log "Assembler: $ASSEMBLER"
log "Trim: $TRIM"
log "Genus global: ${GENUS:-NA}"
log "Species global: ${SPECIES:-NA}"
log "Gram global: $GRAM"
log "Prokka genetic code: $GENETIC_CODE"
log "Prokka use genus: $USE_PROKKA_GENUS"
log "Bakta complete: $BAKTA_COMPLETE"
log "Gambit DB: $GAMBIT_DB"
log "Bakta DB: $BAKTA_DB"
log "PlasmidFinder DB: $PLASMIDFINDER_DB"
[[ "$ASSEMBLER" == "shovill" ]] && log "Shovill gsize global: ${SHOVILL_GSIZE_DEFAULT:-NA}"



# =========================
# FUNCTIONS
# =========================
create_directory() {
    mkdir -p "$1"
}

register_sample() {
    local sample="$1"
    local r1="$2"
    local r2="$3"
    local genus="$4"
    local species="$5"
    local gram="$6"
    local gsize="$7"
    local source="$8"

    safe_sample_name "$sample"
    check_file "$r1"
    check_file "$r2"

    [[ -z "${SAMPLE_R1[$sample]+x}" ]] || die "Amostra duplicada detectada: $sample"

    genus="$(normalize_na "$genus")"
    species="$(normalize_na "$species")"
    gram="$(normalize_na "$gram")"
    gsize="$(normalize_na "$gsize")"

    is_gram "$gram" || die "Valor inválido para gram na amostra '$sample': '$gram'"

    SAMPLE_ORDER+=("$sample")
    SAMPLE_R1["$sample"]="$r1"
    SAMPLE_R2["$sample"]="$r2"
    SAMPLE_GENUS["$sample"]="$genus"
    SAMPLE_SPECIES["$sample"]="$species"
    SAMPLE_GRAMPHENO["$sample"]="$gram"
    SAMPLE_GSIZE["$sample"]="$gsize"
    SAMPLE_SOURCE["$sample"]="$source"
}

load_samplesheet() {
    local parsed_tsv
    parsed_tsv="$(mktemp)"

    python3 - "$SAMPLESHEET" > "$parsed_tsv" <<'PY'
import csv
import sys

samplesheet = sys.argv[1]

required = ["sample", "r1", "r2"]
allowed = ["sample", "r1", "r2", "genus", "species", "gram", "gsize"]

def norm(v):
    if v is None:
        return ""
    v = v.strip()
    if v.lower() in {"", "na", "n/a", "null", "none", "."}:
        return ""
    return v

with open(samplesheet, newline="", encoding="utf-8-sig") as fh:
    reader = csv.DictReader(fh)
    if reader.fieldnames is None:
        raise SystemExit("Erro: samplesheet sem cabeçalho")

    fieldnames = [f.strip() for f in reader.fieldnames]
    lowered = [f.lower() for f in fieldnames]

    for req in required:
        if req not in lowered:
            raise SystemExit(f"Erro: coluna obrigatória ausente no samplesheet: {req}")

    unknown = [f for f in lowered if f not in allowed]
    if unknown:
        raise SystemExit("Erro: colunas não reconhecidas no samplesheet: " + ", ".join(unknown))

    name_map = {}
    for original, low in zip(fieldnames, lowered):
        name_map[low] = original

    line_no = 1
    for row in reader:
        line_no += 1
        if row is None:
            continue

        values = {k.lower().strip(): (v if v is not None else "") for k, v in row.items() if k is not None}

        if all((str(v).strip() == "") for v in values.values()):
            continue

        sample = norm(values.get("sample", ""))
        r1 = norm(values.get("r1", ""))
        r2 = norm(values.get("r2", ""))
        genus = norm(values.get("genus", ""))
        species = norm(values.get("species", ""))
        gram = norm(values.get("gram", ""))
        gsize = norm(values.get("gsize", ""))

        if not sample:
            raise SystemExit(f"Erro: linha {line_no} do CSV sem campo sample")
        if not r1:
            raise SystemExit(f"Erro: linha {line_no} do CSV sem campo r1")
        if not r2:
            raise SystemExit(f"Erro: linha {line_no} do CSV sem campo r2")

        print("\t".join([sample, r1, r2, genus, species, gram, gsize]))
PY

    local line_no=0
    local sample r1 r2 genus species gram gsize

    while IFS=$'\t' read -r sample r1 r2 genus species gram gsize; do
        line_no=$((line_no + 1))
        [[ -n "$sample" ]] || continue
        register_sample "$sample" "$r1" "$r2" "$genus" "$species" "$gram" "$gsize" "samplesheet"
    done < "$parsed_tsv"

    rm -f "$parsed_tsv"

    [[ "${#SAMPLE_ORDER[@]}" -gt 0 ]] || die "Nenhuma amostra válida foi carregada do samplesheet"
}



discover_samples_from_input() {
    local found=0
    local r1 sample r2

    for r1 in "$INPUT_DIR"/*_R1_001.fastq.gz "$INPUT_DIR"/*_1.fastq.gz; do
        [[ -f "$r1" ]] || continue

        if [[ "$r1" == *_R1_001.fastq.gz ]]; then
            sample="$(basename "$r1" _R1_001.fastq.gz)"
            r2="$INPUT_DIR/${sample}_R2_001.fastq.gz"
        else
            sample="$(basename "$r1" _1.fastq.gz)"
            r2="$INPUT_DIR/${sample}_2.fastq.gz"
        fi

        if [[ ! -f "$r2" ]]; then
            log "Aviso: par R2 ausente para $sample. Amostra ignorada."
            continue
        fi

        found=1
        register_sample "$sample" "$r1" "$r2" "$GENUS" "$SPECIES" "$GRAM" "$SHOVILL_GSIZE_DEFAULT" "input-dir"
    done

    [[ "$found" -eq 1 ]] || die "Nenhuma amostra paired-end válida encontrada em $INPUT_DIR"
}



validate_shovill_gsizes() {
    local sample gsize

    if [[ "$ASSEMBLER" != "shovill" ]]; then
        return 0
    fi

    for sample in "${SAMPLE_ORDER[@]}"; do
        gsize="${SAMPLE_GSIZE[$sample]}"
        if [[ -z "$gsize" ]]; then
            if [[ -n "$SHOVILL_GSIZE_DEFAULT" ]]; then
                SAMPLE_GSIZE["$sample"]="$SHOVILL_GSIZE_DEFAULT"
            else
                die "Amostra '$sample' sem gsize. Forneça a coluna gsize no CSV ou use --gsize."
            fi
        fi
    done
}



run_fastqc() {
    local outdir="$1"
    local r1="$2"
    local r2="$3"

    echo

    log "FastQC: $(basename "$r1"), $(basename "$r2")"
    fastqc -t "$THREADS_FASTQC" -o "$outdir" "$r1" "$r2"
}



run_trimmomatic() {
    local r1="$1"
    local r2="$2"
    local p1="$3"
    local u1="$4"
    local p2="$5"
    local u2="$6"
    local summary="$7"

    echo

    log "Trimmomatic: $(basename "$r1")"
    _JAVA_OPTIONS="-Xmx32g" trimmomatic PE -phred33 \
        "$r1" "$r2" \
        "$p1" "$u1" \
        "$p2" "$u2" \
        ILLUMINACLIP:"$ADAPTERS_FILE":2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:50 \
        -threads "$THREADS_TRIMMO" \
        -summary "$summary"
}



run_spades() {
    local r1="$1"
    local r2="$2"
    local outdir="$3"

    echo

    log "SPAdes"
    spades.py \
        -1 "$r1" \
        -2 "$r2" \
        -o "$outdir" \
        --isolate \
        --cov-cutoff auto \
        -t "$THREADS_SPADES" \
        --memory "$MEMORY_SPADES" \
        -k 21,33,55,77,99,127
}



run_shovill() {
    local r1="$1"
    local r2="$2"
    local outdir="$3"
    local sample="$4"
    local gsize="$5"

    echo

    [[ -n "$gsize" ]] || die "gsize ausente para amostra '$sample' usando Shovill"

    echo

    log "Shovill"
    shovill \
        --R1 "$r1" \
        --R2 "$r2" \
        --outdir "$outdir" \
        --assembler spades \
        --minlen 300 \
        --mincov 10 \
        --force \
        --keepfiles \
        --noreadcorr \
        --namefmt "${sample}_%05d" \
        --depth 100 \
        --gsize "$gsize" \
        --cpus "$THREADS_SHOVILL" \
        --ram "$RAM_SHOVILL"
}



run_unicycler() {
    local r1="$1"
    local r2="$2"
    local outdir="$3"

    echo

    log "Unicycler"
    unicycler \
        -1 "$r1" \
        -2 "$r2" \
        -o "$outdir" \
        --threads "$THREADS_UNICYCLER" \
        --min_fasta_length 300 \
        --keep 2 \
        --mode normal
}



rename_and_copy_assembly() {
    local assembler="$1"
    local sample_dir="$2"
    local sample="$3"
    local out_final="$OUTPUT_DIR/Assembly_final/${sample}_${assembler}.fasta"

    case "$assembler" in
        spades)
            [[ -f "$sample_dir/${sample}_spades/scaffolds.fasta" ]] || die "scaffolds.fasta não encontrado para $sample"
            cp "$sample_dir/${sample}_spades/scaffolds.fasta" "$out_final"
            ;;
        shovill)
            [[ -f "$sample_dir/${sample}_shovill/contigs.fa" ]] || die "contigs.fa não encontrado para $sample"
            cp "$sample_dir/${sample}_shovill/contigs.fa" "$out_final"
            ;;
        unicycler)
            [[ -f "$sample_dir/${sample}_unicycler/assembly.fasta" ]] || die "assembly.fasta não encontrado para $sample"
            cp "$sample_dir/${sample}_unicycler/assembly.fasta" "$out_final"
            ;;
        *)
            die "Assembler não reconhecido: $assembler"
            ;;
    esac
}

echo

make_prokka_fasta() {
    local fasta="$1"
    local sample="$2"
    local fasta_clean="${fasta%.fasta}_prokka_clean.fasta"

    awk -v prefix="${sample}" \
        'BEGIN{n=0} /^>/{n++; printf ">%s_contig_%05d\n", prefix, n; next} {print}' \
        "$fasta" > "$fasta_clean"

    log "Prokka: contigs renomeados em $(basename "$fasta_clean") (IDs ≤37 chars)"
    echo "$fasta_clean"
}



run_quast() {
    local sample="$1"
    local sample_dir="$2"

    echo

    quast.py \
        "$OUTPUT_DIR/Assembly_final/${sample}_${ASSEMBLER}.fasta" \
        -o "$sample_dir/${sample}_quast" \
        --circos \
        --glimmer \
        --rna-finding \
        --threads "$THREADS_QUAST"
}



run_assembly_scan() {
    local sample="$1"
    local sample_dir="$2"

    echo

    bash -c "assembly-scan \
        '$OUTPUT_DIR/Assembly_final/${sample}_${ASSEMBLER}.fasta' \
        --prefix '$sample' \
        > '$sample_dir/${sample}_${ASSEMBLER}.tsv'"
}



run_prokka() {
    local fasta="$1"
    local outdir="$2"
    local sample="$3"
    local sample_genus="$4"
    local sample_species="$5"

    local fasta_clean
    fasta_clean="$(make_prokka_fasta "$fasta" "$sample")"

    local cmd=(
        prokka
        --outdir "$outdir"
        --prefix prokka
        --kingdom Bacteria
        --gcode "$GENETIC_CODE"
        --rfam
        --force
        --cpus "$THREADS_PROKKA"
    )

    if [[ -n "$sample_genus" ]]; then
        cmd+=(--genus "$sample_genus")
        [[ "$USE_PROKKA_GENUS" == "yes" ]] && cmd+=(--usegenus)
    fi

    if [[ -n "$sample_species" ]]; then
        cmd+=(--species "$sample_species")
    fi

    echo

    log "Prokka | sample=$sample | genus=${sample_genus:-NA} | species=${sample_species:-NA}"
    "${cmd[@]}" "$fasta_clean"

    echo

    rm -f "$fasta_clean"
    log "Prokka: FASTA temporário removido"
}



run_bakta() {
    local fasta="$1"
    local outdir="$2"
    local sample="$3"
    local sample_genus="$4"
    local sample_species="$5"
    local sample_gram="$6"

    local cmd=(
        bakta
        --db "$BAKTA_DB"
        --output "$outdir"
        --prefix "$sample"
        --force
        --verbose
        --threads "$THREADS_BAKTA"
    )

    [[ -n "$sample_gram" && "$sample_gram" != "?" ]] && cmd+=(--gram "$sample_gram")
    [[ -n "$sample_genus" ]] && cmd+=(--genus "$sample_genus")
    [[ -n "$sample_species" ]] && cmd+=(--species "$sample_species")
    [[ "$BAKTA_COMPLETE" == "yes" ]] && cmd+=(--complete)

    echo

    log "Bakta | sample=$sample | genus=${sample_genus:-NA} | species=${sample_species:-NA} | gram=${sample_gram:-NA}"
    "${cmd[@]}" "$fasta"
}



run_amrfinder() {
    local sample="$1"
    local sample_dir="$2"

    echo

    amrfinder \
        --nucleotide "$OUTPUT_DIR/Assembly_final/${sample}_${ASSEMBLER}.fasta" \
        --plus \
        --threads "$THREADS_AMRFINDER" \
        -o "$sample_dir/${sample}_amrfinder/${sample}.tsv" \
        --nucleotide_output "$sample_dir/${sample}_amrfinder/${sample}_AMR_genes.fasta" \
        --coverage_min 0.75
}



run_plasmidfinder() {
    local sample="$1"
    local sample_dir="$2"

    echo

    plasmidfinder.py \
        -i "$OUTPUT_DIR/Assembly_final/${sample}_${ASSEMBLER}.fasta" \
        -o "$sample_dir/${sample}_plasmidfinder" \
        -p "$PLASMIDFINDER_DB" \
        -mp blastn
}



run_pre_annotation_sample_pipeline() {
    local assembler="$1"
    local r1="$2"
    local r2="$3"
    local sample_dir="$4"
    local sample="$5"
    local sample_gsize="$6"

    echo
    create_directory "$sample_dir"
    create_directory "$sample_dir/${sample}_fastqc"

    if [[ "$TRIM" == "yes" ]]; then
        create_directory "$sample_dir/${sample}_trimmomatic"
        create_directory "$sample_dir/${sample}_fastqc_trimmed"

        run_step "${sample}_fastqc_raw" "FastQC (raw): $sample" \
            run_fastqc "$sample_dir/${sample}_fastqc" "$r1" "$r2"

        run_step "${sample}_trimmomatic" "Trimmomatic: $sample" \
            run_trimmomatic \
                "$r1" "$r2" \
                "$sample_dir/${sample}_trimmomatic/${sample}_1_paired.fastq.gz" \
                "$sample_dir/${sample}_trimmomatic/${sample}_1_unpaired.fastq.gz" \
                "$sample_dir/${sample}_trimmomatic/${sample}_2_paired.fastq.gz" \
                "$sample_dir/${sample}_trimmomatic/${sample}_2_unpaired.fastq.gz" \
                "$sample_dir/${sample}_trimmomatic/statsSummaryFile.txt"

        r1="$sample_dir/${sample}_trimmomatic/${sample}_1_paired.fastq.gz"
        r2="$sample_dir/${sample}_trimmomatic/${sample}_2_paired.fastq.gz"

        run_step "${sample}_fastqc_trimmed" "FastQC (trimmed): $sample" \
            run_fastqc "$sample_dir/${sample}_fastqc_trimmed" "$r1" "$r2"
    else
        run_step "${sample}_fastqc_raw" "FastQC: $sample" \
            run_fastqc "$sample_dir/${sample}_fastqc" "$r1" "$r2"
    fi

    case "$assembler" in
        spades)
            run_step "${sample}_spades" "SPAdes: $sample" \
                run_spades "$r1" "$r2" "$sample_dir/${sample}_spades"
            ;;
        shovill)
            run_step "${sample}_shovill" "Shovill: $sample" \
                run_shovill "$r1" "$r2" "$sample_dir/${sample}_shovill" "$sample" "$sample_gsize"
            ;;
        unicycler)
            run_step "${sample}_unicycler" "Unicycler: $sample" \
                run_unicycler "$r1" "$r2" "$sample_dir/${sample}_unicycler"
            ;;
    esac

    run_step "${sample}_rename_assembly" "Renomear e copiar assembly: $sample" \
        rename_and_copy_assembly "$assembler" "$sample_dir" "$sample"

    run_step "${sample}_quast" "QUAST: $sample" \
        run_quast "$sample" "$sample_dir"

    run_step "${sample}_assembly_scan" "Assembly-scan: $sample" \
        run_assembly_scan "$sample" "$sample_dir"
}



infer_taxonomy_from_predicted_name() {
    local predicted="$1"
    predicted="$(trim_spaces "$predicted")"
    predicted="${predicted//$'\r'/}"

    awk -v s="$predicted" '
        function is_bad_token(x) {
            x_l=tolower(x)
            return (x_l=="" || x_l=="sp." || x_l=="sp" || x_l=="spp." || x_l=="spp" || x_l=="bacterium" || x_l=="bacteria" || x_l=="cf." || x_l=="cf" || x_l=="aff." || x_l=="aff" || x_l=="complex" || x_l=="group" || x_l=="subsp." || x_l=="subsp")
        }
        BEGIN {
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", s)
            n=split(s, a, /[[:space:]]+/)
            genus=""
            species=""

            if (n >= 1) genus=a[1]
            if (n >= 2 && !is_bad_token(a[2])) species=a[2]

            if (tolower(genus)=="uncultured" || tolower(genus)=="unknown" || tolower(genus)=="unclassified") genus=""
            if (tolower(species)=="unknown" || tolower(species)=="unclassified") species=""

            print genus "\t" species
        }
    '
}



infer_taxonomy_for_sample_from_gambit() {
    local sample="$1"
    local gambit_csv="$OUTPUT_DIR/taxonomic_identification_gambit.tsv"

    [[ -f "$gambit_csv" ]] || return 0

    echo
    python3 - "$gambit_csv" "$sample" "$ASSEMBLER" <<'PY'
import csv
import os
import sys

gambit_csv = sys.argv[1]
sample = sys.argv[2]
assembler = sys.argv[3]

with open(gambit_csv, newline="", encoding="utf-8-sig") as fh:
    reader = csv.DictReader(fh)
    if reader.fieldnames is None:
        sys.exit(0)

    lower_map = {name.lower().strip(): name for name in reader.fieldnames}
    if "query" not in lower_map or "predicted.name" not in lower_map:
        sys.exit(0)

    query_key = lower_map["query"]
    pred_key = lower_map["predicted.name"]

    for row in reader:
        q = (row.get(query_key) or "").strip()
        p = (row.get(pred_key) or "").strip()

        q = os.path.basename(q)
        for ext in (".fasta", ".fa", ".fna"):
            if q.endswith(ext):
                q = q[:-len(ext)]

        if q == f"{sample}_{assembler}" or q == sample or q.startswith(f"{sample}_"):
            print(p)
            break
PY
}



resolve_taxonomy_for_sample() {
    local sample="$1"
    local metadata_genus="${SAMPLE_GENUS[$sample]}"
    local metadata_species="${SAMPLE_SPECIES[$sample]}"
    local global_genus="${GENUS:-}"
    local global_species="${SPECIES:-}"

    local final_genus=""
    local final_species=""
    local genus_source=""
    local species_source=""

    if [[ -n "$metadata_genus" ]]; then
        final_genus="$metadata_genus"
        genus_source="samplesheet"
    elif [[ -n "$global_genus" ]]; then
        final_genus="$global_genus"
        genus_source="global"
    fi

    if [[ -n "$metadata_species" ]]; then
        final_species="$metadata_species"
        species_source="samplesheet"
    elif [[ -n "$global_species" ]]; then
        final_species="$global_species"
        species_source="global"
    fi

    if [[ "$RUN_GAMBIT" == "yes" ]]; then
        if [[ -z "$final_genus" || -z "$final_species" ]]; then
            local predicted_name=""
            predicted_name="$(infer_taxonomy_for_sample_from_gambit "$sample" || true)"

            if [[ -n "$predicted_name" ]]; then
                local inferred inferred_genus inferred_species
                inferred="$(infer_taxonomy_from_predicted_name "$predicted_name")"
                inferred_genus="$(cut -f1 <<< "$inferred")"
                inferred_species="$(cut -f2 <<< "$inferred")"

                if [[ -z "$final_genus" && -n "$inferred_genus" ]]; then
                    final_genus="$inferred_genus"
                    genus_source="gambit_predicted.name"
                fi

                if [[ -z "$final_species" && -n "$inferred_species" ]]; then
                    final_species="$inferred_species"
                    species_source="gambit_predicted.name"
                fi
            fi
        fi
    fi

    SAMPLE_FINAL_GENUS["$sample"]="$final_genus"
    SAMPLE_FINAL_SPECIES["$sample"]="$final_species"
    SAMPLE_GENUS_SOURCE["$sample"]="${genus_source:-NA}"
    SAMPLE_SPECIES_SOURCE["$sample"]="${species_source:-NA}"
}



write_taxonomy_resolution_table() {
    local outfile="$OUTPUT_DIR/sample_taxonomy_resolution.tsv"
    {
        echo -e "sample\tgenus_input\tspecies_input\tgram_input\tgenus_final\tspecies_final\tgenus_source\tspecies_source"
        for sample in "${SAMPLE_ORDER[@]}"; do
            echo -e "${sample}\t${SAMPLE_GENUS[$sample]:-}\t${SAMPLE_SPECIES[$sample]:-}\t${SAMPLE_GRAMPHENO[$sample]:-}\t${SAMPLE_FINAL_GENUS[$sample]:-}\t${SAMPLE_FINAL_SPECIES[$sample]:-}\t${SAMPLE_GENUS_SOURCE[$sample]:-NA}\t${SAMPLE_SPECIES_SOURCE[$sample]:-NA}"
        done
    } > "$outfile"

    echo
    log "Resumo taxonômico salvo em: $outfile"
}



run_annotation_steps_for_sample() {
    local sample="$1"
    local sample_dir="$OUTPUT_DIR/$sample"
    local fasta="$OUTPUT_DIR/Assembly_final/${sample}_${ASSEMBLER}.fasta"

    local sample_genus="${SAMPLE_FINAL_GENUS[$sample]}"
    local sample_species="${SAMPLE_FINAL_SPECIES[$sample]}"
    local sample_gram="${SAMPLE_GRAMPHENO[$sample]}"

    if [[ -z "$sample_gram" && -n "$GRAM" && "$GRAM" != "?" ]]; then
        sample_gram="$GRAM"
    fi

    create_directory "$sample_dir/${sample}_prokka"
    run_step "${sample}_prokka" "Prokka: $sample" \
        run_prokka "$fasta" "$sample_dir/${sample}_prokka" "$sample" "$sample_genus" "$sample_species"

    create_directory "$sample_dir/${sample}_bakta"
    run_step "${sample}_bakta" "Bakta: $sample" \
        run_bakta "$fasta" "$sample_dir/${sample}_bakta" "$sample" "$sample_genus" "$sample_species" "$sample_gram"

    create_directory "$sample_dir/${sample}_amrfinder"
    run_step "${sample}_amrfinder" "AMRFinder: $sample" \
        run_amrfinder "$sample" "$sample_dir"

    create_directory "$sample_dir/${sample}_plasmidfinder"
    run_step "${sample}_plasmidfinder" "PlasmidFinder: $sample" \
        run_plasmidfinder "$sample" "$sample_dir"
}



run_gambit_on_assemblies() {

    echo
    gambit -d "$GAMBIT_DB" query \
        -o "$OUTPUT_DIR/taxonomic_identification_gambit.tsv" \
        "$OUTPUT_DIR/Assembly_final/"*.fasta
}

echo

run_multiqc() {
    if [[ "$TRIM" == "yes" ]]; then

    echo
        multiqc -o "$OUTPUT_DIR/multiqc_output" \
            "$OUTPUT_DIR"/*/*_fastqc \
            "$OUTPUT_DIR"/*/*_fastqc_trimmed
    else
        multiqc -o "$OUTPUT_DIR/multiqc_output" \
            "$OUTPUT_DIR"/*/*_fastqc
    fi
}

echo

run_gtdbtk() {

    echo
    log "Trocando para ambiente GTDB-Tk: $GTDBTK_ENV"
    conda deactivate || true
    eval "$(conda shell.bash hook)" || die "Falha ao reinicializar o conda shell hook"
    conda activate "$GTDBTK_ENV" || die "Não foi possível ativar o ambiente conda: $GTDBTK_ENV"

    check_command gtdbtk

    local gtdbtk_scratch="$OUTPUT_DIR/.gtdbtk_scratch"
    mkdir -p "$gtdbtk_scratch"

    echo
    log "GTDB-Tk"
    gtdbtk classify_wf \
        --skip_ani_screen \
        --genome_dir "$OUTPUT_DIR/Assembly_final/" \
        --out_dir "$OUTPUT_DIR/taxonomic_identification_gtdbtk" \
        --extension fasta \
        --cpus "$THREADS_GTDBTK" \
        --pplacer_cpus 1 \
        --scratch_dir "$gtdbtk_scratch"

    rm -rf "$gtdbtk_scratch"
    conda deactivate || true

    eval "$(conda shell.bash hook)" || die "Falha ao reativar conda shell hook"
    conda activate "$CONDA_ENV" || die "Não foi possível reativar o ambiente conda principal: $CONDA_ENV"
}



# =========================
# SAMPLE LOADING
# =========================
if [[ -n "$SAMPLESHEET" ]]; then
    load_samplesheet
else
    discover_samples_from_input
fi

validate_shovill_gsizes

log "Total de amostras carregadas: ${#SAMPLE_ORDER[@]}"


# =========================
# PHASE 1: QC + ASSEMBLY + QUAST + ASSEMBLY-SCAN
# =========================
for sample in "${SAMPLE_ORDER[@]}"; do
    sample_dir="$OUTPUT_DIR/$sample"
    log "Processando amostra (fase 1): $sample"
    run_pre_annotation_sample_pipeline \
        "$ASSEMBLER" \
        "${SAMPLE_R1[$sample]}" \
        "${SAMPLE_R2[$sample]}" \
        "$sample_dir" \
        "$sample" \
        "${SAMPLE_GSIZE[$sample]}"
done



# =========================
# PHASE 2: GAMBIT + RESOLUÇÃO TAXONÔMICA
# =========================

echo
if [[ "$RUN_GAMBIT" == "yes" ]]; then
    run_step "global_gambit" "Gambit" run_gambit_on_assemblies
else
    log "Gambit desabilitado: inferência de genus/species não será executada"
fi

echo
for sample in "${SAMPLE_ORDER[@]}"; do
    resolve_taxonomy_for_sample "$sample"
    log "Taxonomia final | $sample | genus=${SAMPLE_FINAL_GENUS[$sample]:-NA} (${SAMPLE_GENUS_SOURCE[$sample]:-NA}) | species=${SAMPLE_FINAL_SPECIES[$sample]:-NA} (${SAMPLE_SPECIES_SOURCE[$sample]:-NA}) | gram=${SAMPLE_GRAMPHENO[$sample]:-NA}"
done

echo
run_step "global_taxonomy_resolution" "Gerar tabela de resolução taxonômica" write_taxonomy_resolution_table



# =========================
# PHASE 3: ANNOTATION + AMR + PLASMIDS
# =========================

echo
for sample in "${SAMPLE_ORDER[@]}"; do
    log "Processando amostra (fase 3): $sample"
    run_annotation_steps_for_sample "$sample"
done



# =========================
# POST-PROCESSING
# =========================

echo
if [[ "$RUN_MULTIQC" == "yes" ]]; then
    run_step "global_multiqc" "MultiQC" run_multiqc
fi

if [[ "$RUN_GTDBTK" == "yes" ]]; then
    if step_done "global_gtdbtk"; then
        log "SKIP (já concluído): GTDB-Tk"
    else
        run_gtdbtk
        mark_done "global_gtdbtk"
    fi
fi

echo
echo "########################################################################"
echo "##### Script desenvolvido por: Jean Phellipe Marques do Nascimento #####"
echo "############# Laboratório de Vigilância Genômica - LACEN-AL ############"
echo "########################################################################"
echo

log "Pipeline finalizado com sucesso"
