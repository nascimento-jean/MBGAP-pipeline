# MBGAP — Mykaella's Bacterial Genome Analysis Pipeline

**Versão 3.0**  
Desenvolvido por **Jean Phellipe Marques do Nascimento**  
Laboratório de Vigilância Genômica — LACEN-AL

---

## O que é o MBGAP?

O MBGAP é uma pipeline automatizada para análise de genomas bacterianos a partir de dados de sequenciamento de nova geração (NGS) no formato **paired-end Illumina**. Com um único comando, ele executa de forma integrada:

1. **Controle de qualidade** das reads brutas (FastQC)
2. **Trimagem** de adaptadores e baixa qualidade (Trimmomatic) — opcional
3. **Montagem** do genoma (*assembly*) — SPAdes, Shovill ou Unicycler
4. **Avaliação da montagem** (QUAST + Assembly-scan)
5. **Identificação taxonômica** (GAMBIT e/ou GTDB-Tk)
6. **Anotação genômica** (Prokka + Bakta)
7. **Detecção de genes de resistência antimicrobiana** (AMRFinder)
8. **Identificação de plasmídeos** (PlasmidFinder)
9. **Relatório consolidado de qualidade** (MultiQC)

---

## Antes de começar — Pré-requisitos

### Sistema operacional

O MBGAP foi desenvolvido para rodar em **Linux** (Ubuntu/Debian recomendado). Não é compatível com Windows diretamente (você pode usar WSL2 no Windows, se necessário).

### Conda

Todas as ferramentas são gerenciadas via **Conda**. Se você ainda não tem o Conda instalado:

1. Baixe o instalador do Miniconda:
   ```bash
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   ```
2. Reinicie o terminal ou execute `source ~/.bashrc`.

### Ambientes Conda necessários

O MBGAP utiliza **dois ambientes Conda**:

| Ambiente | Ferramentas incluídas |
|---|---|
| `bioinfo` (padrão) | FastQC, Trimmomatic, SPAdes, Shovill, Unicycler, QUAST, Prokka, Bakta, AMRFinder, PlasmidFinder, GAMBIT, assembly-scan, MultiQC |
| `gtdbtk` (padrão) | GTDB-Tk |

Observação: Para criação dos ambientes Conda necessários, pode ser utilizado o repositório: https://github.com/nascimento-jean/Criacao_Env_Bioinfo

> O GTDB-Tk requer um ambiente separado pois possui dependências conflitantes com as demais ferramentas.

Você pode usar nomes diferentes para os ambientes — basta informá-los com `--conda-env` e `--gtdbtk-env`.

**Sugestão de instalação do ambiente principal** (adapte conforme necessário):
```bash
conda create -n bioinfo -c conda-forge -c bioconda \
  fastqc trimmomatic spades shovill unicycler quast \
  prokka bakta amrfinder plasmidfinder gambit assembly-scan multiqc
```

### Bancos de dados necessários

Você precisa baixar e configurar **três bancos de dados** antes de usar o MBGAP:

| Banco | Ferramenta | Como obter |
|---|---|---|
| GAMBIT DB | GAMBIT | [github.com/jlumpe/gambit](https://github.com/jlumpe/gambit) |
| Bakta DB | Bakta | `bakta_db download --output /caminho/bakta_db` |
| PlasmidFinder DB | PlasmidFinder | [bitbucket.org/genomicepidemiology/plasmidfinder_db](https://bitbucket.org/genomicepidemiology/plasmidfinder_db) |

---

## Instalação do MBGAP

1. Baixe o script:
   ```bash
   # Via git clone
   git clone https://github.com/SEU_USUARIO/MBGAP.git
   cd MBGAP
   ```
   Ou faça o download direto do arquivo `MBGAP_v3_0.sh`.

2. Torne o script executável:
   ```bash
   chmod +x MBGAP_v3_0.sh
   ```

3. Pronto! Não há instalação adicional além dos ambientes Conda descritos acima.

---

## Como usar

### Modos de entrada

Você pode fornecer os dados de duas formas:

**Opção A — Pasta com arquivos FASTQ**

Se seus arquivos seguem o padrão `AMOSTRA_R1_001.fastq.gz` / `AMOSTRA_R2_001.fastq.gz` (ou `AMOSTRA_1.fastq.gz` / `AMOSTRA_2.fastq.gz`):

```bash
bash MBGAP_v3_0.sh \
  --input /caminho/para/fastqs/ \
  --output /caminho/para/resultados/ \
  --assembler spades \
  --gambit-db /caminho/gambit_db/ \
  --bakta-db /caminho/bakta_db/ \
  --plasmidfinder-db /caminho/plasmidfinder_db/
```

**Opção B — Samplesheet CSV (recomendado para múltiplas amostras com metadados)**

Crie um arquivo CSV com as informações de cada amostra:

```csv
sample,r1,r2,genus,species,gram,gsize
ISO001,/dados/ISO001_R1.fastq.gz,/dados/ISO001_R2.fastq.gz,Escherichia,coli,-,5.2M
ISO002,/dados/ISO002_R1.fastq.gz,/dados/ISO002_R2.fastq.gz,Klebsiella,pneumoniae,-,5.5M
ISO003,/dados/ISO003_R1.fastq.gz,/dados/ISO003_R2.fastq.gz,,,,4.9M
```
Observação: A samplesheet poderá ser criada utilizando o arquivo: create_samplesheet.py

> **Campos obrigatórios:** `sample`, `r1`, `r2`  
> **Campos opcionais:** `genus`, `species`, `gram` (`+`, `-` ou `?`), `gsize` (necessário para Shovill)  
> Deixe campos opcionais em branco quando não souber — o GAMBIT pode inferir automaticamente.

```bash
bash MBGAP_v3_0.sh \
  --samplesheet /caminho/amostras.csv \
  --output /caminho/para/resultados/ \
  --assembler spades \
  --gambit-db /caminho/gambit_db/ \
  --bakta-db /caminho/bakta_db/ \
  --plasmidfinder-db /caminho/plasmidfinder_db/
```

---

## Parâmetros completos

### Entrada e saída

| Parâmetro | Descrição |
|---|---|
| `--input DIR` | Pasta contendo arquivos FASTQ paired-end |
| `--samplesheet FILE` | CSV com metadados por amostra (substitui `--input`) |
| `--output DIR` | **Obrigatório.** Pasta onde os resultados serão salvos |

### Bancos de dados (todos obrigatórios)

| Parâmetro | Descrição |
|---|---|
| `--gambit-db DIR` | Caminho para o banco de dados do GAMBIT |
| `--bakta-db DIR` | Caminho para o banco de dados do Bakta |
| `--plasmidfinder-db DIR` | Caminho para o banco de dados do PlasmidFinder |

### Montagem (*assembly*)

| Parâmetro | Valores | Padrão | Descrição |
|---|---|---|---|
| `--assembler` | `spades`, `shovill`, `unicycler` | `spades` | Programa de montagem a ser usado |
| `--gsize` | ex: `5.0M`, `4500k` | — | Tamanho estimado do genoma (obrigatório para Shovill quando não informado no CSV) |
| `--memory-spades` | número (GB) | `50` | Memória máxima para o SPAdes |
| `--ram-shovill` | número (GB) | `50` | RAM máxima para o Shovill |

### Qualidade e trimagem

| Parâmetro | Valores | Padrão | Descrição |
|---|---|---|---|
| `--trim` | `yes`, `no` | `no` | Ativar trimagem com Trimmomatic |
| `--adapters FILE` | caminho para o arquivo | — | Arquivo de adaptadores do Trimmomatic (**obrigatório** se `--trim yes`) |

### Informações taxonômicas (opcionais globais)

| Parâmetro | Valores | Padrão | Descrição |
|---|---|---|---|
| `--genus` | texto | — | Gênero bacteriano global para todas as amostras |
| `--species` | texto | — | Espécie bacteriana global para todas as amostras |
| `--gram` | `+`, `-`, `?` | `?` | Tipo de coloração de Gram |
| `--genetic-code` | número | `11` | Código genético para o Prokka (11 = bactérias) |
| `--use-prokka-genus` | `yes`, `no` | `no` | Ativar `--usegenus` no Prokka |
| `--bakta-complete` | `yes`, `no` | `no` | Informar ao Bakta que o genoma é completo |

### Ferramentas opcionais

| Parâmetro | Valores | Padrão | Descrição |
|---|---|---|---|
| `--run-gambit` | `yes`, `no` | `yes` | Executar identificação por GAMBIT |
| `--run-gtdbtk` | `yes`, `no` | `yes` | Executar classificação por GTDB-Tk |
| `--run-multiqc` | `yes`, `no` | `yes` | Gerar relatório MultiQC |

### Ambientes Conda

| Parâmetro | Padrão | Descrição |
|---|---|---|
| `--conda-env` | `bioinfo` | Nome do ambiente Conda principal |
| `--gtdbtk-env` | `gtdbtk` | Nome do ambiente Conda do GTDB-Tk |

### Threads (paralelismo)

| Parâmetro | Padrão | Descrição |
|---|---|---|
| `--threads INT` | — | Define todas as threads de uma vez |
| `--threads-fastqc INT` | `8` | Threads para o FastQC |
| `--threads-trimmomatic INT` | `8` | Threads para o Trimmomatic |
| `--threads-spades INT` | `10` | Threads para o SPAdes |
| `--threads-shovill INT` | `10` | Threads para o Shovill |
| `--threads-unicycler INT` | `10` | Threads para o Unicycler |
| `--threads-quast INT` | `8` | Threads para o QUAST |
| `--threads-prokka INT` | `8` | Threads para o Prokka |
| `--threads-bakta INT` | `12` | Threads para o Bakta |
| `--threads-amrfinder INT` | `8` | Threads para o AMRFinder |
| `--threads-gtdbtk INT` | `10` | Threads para o GTDB-Tk |

---

## Estrutura dos resultados

Após a execução, a pasta de saída (`--output`) terá a seguinte estrutura:

```
resultados/
├── pipeline_log.txt                        # Log completo da execução
├── sample_taxonomy_resolution.tsv          # Resumo da identificação taxonômica por amostra
├── taxonomic_identification_gambit.tsv     # Resultado do GAMBIT
├── taxonomic_identification_gtdbtk/        # Resultados do GTDB-Tk
├── multiqc_output/                         # Relatório HTML consolidado (MultiQC)
├── Assembly_final/                         # FASTAs finais de todas as amostras
│   ├── ISO001_spades.fasta
│   └── ISO002_spades.fasta
├── ISO001/                                 # Pasta por amostra
│   ├── ISO001_fastqc/                      # FastQC das reads brutas
│   ├── ISO001_fastqc_trimmed/              # FastQC após trimagem (se --trim yes)
│   ├── ISO001_trimmomatic/                 # Reads trimadas
│   ├── ISO001_spades/                      # Arquivos brutos do SPAdes
│   ├── ISO001_quast/                       # Estatísticas de montagem
│   ├── ISO001_spades.tsv                   # Resumo do assembly-scan
│   ├── ISO001_prokka/                      # Anotação Prokka
│   ├── ISO001_bakta/                       # Anotação Bakta
│   ├── ISO001_amrfinder/                   # Genes de resistência (AMRFinder)
│   └── ISO001_plasmidfinder/               # Plasmídeos identificados
└── ISO002/
    └── ...
```

---

## Exemplos de uso

### Exemplo 1 — Análise simples (sem trimagem, SPAdes, todos os padrões)

```bash
bash MBGAP_v3_0.sh \
  --input /dados/fastqs/ \
  --output /resultados/analise_01/ \
  --assembler spades \
  --gambit-db /databases/gambit/ \
  --bakta-db /databases/bakta/ \
  --plasmidfinder-db /databases/plasmidfinder/
```

### Exemplo 2 — Com trimagem e genus/species globais

```bash
bash MBGAP_v3_0.sh \
  --input /dados/fastqs/ \
  --output /resultados/analise_02/ \
  --assembler spades \
  --trim yes \
  --adapters /databases/TruSeq3-PE.fa \
  --genus Staphylococcus \
  --species aureus \
  --gram + \
  --gambit-db /databases/gambit/ \
  --bakta-db /databases/bakta/ \
  --plasmidfinder-db /databases/plasmidfinder/
```

### Exemplo 3 — Samplesheet com Shovill e 20 threads

```bash
bash MBGAP_v3_0.sh \
  --samplesheet /dados/minhas_amostras.csv \
  --output /resultados/analise_03/ \
  --assembler shovill \
  --threads 20 \
  --gambit-db /databases/gambit/ \
  --bakta-db /databases/bakta/ \
  --plasmidfinder-db /databases/plasmidfinder/
```

### Exemplo 4 — Sem GTDB-Tk (mais rápido, sem o banco GTDB)

```bash
bash MBGAP_v3_0.sh \
  --samplesheet /dados/amostras.csv \
  --output /resultados/analise_04/ \
  --assembler spades \
  --run-gtdbtk no \
  --gambit-db /databases/gambit/ \
  --bakta-db /databases/bakta/ \
  --plasmidfinder-db /databases/plasmidfinder/
```

---

## Retomada automática (*checkpoints*)

O MBGAP salva o progresso de cada etapa. Se a execução for interrompida (queda de energia, erro, etc.), basta rodar o mesmo comando novamente — as etapas já concluídas serão puladas automaticamente.

Os arquivos de checkpoint ficam em `OUTPUT_DIR/.checkpoints/`.

---

## Resolução de problemas comuns

**Erro: `Dependência não encontrada no PATH: fastqc`**  
→ O ambiente Conda não está ativado ou a ferramenta não foi instalada. Verifique com `conda activate bioinfo` e depois `which fastqc`.

**Erro: `Não foi possível ativar o ambiente conda: bioinfo`**  
→ O nome do seu ambiente é diferente. Verifique com `conda env list` e use `--conda-env NOME_DO_SEU_AMBIENTE`.

**Erro: `gsize ausente para amostra X usando Shovill`**  
→ Ao usar `--assembler shovill`, é obrigatório informar o tamanho estimado do genoma. Adicione a coluna `gsize` no CSV ou use `--gsize 5.0M`.

**Erro: `--adapters é obrigatório quando --trim yes`**  
→ Forneça o caminho para o arquivo de adaptadores Illumina com `--adapters /caminho/TruSeq3-PE.fa`.

**Erro: `Amostra duplicada detectada`**  
→ Dois arquivos FASTQ geraram o mesmo nome de amostra ou há linhas duplicadas no CSV.

**O pipeline está lento**  
→ Aumente o número de threads com `--threads 16` (ou mais, conforme seu servidor) e, se necessário, aumente `--memory-spades` e `--ram-shovill`.

---

## Citação / Créditos

Se você utilizar o MBGAP em publicações científicas, por favor cite o repositório e as ferramentas individuais utilizadas:

- **SPAdes:** Bankevich et al., 2012 — *J Comput Biol*
- **Shovill:** [github.com/tseemann/shovill](https://github.com/tseemann/shovill)
- **Unicycler:** Wick et al., 2017 — *PLOS Computational Biology*
- **FastQC:** Andrews, S. — [bioinformatics.babraham.ac.uk](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- **Trimmomatic:** Bolger et al., 2014 — *Bioinformatics*
- **QUAST:** Gurevich et al., 2013 — *Bioinformatics*
- **Prokka:** Seemann, 2014 — *Bioinformatics*
- **Bakta:** Schwengers et al., 2021 — *Microbial Genomics*
- **AMRFinder:** Feldgarden et al., 2021 — *Scientific Reports*
- **PlasmidFinder:** Carattoli et al., 2014 — *Antimicrobial Agents and Chemotherapy*
- **GAMBIT:** Lumpe & Bhatt — [github.com/jlumpe/gambit](https://github.com/jlumpe/gambit)
- **GTDB-Tk:** Chaumeil et al., 2019 — *Bioinformatics*
- **MultiQC:** Ewels et al., 2016 — *Bioinformatics*

---

## Licença

Este projeto está disponível para uso acadêmico e científico.  
Desenvolvido no **Laboratório de Vigilância Genômica — LACEN-AL**.

---

*Dúvidas ou sugestões? Abra uma [issue](../../issues) neste repositório.*
