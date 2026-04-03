#!/usr/bin/env python3
"""
create_samplesheet.py

Gera uma samplesheet CSV para o MBGAP a partir de um diretório de reads paired-end.

Suporta os padrões:
- *_R1_001.fastq.gz  <-> *_R2_001.fastq.gz
- *_1.fastq.gz       <-> *_2.fastq.gz

Saída:
sample,r1,r2,genus,species,gram,gsize

As colunas genus, species, gram e gsize são deixadas vazias para preenchimento manual.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


SUPPORTED_PATTERNS = (
    ("_R1_001.fastq.gz", "_R2_001.fastq.gz"),
    ("_1.fastq.gz", "_2.fastq.gz"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Criar samplesheet CSV a partir de um diretório de reads paired-end."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Diretório contendo os arquivos FASTQ(.gz).",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Arquivo CSV de saída, por exemplo: teste.csv",
    )
    return parser.parse_args()


def find_samples(input_dir: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    seen_samples: set[str] = set()

    for r1_suffix, r2_suffix in SUPPORTED_PATTERNS:
        for r1_path in sorted(input_dir.glob(f"*{r1_suffix}")):
            sample = r1_path.name[: -len(r1_suffix)]

            if sample in seen_samples:
                continue

            r2_path = input_dir / f"{sample}{r2_suffix}"
            if not r2_path.exists():
                print(
                    f"Aviso: par R2 ausente para a amostra '{sample}'. Arquivo esperado: {r2_path}",
                    file=sys.stderr,
                )
                continue

            seen_samples.add(sample)
            rows.append(
                {
                    "sample": sample,
                    "r1": str(r1_path.resolve()),
                    "r2": str(r2_path.resolve()),
                    "genus": "",
                    "species": "",
                    "gram": "",
                    "gsize": "",
                }
            )

    return sorted(rows, key=lambda x: x["sample"])


def write_csv(rows: list[dict[str, str]], output_file: Path) -> None:
    output_file.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = ["sample", "r1", "r2", "genus", "species", "gram", "gsize"]

    with output_file.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()

    input_dir = Path(args.input).expanduser().resolve()
    output_file = Path(args.output).expanduser().resolve()

    if not input_dir.exists():
        print(f"Erro: diretório não encontrado: {input_dir}", file=sys.stderr)
        sys.exit(1)

    if not input_dir.is_dir():
        print(f"Erro: --input não é um diretório: {input_dir}", file=sys.stderr)
        sys.exit(1)

    rows = find_samples(input_dir)

    if not rows:
        print(
            "Erro: nenhuma amostra paired-end válida foi encontrada.\n"
            "Padrões suportados:\n"
            "  *_R1_001.fastq.gz <-> *_R2_001.fastq.gz\n"
            "  *_1.fastq.gz      <-> *_2.fastq.gz",
            file=sys.stderr,
        )
        sys.exit(1)

    write_csv(rows, output_file)

    print(f"Samplesheet criada com sucesso: {output_file}")
    print(f"Total de amostras: {len(rows)}")


if __name__ == "__main__":
    main()
