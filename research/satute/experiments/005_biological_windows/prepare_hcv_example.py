#!/usr/bin/env python3
"""Prepare a codon-preserving HCV genotype reference alignment.

The input is a concatenated GenBank flat file containing the ten accession
versions listed in ``REFERENCES``.  Polyprotein CDS translations are aligned
with MAFFT, nucleotide CDSs are threaded through that amino-acid alignment,
and H77 mat_peptide coordinates define homologous region alignments.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path


REFERENCES = {
    "NC_004102.1": ("HCV_1a_H77", "1a"),
    "AJ238799.1": ("HCV_1b_Con1", "1b"),
    "AB047639.1": ("HCV_2a_JFH1", "2a"),
    "D10988.1": ("HCV_2b_HCJ8", "2b"),
    "D17763.1": ("HCV_3a_NZL1", "3a"),
    "D49374.1": ("HCV_3b_TrKj", "3b"),
    "Y11604.1": ("HCV_4a_ED43", "4a"),
    "AF064490.1": ("HCV_5a_SA13", "5a"),
    "Y12083.1": ("HCV_6a_EUHK2", "6a"),
    "EF108306.2": ("HCV_7a_QC69", "7a"),
}

# Zero-based, half-open amino-acid intervals in the H77 polyprotein.  These
# derive from NC_004102.1 mat_peptide coordinates and exclude the stop codon.
H77_REGIONS = {
    "core": (0, 191),
    "e1": (191, 383),
    "e2": (383, 746),
    "p7": (746, 809),
    "ns2": (809, 1026),
    "ns3": (1026, 1657),
    "ns4a": (1657, 1711),
    "ns4b": (1711, 1972),
    "ns5a": (1972, 2420),
    "ns5b": (2420, 3011),
}


@dataclass(frozen=True)
class Record:
    accession: str
    label: str
    subtype: str
    cds_start: int
    cds_end: int
    cds: str
    translation: str


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fasta_records(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = ""
    chunks: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        if raw.startswith(">"):
            if name:
                records[name] = "".join(chunks)
            name = raw[1:].split()[0]
            chunks = []
        else:
            chunks.append(raw.strip())
    if name:
        records[name] = "".join(chunks)
    return records


def write_fasta(records: list[tuple[str, str]], path: Path, width: int = 80) -> None:
    lines: list[str] = []
    for name, sequence in records:
        lines.append(f">{name}")
        lines.extend(sequence[i : i + width] for i in range(0, len(sequence), width))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def first_cds_feature(record_text: str) -> tuple[int, int, str]:
    match = re.search(r"^     CDS\s+(\d+)\.\.(\d+)\s*$", record_text, re.MULTILINE)
    if not match:
        raise ValueError("A simple primary CDS feature was not found")
    start, end = int(match.group(1)), int(match.group(2))
    feature_start = match.start()
    next_feature = re.search(r"^     \S", record_text[match.end() :], re.MULTILINE)
    feature_end = match.end() + next_feature.start() if next_feature else len(record_text)
    feature = record_text[feature_start:feature_end]
    translation_match = re.search(r'/translation="([A-Z\s*]+)"', feature)
    if not translation_match:
        raise ValueError("The primary CDS has no translation qualifier")
    translation = re.sub(r"\s+", "", translation_match.group(1)).rstrip("*")
    return start, end, translation


def parse_genbank(path: Path) -> list[Record]:
    parsed: dict[str, Record] = {}
    for raw_record in path.read_text(encoding="utf-8").split("\n//"):
        version_match = re.search(r"^VERSION\s+(\S+)", raw_record, re.MULTILINE)
        if not version_match:
            continue
        accession = version_match.group(1)
        if accession not in REFERENCES:
            continue
        origin_match = re.search(r"^ORIGIN\s*\n(.*)$", raw_record, re.MULTILINE | re.DOTALL)
        if not origin_match:
            raise ValueError(f"No ORIGIN sequence in {accession}")
        genome = re.sub(r"[^A-Za-z]", "", origin_match.group(1)).upper()
        start, end, translation = first_cds_feature(raw_record)
        cds = genome[start - 1 : end]
        if len(cds) % 3:
            raise ValueError(f"CDS length is not divisible by three for {accession}")
        if len(cds) // 3 not in {len(translation), len(translation) + 1}:
            raise ValueError(
                f"CDS/translation length mismatch for {accession}: "
                f"{len(cds) // 3} codons versus {len(translation)} residues"
            )
        cds = cds[: 3 * len(translation)]
        label, subtype = REFERENCES[accession]
        parsed[accession] = Record(accession, label, subtype, start, end, cds, translation)

    missing = sorted(set(REFERENCES) - set(parsed))
    if missing:
        raise ValueError(f"Missing expected accession versions: {missing}")
    return [parsed[accession] for accession in REFERENCES]


def aligned_codons(record: Record, aligned_protein: str) -> str:
    codons = [record.cds[i : i + 3] for i in range(0, len(record.cds), 3)]
    output: list[str] = []
    index = 0
    for residue in aligned_protein:
        if residue == "-":
            output.append("---")
        else:
            if index >= len(codons):
                raise ValueError(f"Protein alignment overruns CDS for {record.accession}")
            output.append(codons[index])
            index += 1
    if index != len(codons):
        raise ValueError(f"Protein alignment used {index}/{len(codons)} codons for {record.accession}")
    return "".join(output)


def h77_region_columns(aligned_h77: str) -> dict[str, list[int]]:
    columns = {region: [] for region in H77_REGIONS}
    residue_index = -1
    for column, residue in enumerate(aligned_h77):
        if residue != "-":
            residue_index += 1
        for region, (start, end) in H77_REGIONS.items():
            if residue != "-" and start <= residue_index < end:
                columns[region].append(column)
    for region, (start, end) in H77_REGIONS.items():
        if len(columns[region]) != end - start:
            raise ValueError(f"Unexpected H77 column count for {region}")
    return columns


def write_partition_windows(
    records: list[tuple[str, str]],
    partition: str,
    aligned_h77: str,
    output_root: Path,
    window_size: int,
    window_step: int,
) -> list[dict[str, str | int]]:
    sequence_length = len(records[0][1])
    if any(len(sequence) != sequence_length for _label, sequence in records):
        raise ValueError(f"Inconsistent sequence lengths for {partition}")
    divisor = 2 if partition == "pos12" else 1
    if window_size % divisor or window_step % divisor:
        raise ValueError(f"Window size and step must be divisible by {divisor} for {partition}")
    output_root.mkdir(parents=True, exist_ok=True)
    h77_position_by_column: list[int] = []
    position = 0
    for residue in aligned_h77:
        if residue != "-":
            position += 1
        h77_position_by_column.append(position)

    rows: list[dict[str, str | int]] = []
    for start in range(0, sequence_length - window_size + 1, window_step):
        end = start + window_size
        column_start = start // divisor
        column_end = end // divisor
        midpoint_column = (column_start + column_end - 1) // 2
        filename = f"hcv_{partition}_w{start + 1:04d}_{end:04d}.fasta"
        write_fasta(
            [(label, sequence[start:end]) for label, sequence in records],
            output_root / filename,
        )
        rows.append(
            {
                "partition": partition,
                "start_site": start + 1,
                "end_site": end,
                "aligned_codon_column_start": column_start + 1,
                "aligned_codon_column_end": column_end,
                "h77_polyprotein_midpoint": h77_position_by_column[midpoint_column],
                "filename": filename,
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("genbank", type=Path)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--mafft", default="mafft")
    parser.add_argument("--write-windows", action="store_true")
    parser.add_argument("--window-size", type=int, default=36)
    parser.add_argument("--window-step", type=int, default=12)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    records = parse_genbank(args.genbank)
    proteins = args.outdir / "hcv_polyprotein_unaligned.fasta"
    aligned_proteins = args.outdir / "hcv_polyprotein_aa.fasta"
    aligned_cds = args.outdir / "hcv_polyprotein_codon.fasta"
    aligned_cds_pos12 = args.outdir / "hcv_polyprotein_pos12.fasta"
    aligned_cds_pos3 = args.outdir / "hcv_polyprotein_pos3.fasta"
    write_fasta([(record.label, record.translation) for record in records], proteins)

    with aligned_proteins.open("w", encoding="utf-8") as handle:
        subprocess.run(
            [args.mafft, "--auto", "--reorder", str(proteins)],
            stdout=handle,
            check=True,
        )
    protein_alignment = fasta_records(aligned_proteins)
    expected_labels = {record.label for record in records}
    if set(protein_alignment) != expected_labels:
        raise ValueError("MAFFT output labels differ from the input labels")

    record_by_label = {record.label: record for record in records}
    codon_alignment = {
        label: aligned_codons(record_by_label[label], sequence)
        for label, sequence in protein_alignment.items()
    }
    lengths = {len(sequence) for sequence in codon_alignment.values()}
    if len(lengths) != 1:
        raise ValueError("Codon alignment lengths are inconsistent")
    write_fasta(list(codon_alignment.items()), aligned_cds)
    write_fasta(
        [
            (label, "".join(sequence[i : i + 2] for i in range(0, len(sequence), 3)))
            for label, sequence in codon_alignment.items()
        ],
        aligned_cds_pos12,
    )
    write_fasta(
        [
            (label, "".join(sequence[i + 2] for i in range(0, len(sequence), 3)))
            for label, sequence in codon_alignment.items()
        ],
        aligned_cds_pos3,
    )
    pos12_records = fasta_records(aligned_cds_pos12)
    pos3_records = fasta_records(aligned_cds_pos3)

    columns = h77_region_columns(protein_alignment["HCV_1a_H77"])
    region_outputs: dict[str, dict[str, str | int]] = {}
    for region, selected in columns.items():
        aa_path = args.outdir / f"hcv_{region}_aa.fasta"
        nt_path = args.outdir / f"hcv_{region}_codon.fasta"
        pos12_path = args.outdir / f"hcv_{region}_pos12.fasta"
        pos3_path = args.outdir / f"hcv_{region}_pos3.fasta"
        aa_records = [
            (label, "".join(sequence[column] for column in selected))
            for label, sequence in protein_alignment.items()
        ]
        nt_records = [
            (
                label,
                "".join(codon_alignment[label][3 * column : 3 * column + 3] for column in selected),
            )
            for label in protein_alignment
        ]
        write_fasta(aa_records, aa_path)
        write_fasta(nt_records, nt_path)
        write_fasta(
            [
                (label, "".join(sequence[i : i + 2] for i in range(0, len(sequence), 3)))
                for label, sequence in nt_records
            ],
            pos12_path,
        )
        write_fasta(
            [
                (label, "".join(sequence[i + 2] for i in range(0, len(sequence), 3)))
                for label, sequence in nt_records
            ],
            pos3_path,
        )
        region_outputs[region] = {
            "h77_residues": len(selected),
            "amino_acid_alignment": aa_path.name,
            "codon_alignment": nt_path.name,
            "codon_sites": 3 * len(selected),
            "codon_positions_1_2_alignment": pos12_path.name,
            "codon_positions_1_2_sites": 2 * len(selected),
            "codon_position_3_alignment": pos3_path.name,
            "codon_position_3_sites": len(selected),
        }

    manifest = {
        "source": "NCBI Nucleotide GenBank flat files fetched by accession",
        "source_genbank": str(args.genbank),
        "source_genbank_sha256": sha256(args.genbank),
        "coordinate_reference": "NC_004102.1 (H77) mat_peptide features",
        "accessions": [
            {
                "accession": record.accession,
                "label": record.label,
                "subtype": record.subtype,
                "cds_start": record.cds_start,
                "cds_end": record.cds_end,
                "polyprotein_residues": len(record.translation),
            }
            for record in records
        ],
        "mafft_command": [args.mafft, "--auto", "--reorder", str(proteins)],
        "polyprotein_amino_acid_alignment": aligned_proteins.name,
        "polyprotein_codon_alignment": aligned_cds.name,
        "polyprotein_codon_sites": next(iter(lengths)),
        "polyprotein_codon_positions_1_2_alignment": aligned_cds_pos12.name,
        "polyprotein_codon_position_3_alignment": aligned_cds_pos3.name,
        "regions": region_outputs,
    }
    if args.write_windows:
        window_root = args.outdir / "windows"
        window_rows = []
        window_rows.extend(
            write_partition_windows(
                list(pos12_records.items()),
                "pos12",
                protein_alignment["HCV_1a_H77"],
                window_root / "pos12",
                args.window_size,
                args.window_step,
            )
        )
        window_rows.extend(
            write_partition_windows(
                list(pos3_records.items()),
                "pos3",
                protein_alignment["HCV_1a_H77"],
                window_root / "pos3",
                args.window_size,
                args.window_step,
            )
        )
        fieldnames = list(window_rows[0])
        with (args.outdir / "window_manifest.tsv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(window_rows)
        manifest["windows"] = {
            "size_sites": args.window_size,
            "step_sites": args.window_step,
            "manifest": "window_manifest.tsv",
            "counts": {
                partition: sum(1 for row in window_rows if row["partition"] == partition)
                for partition in ("pos12", "pos3")
            },
        }
    (args.outdir / "preparation_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Wrote {args.outdir / 'preparation_manifest.json'}")


if __name__ == "__main__":
    main()
