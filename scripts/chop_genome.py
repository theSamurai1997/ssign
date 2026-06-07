#!/usr/bin/env python3
"""Chop a single-contig GenBank record into N equal-size artificial contigs.

Controlled test for the "fragmentation reduces substrate detection" hypothesis.
Take the complete K-12 chromosome and slice it into N pieces of equal length.
Features fully contained in a chunk get their coordinates shifted to be
chunk-local; features that span chunk boundaries are dropped (mirrors what
a real WGS assembler would do at a contig break).

Usage:
    python chop_genome.py <input.gbff> <output.gbff> <n_contigs>
"""

from __future__ import annotations

import sys

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


def chop(input_path: str, output_path: str, n_contigs: int) -> None:
    records = list(SeqIO.parse(input_path, "genbank"))
    if len(records) != 1:
        raise SystemExit(f"Expected exactly 1 GenBank record, got {len(records)}")
    src = records[0]
    seq = src.seq
    total = len(seq)
    chunk_size = total // n_contigs

    chunks: list[SeqRecord] = []
    cds_total = sum(1 for f in src.features if f.type == "CDS")
    cds_kept = 0
    cds_dropped = 0

    for i in range(n_contigs):
        c_start = i * chunk_size
        c_end = total if i == n_contigs - 1 else (i + 1) * chunk_size

        chunk_seq = seq[c_start:c_end]
        chunk_features: list[SeqFeature] = []
        for f in src.features:
            if f.type != "CDS":
                continue
            try:
                f_start = int(f.location.start)
                f_end = int(f.location.end)
            except (TypeError, ValueError):
                continue
            if f_start >= c_start and f_end <= c_end:
                shifted = FeatureLocation(
                    f_start - c_start,
                    f_end - c_start,
                    strand=f.location.strand,
                )
                chunk_features.append(
                    SeqFeature(
                        location=shifted,
                        type=f.type,
                        qualifiers=dict(f.qualifiers),
                    )
                )

        chunk_id = f"chunk{i + 1:04d}"
        chunk_rec = SeqRecord(
            seq=chunk_seq,
            id=chunk_id,
            name=chunk_id,
            description=f"{src.description} | chunk {i + 1}/{n_contigs}",
            annotations={"molecule_type": "DNA"},
            features=chunk_features,
        )
        chunks.append(chunk_rec)
        cds_kept += len(chunk_features)

    cds_dropped = cds_total - cds_kept

    SeqIO.write(chunks, output_path, "genbank")
    print(
        f"chopped {input_path} -> {output_path}: "
        f"{total:,} bp into {n_contigs} chunks of ~{chunk_size:,} bp each. "
        f"CDS kept {cds_kept}/{cds_total} (dropped {cds_dropped} cross-boundary)."
    )


def main() -> int:
    if len(sys.argv) != 4:
        print(__doc__, file=sys.stderr)
        return 2
    chop(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
