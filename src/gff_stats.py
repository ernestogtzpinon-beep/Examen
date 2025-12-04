#!/usr/bin/env python3
"""Parse a GFF file and produce JSON stats.

Features:
- Uses argparse for CLI
- Uses dictionaries and comprehensions
- Basic file handling
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from typing import Dict, Tuple


def compute_stats_from_gff(path: str, filter_type: str | None = None) -> Dict:
    """Compute statistics from a GFF file.

    Args:
        path: Path to GFF file
        filter_type: If provided, only compute stats for this feature type

    Returns a dict with keys: total_features, by_type, avg_length, strand_distribution
    """
    total = 0
    counts = defaultdict(int)
    length_sums = defaultdict(int)
    strand_counts = defaultdict(int)

    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                # skip malformed lines
                continue
            feature_type = cols[2]
            
            # Apply filter if specified
            if filter_type and feature_type != filter_type:
                continue
            
            try:
                start = int(cols[3])
                end = int(cols[4])
            except ValueError:
                continue
            length = end - start + 1
            strand = cols[6]

            total += 1
            counts[feature_type] += 1
            length_sums[feature_type] += length
            if strand in ("+", "-"):
                strand_counts[strand] += 1
            else:
                strand_counts.setdefault(".", 0)

    # average lengths as floats with one decimal
    avg_length = {
        k: round(length_sums[k] / counts[k], 1) for k in counts
    } if counts else {}

    # strand distribution as percentages (integers summing to ~100)
    total_strands = sum(v for k, v in strand_counts.items() if k in ("+", "-"))
    if total_strands:
        strand_distribution = {
            "+": int(round(100 * strand_counts.get("+", 0) / total_strands)),
            "-": int(round(100 * strand_counts.get("-", 0) / total_strands)),
        }
    else:
        strand_distribution = {"+": 0, "-": 0}

    result = {
        "total_features": total,
        "by_type": dict(counts),
        "avg_length": avg_length,
        "strand_distribution": strand_distribution,
    }
    
    # Add filter info if applied
    if filter_type:
        result["filter_type"] = filter_type

    return result


def cli(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Compute stats from GFF file and write JSON output"
    )
    parser.add_argument(
        "--gff",
        default="input.gff",
        help="Input GFF file (default: input.gff)"
    )
    parser.add_argument(
        "--out",
        default="output.json",
        help="Output JSON file (default: output.json)"
    )
    parser.add_argument(
        "--filter-type",
        default=None,
        help="Filter statistics by feature type (e.g., gene, CDS, mRNA)"
    )
    args = parser.parse_args(argv)

    stats = compute_stats_from_gff(args.gff, filter_type=args.filter_type)
    with open(args.out, "w", encoding="utf-8") as outfh:
        json.dump(stats, outfh, indent=2, ensure_ascii=False)

    return 0


if __name__ == "__main__":
    raise SystemExit(cli())
