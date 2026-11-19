#!/usr/bin/env python3
"""Group AFDB v4 proteome tarballs by NCBI TaxID and preview the .bpdbpack layout."""

import argparse
import os
import re
import sys
from collections import defaultdict

TAR_RE = re.compile(r"proteome-tax_id-(\d+)-(\d+)_v4\.tar$")
MIB = 1024 * 1024
GIB = 1024 * MIB
TIB = 1024 * GIB

# Empirical ratio from v4 E. coli & Human runs: bPDB is ~10.5% of raw GCS tar size
EST_BPDB_RATIO = 0.105


def fmt_bytes(n: int) -> str:
    if n >= TIB:
        return f"{n / TIB:6.2f} TiB"
    if n >= GIB:
        return f"{n / GIB:6.2f} GiB"
    if n >= MIB:
        return f"{n / MIB:6.2f} MiB"
    return f"{n / 1024:6.1f} KiB"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("manifest", help="Path to proteomes_manifest.txt from `gcloud storage ls -l`")
    ap.add_argument("--threshold-mib", type=float, default=100.0,
                    help="Min total tar size (MiB) for a standalone per-TaxID .bpdbpack (default: 100)")
    ap.add_argument("--small-batch-gib", type=float, default=10.0,
                    help="Target raw tar size (GiB) per clustered small-taxa pack (~1 GiB bPDB; default: 10)")
    ap.add_argument("--out-dir", default=None,
                    help="Optional directory to write large_proteomes.tsv and small_batches.tsv")
    args = ap.parse_args()

    # tax_id -> list of (shard_idx, size_bytes, filename)
    taxa: dict[int, list[tuple[int, int, str]]] = defaultdict(list)
    total_tar_bytes = 0
    total_shards = 0

    with open(args.manifest, "r", encoding="utf-8") as f:
        for line in f:
            if not line.rstrip().endswith("_v4.tar"):
                continue
            parts = line.split()
            if len(parts) < 3 or not parts[0].isdigit():
                continue
            size = int(parts[0])
            uri = parts[2]
            fname = uri.rsplit("/", 1)[-1]
            m = TAR_RE.match(fname)
            if not m:
                continue
            tax_id = int(m.group(1))
            shard_idx = int(m.group(2))
            taxa[tax_id].append((shard_idx, size, fname))
            total_tar_bytes += size
            total_shards += 1

    if not taxa:
        print("Error: no *_v4.tar entries found in manifest.", file=sys.stderr)
        return 1

    # Summarize per TaxID
    records = []  # (tax_id, total_bytes, shards_sorted)
    multi_shard_taxa = 0
    shard_gaps = 0
    for tax_id, shards in taxa.items():
        shards.sort(key=lambda x: x[0])
        sz = sum(s[1] for s in shards)
        if len(shards) > 1:
            multi_shard_taxa += 1
        if [s[0] for s in shards] != list(range(len(shards))):
            shard_gaps += 1
        records.append((tax_id, sz, shards))

    print("=" * 72)
    print(f"AFDB v4 Manifest Overview ({args.manifest})")
    print("=" * 72)
    print(f"Unique TaxIDs:          {len(records):>10,}")
    print(f"Total v4 tar shards:    {total_shards:>10,}  (multi-shard TaxIDs: {multi_shard_taxa:,})")
    print(f"Total raw tar volume:   {fmt_bytes(total_tar_bytes):>10}  (est. bPDB: ~{fmt_bytes(int(total_tar_bytes * EST_BPDB_RATIO)).strip()})")
    if shard_gaps:
        print(f"TaxIDs w/ shard gaps:   {shard_gaps:>10,}  (non-contiguous shard indices)")

    # Size distribution by TaxID (not individual shard)
    buckets = [
        ("< 1 MiB",          0,            1 * MIB),
        ("1 - 10 MiB",       1 * MIB,     10 * MIB),
        ("10 - 100 MiB",    10 * MIB,    100 * MIB),
        ("100 - 500 MiB",  100 * MIB,    500 * MIB),
        ("500 MiB - 1 GiB", 500 * MIB,     1 * GIB),
        ("1 - 5 GiB",        1 * GIB,      5 * GIB),
        ("5 - 15 GiB",       5 * GIB,     15 * GIB),
        (">= 15 GiB",       15 * GIB,     1 << 62),
    ]
    print("\nGrouped by TaxID (sum of all shards per organism):")
    print(f"  {'BUCKET':<18} {'TAXIDS':>9} {'SHARDS':>9} {'RAW TAR':>12} {'EST bPDB':>12} {'SHARE':>7}")
    print(f"  {'-'*18} {'-'*9} {'-'*9} {'-'*12} {'-'*12} {'-'*7}")
    for label, lo, hi in buckets:
        sub = [r for r in records if lo <= r[1] < hi]
        n_tax = len(sub)
        n_sh = sum(len(r[2]) for r in sub)
        b_tar = sum(r[1] for r in sub)
        b_bpdb = int(b_tar * EST_BPDB_RATIO)
        pct = 100.0 * b_tar / total_tar_bytes
        print(f"  {label:<18} {n_tax:>9,} {n_sh:>9,} {fmt_bytes(b_tar):>12} {fmt_bytes(b_bpdb):>12} {pct:>6.2f}%")

    # Top 15 largest proteomes (critical for sizing peak RAM per single-file proteome)
    top = sorted(records, key=lambda r: r[1], reverse=True)[:15]
    print("\nTop 15 Largest Proteomes (peak single-pack memory check):")
    print(f"  {'TAXID':>10} {'SHARDS':>7} {'RAW TAR':>12} {'EST bPDB':>12}")
    print(f"  {'-'*10} {'-'*7} {'-'*12} {'-'*12}")
    for tax_id, sz, shards in top:
        print(f"  {tax_id:>10} {len(shards):>7} {fmt_bytes(sz):>12} {fmt_bytes(int(sz * EST_BPDB_RATIO)):>12}")

    # Simulate two-tier packing plan
    thresh = int(args.threshold_mib * MIB)
    batch_target = int(args.small_batch_gib * GIB)

    large = sorted([r for r in records if r[1] >= thresh], key=lambda r: r[0])
    small = sorted([r for r in records if r[1] < thresh], key=lambda r: r[0])

    small_batches: list[list[tuple[int, int, list[tuple[int, int, str]]]]] = []
    cur_batch = []
    cur_bytes = 0
    for r in small:
        if cur_batch and cur_bytes + r[1] > batch_target:
            small_batches.append(cur_batch)
            cur_batch = []
            cur_bytes = 0
        cur_batch.append(r)
        cur_bytes += r[1]
    if cur_batch:
        small_batches.append(cur_batch)

    large_tar = sum(r[1] for r in large)
    small_tar = sum(r[1] for r in small)
    large_shards = sum(len(r[2]) for r in large)
    small_shards = sum(len(r[2]) for r in small)

    print("\n" + "=" * 72)
    print(f"Two-Tier .bpdbpack Plan (threshold = {args.threshold_mib:.0f} MiB, small batch = {args.small_batch_gib:.1f} GiB tar)")
    print("=" * 72)
    print(f"  Large proteomes (1 pack/TaxID): {len(large):>8,} packs  | {large_shards:>7,} tarballs | {fmt_bytes(large_tar)} tar -> ~{fmt_bytes(int(large_tar * EST_BPDB_RATIO)).strip()} bPDB")
    print(f"  Small taxa (clustered packs):   {len(small_batches):>8,} packs  | {small_shards:>7,} tarballs | {fmt_bytes(small_tar)} tar -> ~{fmt_bytes(int(small_tar * EST_BPDB_RATIO)).strip()} bPDB")
    print(f"  TOTAL .bpdbpack files on SSD:   {len(large) + len(small_batches):>8,} files")

    if args.out_dir:
        os.makedirs(args.out_dir, exist_ok=True)
        large_path = os.path.join(args.out_dir, "large_proteomes.tsv")
        small_path = os.path.join(args.out_dir, "small_batches.tsv")

        with open(large_path, "w", encoding="utf-8") as out:
            out.write("# pack_id\ttax_id\ttotal_tar_bytes\tshard_count\ttar_files\n")
            for tax_id, sz, shards in large:
                pack_id = f"proteomes/{tax_id % 100:02d}/tax_{tax_id}"
                files = ",".join(s[2] for s in shards)
                out.write(f"{pack_id}\t{tax_id}\t{sz}\t{len(shards)}\t{files}\n")

        with open(small_path, "w", encoding="utf-8") as out:
            out.write("# pack_id\ttax_count\ttotal_tar_bytes\tshard_count\ttar_files\n")
            for bi, batch in enumerate(small_batches):
                pack_id = f"small/small_{bi:04d}"
                btar = sum(r[1] for r in batch)
                bfiles = [s[2] for r in batch for s in r[2]]
                out.write(f"{pack_id}\t{len(batch)}\t{btar}\t{len(bfiles)}\t{','.join(bfiles)}\n")

        print(f"\nWrote plan files to {args.out_dir}/:")
        print(f"  {large_path}")
        print(f"  {small_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
