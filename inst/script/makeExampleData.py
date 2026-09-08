#!/usr/bin/env python3
"""Generate the example peak files shipped in inst/extdata.

Three replicates over a stretch of chr1 and chr2. Most positions are
shared, a portion is deliberately weak in one replicate so the rescue
step has something to do, and the third replicate carries extra noise
peaks so the discard step does too. A BED3 copy of the first replicate is
written as well, to exercise the presence-only path.

Run with: python3 makeExampleData.py ../extdata
"""

import random
import sys
from pathlib import Path

random.seed(20240517)

N_SHARED = 220
N_WEAK = 60
N_NOISE = 90
CHROMS = {"chr1": 12_000_000, "chr2": 9_000_000}


def make_positions(n, used, chrom_lengths):
    """Draw non-overlapping peak starts, keeping a gap between them."""
    out = []
    while len(out) < n:
        chrom = random.choice(list(chrom_lengths))
        start = random.randint(1, chrom_lengths[chrom] - 2000)
        block = start // 5000
        if (chrom, block) in used:
            continue
        used.add((chrom, block))
        out.append((chrom, start))
    return out


def jitter(start, width, amount=90):
    """Shift a peak slightly, the way independent calls never agree exactly."""
    shift = random.randint(-amount, amount)
    return max(1, start + shift), width + random.randint(-40, 40)


def write_narrowpeak(path, records):
    with open(path, "w") as handle:
        for i, (chrom, start, width, neglog10p) in enumerate(records):
            end = start + width
            summit = width // 2 + random.randint(-30, 30)
            score = min(1000, int(neglog10p * 12))
            handle.write(
                f"{chrom}\t{start}\t{end}\tpeak_{i + 1}\t{score}\t.\t"
                f"{neglog10p / 3:.5f}\t{neglog10p:.5f}\t"
                f"{max(0.0, neglog10p - 1):.5f}\t{summit}\n"
            )


def write_bed3(path, records):
    with open(path, "w") as handle:
        for chrom, start, width, _ in records:
            handle.write(f"{chrom}\t{start}\t{start + width}\n")


def build():
    used = set()
    shared = make_positions(N_SHARED, used, CHROMS)
    weak = make_positions(N_WEAK, used, CHROMS)

    replicates = []
    # the third replicate is the shallow one: weaker everywhere and noisier
    for depth in (1.0, 0.95, 0.55):
        records = []

        for chrom, start, in shared:
            width = random.randint(280, 900)
            s, w = jitter(start, width)
            strength = random.gauss(11, 3) * depth
            records.append((chrom, s, max(w, 150), max(strength, 4.2)))

        # positions that only clear the weak threshold anywhere
        for chrom, start in weak:
            width = random.randint(250, 600)
            s, w = jitter(start, width)
            strength = random.gauss(5.2, 0.9) * depth
            records.append((chrom, s, max(w, 150), max(strength, 4.05)))

        # replicate-specific noise, no support expected
        noise_used = set()
        for chrom, start in make_positions(N_NOISE, noise_used, CHROMS):
            width = random.randint(200, 450)
            strength = max(random.gauss(4.9, 0.6), 4.02)
            records.append((chrom, start, width, strength))

        records.sort(key=lambda r: (r[0], r[1]))
        replicates.append(records)

    return replicates


def main():
    target = Path(sys.argv[1] if len(sys.argv) > 1 else "../extdata")
    target.mkdir(parents=True, exist_ok=True)

    replicates = build()
    for i, records in enumerate(replicates, start=1):
        write_narrowpeak(target / f"rep{i}.narrowPeak", records)

    write_bed3(target / "rep1_minimal.bed", replicates[0])
    print(f"wrote {len(replicates)} narrowPeak files to {target}")


if __name__ == "__main__":
    main()
