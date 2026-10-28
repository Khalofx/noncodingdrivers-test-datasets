#!/usr/bin/env python3
"""Synthetic SNV cohort with planted non-coding drivers, written as a MAF.

Background: random SNVs across chr1-22/X with a fixed substitution spectrum.
Planted: six regulatory elements with known burden, position choice and spectrum.
Everything is reproducible from --seed. The truth file is the answer key.
"""
import argparse
import bisect
import csv
import gzip
import json
import random

import pysam

COMP = str.maketrans("ACGT", "TGCA")
# Pyrimidine-context spectrum: C>T dominant, T>C next, as in many cancer genomes.
C_ALTS = [("A", 0.25), ("G", 0.15), ("T", 0.60)]
T_ALTS = [("A", 0.30), ("C", 0.45), ("G", 0.25)]
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
TERT_HOTSPOTS = [(1295113, "G", "A"), (1295135, "G", "A")]  # C228T, C250T on the minus strand


def pick(rng, weighted):
    r, acc = rng.random(), 0.0
    for value, w in weighted:
        acc += w
        if r < acc:
            return value
    return weighted[-1][0]


def mutate(rng, ref):
    ref = ref.upper()
    if ref in "CT":
        return pick(rng, C_ALTS if ref == "C" else T_ALTS)
    # purine: draw in pyrimidine context, then complement
    alt = pick(rng, C_ALTS if ref == "G" else T_ALTS)
    return alt.translate(COMP)


def load_elements(path):
    out = []
    with open(path) as fh:
        for line in fh:
            c, s, e, name = line.rstrip("\n").split("\t")[:4]
            out.append((c, int(s), int(e), name))
    return out


def cadd_best(tbx, chrom, start, end, ref_of):
    """All (pos, alt, phred) within [start, end), 1-based pos."""
    rows = []
    for r in tbx.fetch(chrom.removeprefix("chr"), start, end):
        f = r.split("\t")
        pos, ref, alt, phred = int(f[1]), f[2], f[3], float(f[5])
        if ref_of(pos) == ref:
            rows.append((pos, alt, phred))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--elements", required=True)
    ap.add_argument("--cadd", required=True)
    ap.add_argument("--out-maf", required=True)
    ap.add_argument("--out-truth", required=True)
    ap.add_argument("--donors", type=int, default=60)
    ap.add_argument("--seed", type=int, default=20260923)
    a = ap.parse_args()

    rng = random.Random(a.seed)
    fa = pysam.FastaFile(a.fasta)
    tbx = pysam.TabixFile(a.cadd)
    base = lambda c, p: fa.fetch(c, p - 1, p).upper()
    donors = [f"SYN{i:03d}" for i in range(1, a.donors + 1)]
    muts = []  # (donor, chrom, pos, ref, alt, source)
    used = set()

    def add(d, c, p, ref, alt, src):
        if (d, c, p) in used:
            return False
        used.add((d, c, p))
        muts.append((d, c, p, ref, alt, src))
        return True

    # ---- planted elements -------------------------------------------------
    elements = load_elements(a.elements)
    tert = [e for e in elements if e[0] == "chr5" and e[1] < 1295113 <= e[2] and e[3].startswith("promoter")]
    assert tert, "No promoter element covers the TERT hotspot"
    chosen = {"E1": tert[0]}
    pool = [e for e in elements if e[3].startswith("enhancer_") and e[0] in CHROMS[:22] and 600 <= e[2] - e[1] <= 1600]
    rng.shuffle(pool)
    taken = [tert[0]]

    def far(e):
        return all(e[0] != t[0] or abs(e[1] - t[1]) > 5_000_000 for t in taken)

    def at_frac(e):
        s = fa.fetch(e[0], e[1], e[2]).upper()
        return (s.count("A") + s.count("T")) / max(1, len(s) - s.count("N")), s.count("N")

    for key in ["E2", "E3", "E4", "E5", "E6"]:
        for e in pool:
            if not far(e):
                continue
            atf, n = at_frac(e)
            if n:
                continue
            if key == "E6" and not atf <= 0.42:
                continue
            chosen[key] = e
            taken.append(e)
            break

    truth = []
    order = donors[:]
    rng.shuffle(order)
    # E1: TERT hotspots in 18 donors (12 at C228T, 6 at C250T)
    for i, d in enumerate(order[:18]):
        p, ref, alt = TERT_HOTSPOTS[0 if i < 12 else 1]
        assert base("chr5", p) == ref, "TERT reference base mismatch"
        add(d, "chr5", p, ref, alt, "E1")
    plan = {"E2": (14, "high_cadd"), "E3": (10, "random"), "E4": (4, "high_cadd"),
            "E5": (3, "random"), "E6": (12, "at_only")}
    for key, (n, mode) in plan.items():
        c, s, e, name = chosen[key]
        ds = rng.sample(donors, n)
        if mode == "high_cadd":
            cands = sorted(cadd_best(tbx, c, s, e, lambda p: base(c, p)), key=lambda x: -x[2])
            top = cands[: max(n * 3, 30)]
            picks = rng.sample(top, n)
        else:
            positions = [p for p in range(s + 1, e + 1) if base(c, p) in ("AT" if mode == "at_only" else "ACGT")]
            picks = [(p, None, None) for p in rng.sample(positions, n)]
        for d, (p, alt, _) in zip(ds, picks):
            ref = base(c, p)
            add(d, c, p, ref, alt or mutate(rng, ref), key)
    for key, (c, s, e, name) in chosen.items():
        rows = [m for m in muts if m[5] == key]
        atf, _ = at_frac((c, s, e, name))
        phred = {(p, alt): ph for p, alt, ph in cadd_best(tbx, c, s, e, lambda p: base(c, p))}
        truth.append({
            "planted_id": key, "element_id": name, "chrom": c, "start0": s, "end": e,
            "length_bp": e - s, "donors_mutated": len({m[0] for m in rows}), "mutations": len(rows),
            "element_at_fraction": round(atf, 3),
            "mutation_at_fraction": round(sum(m[3] in "AT" for m in rows) / len(rows), 3),
            "median_cadd_phred": sorted(phred.get((m[2], m[4]), float("nan")) for m in rows)[len(rows) // 2],
            "design": {"E1": "TERT promoter hotspots C228T (12 donors) + C250T (6 donors)",
                       "E2": "14 donors, positions drawn from the top-CADD SNVs in the element",
                       "E3": "10 donors, random positions and spectrum-drawn alleles",
                       "E4": "4 donors, top-CADD SNVs (impact without much burden)",
                       "E5": "3 donors, random positions (at or below detection limit)",
                       "E6": "12 donors, A/T reference bases only in a GC-rich element (spectrum decoy)"}[key],
        })

    # ---- background --------------------------------------------------------
    lens = [fa.get_reference_length(c) for c in CHROMS]
    cum = [0]
    for L in lens:
        cum.append(cum[-1] + L)
    planted_spans = [(t["chrom"], t["start0"] - 50_000, t["end"] + 50_000) for t in truth]
    per_donor = {d: max(1500, int(rng.lognormvariate(8.25, 0.35))) for d in donors}
    for d in donors:
        k = 0
        while k < per_donor[d]:
            g = rng.randrange(cum[-1])
            i = bisect.bisect_right(cum, g) - 1
            c, p = CHROMS[i], g - cum[i] + 1
            if any(c == pc and ps <= p <= pe for pc, ps, pe in planted_spans):
                continue  # keep planted windows clean so the truth is exact
            ref = base(c, p)
            if ref not in "ACGT":
                continue
            if add(d, c, p, ref, mutate(rng, ref), "background"):
                k += 1

    muts.sort(key=lambda m: (CHROMS.index(m[1]), m[2], m[0]))
    with gzip.open(a.out_maf, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Variant_Classification",
                    "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
                    "Tumor_Sample_Barcode", "Project_Code"])
        for d, c, p, ref, alt, _ in muts:
            w.writerow(["Unknown", c, p, p, "IGR", "SNP", ref, ref, alt, d, "SYNTH-NC"])
    json.dump({"seed": a.seed, "donors": a.donors, "total_mutations": len(muts),
               "background_mutations": sum(m[5] == "background" for m in muts),
               "background_per_donor_range": [min(per_donor.values()), max(per_donor.values())],
               "planted": truth}, open(a.out_truth, "w"), indent=2)
    print(json.dumps({"total": len(muts), "planted": {t["planted_id"]: (t["element_id"], t["mutations"]) for t in truth}}, indent=1))


if __name__ == "__main__":
    main()
