#!/usr/bin/env python3
"""Synthetic SNV cohort with planted non-coding drivers, written as a MAF.

Background: random SNVs across chr1-22/X with a fixed substitution spectrum.
Planted: six regulatory elements with known burden, position choice and spectrum.
Everything is reproducible from --seed. The truth file is the answer key.

--catalogue-style regulatory_build (default) is the SYNTH-NC v1 generator, unchanged: elements are
single catalogue rows of the Regulatory Build catalogue. --catalogue-style pcawg_style (SYNTH-NC v2)
reads the PCAWG-style catalogue: an element is all rows sharing an id; E1 is the TERT promoter core
(promCore__TERT__*); E2-E6 are enhancers on chr1-22 of 600-1600 bp in total that overlap no other
catalogue element; planted positions are drawn inside the element's fragments, and the background is
kept 50 kb away from each element's span.
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
    ap.add_argument("--catalogue-style", choices=["regulatory_build", "pcawg_style"], default="regulatory_build")
    a = ap.parse_args()
    if a.catalogue_style == "pcawg_style":
        return main_pcawg_style(a)

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


def load_grouped(path):
    """PCAWG-style catalogue: id -> sorted fragments [(chrom, start0, end)]; also per-chromosome sorted rows."""
    groups, rows = {}, {}
    with open(path) as fh:
        for line in fh:
            c, s, e, name = line.rstrip("\n").split("\t")[:4]
            s, e = int(s), int(e)
            groups.setdefault(name, []).append((c, s, e))
            rows.setdefault(c, []).append((s, e, name))
    for v in groups.values():
        v.sort()
    for v in rows.values():
        v.sort()
    return groups, rows


def overlaps_other(rows, frags, name):
    """True if any fragment of `name` overlaps a row of another element (rows sorted by start)."""
    for c, s, e in frags:
        lst = rows.get(c, [])
        i = bisect.bisect_left(lst, (s - 1_000_000, 0, ""))
        while i < len(lst) and lst[i][0] < e:
            rs, re_, rn = lst[i]
            if rn != name and rs < e and s < re_:
                return True
            i += 1
    return False


def main_pcawg_style(a):
    rng = random.Random(a.seed)
    fa = pysam.FastaFile(a.fasta)
    tbx = pysam.TabixFile(a.cadd)
    base = lambda c, p: fa.fetch(c, p - 1, p).upper()
    donors = [f"SYN{i:03d}" for i in range(1, a.donors + 1)]
    muts = []
    used = set()

    def add(d, c, p, ref, alt, src):
        if (d, c, p) in used:
            return False
        used.add((d, c, p))
        muts.append((d, c, p, ref, alt, src))
        return True

    groups, rows = load_grouped(a.elements)
    tert = sorted(n for n, fr in groups.items() if n.startswith("promCore__TERT__")
                  and any(c == "chr5" and s < 1295113 <= e and s < 1295135 <= e for c, s, e in fr))
    assert tert, "No promCore__TERT element covers both TERT hotspots"
    chosen = {"E1": tert[0]}
    pool = []
    for n in sorted(groups):
        fr = groups[n]
        if not n.startswith("enhancer_") or len({f[0] for f in fr}) != 1 or fr[0][0] not in CHROMS[:22]:
            continue
        if not 600 <= sum(e - s for _, s, e in fr) <= 1600:
            continue
        pool.append(n)
    rng.shuffle(pool)
    span = lambda n: (groups[n][0][0], min(s for _, s, _ in groups[n]), max(e for _, _, e in groups[n]))
    taken = [span(tert[0])]

    def far(n):
        c, s, _ = span(n)
        return all(c != t[0] or abs(s - t[1]) > 5_000_000 for t in taken)

    def seq(n):
        return "".join(fa.fetch(c, s, e).upper() for c, s, e in groups[n])

    def at_frac(n):
        s = seq(n)
        return (s.count("A") + s.count("T")) / max(1, len(s) - s.count("N")), s.count("N")

    for key in ["E2", "E3", "E4", "E5", "E6"]:
        for n in pool:
            if not far(n) or overlaps_other(rows, groups[n], n):
                continue
            atf, nn = at_frac(n)
            if nn:
                continue
            if key == "E6" and not atf <= 0.42:
                continue
            chosen[key] = n
            taken.append(span(n))
            break

    def cadd_all(n):
        out = []
        for c, s, e in groups[n]:
            out += cadd_best(tbx, c, s, e, lambda p, c=c: base(c, p))
        return out

    truth = []
    order = donors[:]
    rng.shuffle(order)
    for i, d in enumerate(order[:18]):
        p, ref, alt = TERT_HOTSPOTS[0 if i < 12 else 1]
        assert base("chr5", p) == ref, "TERT reference base mismatch"
        add(d, "chr5", p, ref, alt, "E1")
    plan = {"E2": (14, "high_cadd"), "E3": (10, "random"), "E4": (4, "high_cadd"),
            "E5": (3, "random"), "E6": (12, "at_only")}
    for key, (k, mode) in plan.items():
        n = chosen[key]
        c = groups[n][0][0]
        ds = rng.sample(donors, k)
        if mode == "high_cadd":
            cands = sorted(cadd_all(n), key=lambda x: -x[2])
            top = cands[: max(k * 3, 30)]
            picks = rng.sample(top, k)
        else:
            positions = [p for _, s, e in groups[n] for p in range(s + 1, e + 1)
                         if base(c, p) in ("AT" if mode == "at_only" else "ACGT")]
            picks = [(p, None, None) for p in rng.sample(positions, k)]
        for d, (p, alt, _) in zip(ds, picks):
            ref = base(c, p)
            add(d, c, p, ref, alt or mutate(rng, ref), key)
    designs = {"E1": "TERT promoter hotspots C228T (12 donors) + C250T (6 donors)",
               "E2": "14 donors, positions drawn from the top-CADD SNVs in the element",
               "E3": "10 donors, random positions and spectrum-drawn alleles",
               "E4": "4 donors, top-CADD SNVs (impact without much burden)",
               "E5": "3 donors, random positions (at or below detection limit)",
               "E6": "12 donors, A/T reference bases only in a GC-rich element (spectrum decoy)"}
    for key, n in chosen.items():
        rows_k = [m for m in muts if m[5] == key]
        c, s0, e0 = span(n)
        atf, _ = at_frac(n)
        phred = {(p, alt): ph for p, alt, ph in cadd_all(n)}
        truth.append({
            "planted_id": key, "element_id": n, "chrom": c, "start0": s0, "end": e0,
            "length_bp": sum(e - s for _, s, e in groups[n]), "fragments": [[s, e] for _, s, e in groups[n]],
            "donors_mutated": len({m[0] for m in rows_k}), "mutations": len(rows_k),
            "element_at_fraction": round(atf, 3),
            "mutation_at_fraction": round(sum(m[3] in "AT" for m in rows_k) / len(rows_k), 3),
            "median_cadd_phred": sorted(phred.get((m[2], m[4]), float("nan")) for m in rows_k)[len(rows_k) // 2],
            "design": designs[key],
        })

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
                continue
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
    json.dump({"version": 2, "catalogue_style": "pcawg_style", "seed": a.seed, "donors": a.donors,
               "total_mutations": len(muts),
               "background_mutations": sum(m[5] == "background" for m in muts),
               "background_per_donor_range": [min(per_donor.values()), max(per_donor.values())],
               "planted": truth}, open(a.out_truth, "w"), indent=2)
    print(json.dumps({"total": len(muts), "planted": {t["planted_id"]: (t["element_id"], t["mutations"]) for t in truth}}, indent=1))


if __name__ == "__main__":
    main()
