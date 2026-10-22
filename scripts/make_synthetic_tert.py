#!/usr/bin/env python3
"""Deterministic synthetic ActiveDriverWGS fixture on real GRCh38 bases (no patient data).
Usage: make_synthetic_tert.py <hg38.fa> <hg38.fa.fai> <outdir>"""
import random, sys, pathlib
fa, fai, out = sys.argv[1], sys.argv[2], pathlib.Path(sys.argv[3]); out.mkdir(parents=True, exist_ok=True)
idx = {}
for l in open(fai):
    f = l.rstrip("\n").split("\t"); idx[f[0]] = tuple(map(int, f[1:5]))
def fetch(chrom, start1, end1):
    L, off, lb, lB = idx[chrom]; s0 = start1 - 1; n = end1 - start1 + 1
    with open(fa, "rb") as fh:
        fh.seek(off + (s0 // lb) * lB + s0 % lb)
        raw = fh.read((n // lb + 2) * lB).replace(b"\n", b"").replace(b"\r", b"")
    return raw[:n].decode().upper()
W0, W1 = 1_200_001, 1_400_000
win = fetch("chr5", W0, W1)
base = lambda p: win[p - W0]
assert base(1295113) == "G", f"chr5:1295113 is {base(1295113)}, expected G (TERT C228T)"
rng = random.Random(20261020)
rows = []
donors = [f"D{i:02d}" for i in range(1, 31)]
for d in donors:
    k = 0
    while k < 15:
        p = rng.randint(W0 + 1000, W1 - 1000); r = base(p)
        if r not in "ACGT": continue
        rows.append((d, "chr5", p, r, rng.choice([b for b in "ACGT" if b != r]), "SNP")); k += 1
for d in donors[:10]:
    rows.append((d, "chr5", 1295113, "G", "A", "SNP"))
# rows the module must drop before ActiveDriverWGS: an indel, an MNV and a chrM SNV
p = 1_310_000; rows.append(("D11", "chr5", p, base(p) + base(p + 1), base(p), "DEL"))
p = 1_320_000; rows.append(("D12", "chr5", p, base(p) + base(p + 1), "".join("A" if b != "A" else "C" for b in base(p) + base(p + 1)), "DNP"))
if "chrM" in idx:
    r = fetch("chrM", 100, 100); rows.append(("D13", "chrM", 100, r, "A" if r != "A" else "G", "SNP"))
rows = sorted(set(rows), key=lambda x: (x[1], x[2], x[0]))
with open(out / "mutations.tsv", "w") as fh:
    fh.write("sample\tchrom\tpos\tref\talt\tvtype\n")
    for r in rows: fh.write("\t".join(map(str, r)) + "\n")
els = [("chr5", 1220000, 1220600, "elem_A"), ("chr5", 1250000, 1250800, "elem_B"),
       ("chr5", 1294989, 1295601, "TERT_promoter_fixture"), ("chr5", 1330000, 1330500, "elem_C"),
       ("chr5", 1360000, 1360700, "elem_D"), ("chr5", 1380000, 1380400, "elem_E"),
       ("chrUn_KI270742v1", 100, 600, "elem_unplaced")]       # dropped by the BED4 adapter
with open(out / "elements.bed", "w") as fh:
    for e in els: fh.write("\t".join(map(str, e)) + "\n")
print(f"fixture: {len(rows)} mutation rows, {len(donors)} donors, {len(els)} elements -> {out}")
