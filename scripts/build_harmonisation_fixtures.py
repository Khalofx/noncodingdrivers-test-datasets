#!/usr/bin/env python3
"""Build the input-harmonisation fixtures from the SYNTH-TERT canonical table (stdlib only).

Writes into <outdir> (data/synthetic_tert/harmonisation of noncodingdrivers-test-datasets):
  synth_tert.hg19.maf            all 30 donors, hg19 = hg38 + 115 (one ungapped chain block),
                                 PCAWG style: Ensembl contig names, Start_position spelling, CRLF;
                                 plus the chrM row as 'MT' (dropped by maf_split.awk) and one
                                 row outside the fixture chain (rejected by Picard as NoTarget)
  synth_tert_tiers.hg19.maf      the MAF donors of the mixed-tier samplesheet
  mutations_tiers.tsv            canonical table of the mutations donors of the mixed-tier sheet
  vcf/*.vcf                      per-sample VCFs (bgzipped and indexed by the caller script)
  expected.snvs.tsv              expected SNV rows of both MAF and mixed-tier runs (sample chrom pos ref alt)
"""
import collections, os, subprocess, sys

src, outdir, fa38, fa19, ttv = sys.argv[1:6]
OFFSET = 115
CHAIN_T0 = 1100000          # hg19 start (0-based) of the fixture chain block
os.makedirs(f"{outdir}/vcf", exist_ok=True)


class Fasta:
    def __init__(self, path):
        self.fh = open(path, "rb"); self.idx = {}
        for line in open(path + ".fai"):
            n, l, o, lb, lw = line.split("\t")[:5]; self.idx[n] = (int(l), int(o), int(lb), int(lw))
    def get(self, c, s, e):
        l, o, lb, lw = self.idx[c]
        a = o + (s - 1) // lb * lw + (s - 1) % lb; b = o + (e - 1) // lb * lw + (e - 1) % lb
        self.fh.seek(a); return self.fh.read(b - a + 1).replace(b"\n", b"").decode().upper()


g38, g19 = Fasta(fa38), Fasta(fa19)
rows = [l.rstrip("\n").split("\t") for l in open(src)]
hdr, rows = rows[0], rows[1:]
assert hdr == ["sample", "chrom", "pos", "ref", "alt", "vtype"], hdr
rows = [dict(zip(hdr, r)) for r in rows]
donors = sorted({r["sample"] for r in rows})
print(f"{len(rows)} rows, {len(donors)} donors, vtypes {dict(collections.Counter(r['vtype'] for r in rows))}")
print("non-chr5 rows:", [(r["sample"], r["chrom"], r["pos"], r["ref"], r["alt"], r["vtype"]) for r in rows if r["chrom"] != "chr5"])

# hg38 REF check and hg19 REF identity at the lifted position
bad38 = bad19 = 0
for r in rows:
    if r["chrom"] != "chr5" or r["ref"] == "-":
        continue
    p, ref = int(r["pos"]), r["ref"]
    if g38.get("chr5", p, p + len(ref) - 1) != ref: bad38 += 1
    if g19.get("chr5", p + OFFSET, p + OFFSET + len(ref) - 1) != ref: bad19 += 1
print(f"REF vs hg38 fixture mismatches: {bad38}; REF vs hg19 fixture at pos+115 mismatches: {bad19}")


def maf_line(r, chrom, start):
    ref, alt, vt = r["ref"], r["alt"], r["vtype"]
    end = start + 1 if vt == "INS" else start + len(ref) - 1
    return [r["sample"], chrom, str(start), str(end), "+", vt, ref, ref, alt, r["sample"], r["sample"] + "_N", "SYNTH-TERT", "30", "10"]


MAF_HDR = ["Hugo_Symbol", "Chromosome", "Start_position", "End_position", "Strand", "Variant_Type", "Reference_Allele",
           "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode",
           "Project_Code", "t_ref_count", "t_alt_count"]


def write_maf(path, donor_set, extra_outside=True):
    lines = []
    for r in rows:
        if r["sample"] not in donor_set:
            continue
        if r["chrom"] == "chr5":
            lines.append(maf_line(r, "5", int(r["pos"]) + OFFSET))
        elif r["chrom"] == "chrM":
            lines.append(maf_line(r, "MT", int(r["pos"])))
        else:
            sys.exit(f"unexpected contig {r['chrom']}")
    n_outside = 0
    if extra_outside:
        # one SNV 50 kb before the chain block: lifts nowhere (Picard NoTarget)
        d = sorted(donor_set)[0]; p = CHAIN_T0 - 50000; ref = g19.get("chr5", p, p)
        alt = {"A": "C", "C": "T", "G": "A", "T": "G"}[ref]
        lines.append(["OUTSIDE", "5", str(p), str(p), "+", "SNP", ref, ref, alt, d, d + "_N", "SYNTH-TERT", "30", "10"])
        n_outside = 1
    lines[0:0] = [MAF_HDR]
    with open(path, "w", newline="") as fh:
        fh.write("#version 2.4\r\n")
        for l in lines: fh.write("\t".join(l) + "\r\n")
    return len(lines) - 1, n_outside


def snv_rows(donor_set, chroms):
    out = sorted((r["sample"], r["chrom"], int(r["pos"]), r["ref"], r["alt"]) for r in rows
                 if r["sample"] in donor_set and r["vtype"] == "SNP" and r["chrom"] in chroms)
    return out


def write_expected(path, recs):
    with open(path, "w") as fh:
        fh.write("sample\tchrom\tpos\tref\talt\n")
        for s, c, p, ref, alt in sorted(recs, key=lambda x: (x[0], x[1], x[2], x[3], x[4])):
            fh.write(f"{s}\t{c}\t{p}\t{ref}\t{alt}\n")


# 1. whole-cohort hg19 MAF
n, n_out = write_maf(f"{outdir}/synth_tert.hg19.maf", set(donors))
all_snvs = snv_rows(set(donors), {"chr5"})
write_expected(f"{outdir}/expected.snvs.tsv", all_snvs)
print(f"synth_tert.hg19.maf: {n} rows ({n_out} outside the chain)")

# 2. mixed tiers: vcf (UCSC names), vcf (Ensembl names), sarek, maf (hg19), mutations
vcf_ucsc, vcf_ens, sarek = donors[0], donors[1], donors[2]
maf_donors = set(donors[3:15]); mut_donors = set(donors[15:])
n2, _ = write_maf(f"{outdir}/synth_tert_tiers.hg19.maf", maf_donors)
with open(f"{outdir}/mutations_tiers.tsv", "w") as fh:
    fh.write("\t".join(hdr) + "\n")
    for r in rows:
        if r["sample"] in mut_donors: fh.write("\t".join(r[k] for k in hdr) + "\n")
exp = snv_rows(mut_donors, {"chr5", "chrM"}) + snv_rows({vcf_ucsc, vcf_ens, sarek} | maf_donors, {"chr5"})
assert sorted(exp) == sorted(all_snvs), "mixed-tier SNVs differ from the whole-cohort set"
print(f"tiers: vcf {vcf_ucsc}, vcf(ensembl) {vcf_ens}, sarek {sarek}, maf {len(maf_donors)} donors ({n2} rows), mutations {len(mut_donors)} donors")

# per-sample GRCh38 VCFs through the gate tool (indels anchored from the hg38 fixture FASTA)
tab = f"{outdir}/vcf/_table.tsv"
with open(tab, "w") as fh:
    fh.write("\t".join(hdr) + "\n")
    for r in rows:
        if r["chrom"] == "chr5" and r["sample"] in {vcf_ucsc, vcf_ens, sarek}:
            fh.write("\t".join(r[k] for k in hdr) + "\n")
    # all-donor SNVs as one sample: the 'GRCh38 calls declared hg19' mismatch fixture
    seen = set()
    for r in rows:
        k = (r["chrom"], r["pos"], r["ref"], r["alt"])
        if r["chrom"] == "chr5" and r["vtype"] == "SNP" and k not in seen:
            seen.add(k); fh.write("\t".join(["synth_tert_all", r["chrom"], r["pos"], r["ref"], r["alt"], "SNP"]) + "\n")
subprocess.run([sys.executable, ttv, tab, f"{outdir}/vcf", "--fasta", fa38], check=True)
os.remove(tab)
# Ensembl contig names for one sample
p = f"{outdir}/vcf/{vcf_ens}.vcf"; txt = open(p).read().replace("##contig=<ID=chr5,", "##contig=<ID=5,")
txt = "\n".join(l[3:] if l.startswith("chr5\t") else l for l in txt.split("\n"))
open(f"{outdir}/vcf/{vcf_ens}.ensembl.vcf", "w").write(txt); os.remove(p)
# sarek-style tumour/normal VCF: tumour column <donor>_T, normal <donor>_N (0/0)
p = f"{outdir}/vcf/{sarek}.vcf"; out = []
for l in open(p).read().rstrip("\n").split("\n"):
    if l.startswith("#CHROM"): out.append("\t".join(l.split("\t")[:9] + [f"{sarek}_N", f"{sarek}_T"]))
    elif l.startswith("#"): out.append(l)
    else: f = l.split("\t"); out.append("\t".join(f[:9] + ["0/0", f[9]]))
open(f"{outdir}/vcf/{sarek}.sarek.vcf", "w").write("\n".join(out) + "\n"); os.remove(p)
print("tumour column of sarek fixture:", f"{sarek}_T")
open(f"{outdir}/_donors.txt", "w").write(f"{vcf_ucsc} {vcf_ens} {sarek}\n")
