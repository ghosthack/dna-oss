"""APOE (ε2/ε3/ε4) genotypes for every family member.

Usage: python3 -I scripts/apoe.py
Writes results/apoe.tsv. Uses measured chip calls when present, otherwise the 23andMe
imputed haplotype dosages (GRCh38), which also phase double heterozygotes.
"""
import glob
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dnalib import BCFTOOLS, RESULTS, WORK, has_genotypes, load_person, people

# (rs429358 allele, rs7412 allele) on one chromosome -> isoform
ISOFORM = {("T", "T"): "ε2", ("T", "C"): "ε3", ("C", "C"): "ε4", ("C", "T"): "ε1"}
GRCH38 = {"rs429358": ("chr19", 44908684, "T", "C"), "rs7412": ("chr19", 44908822, "C", "T")}
# Alzheimer odds vs. the common ε3/ε3 type (Farrer 1997 meta-analysis, white clinic/autopsy samples: ε2/ε2 and
# ε2/ε3 ~0.6, ε2/ε4 2.6, ε3/ε4 3.2, ε4/ε4 14.9; population-based studies give somewhat lower odds).
RISK = {
    "ε2/ε2": "Lower than ε3/ε3 (odds ~0.6×; ε2 is protective)",
    "ε2/ε3": "Lower than ε3/ε3 (odds ~0.6×; ε2 is protective)",
    "ε3/ε3": "Average (most common type; the comparison group)",
    "ε2/ε4": "Higher: odds ~2.5× vs ε3/ε3 (close to ε3/ε4)",
    "ε3/ε4": "Higher: odds ~3× vs ε3/ε3", "ε4/ε4": "Much higher: odds ~12–15× vs ε3/ε3",
}


def imputed_haplotypes(person):
    """[(rs429358, rs7412) per haplotype] from imputed HDS, or None."""
    bcf = glob.glob(os.path.join(WORK, "imputed", person, "*", "chr19.bcf"))
    if not bcf:
        return None
    alleles = {}
    for rs, (chrom, pos, ref, alt) in GRCH38.items():
        out = subprocess.run([BCFTOOLS, "query", "-r", f"{chrom}:{pos}", "-f", "%POS\t%REF\t%ALT\t[%HDS]\n", bcf[0]],
                             capture_output=True, text=True).stdout
        rec = [l.split("\t") for l in out.splitlines() if l.split("\t")[1:3] == [ref, alt]]
        if not rec:
            return None
        try:
            hds = [float(x) for x in rec[0][3].strip().split(",")]
        except ValueError:  # "." = no dosage at this site
            return None
        if any(0.1 < h < 0.9 for h in hds):  # uncertain imputation
            return None
        alleles[rs] = [alt if h >= 0.5 else ref for h in hds]
    return list(zip(alleles["rs429358"], alleles["rs7412"]))


def call(person):
    data, _ = load_person(person)
    g1 = data.get("rs429358", (0, 0, None))[2]
    g2 = data.get("rs7412", (0, 0, None))[2]
    if g1 and g2:
        het1, het2 = g1[0] != g1[1], g2[0] != g2[1]
        if not (het1 and het2):
            a = [g1[0], g1[1]] if het1 else [g1[0]] * 2
            b = [g2[0], g2[1]] if het2 else [g2[0]] * 2
            return "/".join(sorted(ISOFORM[x] for x in zip(a, b))), "measured"
        haps = imputed_haplotypes(person)  # double heterozygote: phase with imputed haplotypes
        if haps:
            return "/".join(sorted(ISOFORM[h] for h in haps)), "measured, phased by imputation"
        return "ε2/ε4", "measured (assumed ε2/ε4; rare ε1/ε3 not excluded)"
    haps = imputed_haplotypes(person)
    if haps:
        return "/".join(sorted(ISOFORM[h] for h in haps)), "imputed (chip no-call)"
    return None, "not determinable"


def main():
    with open(os.path.join(RESULTS, "apoe.tsv"), "w") as fh:
        fh.write("person\tapoe\tsource\trisk\n")
        for p in people():
            if not has_genotypes(p):
                continue
            geno, source = call(p)
            risk = RISK.get(geno, "")
            fh.write(f"{p}\t{geno or ''}\t{source}\t{risk}\n")
            print(f"{p:17s} {geno or '?':7s} {source:40s} {risk}")


if __name__ == "__main__":
    main()
