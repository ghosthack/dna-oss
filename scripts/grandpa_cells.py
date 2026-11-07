"""Grandpa* column for the landing-page trait tables, from the virtual-grandpa haplotype (scripts/grandpa.py).

We only know the ONE copy grandpa passed to mom, so each cell names that copy and leaves his other copy as "?".
Effect alleles are the plus-strand alleles used by wellness.py (same studies, same directions).

Usage: python3 -I scripts/grandpa_cells.py  → results/grandpa_cells.json
"""
import csv
import gzip
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dnalib import FAMILY, RAW, RESULTS, WORK
from pgx import SNPS as PGX_SNPS, ancestry_json, mostly_european
from pigmentation import MC1R_STRONG, MC1R_WEAK
from wellness import CARRIER_PANELS

HAP = os.path.join(WORK, "grandpa", "grandpa_haplotype.tsv.gz")

# report → (rsids, function(alleles dict) → label of his transmitted copy)
def one(rsid, effect, yes, no):
    return [rsid], lambda a: yes if a[rsid] == effect else no


def caffeine(a):
    n = (a["rs2472297"] == "T") + (a["rs4410790"] == "C")
    return f"{n} of 2"


def amd(a):
    risk = ["CFH risk"] * (a["rs1061170"] == "C") + ["ARMS2 risk"] * (a["rs10490924"] == "T")
    return " + ".join(risk) or "No risk"


def hfe(a):
    return "C282Y" if a["rs1800562"] == "A" else "H63D" if a["rs1799945"] == "G" else "Normal"


def clotting(a):
    if "rs1799963" not in a:  # prothrombin has no regular chip probe on some chips, so it can be missing
        return "Factor V Leiden" if a["rs6025"] == "T" else "No FVL"
    found = ["Factor V Leiden"] * (a["rs6025"] == "T") + ["Prothrombin variant"] * (a["rs1799963"] == "A")
    return " + ".join(found) or "None"


def alpha1(a):
    return "Z" if a["rs28929474"] == "T" else "S" if a["rs17580"] == "A" else "M"


def mthfr(a):
    return "677T" if a["rs1801133"] == "A" else "1298C" if a["rs1801131"] == "G" else "Typical"


def celiac(a):
    # Same tag rules as wellness.celiac(): DQ2.5 = rs2187668 T, DQ8 = rs7454108 C; DQ2.2 possible only with
    # rs2395182 T + rs7775228 C (rs2395182 G or rs7775228 T excludes it); rs4713586 G = maybe DQ4 (caveat only).
    found = ["DQ2.5"] * (a["rs2187668"] == "T") + ["DQ8"] * (a["rs7454108"] == "C")
    if a["rs2395182"] == "T" and a["rs7775228"] == "C":
        dq22 = "possible DQ2.2" + (" (or DQ4)" if a.get("rs4713586") == "G" else "")
        return " + ".join(found + [dq22]) if found else dq22[0].upper() + dq22[1:]
    return " + ".join(found) or "No DQ2/8"


def bitter(a):
    hap = (a["rs713598"], a["rs1726866"], a["rs10246939"])
    return {("G", "G", "C"): "PAV (taster)", ("C", "A", "T"): "AVI"}.get(hap, "Rare")


CELLS = {
    "Alcohol flush": one("rs671", "A", "Flush variant", "Normal"),
    "Caffeine consumption": (["rs2472297", "rs4410790"], caffeine),
    "Deep sleep": one("rs73598374", "T", "Deep-sleep variant", "Typical"),
    "Lactose intolerance": one("rs4988235", "A", "Tolerant-type", "Intolerant-type"),
    "Muscle composition": one("rs1815739", "T", "Endurance-type", "Power-type"),
    "Saturated fat & weight": one("rs5082", "G", "G", "Not sensitive"),
    "Sleep movement": one("rs3923809", "A", "More", "Less"),
    "Macular degeneration": (["rs1061170", "rs10490924"], amd),
    "Hemochromatosis (iron)": (["rs1800562", "rs1799945"], hfe),
    "Inherited clotting": (["rs6025"], clotting),
    "Alpha-1 antitrypsin": (["rs28929474", "rs17580"], alpha1),
    "Folate / homocysteine": (["rs1801133", "rs1801131"], mthfr),
    "Vitamin D": one("rs2282679", "G", "Lower-type", "Typical"),
    "Vitamin B12": one("rs602662", "G", "Lower-type", "Higher-type"),
    "Vitamin A (from β-carotene)": one("rs7501331", "T", "Reduced", "Typical"),
    "Celiac genes (HLA-DQ2/DQ8)": (["rs2187668", "rs7454108", "rs2395182", "rs7775228"], celiac),
    "Bitter taste": (["rs713598", "rs1726866", "rs10246939"], bitter),
    "Cilantro soapiness": one("rs72921001", "A", "Less soapy", "Soapy-prone"),
    "Asparagus smell": one("rs4481887", "A", "Higher-odds copy", "Lower-odds copy"),
    "Omega-3 conversion": one("rs174537", "T", "Slow", "Efficient"),
    "Lipoprotein(a)": (["rs10455872", "rs3798220"],
                       lambda a: "Risk variant" if a["rs10455872"] == "G" or a["rs3798220"] == "C" else "No risk variant"),
    "HLA-B27 (spine arthritis)": one("rs4349859", "A", "B27 likely", "No B27"),
    "Earwax type": one("rs17822931", "T", "Dry-type", "Wet!"),  # wet is dominant: a wet copy means wet earwax
    "Parkinson's genes": (["rs34637584", "i4000415"],
                          lambda a: " + ".join(["LRRK2 G2019S"] * (a["rs34637584"] == "A") + ["GBA N370S"] * (a["i4000415"] == "C"))
                          or "Not detected"),
    "Hereditary amyloidosis": (["rs28933979", "rs76992529", "rs121918070"],
                               lambda a: "Detected" if a["rs28933979"] == "A" or a["rs76992529"] == "A"
                               or a["rs121918070"] == "G" else "Not detected"),
    "Kidney disease (APOL1)": (["rs73885319"], lambda a: "G1" if a["rs73885319"] == "G" else "Not detected (G1 only)"),
    "Hereditary prostate cancer": one("rs138213197", "T", "G84E", "Not detected"),
    "Colon polyp genes (MUTYH)": (["rs34612342", "rs36053993"],
                                  lambda a: "Carrier (Y179C)" if a["rs34612342"] == "C" else
                                  "Carrier (G396D)" if a["rs36053993"] == "T" else "None"),
}
OPTIONAL = {"rs4713586"}  # caveat-only SNP (as in wellness.py): read when reconstructable, never required
CARRIER_REPORTS = {"HBB": "Beta-thalassemia & sickle cell", "CFTR": "Cystic fibrosis", "GJB2": "Inherited hearing loss",
                   "MEFV": "Familial Mediterranean Fever", "PAH": "Phenylketonuria (PKU)", "GBA": "Gaucher disease",
                   "ALDOB": "Hereditary fructose intolerance", "ACADM": "MCAD deficiency (fat metabolism)", "GAA": "Pompe disease",
                   "PMM2": "PMM2-CDG (metabolic disorder)", "SLC26A4": "Pendred syndrome (hearing)", "PKLR": "Pyruvate kinase deficiency"}


# Genes whose copy is described by its variant rather than a star allele.
COPY_WORDS = {"VKORC1": ("Sensitive", "Normal"), "SLCO1B1": ("Reduced", "Normal"), "ABCG2": ("Reduced", "Normal"),
              "CYP2B6": ("Reduced", "Normal")}


def med_cell(gene, alleles):
    snps = PGX_SNPS[gene]
    have = {rs: alleles[rs] for rs in snps if rs in alleles}
    if not have:
        return None
    hits = [label for rs, (alt, label) in snps.items() if have.get(rs) == alt]
    partial = len(have) < len(snps)
    if gene.startswith("HLA"):
        return "Possible" if hits else "Unlikely"
    if gene in COPY_WORDS:
        return COPY_WORDS[gene][0] if hits else COPY_WORDS[gene][1]
    if gene == "TPMT" and {"*3B", "*3C"} <= set(hits):
        hits = ["*3A"]
    if gene == "CYP2D6" and "*4" in hits and "*10" in hits:
        hits.remove("*10")  # 100C>T (rs1065852) is part of the *4 haplotype, as in pgx.py
    if gene == "G6PD":  # X-linked: grandpa's only X went to mom, so this is his full genotype
        # A- = 202A on 376G; 202A alone counts when 376 wasn't resolved (same rule as pgx.py)
        deficient = ((have.get("rs1050828") == "T" and have.get("rs1050829", "C") == "C")
                     or have.get("rs5030868") == "A")
        return "Deficient" if deficient else "Normal"
    label = "+".join(hits) if hits else "*1"
    return label + ("?" if partial and not hits else "")


def why_missing(rsids, status):
    """Why these rsIDs have no grandpa allele, from the haplotype file's method column ({rsid: (method, grandma)})."""
    by = {}
    for rs in rsids:
        method, gma = status.get(rs, (None, None))
        why = ("not on mom's chip" if method is None else
               "mom/grandma genotype conflict" if method == "mendel_error" else
               "not a single-base call" if method not in ("unresolved", "kid_conflict") else
               "ambiguous (mom and grandma both heterozygous)" if gma else "ambiguous (mom heterozygous, no grandma call)")
        by.setdefault(why, []).append(rs)
    return "; ".join(f"{', '.join(rs)} {why}" for why, rs in by.items())


def red_haired_daughter():
    """True when grandpa is known to have a red-haired daughter (mom's sister): family.json
    "maternal_grandfather": {"red_haired_daughter": true}, or a maternal aunt with red hair in raw/family_notes.json."""
    if (FAMILY.get("maternal_grandfather") or {}).get("red_haired_daughter"):
        return True
    path = os.path.join(RAW, "family_notes.json")
    aunt = (json.load(open(path, encoding="utf-8")) if os.path.exists(path) else {}).get("relatives", {}).get("aunt_maternal")
    return bool(aunt) and bool(re.search(r"red|auburn|ginger|copper", json.dumps(aunt, ensure_ascii=False), re.I))


def main():
    if not os.path.exists(HAP) or not os.path.exists(os.path.join(RESULTS, "grandpa.json")):
        print("- grandpa cells skipped: no virtual-grandpa haplotype (run grandpa.py; needs mom and grandma)")
        return
    wanted = {r for rsids, _ in CELLS.values() for r in rsids} | OPTIONAL
    wanted |= {rs for panel in CARRIER_PANELS.values() for rs, _, _ in panel}
    wanted |= {rs for snps in PGX_SNPS.values() for rs in snps} | set(MC1R_STRONG) | set(MC1R_WEAK) | {"rs12913832"}
    alleles, method, status = {}, {}, {}
    with gzip.open(HAP, "rt") as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[2] in wanted:
                status[f[2]] = (f[6], f[4])  # unresolved sites too, to say why a cell can't be filled
            if f[2] in wanted and f[5] in "ACGT" and len(f[5]) == 1:
                alleles[f[2]], method[f[2]] = f[5], f[6]
    # Imputed-phase calls from grandpa.py (its MARKERS) for SNPs missing from the chip haplotype.
    with open(os.path.join(RESULTS, "grandpa.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            a = r["grandpa_transmitted_allele"]
            if r["rsid"] in wanted and r["rsid"] not in alleles and a in ("A", "C", "G", "T"):
                alleles[r["rsid"]], method[r["rsid"]] = a, "imputed phase"
    out = {}
    for report, (rsids, fn) in CELLS.items():
        missing = [r for r in rsids if r not in alleles]
        if missing:
            out[report] = {"cell": "?", "detail": f"Not reconstructable: {why_missing(missing, status)}"}
            continue
        how = ", ".join(sorted({method[r] for r in rsids}))
        if report == "Kidney disease (APOL1)":
            how += "; G2 deletion (rs71785313) not reconstructable"
        if report == "Inherited clotting" and "rs1799963" not in alleles:
            how += "; prothrombin G20210A not reconstructable"
        if report == "Saturated fat & weight" and alleles["rs5082"] != "G":
            out[report] = {"cell": "Not sensitive!",
                           "detail": f"Sensitivity needs G/G; his {alleles['rs5082']} copy rules it out"}
            continue
        out[report] = {"cell": fn(alleles), "detail": f"His copy passed to mom ({len(rsids)} SNP{'s' * (len(rsids) > 1)}; {how}); other copy unknown"}
    for gene, report in CARRIER_REPORTS.items():  # a variant on his transmitted copy means at least one copy
        have = [(rs, alt, name) for rs, alt, name in CARRIER_PANELS[gene] if rs in alleles]
        hits = [(rs, name) for rs, alt, name in have if alleles[rs] == alt]
        # No "!": his other copy is unknown (could also carry a variant), and rare-variant chip calls can be false.
        out[report] = ({"cell": f"Carrier ({hits[0][1]})",
                        "detail": f"Variant on the copy he passed to mom ({method[hits[0][0]]}): at least one copy; "
                                  "other copy unknown. Chip calls of rare variants can be false: confirm with a clinical test"}
                       if hits else {"cell": "None" if have else "?",
                                     "detail": f"No variant on his copy passed to mom ({len(have)} of "
                                               f"{len(CARRIER_PANELS[gene])} variant{'s' * (len(CARRIER_PANELS[gene]) != 1)} reconstructable); other copy unknown"})
    out["BRCA founder variants"] = {"cell": "?", "detail": "Indel probes; not reconstructable from mom's chip haplotype"}
    out["Sleep paralysis (23andMe score)"] = {"cell": "?", "detail": "23andMe-only polygenic score; can't be computed"}
    g = json.load(open(os.path.join(RESULTS, "grandpa.json")))
    apoe = g.get("apoe_detail", {}).get("grandpa_transmitted")
    no_call = ": mom's chip no-call" if g.get("apoe_detail", {}).get("mom_chip_missing") else ""
    out["Alzheimer's (APOE)"] = ({"cell": apoe, "detail": f"His copy passed to mom (imputed only{no_call}); other copy unknown"}
                                 if apoe else {"cell": "?", "detail": "Not reconstructable"})
    meds = {}
    european = mostly_european(ancestry_json().get("grandpa_inferred"))
    for gene in PGX_SNPS:
        c = med_cell(gene, alleles)
        if c is None:
            on_chip = [rs for rs in PGX_SNPS[gene] if rs in status]
            meds[gene] = {"cell": "?", "detail": f"Not reconstructable: {why_missing(on_chip, status)}" if on_chip
                          else "Not on mom's chip"}
        elif gene == "G6PD":
            meds[gene] = {"cell": c, "detail": "X-linked: his only copy, fully passed to mom"}
        elif gene == "HLA-B*58:01" and c == "Possible":  # same caveat as pgx.py's B*58:01 phenotype
            why = ("low predictive value in Europeans, most likely a false positive" if european else
                   "validated in East Asians; its accuracy varies by ancestry")
            meds[gene] = {"cell": c + " / ?", "detail": f"His copy passed to mom carries the tag SNP ({why}); other copy "
                          "unknown. HLA-B typing needed before allopurinol"}
        else:
            meds[gene] = {"cell": c + " / ?", "detail": "His copy passed to mom; other copy unknown"}
    out["_meds"] = meds
    # Only sites reconstructed from mom's chip haplotype (grandpa.tsv has no MC1R rows; indels such as
    # N29insA rs312262906 are dropped by the single-base filter above), so report coverage.
    have_r = [rs for rs in MC1R_STRONG if rs in alleles]
    red = [MC1R_STRONG[rs] for rs in have_r if alleles[rs] == {
        "rs1805007": "T", "rs1805008": "T", "rs1805009": "C", "rs11547464": "A", "rs1805006": "A",
        "rs1110400": "C", "rs312262906": "A", "rs201326893": "A"}[rs]]
    have_w = [rs for rs in MC1R_WEAK if rs in alleles]
    weak = [MC1R_WEAK[rs] for rs in have_w if alleles[rs] == {"rs1805005": "T", "rs2228479": "A", "rs885479": "A"}[rs]]
    checked = (f"{len(have_r)} of {len(MC1R_STRONG)} strong (R) and {len(have_w)} of {len(MC1R_WEAK)} weaker (r) "
               "variants reconstructable; imputed data not checked for MC1R")
    # A red-haired daughter (mom's sister) constrains his untransmitted copy; without one it stays unknown.
    aunt = red_haired_daughter()
    if red:
        red_cell = {"cell": "+".join(red) + " / ?",
                    "detail": f"Red variants on the copy he passed to mom ({checked}); other copy unknown"}
    elif weak:
        # A red-haired aunt may have inherited this same copy (with an R from grandma), so his other copy is unconstrained.
        red_cell = {"cell": "+".join(weak) + " (r) / ?",
                    "detail": f"Weaker r variant on the copy he passed to mom ({checked}); "
                              + ("the aunt may have inherited this copy, so nothing can be said about his other copy"
                                 if aunt else "other copy unknown")}
    elif have_r and aunt:
        # Assumes the aunt is a full sister who got an R from grandma; R/r (weak V60L, V92M, R163Q) can also give red
        # hair. Only valid if his transmitted copy has no R or r at the sites checked (else the aunt could share it).
        red_cell = {"cell": "None found / likely variant (R or r)",
                    "detail": f"No red variant (R or r) found on the copy he passed to mom ({checked}); his other copy "
                              "probably carried an MC1R variant, strong R or weaker r (inferred from a red-haired "
                              "daughter, see note below)"}
    elif have_r:
        red_cell = {"cell": "None found / ?",
                    "detail": f"No red variant (R or r) found on the copy he passed to mom ({checked}); other copy unknown"}
    else:
        red_cell = {"cell": "?", "detail": f"Not reconstructable ({checked})"}
    herc2 = alleles.get("rs12913832")
    out["_pigmentation"] = {
        "Eye color (predicted)": {"cell": {"A": "Brown-type / ?", "G": "Blue-type / ?"}.get(herc2, "?"),
                                  "detail": "HERC2 rs12913832 on the copy he passed to mom; the eye model needs both copies"},
        "Red hair (MC1R)": red_cell,
    }
    json.dump(out, open(os.path.join(RESULTS, "grandpa_cells.json"), "w"), ensure_ascii=False, indent=1)
    print(os.path.join(RESULTS, "grandpa_cells.json"))


if __name__ == "__main__":
    main()
