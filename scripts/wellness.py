"""Re-creates 23andMe's single-marker Wellness reports for every family member.

Usage: python3 -I scripts/wellness.py
Writes results/wellness.tsv. Marker logic follows the 23andMe report text
(as archived in raw/<person>/23andme_reports/); genotypes are plus strand, build 37, as in the raw files.
Chip no-calls are filled from the phased download where possible (marked * in genotypes).
Genetic Weight and Sleep Paralysis are 23andMe polygenic models and can't be reproduced (see report_archive.py).
"""
import glob
import os
import sys
import zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dnalib import RAW, RESULTS, has_genotypes, imputed_genotypes, load_person, people, sex as sex_of


def count(gt, allele):
    return gt.count(allele) if gt else None


def alcohol(g):
    n = count(g["rs671"], "A")
    return [ "Unlikely to flush after drinking alcohol", "Likely to flush (one ALDH2*2 copy)",
             "Very likely to flush (two ALDH2*2 copies)"][n]


def caffeine(g):
    n = count(g["rs2472297"], "T") + count(g["rs4410790"], "C")
    return (f"{n} variant{'s' if n != 1 else ''} associated with consuming more caffeine — "
            + ("likely to consume a little less than average" if n == 0 else
               "about average" if n <= 2 else "likely to consume a little more than average"))


def deep_sleep(g):
    n = count(g["rs73598374"], "T")
    return "Not likely to be an especially deep sleeper" if n == 0 else "Likely to be an especially deep sleeper"


def lactose(g):
    n = count(g["rs4988235"], "A")  # A keeps LCT on in adulthood
    return "Likely lactose intolerant" if n == 0 else "Not likely to be lactose intolerant"


def muscle(g):
    n = count(g["rs1815739"], "T")  # T = non-functional ACTN3 (577X)
    return ["Common in elite power athletes (two working copies)",
            "Common in elite power athletes (one working copy)",
            "Common in elite endurance athletes (no working copy)"][n]


def sat_fat(g):
    n = count(g["rs5082"], "G")
    return ("Likely to weigh more on a high saturated-fat diet (G/G)" if n == 2 else
            "Not likely to weigh more on a high saturated-fat diet")


def sleep_move(g):
    n = count(g["rs3923809"], "A")
    return ["Likely to move less than average during sleep", "About average movement during sleep",
            "Likely to move more than average during sleep"][n]


def bitter(g):
    # Taster haplotype PAV = rs713598 G (Pro49), rs1726866 G (Ala262), rs10246939 C (Val296); non-taster AVI = C, A, T.
    # Allele-to-amino-acid mapping verified with Ensembl VEP (TAS2R38 is on the minus strand).
    n = [count(g["rs713598"], "G"), count(g["rs1726866"], "G"), count(g["rs10246939"], "C")]
    if len(set(n)) > 1:
        return f"Rare TAS2R38 haplotype (taster alleles {n}) — intermediate bitter sensitivity likely"
    # Brassica link only (Sandell & Breslin 2006); hop bitterness in beer acts on other receptors (TAS2R1/14/40).
    return ["Non-taster (AVI/AVI): broccoli, Brussels sprouts and other cabbage-family greens likely taste mild",
            "Taster (PAV/AVI): medium bitter sensitivity",
            "Strong taster (PAV/PAV): broccoli, Brussels sprouts and other cabbage-family greens likely taste more "
            "bitter"][n[0]]


def cilantro(g):
    n = count(g["rs72921001"], "A")  # Eriksson 2012: OR 0.81 per A allele for soapy taste
    return ["More likely to find cilantro soapy (C/C)", "Average odds of cilantro tasting soapy",
            "Less likely to find cilantro soapy (A/A)"][n]


def asparagus(g):
    # Eriksson 2010 (Table 9, plus strand G/A): ~57% of G/G, ~71% of A/G and ~74% of A/A report smelling it, so
    # every genotype is more likely than not to smell it; A only raises the odds (23andMe uses this one SNP too).
    n = count(g["rs4481887"], "A")
    return ["Likely can smell asparagus pee (G/G: ~57% of people with this genotype smell it)",
            "Likely can smell asparagus pee (one A: ~71% of people with this genotype smell it)",
            "Likely can smell asparagus pee (A/A: ~74% of people with this genotype smell it)"][n]


def earwax(g):
    # ABCC11 G180R (gene on −): plus-strand T = dry earwax, recessive (Yoshiura 2006); T/T is common in East Asia, rare
    # in Europe. Dry-type also means less underarm odor.
    n = count(g["rs17822931"], "T")
    return ["Wet earwax (C/C)", "Wet earwax (one dry-type copy)", "Dry earwax (T/T)"][n]


def omega3(g):
    n = count(g["rs174537"], "T")  # T = lower FADS1 activity (lower AA, slower ALA → EPA/DHA)
    return ["Efficient omega-3 converter (FADS1 G/G)",
            "Moderate omega-3 conversion (FADS1 G/T)",
            "Slow omega-3 converter (FADS1 T/T): fish/algae EPA-DHA more useful than flax/walnut ALA"][n]


def celiac(g):
    # Monsuur 2008 tag SNPs (plus strand): DQ2.5 = rs2187668 T; DQ8 = rs7454108 C (paper's G/A assay is the
    # minus strand). DQ2.2 = rs2395182 T (major, same strand as the paper) + rs7775228 C (minor; paper's G is the
    # minus strand); rs2395182 G/G or no rs7775228 C excludes DQ2.2. rs4713586 G (A/G on plus, as in the paper)
    # tags DQ4 rather than DQ2.2; genotypes are unphased, so it only adds a caveat. DQ7.5 (rs4639334) isn't on the chip.
    dq25 = count(g["rs2187668"], "T")
    dq8 = count(g["rs7454108"], "C")
    dq22_possible = count(g["rs2395182"], "T") > 0 and count(g["rs7775228"], "C") > 0
    dq22 = "DQ2.2" + (" (or DQ4)" if g.get("rs4713586") and count(g["rs4713586"], "G") else "")
    found = [f"DQ2.5×{dq25}"] * bool(dq25) + [f"DQ8×{dq8}"] * bool(dq8)
    if dq22_possible and not found:
        return f"Possible {dq22} — tag SNPs can't confirm it; no DQ2.5/DQ8"
    if not found:
        return ("No DQ2/DQ8: celiac very unlikely (about 95–99% of people with celiac carry DQ2 or DQ8; "
                "DQ7.5 not tested)")
    return (f"Carries {', '.join(found)}{' + possible ' + dq22 if dq22_possible else ''}: celiac possible but uncommon "
            "(~30–40% of people carry DQ2/DQ8; few develop celiac)")


def amd(g):
    # CFH Y402H rs1061170 C = risk; ARMS2 A69S rs10490924 T = risk (both plus strand, the two major AMD loci).
    c, a = count(g["rs1061170"], "C"), count(g["rs10490924"], "T")
    total = c + a
    # Published ORs are vs. zero risk alleles; ~1.1 risk alleles is typical in Europeans (CFH C ~0.36, ARMS2 T ~0.20).
    level = {0: "Lower than average", 1: "About average (one major risk allele; about 1 is typical in Europeans)",
             2: "Somewhat raised (two major risk alleles)"}.get(total, "Raised (several major risk alleles)")
    return f"{level} — CFH {c}, ARMS2 {a} risk copies; not smoking matters most"


def mthfr(g):
    # C677T rs1801133: plus A = 677T. A1298C rs1801131: plus G = 1298C. MTHFR is on the minus strand.
    t677, c1298 = count(g["rs1801133"], "A"), count(g["rs1801131"], "G")
    if t677 == 2:
        return "677TT: ~30% enzyme activity — homocysteine runs higher when folate is low; eat folate-rich foods"
    if t677 == 1 and c1298 >= 1:
        return (f"677CT + 1298{'AC' if c1298 == 1 else 'CC'}: mildly reduced activity — adequate folate keeps "
                "homocysteine normal")
    if c1298 == 2:  # van der Put 1998: ~60% activity, similar to the 677CT + 1298AC compound heterozygote
        return "1298CC: mildly reduced activity (~60%) — homocysteine usually normal; adequate folate is enough"
    if t677 == 1 or c1298 == 1:
        return "One reduced copy: minor effect, no action needed"
    return "Typical MTHFR activity"


def vitamin_d(g):
    n = count(g["rs2282679"], "G")  # GC (vitamin D binding protein); plus G = lower 25(OH)D (Wang 2010)
    return ["Typical vitamin D levels", "Slightly lower vitamin D tendency (one GC copy)",
            "Lower vitamin D tendency (two GC copies) — worth checking levels in winter"][n]


def vitamin_b12(g):
    n = count(g["rs602662"], "G")  # FUT2; A = higher B12 (non-secretor-linked), G = lower (Hazra 2008)
    return ["Higher B12 tendency (FUT2 A/A)", "Slightly lower B12 tendency (FUT2 A/G)",
            "Lower B12 tendency (FUT2 G/G) — matters most with little meat/dairy"][n]


def vitamin_a(g):
    n = count(g["rs7501331"], "T")  # BCMO1 A379V; T = reduced β-carotene → vitamin A conversion (Leung 2009)
    return ["Typical β-carotene → vitamin A conversion", "Somewhat reduced β-carotene conversion",
            "Reduced β-carotene conversion — animal sources of vitamin A are more reliable"][n]


def hemochromatosis(g):
    # HFE (plus strand, gene on +): rs1800562 A = C282Y, rs1799945 G = H63D (Ensembl VEP).
    c282y, h63d = count(g["rs1800562"], "A"), count(g["rs1799945"], "G")
    if c282y == 2:
        return "C282Y/C282Y: at risk of iron overload — ask a doctor for ferritin + transferrin saturation"
    if c282y == 1 and h63d >= 1:
        return "C282Y/H63D: mildly increased iron-overload risk — a ferritin check is reasonable"
    if h63d == 2:
        return "H63D/H63D: slightly raised iron levels possible; overload is rare"
    if c282y == 1:
        return "C282Y carrier: no disease risk from this alone"
    if h63d == 1:
        return "H63D carrier: no disease risk from this alone"
    return "No HFE variants"


def clotting(g, sex=None):
    # F5 Leiden rs6025 T (R534Q, gene on −); F2 G20210A rs1799963 A (gene on +, 3'UTR).
    # Either SNP can be missing (rs1799963 has no regular probe on some chips): report what was called.
    fvl, pt = count(g["rs6025"], "T"), count(g["rs1799963"], "A")
    found = [f"Factor V Leiden ×{fvl}"] * bool(fvl) + [f"Prothrombin G20210A ×{pt}"] * bool(pt)
    untested = ["Factor V Leiden"] * (fvl is None) + ["prothrombin G20210A"] * (pt is None)
    caveat = f" ({' and '.join(untested)} not tested)" if untested else ""
    if not found:
        if fvl is None and pt is None:
            return "Not determinable (no-call)"
        return ("No Factor V Leiden or prothrombin variant" if not untested else
                ("No Factor V Leiden variant (prothrombin G20210A not tested)" if pt is None
                 else "No prothrombin G20210A variant (F5 Leiden probe not tested)"))
    when = {"F": "estrogen pills/HRT, surgery, pregnancy, long flights",
            "M": "surgery, immobilisation, long flights"}.get(sex, "surgery, long flights and, in women, estrogen "
                                                                  "pills/HRT and pregnancy")
    return ", ".join(found) + ": raised clot risk — matters for " + when + caveat


def alpha1(g):
    # SERPINA1 (gene on −): PiZ rs28929474 T (E366K, legacy E342K); PiS rs17580 A (E288V, legacy E264V).
    z, s = count(g["rs28929474"], "T"), count(g["rs17580"], "A")
    alleles = ["Z"] * z + ["S"] * s
    alleles += ["M"] * (2 - len(alleles))
    geno = "Pi" + "".join(sorted(alleles, key="MSZ".index))
    if geno in ("PiZZ", "PiSZ"):
        return f"{geno}: alpha-1 antitrypsin deficiency risk — lung/liver; never smoke, tell a doctor"
    if geno in ("PiMZ", "PiSS", "PiMS"):
        return f"{geno}: carrier/mildly low levels — avoid smoking and dusty air; no disease by itself"
    return f"{geno}: normal alpha-1 antitrypsin"


def lpa(g):
    # LPA (gene on −): rs10455872 G and rs3798220 C (I4399M) are the plus-strand risk alleles (Clarke 2009, Ensembl VEP).
    n = count(g["rs10455872"], "G") + count(g["rs3798220"], "C")
    if n:
        return (f"LPA risk variant ×{n}: likely high Lp(a), an inherited heart-risk factor normal cholesterol tests "
                "miss — ask for a one-time Lp(a) blood test")
    return ("No LPA risk variant: lowers the odds of high Lp(a) but doesn't rule it out (many people with high Lp(a) "
            "carry neither variant; levels are mainly set by the LPA kringle IV-2 repeat number) — a one-time blood "
            "test is still recommended")


def hla_b27(g):
    # rs4349859 A tags HLA-B*27 (Evans 2011: 98% sensitivity, 99% specificity in Europeans; plus strand G/A per
    # Ensembl GRCh37); rs116488202 T is a near-perfect proxy (Cortes 2013) used as a cross-check.
    n, check = count(g["rs4349859"], "A"), count(g["rs116488202"], "T")
    if check is not None and bool(n) != bool(check):
        return "Inconclusive: the two HLA-B27 tag SNPs disagree — clinical HLA-B27 typing if it matters"
    if n:
        return (f"HLA-B27 likely (×{n}): carried by ~8% of Europeans, most never develop spine arthritis; it matters "
                "if there is chronic inflammatory back pain or eye inflammation (uveitis)")
    return ("HLA-B27 not detected: ankylosing spondylitis much less likely (the tag misses some non-European B27 "
            "types; clinical typing if symptoms)")


def parkinsons(g):
    # LRRK2 G2019S rs34637584 A (penetrance ~25–40% by 80: Marder 2015, Lee 2017; clinic-based Healy 2008 higher);
    # GBA N370S at i4000415 C (gene on −, plus C; PD OR ~2–3 for this mild variant, Gan-Or 2015; ~5 for GBA variants
    # overall, Sidransky 2009). 23andMe's Parkinson's report tests exactly these two. No-calls are skipped.
    found = [f"LRRK2 G2019S×{count(g['rs34637584'], 'A')}"] * bool(count(g["rs34637584"], "A")) + \
            [f"GBA N370S×{count(g['i4000415'], 'C')}"] * bool(count(g["i4000415"], "C"))
    if found:
        return (f"Variant detected ({', '.join(found)}): higher Parkinson's risk, but most carriers never develop it — "
                "confirm with a clinical test")
    return "Variants not detected (LRRK2 G2019S, GBA N370S); most Parkinson's has no single-gene cause"


def ttr(g):
    # TTR (gene on +): V30M rs28933979 A, V122I rs76992529 A, T60A rs121918070 G — the three variants 23andMe tests.
    # V30M: founder in Portugal, Sweden, Japan (also Majorca/Cyprus); T60A: Irish founder; V122I mostly African ancestry.
    # The Italian founder variants (F64L, E89Q, I68L) are not on 23andMe's panel.
    found = [name for rs, alt, name in (("rs28933979", "A", "V30M"), ("rs76992529", "A", "V122I"),
                                        ("rs121918070", "G", "T60A")) if count(g[rs], alt)]
    if found:
        return (f"Variant detected ({', '.join(found)}): risk of hereditary amyloidosis (nerve/heart), incomplete "
                "penetrance — confirm with a clinical test; treatments exist")
    return "Variants not detected (V30M, V122I, T60A; the Italian founder variants F64L, E89Q, I68L are not tested)"


def apol1(g):
    # APOL1 G1 tag rs73885319 G (S342G; partner SNP rs60910145 not on the chip) and G2 rs71785313 deletion (D).
    # Risk is recessive: two risk alleles raise kidney-disease odds ~7–10× (Genovese 2010); almost only in
    # African ancestry.
    n = (count(g["rs73885319"], "G") or 0) + (count(g["rs71785313"], "D") or 0)
    if n >= 2:
        return "Two APOL1 risk variants: higher kidney-disease risk — check blood pressure, urine protein and kidney function"
    if n == 1:
        return "One APOL1 risk variant: no meaningful added kidney risk on its own"
    if g["rs71785313"] is None:  # G2 deletion probe absent (v4 chip)
        return "Variants not detected (G1 only; G2 not tested on this chip)"
    return "Variants not detected (G1, G2)"


def hoxb13(g, sex=None):
    # HOXB13 G84E rs138213197 (gene on −; plus T): population OR ~3–5 (Karlsson 2014; meta-analyses), ~20× enrichment
    # in familial/early-onset cases (Ewing 2012); lifetime risk for male carriers roughly 33–60%. Women can pass it on.
    n = count(g["rs138213197"], "T")
    if n and sex == "F":
        return (f"G84E detected (×{n}): no prostate cancer risk for her; relevant to male relatives (father, "
                "brothers, sons), who can be offered testing and an earlier PSA discussion; confirm with a clinical test")
    if n:
        return (f"G84E detected (×{n}): " + ("men have" if sex == "M" else "male carriers have")
                + " ~3–5× prostate cancer odds — earlier PSA discussion; confirm with a clinical test")
    return "G84E not detected"


# Carrier panels: (probe, plus-strand variant allele, name). Identities and alleles verified against Ensembl
# GRCh37/VEP and ClinVar; i-probes are 23andMe internal IDs (F508del and 35delG are I/D-coded,
# i3003137 is the A/T HbS probe at rs334). Other HBB/CFTR/GJB2/MEFV probes on the chip were unconfirmed and are skipped:
# the next most common Southern European CF variants (G542X, N1303K, W1282X, R1162X, 1717-1G>A, 2183AA>G) sit at
# internal i-probes whose reported alleles couldn't be checked, so CF coverage is F508del plus rarer variants.
CARRIER_PANELS = {
    "HBB": [("rs11549407", "A", "β-thal codon 39"), ("rs35724775", "G", "β-thal IVS-I-6"),
            ("rs35004220", "T", "β-thal IVS-I-110"), ("rs33971440", "T", "β-thal IVS-I-1"),
            ("rs34690599", "C", "β-thal IVS-II-745"), ("i3003137", "A", "HbS (sickle)"),
            ("rs33930165", "T", "HbC")],
    "CFTR": [("i3000001", "D", "F508del"), ("i4000292", "D", "I507del"), ("rs78655421", "A", "R117H"),
             ("rs121909005", "G", "S549R"), ("i4000299", "T", "V520F"),
             ("rs80224560", "A", "2789+5G>A")],  # c.2657+5G>A; CFTR is on the + strand (Ensembl VEP, ClinVar)
    "GJB2": [("i4000434", "D", "35delG"), ("rs104894396", "T", "W24X"), ("rs111033217", "G", "K15T"),
             ("rs35887622", "G", "M34T (mild)"), ("rs72474224", "T", "V37I (mild)")],
    "MEFV": [("rs61752717", "C", "M694V"), ("rs28940580", "G", "M680I"), ("rs28940579", "G", "V726A")],
    # Extra 23andMe carrier conditions, several common in Southern European populations (checked: rsID ↔ probe and
    # allele sets vs Ensembl GRCh37; plus-strand alleles from ClinVar/VEP; minus-strand genes complemented).
    "PAH": [("rs5030855", "T", "IVS10-11G>A"), ("rs5030858", "A", "R408W"), ("rs5030849", "T", "R261Q"),
            ("i3003399", "A", "R252W"), ("i3003401", "A", "P281L"), ("i4000472", "G", "I65T"),
            ("rs5030841", "G", "L48S"), ("i3003397", "T", "R158Q"), ("i3003405", "C", "Y414C"),
            ("i4000478", "T", "E280K"), ("rs5030861", "T", "IVS12+1G>A"), ("i3003403", "C", "E390G"),
            ("i4000467", "A", "G272X")],
    "GBA": [("i4000415", "C", "N370S")],
    "ALDOB": [("rs1800546", "G", "A149P"), ("rs76917243", "T", "A174D"), ("i5012664", "C", "N334K")],
    "ACADM": [("i5012759", "G", "K304E")],
    "GAA": [("rs386834236", "G", "IVS1-13T>G")],
    "PMM2": [("i5012680", "A", "R141H"), ("i5012679", "A", "F119L")],
    "SLC26A4": [("i5012616", "C", "L236P"), ("i5012618", "C", "T416P"), ("rs80338849", "A", "IVS8+1G>A"),
                ("i5000002", "G", "H723R")],
    "PKLR": [("rs116100695", "A", "R486W")],
}
CARRIER_ONE = {  # one copy (matched on variant name or gene): usually healthy, with these exceptions
    "β-thal": "usually healthy, but β-thal carriers have small red cells (thalassemia minor) that can be mistaken "
              "for iron deficiency",
    "HbS": "usually healthy; sickle trait carries rare risks with extreme exertion or dehydration",
    "MEFV": "usually healthy, though one copy (especially M694V) occasionally causes mild FMF-like fevers",
    "GBA": "no Gaucher disease from one copy, but N370S carriers have roughly 2–3× Parkinson's odds (most never develop it)",
}
CARRIER_TWO = {  # what two variant copies (same or different variants) can mean
    "HBB": "two variants: possible thalassemia/sickle disease or a variant combination — see a doctor",
    "CFTR": "two variants: possible cystic fibrosis or a CFTR-related condition — see a doctor",
    "GJB2": "two variants: possible GJB2 hearing loss (mild for M34T/V37I) — see a doctor",
    "MEFV": "two variants: possible Familial Mediterranean Fever (recurrent fevers with pain) — see a doctor",
    "PAH": "two variants: possible PKU/high phenylalanine — see a doctor",
    "GBA": "two variants (N370S ×2): possible Gaucher disease (often mild) — see a doctor",
    "ALDOB": "two variants: possible hereditary fructose intolerance — see a doctor",
    "ACADM": "two variants: possible MCAD deficiency — see a doctor",
    "GAA": "two variants: possible Pompe disease — see a doctor",
    "PMM2": "two variants: possible PMM2-CDG — see a doctor",
    "SLC26A4": "two variants: possible Pendred syndrome/hearing loss — see a doctor",
    "PKLR": "two variants: possible pyruvate kinase deficiency — see a doctor",
}


def carrier(gene):
    def fn(g):
        hits = [f"{name}{' ×2' if count(g[rs], alt) == 2 else ''}"
                for rs, alt, name in CARRIER_PANELS[gene] if g[rs] and count(g[rs], alt)]
        copies = sum(count(g[rs], alt) for rs, alt, _ in CARRIER_PANELS[gene] if g[rs])
        tested = sum(1 for rs, _, _ in CARRIER_PANELS[gene] if g[rs])
        if copies >= 2:
            return f"{', '.join(hits)} — {CARRIER_TWO[gene]}"
        if copies == 1:
            one = next((v for k, v in CARRIER_ONE.items() if k in hits[0] or k == gene),
                       "usually no disease from one copy")
            return (f"Carrier: {hits[0]} — {one}; matters for children only if a partner is also "
                    "a carrier. Chip calls of rare variants can be false: confirm with a clinical test")
        return f"Not a carrier of the {tested} variant{'s' * (tested != 1)} tested (only the variants listed; others not tested)"
    return fn


CONDITIONS = [
    ("Hemochromatosis (iron)", "HFE", ["rs1800562", "rs1799945"], hemochromatosis),
    ("Inherited clotting", "F5 / F2", ["rs6025", "rs1799963"], clotting),
    ("Alpha-1 antitrypsin", "SERPINA1", ["rs28929474", "rs17580"], alpha1),
    ("Lipoprotein(a)", "LPA", ["rs10455872", "rs3798220"], lpa),
    ("HLA-B27 (spine arthritis)", "HLA-B", ["rs4349859", "rs116488202"], hla_b27),
    ("Parkinson's genes", "LRRK2 / GBA", ["rs34637584", "i4000415"], parkinsons),
    ("Hereditary amyloidosis", "TTR", ["rs28933979", "rs76992529", "rs121918070"], ttr),
    ("Kidney disease (APOL1)", "APOL1", ["rs73885319", "rs71785313"], apol1),
    ("Hereditary prostate cancer", "HOXB13", ["rs138213197"], hoxb13),
]

CARRIERS = [
    ("Beta-thalassemia & sickle cell", "HBB", [rs for rs, _, _ in CARRIER_PANELS["HBB"]], carrier("HBB")),
    ("Cystic fibrosis", "CFTR", [rs for rs, _, _ in CARRIER_PANELS["CFTR"]], carrier("CFTR")),
    ("Inherited hearing loss", "GJB2", [rs for rs, _, _ in CARRIER_PANELS["GJB2"]], carrier("GJB2")),
    ("Familial Mediterranean Fever", "MEFV", [rs for rs, _, _ in CARRIER_PANELS["MEFV"]], carrier("MEFV")),
    ("Phenylketonuria (PKU)", "PAH", [rs for rs, _, _ in CARRIER_PANELS["PAH"]], carrier("PAH")),
    ("Gaucher disease", "GBA", [rs for rs, _, _ in CARRIER_PANELS["GBA"]], carrier("GBA")),
    ("Hereditary fructose intolerance", "ALDOB", [rs for rs, _, _ in CARRIER_PANELS["ALDOB"]], carrier("ALDOB")),
    ("MCAD deficiency (fat metabolism)", "ACADM", [rs for rs, _, _ in CARRIER_PANELS["ACADM"]], carrier("ACADM")),
    ("Pompe disease", "GAA", [rs for rs, _, _ in CARRIER_PANELS["GAA"]], carrier("GAA")),
    ("PMM2-CDG (metabolic disorder)", "PMM2", [rs for rs, _, _ in CARRIER_PANELS["PMM2"]], carrier("PMM2")),
    ("Pendred syndrome (hearing)", "SLC26A4", [rs for rs, _, _ in CARRIER_PANELS["SLC26A4"]], carrier("SLC26A4")),
    ("Pyruvate kinase deficiency", "PKLR", [rs for rs, _, _ in CARRIER_PANELS["PKLR"]], carrier("PKLR")),
]

NUTRITION = [
    ("Macular degeneration", "CFH / ARMS2", ["rs1061170", "rs10490924"], amd),
    ("Folate / homocysteine", "MTHFR", ["rs1801133", "rs1801131"], mthfr),
    ("Vitamin D", "GC", ["rs2282679"], vitamin_d),
    ("Vitamin B12", "FUT2", ["rs602662"], vitamin_b12),
    ("Vitamin A (from β-carotene)", "BCMO1", ["rs7501331"], vitamin_a),
]

TASTE = [
    ("Celiac genes (HLA-DQ2/DQ8)", "HLA-DQA1/DQB1", ["rs2187668", "rs7454108", "rs2395182", "rs7775228", "rs4713586"],
     celiac),
    ("Bitter taste", "TAS2R38", ["rs713598", "rs1726866", "rs10246939"], bitter),
    ("Cilantro soapiness", "OR6A2", ["rs72921001"], cilantro),
    ("Asparagus smell", "OR2M7 region", ["rs4481887"], asparagus),
    ("Omega-3 conversion", "FADS1", ["rs174537"], omega3),
    ("Earwax type", "ABCC11", ["rs17822931"], earwax),
]

REPORTS = [
    ("Alcohol flush", "ALDH2", ["rs671"], alcohol),
    ("Caffeine consumption", "CYP1A2 / AHR", ["rs2472297", "rs4410790"], caffeine),
    ("Deep sleep", "ADA", ["rs73598374"], deep_sleep),
    ("Lactose intolerance", "LCT", ["rs4988235"], lactose),
    ("Muscle composition", "ACTN3", ["rs1815739"], muscle),
    ("Saturated fat & weight", "APOA2", ["rs5082"], sat_fat),
    ("Sleep movement", "BTBD9", ["rs3923809"], sleep_move),
]


def phased_calls(person, rsids):
    """Fill chip no-calls from the 23andMe phased download (same build, rsid-keyed)."""
    out = {}
    for path in glob.glob(os.path.join(RAW, person, "phased_genotype_data_*.zip")):
        z = zipfile.ZipFile(path)
        name = next((n for n in z.namelist() if n.startswith("phased_genome_")
                     and "statistical" not in n and "parent" not in n), None)
        if name is None:  # an unexpected layout: skip this optional cross-check rather than stop the run
            print(f"note: no plain phased_genome_* file in {os.path.basename(path)}; skipped")
            continue
        with z.open(name) as fh:
            for line in fh:
                f = line.decode().strip().split("\t")
                if len(f) == 5 and f[0] in rsids and f[3] in "ACGT" and f[4] in "ACGT" and f[3] and f[4]:
                    out[f[0]] = "".join(sorted(f[3] + f[4]))
    return out


# Caveat-only or cross-check SNPs: a no-call doesn't block the report (hla_b27 copes without its cross-check).
OPTIONAL = {"rs4713586", "rs116488202"}
# "Any variant detected" reports: like carrier panels, report on whichever probes were called (their functions name
# what wasn't tested). A missing prothrombin probe must not hide a Factor V Leiden result.
ANY_CALLED = {"Parkinson's genes", "Hereditary amyloidosis", "Kidney disease (APOL1)", "Inherited clotting"}


def report_ok(category, name, g):
    """Whether report `name` can be scored from genotypes g ({rsid: genotype or None})."""
    return (all(v for rs, v in g.items() if rs not in OPTIONAL)
            or ((category == "Carrier status" or name in ANY_CALLED) and any(g.values())))
BY_SEX = {clotting, hoxb13}  # wording depends on family.json sex (None = unknown: neutral wording)


def main():
    rows = []
    sets = [("23andMe wellness", REPORTS), ("Taste & diet", TASTE), ("Nutrition & health", NUTRITION),
            ("Inherited conditions", CONDITIONS), ("Carrier status", CARRIERS)]
    all_snps = {rs for _, reports in sets for _, _, snps, _ in reports for rs in snps}
    for p in people():
        if not has_genotypes(p):
            continue
        data, _ = load_person(p)
        missing = {rs for rs in all_snps if rs not in data}
        filled = phased_calls(p, missing) if missing else {}
        filled.update(imputed_genotypes(p, missing - filled.keys()))  # chip lacks the probe (v4) or uses an internal ID
        for category, name, gene, snps, fn in [(c, *r) for c, reports in sets for r in reports]:
            g = {rs: data[rs][2] if rs in data else filled.get(rs) for rs in snps}
            genos = " ".join(f"{rs}={g[rs] or '--'}{'*' if rs in filled else ''}" for rs in snps)
            # Carrier panels report on whatever probes were called; other reports need every SNP.
            ok = report_ok(category, name, g)
            result = (fn(g, sex_of(p)) if fn in BY_SEX else fn(g)) if ok else "Not determinable (no-call)"
            rows.append((p, category, name, gene, genos, result))
            print(f"{p:17s} {name:24s} {genos:34s} {result}")
    with open(os.path.join(RESULTS, "wellness.tsv"), "w") as fh:
        fh.write("person\tcategory\treport\tgene\tgenotypes\tresult\n")
        fh.writelines("\t".join(r) + "\n" for r in rows)


if __name__ == "__main__":
    main()
