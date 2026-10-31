"""Eye and hair colour prediction with the published IrisPlex / HIrisPlex models.

Run:  python3 -I scripts/pigmentation.py
Out:  results/pigmentation.json, results/pigmentation.tsv

MODEL SOURCE (every coefficient below)
--------------------------------------
Walsh et al. (2013) "The HIrisPlex system for simultaneous prediction of hair and
eye colour from DNA", Forensic Sci. Int. Genet. 7:98-115, doi:10.1016/j.fsigen.2012.07.005.
The multinomial logistic regression parameters were read out of that paper's own
supplementary prediction tool, Supplementary file mmc7.xls ("HIrisPlex hair & eye
colour prediction tool"), downloaded from
https://ars.els-cdn.com/content/image/1-s2.0-S1872497312001810-mmc7.xls
 * sheet 1  = user input sheet, giving the counted ("minor") allele per SNP,
 * sheet 2  = the model sheet: hair constants + betas (Brown/Red/Black, reference
              category Blond), hair-shade constants + betas (Dark, reference Light),
              and the IrisPlex eye model constants + betas (Blue/Intermediate,
              reference Brown; the eye model is Walsh et al. 2011, FSIG 5:170-180).
The parsed constants were verified against the workbook's own cached probabilities
for an all-zero input (hair 0.1029/0.0014/0.0743/0.8214, shade 0.0739/0.9261,
eye 0.9500/0.0296/0.0204) -- softmax of the constants reproduces them exactly.

Hair colour wording follows the HIrisPlex prediction guide (introduced in Walsh et
al. 2013, FSIG 7:98-115; reproduced in the HIrisPlex-S webtool user manual v2.0, 2018,
https://hirisplex.erasmusmc.nl/pdf/hirisplex.erasmusmc.nl.pdf p.15), in the
machine-readable form published as classifHISplex in Maroti et al.,
https://github.com/zmaroti/aHISplex (cmd/classifHISplex/main.go).

SKIN: the HIrisPlex-S 5-category skin model intercepts are not published.  Walsh et
al. (2017) "Global skin colour prediction from DNA", Hum Genet 136:847-863
(PMC5487854) Table 2 prints only per-category beta coefficients for the 36 skin
SNPs, with no intercepts, and neither that paper's supplements nor the HIrisPlex-S
paper (Chaitanya et al. 2018, FSIG 35:123-135) supplements contain the model.  The
webtool computes server-side.  So no skin probabilities are produced here; the
script reports the genotypes of the HIrisPlex-S skin markers it can read instead.

ALLELE ORIENTATION
------------------
Raw data is GRCh37 plus strand; the models count their own effect allele, which for
several SNPs is written on the minus strand.  The plus-strand equivalent of every
effect allele is given in EFFECT below, with the strand flag taken from the
peer-reviewed translation table of Maroti et al. (aHISplex_data/trans_GRCh37.tsv,
"F" = as written, "R" = complement), cross-checked against the plus-strand allele
string from Ensembl (ALLELES below, REST /variation/homo_sapiens, GRCh38 mapping;
same strand as GRCh37 here).  The two palindromic SNPs were resolved independently:
 * rs16891982 (C/G): counted allele C.  1000G phase 3 frequencies (Ensembl) give
   plus-strand C = 0.062 EUR / 0.964 AFR, i.e. plus C is the dark/ancestral allele,
   which is what the model's negative blue-eye beta (-1.5257) requires.  No flip.
 * rs1805009 (G/C): counted allele C.  Ensembl VEP: plus-strand C = MC1R D294H, the
   red-hair variant the model scores (+4.59 red); plus C = 0.008 EUR.  No flip.
Non-palindromic SNPs were matched by letter or complement and every effect allele
was asserted to be a member of the Ensembl plus-strand allele set (see check below).

MISSING DATA: the webtool re-fits its models on whatever subset of SNPs is supplied,
which cannot be reproduced without the training data, so a SNP that cannot be read
at all is simply dropped from the linear predictor (equivalent to the webtool's "NA"
input) and listed in *_snps_missing, with the model's own hard rules applied: no eye
prediction without HERC2 rs12913832, no hair prediction without any MC1R marker.

GENOTYPE SOURCES, in order: measured chip (dnalib.load_person) -> 23andMe phased
download -> 23andMe imputed GRCh38 BCF (FORMAT/HDS haploid dosages).  The source
used is recorded per SNP; imputed dosages are fractional and flagged.
"""
import glob
import json
import os
import subprocess
import sys
import zipfile
from math import exp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dnalib import (BCFTOOLS, RAW, RESULTS, WORK, children, has_genotypes, load_person,  # noqa: E402
                    person, siblings)

# Report order: self, the parents and maternal grandmother, then siblings (missing roles drop out).
PEOPLE = [p for p in (person("self"), person("mother"), person("father"), person("maternal_grandmother"),
                      *siblings()) if p]
# child -> parents checked for Mendelian consistency (missing parents are dropped in family_check).
PARENTS = {c: tuple(p for p in (person("mother"), person("father")) if p) for c in children()}
if person("mother") and person("maternal_grandmother"):
    PARENTS[person("mother")] = (person("maternal_grandmother"),)

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}

# rsid -> (model effect allele, strand flag vs GRCh37 plus, plus-strand effect allele)
# flag from aHISplex trans_GRCh37.tsv (F = same strand, R = complement).
EFFECT = {
    "rs312262906": ("A", "F", "A"),   # MC1R N29insA (insertion of A)
    "rs11547464":  ("A", "F", "A"),
    "rs885479":    ("T", "R", "A"),
    "rs1805008":   ("T", "F", "T"),
    "rs1805005":   ("T", "F", "T"),
    "rs1805006":   ("A", "F", "A"),
    "rs1805007":   ("T", "F", "T"),
    "rs1805009":   ("C", "F", "C"),   # palindromic, resolved by VEP + frequency
    "rs201326893": ("A", "F", "A"),   # MC1R Y152OCH
    "rs2228479":   ("A", "F", "A"),
    "rs1110400":   ("C", "F", "C"),
    "rs28777":     ("C", "F", "C"),
    "rs16891982":  ("C", "F", "C"),   # palindromic, resolved by frequency
    "rs12821256":  ("G", "R", "C"),
    "rs4959270":   ("A", "F", "A"),
    "rs12203592":  ("T", "F", "T"),
    "rs1042602":   ("T", "R", "A"),
    "rs1800407":   ("A", "R", "T"),
    "rs2402130":   ("G", "F", "G"),
    "rs12913832":  ("T", "R", "A"),
    "rs2378249":   ("C", "R", "G"),
    "rs683":       ("G", "R", "C"),
    "rs12896399":  ("G", "F", "G"),   # eye model counts G (webtool counts T)
    "rs1393350":   ("T", "R", "A"),
}

# Ensembl plus-strand allele strings (GRCh38), used to validate EFFECT.
ALLELES = {
    "rs11547464": "GA", "rs885479": "GAC", "rs1805008": "CT", "rs1805005": "GT",
    "rs1805006": "CAG", "rs1805007": "CAGT", "rs1805009": "GAC", "rs201326893": "CA",
    "rs2228479": "GAC", "rs1110400": "TC", "rs28777": "CAG", "rs16891982": "CAG",
    "rs12821256": "TACG", "rs4959270": "CAGT", "rs12203592": "CGT", "rs1042602": "CA",
    "rs1800407": "CT", "rs2402130": "GACT", "rs12913832": "ACG", "rs2378249": "GAC",
    "rs683": "CA", "rs12896399": "GAT", "rs1393350": "GA",
}

# GRCh38 positions (aHISplex trans_GRCh38.tsv) for the imputed-BCF fallback.
GRCH38 = {
    "rs312262906": ("16", 89919342), "rs11547464": ("16", 89919683),
    "rs885479": ("16", 89919746), "rs1805008": ("16", 89919736),
    "rs1805005": ("16", 89919436), "rs1805006": ("16", 89919510),
    "rs1805007": ("16", 89919709), "rs1805009": ("16", 89920138),
    "rs201326893": ("16", 89919714), "rs2228479": ("16", 89919532),
    "rs1110400": ("16", 89919722), "rs28777": ("5", 33958854),
    "rs16891982": ("5", 33951588), "rs12821256": ("12", 88934558),
    "rs4959270": ("6", 457748), "rs12203592": ("6", 396321),
    "rs1042602": ("11", 89178528), "rs1800407": ("15", 27985172),
    "rs2402130": ("14", 92334859), "rs12913832": ("15", 28120472),
    "rs2378249": ("20", 34630286), "rs12896399": ("14", 92307319),
    "rs1393350": ("11", 89277878), "rs683": ("9", 12709305),
    "rs1426654": ("15", 48134287), "rs1800414": ("15", 27951891),
    "rs3212355": ("16", 89917970), "rs1129038": ("15", 28111713),
    "rs1667394": ("15", 28285036), "rs2238289": ("15", 28208069),
    "rs6497292": ("15", 28251049), "rs1470608": ("15", 28042975),
    "rs1545397": ("15", 27942626), "rs12441727": ("15", 28026629),
    "rs17128291": ("14", 92416482), "rs1126809": ("11", 89284793),
    "rs10756819": ("9", 16858086), "rs3114908": ("16", 89317317),
    "rs8051733": ("16", 89957798), "rs6119471": ("20", 34197406),
    "rs6059655": ("20", 34077942),
}

# --- models (Walsh et al. 2013 supplementary tool mmc7.xls, sheet 2) -----------
EYE_CATS = ["blue", "intermediate"]          # reference category: brown
EYE_CONST = [3.8402, 0.372]
EYE_BETA = {                                  # rsid: [blue, intermediate]
    "rs12913832": [-4.8727, -1.99],
    "rs1800407": [1.1536, 1.052],
    "rs12896399": [-0.5253, -0.0077],
    "rs16891982": [-1.5257, -0.7363],
    "rs1393350": [0.4431, 0.259],
    "rs12203592": [0.5951, 0.6872],
}

HAIR_CATS = ["brown", "red", "black"]        # reference category: blond
HAIR_CONST = [-2.0769, -6.3953, -2.4029]
HAIR_BETA = {                                 # rsid: [brown, red, black]
    "rs312262906": [1.2522, 25.508, 2.732],
    "rs11547464": [-0.61155, 2.5381, -16.969],
    "rs885479": [0.2937, -0.20889, 0.39983],
    "rs1805008": [-0.50143, 2.801, -0.86062],
    "rs1805005": [0.21172, 0.93493, -0.0029013],
    "rs1805006": [1.9293, 3.65, -16.088],
    "rs1805007": [-0.32318, 3.4408, -1.3757],
    "rs1805009": [0.60861, 4.5868, 0.060631],
    "rs201326893": [0.25624, 22.107, 3.9824],
    "rs2228479": [-0.054143, 0.62307, 0.17012],
    "rs1110400": [-0.56315, 1.4453, 0.29143],
    "rs28777": [0.52168, 0.70401, 0.82228],
    "rs16891982": [0.75284, -0.41869, 1.1617],
    "rs12821256": [-0.34957, -0.57964, -0.89824],
    "rs4959270": [-0.19171, 0.24861, -0.36359],
    "rs12203592": [1.6475, 0.90233, 1.997],
    "rs1042602": [0.16092, 0.45003, 0.065432],
    "rs1800407": [-0.19111, -0.27606, -0.49601],
    "rs2402130": [0.35821, 0.28313, 0.26536],
    "rs12913832": [1.214, -0.093776, 1.9391],
    "rs2378249": [0.12669, 0.76634, -0.089509],
    "rs683": [0.21172, -0.053427, 0.15796],
}

SHADE_CONST = [-2.528]                        # dark; reference category: light
SHADE_BETA = {
    "rs11547464": [-16.575], "rs885479": [0.18709], "rs1805008": [-0.9331],
    "rs1805005": [-0.030452], "rs1805006": [-15.305], "rs1805007": [-1.7901],
    "rs1805009": [-0.078426], "rs2228479": [0.0049399], "rs1110400": [-0.15892],
    "rs28777": [1.4594], "rs16891982": [0.78071], "rs12821256": [-0.7757],
    "rs4959270": [-0.44286], "rs12203592": [1.8636], "rs1042602": [0.18179],
    "rs1800407": [-0.59583], "rs2402130": [0.33304], "rs12913832": [1.9622],
    "rs2378249": [0.090083], "rs683": [0.18114],
}

MC1R_STRONG = {  # "R" alleles, strong red-hair effect
    "rs1805007": "R151C", "rs1805008": "R160W", "rs1805009": "D294H",
    "rs11547464": "R142H", "rs1805006": "D84E", "rs1110400": "I155T",
    "rs312262906": "N29insA", "rs201326893": "Y152OCH",
}
MC1R_WEAK = {"rs1805005": "V60L", "rs2228479": "V92M", "rs885479": "R163Q"}

# HIrisPlex-S skin markers that are not already in the eye/hair models; reported
# as genotypes only (no published model to run).
SKIN_EXTRA = ["rs1426654", "rs1800414", "rs3212355", "rs1129038", "rs1667394",
              "rs2238289", "rs6497292", "rs1470608", "rs1545397", "rs12441727",
              "rs17128291", "rs1126809", "rs10756819", "rs3114908", "rs8051733",
              "rs6119471", "rs6059655"]

MODEL_SNPS = sorted(set(EYE_BETA) | set(HAIR_BETA) | set(SHADE_BETA))
ALL_SNPS = sorted(set(MODEL_SNPS) | set(SKIN_EXTRA))


def check_orientation():
    """Every plus-strand effect allele must be one of the plus-strand alleles."""
    for rs, (model_allele, flag, plus) in EFFECT.items():
        assert plus == (model_allele if flag == "F" else COMP[model_allele]), rs
        if rs in ALLELES:
            assert plus in ALLELES[rs], f"{rs}: {plus} not in {ALLELES[rs]}"


def phased_calls(person, rsids):
    """Fill chip no-calls from the 23andMe phased download (same build, rsid-keyed)."""
    out = {}
    for path in glob.glob(os.path.join(RAW, person, "phased_genotype_data_*.zip")):
        z = zipfile.ZipFile(path)
        try:
            name = next(n for n in z.namelist() if n.startswith("phased_genome_")
                        and "statistical" not in n and "parent" not in n)
        except StopIteration:
            continue
        with z.open(name) as fh:
            for line in fh:
                f = line.decode().strip().split("\t")
                if len(f) == 5 and f[0] in rsids and f[3] in "ACGT" and f[4] in "ACGT":
                    out[f[0]] = "".join(sorted(f[3] + f[4]))
    return out


def imputed_dosages(person, rsids):
    """Effect-allele dosage from the imputed GRCh38 BCFs (FORMAT/HDS)."""
    out = {}
    if not os.path.exists(BCFTOOLS):
        return out
    by_chrom = {}
    for rs in rsids:
        if rs in GRCH38 and rs in EFFECT:
            ch, pos = GRCH38[rs]
            by_chrom.setdefault(ch, []).append((pos, rs))
    for ch, items in by_chrom.items():
        paths = glob.glob(os.path.join(WORK, "imputed", person, "*", f"chr{ch}.bcf"))
        if not paths:
            continue
        regions = ",".join(f"chr{ch}:{pos}" for pos, _ in items)
        cmd = [BCFTOOLS, "query", "-f", "%POS %REF %ALT [%HDS]\n", "-r", regions, paths[0]]
        try:
            res = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        except (OSError, subprocess.SubprocessError):
            continue
        rows = {}
        for line in res.stdout.splitlines():
            f = line.split()
            if len(f) == 4:
                rows[int(f[0])] = (f[1], f[2], f[3])
        for pos, rs in items:
            if pos not in rows:
                continue
            ref, alt, hds = rows[pos]
            try:
                alt_dose = sum(float(x) for x in hds.split(","))
            except ValueError:
                continue
            plus = EFFECT[rs][2]
            if rs == "rs312262906":          # insertion: ALT = REF + "A"
                eff_is_alt = len(alt) > len(ref)
            elif plus == alt:
                eff_is_alt = True
            elif plus == ref:
                eff_is_alt = False
            else:
                continue                      # orientation mismatch: do not use
            out[rs] = round(alt_dose if eff_is_alt else 2 - alt_dose, 3)
    return out


def dosages(person):
    """Return {rsid: (dosage, source, genotype_or_None)} for all pigmentation SNPs."""
    chip, stats = load_person(person)
    out = {}
    for rs in ALL_SNPS:
        if rs in chip:
            gt = chip[rs][2]
            if rs == "rs312262906":
                # 23andMe codes this indel as D/I; I = the inserted A the model counts.
                out[rs] = (float(gt.count("I")), "chip", gt) if set(gt) <= {"D", "I"} \
                    else (float(gt.count("A")), "chip", gt)
            elif rs in EFFECT:
                out[rs] = (float(gt.count(EFFECT[rs][2])), "chip", gt)
            else:
                out[rs] = (None, "chip", gt)
    missing = [rs for rs in ALL_SNPS if rs not in out]
    for rs, gt in phased_calls(person, set(missing)).items():
        out[rs] = (float(gt.count(EFFECT[rs][2])) if rs in EFFECT else None, "phased", gt)
    missing = [rs for rs in ALL_SNPS if rs not in out]
    for rs, d in imputed_dosages(person, missing).items():
        out[rs] = (d, "imputed", None)
    return out, stats


def softmax(logits):
    """Probabilities for len(logits) categories plus an implicit reference (logit 0)."""
    vals = [0.0] + list(logits)
    m = max(vals)
    e = [exp(v - m) for v in vals]
    s = sum(e)
    return [x / s for x in e]


def predict(const, beta, cats, ref_cat, dose):
    used, missing = [], []
    logits = list(const)
    for rs, b in beta.items():
        d = dose.get(rs, (None,))[0]
        if d is None:
            missing.append(rs)
            continue
        used.append(rs)
        for i, bi in enumerate(b):
            logits[i] += bi * d
    probs = softmax(logits)
    out = {ref_cat: probs[0]}
    for i, c in enumerate(cats):
        out[c] = probs[i + 1]
    return {k: round(v, 4) for k, v in out.items()}, sorted(used), sorted(missing)


def hair_phrase(hair, shade):
    """HIrisPlex prediction guide (Walsh et al. 2013; HIrisPlex-S manual v2.0), classifHISplex rules."""
    bl, br, rd, bk = hair["blond"], hair["brown"], hair["red"], hair["black"]
    light, dark = shade["light"], shade["dark"]
    if bk > bl and bk >= br and bk > rd:
        if bk > 0.7 or dark > 0.5:
            return "black"
        return "dark brown/black"
    if bl >= br and bl > rd and bl > bk:
        if bl > 0.7:
            return "blond" if light > 0.95 else "blond/dark blond"
        return "blond/dark blond" if light > 0.9 else "dark blond/brown"
    if br >= bl and br > rd and br > bk:
        if br > 0.7:
            return "brown" if light > 0.8 else "brown/dark brown"
        return "brown/dark brown" if light > 0.8 else "dark brown/black"
    if rd > bl and rd > br and rd > bk:
        return "red"
    return "undetermined"


def mc1r_counts(dose):
    strong = weak = 0
    detail = {}
    for rs, name in MC1R_STRONG.items():
        d = dose.get(rs, (None,))[0]
        if d:
            strong += d
            detail[name] = round(d, 2)
    for rs, name in MC1R_WEAK.items():
        d = dose.get(rs, (None,))[0]
        if d:
            weak += d
            detail[name] = round(d, 2)
    return round(strong, 2), round(weak, 2), detail


def mendel_ok(child, parents):
    """True if the child's alleles can come from the parents' genotypes."""
    if any(len(g) != 2 for g in [child] + parents):
        return True
    if len(parents) == 2:
        p, q = parents
        combos = {"".join(sorted(a + b)) for a in p for b in q}
        return child in combos
    p = parents[0]
    return bool(set(child) & set(p))


def family_check(genos):
    rows = []
    for child, pars in PARENTS.items():
        pars = [p for p in pars if p in genos]
        if child not in genos or not pars:
            continue
        bad = []
        for rs in ALL_SNPS:
            c = genos[child].get(rs, (None, None, None))[2]
            ps = [genos[p].get(rs, (None, None, None))[2] for p in pars]
            if not c or any(x is None for x in ps):
                continue
            if not mendel_ok(c, ps):
                bad.append(rs)
        tested = sum(1 for rs in ALL_SNPS
                     if genos[child].get(rs, (None, None, None))[2]
                     and all(genos[p].get(rs, (None, None, None))[2] for p in pars))
        rows.append({"child": child, "parents": pars, "snps_tested": tested,
                     "inconsistent": bad})
    return rows


def main():
    check_orientation()
    people = [p for p in PEOPLE if has_genotypes(p)]
    genos, results = {}, {}
    for p in people:
        dose, stats = dosages(p)
        genos[p] = dose
        eye, eye_used, eye_missing = predict(EYE_CONST, EYE_BETA, EYE_CATS, "brown", dose)
        hair, hair_used, hair_missing = predict(HAIR_CONST, HAIR_BETA, HAIR_CATS, "blond", dose)
        shade, shade_used, shade_missing = predict(SHADE_CONST, SHADE_BETA, ["dark"], "light", dose)
        herc2 = dose.get("rs12913832", (None, None, None))
        if herc2[0] is None:
            eye = {"note": "no eye prediction: HERC2 rs12913832 missing (model rule)"}
            phrase_eye = "not predicted"
        else:
            phrase_eye = max(("blue", "intermediate", "brown"), key=lambda c: eye[c])
        mc1r_all = set(MC1R_STRONG) | set(MC1R_WEAK)
        if all(dose.get(rs, (None,))[0] is None for rs in mc1r_all):
            hair = {"note": "no hair prediction: all MC1R markers missing (model rule)"}
            phrase_hair = "not predicted"
        else:
            phrase_hair = hair_phrase(hair, shade)
        strong, weak, mc_detail = mc1r_counts(dose)
        results[p] = {
            "chip_snps": stats.get("total") if stats else None,
            "eye": eye, "eye_call": phrase_eye,
            "eye_snps_used": eye_used, "eye_snps_missing": eye_missing,
            "hair": hair, "hair_shade": shade, "hair_call": phrase_hair,
            "hair_snps_used": hair_used, "hair_snps_missing": hair_missing,
            "shade_snps_missing": shade_missing,
            "skin": {"note": "no HIrisPlex-S skin probabilities: the 5-category skin "
                             "model intercepts are not published (see module docstring)",
                     "genotypes": {rs: dose[rs][2] for rs in SKIN_EXTRA
                                   if rs in dose and dose[rs][2]}},
            "herc2_rs12913832": {"genotype_plus_strand": herc2[2],
                                 "brown_allele_A_dosage": herc2[0],
                                 "source": herc2[1]},
            "mc1r": {"strong_R_allele_count": strong, "weak_r_allele_count": weak,
                     "variants": mc_detail},
            "genotype_sources": {rs: dose[rs][1] for rs in sorted(dose)},
            "imputed_or_phased": sorted(rs for rs in dose if dose[rs][1] != "chip"),
            "confidence_note": "model SNPs not measured on the chip were filled from the "
                               "23andMe phased/imputed files (see imputed_or_phased); those "
                               "dosages are statistical, not measured"
                               if any(dose[rs][1] != "chip" and rs in MODEL_SNPS for rs in dose)
                               else "all model SNPs measured on the chip",
        }

    out = {
        "model": {
            "eye": "IrisPlex (Walsh et al. 2011), coefficients from the Walsh et al. 2013 "
                   "FSIG supplementary tool mmc7.xls",
            "hair_and_shade": "HIrisPlex (Walsh et al. 2013 FSIG 7:98-115), same tool",
            "hair_wording": "HIrisPlex prediction guide (Walsh et al. 2013), as reproduced "
                            "in the HIrisPlex-S manual v2.0 and implemented in aHISplex "
                            "classifHISplex",
            "skin": "not run (model parameters unpublished)",
            "note": "The current HIrisPlex-S webtool uses re-fitted models on a larger "
                    "database, so webtool probabilities will differ somewhat.",
        },
        "orientation": {rs: {"model_allele": v[0], "strand_vs_GRCh37_plus": v[1],
                             "plus_strand_effect_allele": v[2]}
                        for rs, v in sorted(EFFECT.items())},
        "people": results,
        "family_consistency": family_check(genos),
    }
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "pigmentation.json"), "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=False)
        fh.write("\n")

    cols = ["person", "eye_call", "p_blue", "p_intermediate", "p_brown", "hair_call",
            "p_blond", "p_brown_hair", "p_red", "p_black", "p_light", "p_dark",
            "herc2_rs12913832", "mc1r_R", "mc1r_r", "eye_snps_missing", "hair_snps_missing"]
    rows = []
    for p in people:
        r = results[p]
        e, h, s = r["eye"], r["hair"], r["hair_shade"]
        rows.append([p, r["eye_call"],
                     e.get("blue", ""), e.get("intermediate", ""), e.get("brown", ""),
                     r["hair_call"], h.get("blond", ""), h.get("brown", ""),
                     h.get("red", ""), h.get("black", ""), s.get("light", ""), s.get("dark", ""),
                     r["herc2_rs12913832"]["genotype_plus_strand"] or "",
                     r["mc1r"]["strong_R_allele_count"], r["mc1r"]["weak_r_allele_count"],
                     ",".join(r["eye_snps_missing"]) or "-",
                     ",".join(r["hair_snps_missing"]) or "-"])
    with open(os.path.join(RESULTS, "pigmentation.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    for row in rows:
        p, eye_call = row[0], row[1]
        print(f"{p:17s} eye {eye_call:12s} blue={row[2]} int={row[3]} brown={row[4]}")
        print(f"{'':17s} hair {row[5]:20s} blond={row[6]} brown={row[7]} red={row[8]} "
              f"black={row[9]}  light={row[10]} dark={row[11]}")
        print(f"{'':17s} HERC2 rs12913832={row[12]}  MC1R R={row[13]} r={row[14]}")
    for fc in out["family_consistency"]:
        print(f"{fc['child']} vs {'+'.join(fc['parents'])}: {fc['snps_tested']} SNPs, "
              f"{len(fc['inconsistent'])} inconsistent "
              f"{fc['inconsistent'] if fc['inconsistent'] else ''}")


if __name__ == "__main__":
    main()
