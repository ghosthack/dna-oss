"""One-page, printable DNA summaries for clinicians: results/doctor/<person>.html + index.html.

Usage: python3 -I scripts/doctor.py [--date YYYY-MM-DD]   (run after pgx, apoe, wellness and pgs)

Reads existing results only (no genetics is recomputed here). People without any results
(e.g. an empty raw/ folder) are skipped automatically. Pages are self-contained HTML files
that stay on this machine.
"""
import csv
import glob
import html
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dnalib
from dnalib import RAW, RESULTS, children, person as role_person, people, role, siblings, text
from pgx import DRUGS, note

esc = html.escape
OUT = os.path.join(RESULTS, "doctor")
HLA_DRUG = {"HLA-B*57:01": "abacavir", "HLA-A*31:01": "carbamazepine", "HLA-B*58:01": "allopurinol"}
# CPIC wording per allele (B*58:01/allopurinol "contraindicated" is handled in its own branch below).
HLA_ACTION = {"HLA-B*57:01": "abacavir is not recommended (hypersensitivity risk)",
              "HLA-A*31:01": ("if carbamazepine-naive, do not use carbamazepine when alternatives are available; "
                              "if none, consider it with increased clinical monitoring (SJS/TEN/DRESS risk)")}
WARFARIN = ("VKORC1", "CYP2C9", "CYP4F2")


# ---------------------------------------------------------------- inputs
def read_tsv(name):
    path = os.path.join(RESULTS, name)
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_pgs(person):
    path = os.path.join(RESULTS, f"pgs_{person}.json")
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        return json.load(fh)


def load_results_json(name):
    path = os.path.join(RESULTS, name)
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        return json.load(fh)


def sex(person):
    """"M", "F" or None: family.json sex, else the sex inferred from the chip (summary.json qc), else unknown."""
    return dnalib.sex(person) or load_results_json("summary.json").get("qc", {}).get(person, {}).get("sex") or None


def mostly_european(person):
    """True when results/ancestry.json puts this person at >=90% broad European ancestry."""
    broad = load_results_json("ancestry.json").get("people", {}).get(person, {}).get("broad", {})
    return (broad.get("European") or 0) >= 0.9


def label(person):
    return dnalib.label(person) if role(person) else person.replace("-", " ").capitalize()


def patient(person):
    """Page title / index name for a clinician: the family.json name, with the relation in parentheses
    ("Maya (Alex's sister)"); without a configured name, the usual label."""
    own = dnalib.FAMILY["people"].get(person, {}).get("name")
    r = role(person)
    if not own:
        return label(person)
    if r in (None, "self"):
        return own
    rel = (kin("sibling", sex(person)) if r == "sibling" else r.replace("maternal_", "maternal "))
    viewer = dnalib.FAMILY["people"].get(role_person("self") or "", {}).get("name")
    return f"{own} ({viewer}'s {rel})" if viewer else f"{own} ({rel})"


def pron(person, case="subj"):
    """he/she/they ("subj"), his/her/their ("poss"), him/her/them ("obj") from family.json sex."""
    forms = {"M": ("he", "his", "him"), "F": ("she", "her", "her")}.get(sex(person), ("they", "their", "them"))
    return forms[("subj", "poss", "obj").index(case)]


def kin(word, s):
    """'child' -> son/daughter/child, 'sibling' -> brother/sister/sibling, by sex."""
    return {"child": {"M": "son", "F": "daughter"}, "sibling": {"M": "brother", "F": "sister"}}[word].get(s, word)


def people_with_results(pgx, apoe, wellness):
    found = {r["person"] for rows in (pgx, apoe, wellness) for r in rows}
    found |= {os.path.basename(p)[4:-5] for p in glob.glob(os.path.join(RESULTS, "pgs_*.json"))}
    folders = sorted(d for d in os.listdir(RAW) if os.path.isdir(os.path.join(RAW, d)))
    order = people() + [d for d in folders if d not in people()]
    return [p for p in order if p in found]


def data_sources(person):
    has = lambda pat: bool(glob.glob(os.path.join(RAW, person, pat)))
    chips = []
    for f in glob.glob(os.path.join(RAW, person, "genome_*")):
        v = re.search(r"_(v\d)_", os.path.basename(f))
        chips.append(f"23andMe {v.group(1) if v else ''} chip".replace("  ", " "))
    if has("dna-data-*") or has("AncestryDNA*"):
        chips.append("AncestryDNA v2 chip" + (" (merged with 23andMe; discordant calls dropped)" if chips else ""))
    return chips, has("imputed_genotype_data_*")


# ---------------------------------------------------------------- sections
def short_pheno(ph):
    return ph.split(" (")[0].split(" — ")[0]


def show_pheno(ph):
    """Phenotype for the table: keep a CPIC activity score (used for dosing), drop other asides."""
    base = short_pheno(ph)
    if "(heterozygous female)" in ph:  # G6PD: the qualifier is the result
        return ph
    if "(activity score" in ph:
        score = ph.split("(activity score")[1].split(")")[0].split(";")[0].strip()
        return f"{base} (activity score {score})"
    return base


# Affected-drug wording where pgx.DRUGS overstates the guideline coverage.
DRUG_TEXT = {"ABCG2": "rosuvastatin (urate/gout association; allopurinol response research-only, no CPIC guidance)"}


def drugs(gene):
    return DRUG_TEXT.get(gene, DRUGS[gene])


# Current CPIC CYP2D6 activity values (*10 = *41 = 0.25; *41 was 0.5 in the 2019 consensus).
D6_VALUE = {"*1": 1, "*4": 0, "*6": 0, "*10": 0.25, "*41": 0.25}


def d6_pheno(score):
    return ("poor metabolizer" if score == 0 else "intermediate metabolizer" if score < 1.25
            else "normal metabolizer" if score <= 2.25 else "ultrarapid metabolizer")


def is_poor(ph):
    return ph.startswith("Poor")


def rosu_start(slco, abcg2):
    """CPIC 2022 rosuvastatin starting-dose limit from SLCO1B1 + ABCG2 ('' = desired starting dose)."""
    if is_poor(abcg2):
        return "≤10 mg" if is_poor(slco) or slco.startswith("Decreased") else "≤20 mg"
    return "≤20 mg" if is_poor(slco) else ""


def fluva(slco, c9):
    """CPIC 2022 fluvastatin advice from SLCO1B1 + CYP2C9 phenotypes ('' = standard dosing)."""
    c9_im, c9_pm = c9.startswith("Intermediate"), is_poor(c9)
    if is_poor(slco):
        return "use an alternative statin" if c9_im or c9_pm else "≤40 mg/day starting dose"
    if slco.startswith("Decreased"):
        return ("use an alternative statin" if c9_pm else
                "≤20 mg/day starting dose (alternative statin if >20 mg needed)" if c9_im else "")
    return "≤20 mg/day starting dose" if c9_pm else "≤40 mg/day starting dose" if c9_im else ""


def pgx_action(gene, r, by_gene):
    """CPIC-based action text for a non-routine result ('' = routine). Falls back to pgx.note()."""
    ph, dip = r["phenotype"], r["diplotype"]
    other = lambda g: (by_gene.get(g) or {}).get("phenotype", "")
    if gene == "CYP2C19":
        if ph.startswith("Ultrarapid"):
            return ("Faster clearance of several CYP2C19 drugs. CPIC: clopidogrel — standard dose. PPIs — increase "
                    "starting dose by 100%. Citalopram/escitalopram — consider an alternative not mainly cleared by "
                    "CYP2C19. Sertraline — standard starting dose (switch if inadequate response). Tertiary-amine TCAs "
                    "(amitriptyline, clomipramine…) — consider nortriptyline/desipramine. Voriconazole (adults) — choose "
                    "an alternative (isavuconazole, posaconazole, liposomal amphotericin B).")
        if ph.startswith("Rapid"):
            return ("Clopidogrel: standard dose (CPIC strong). PPIs: standard starting dose; consider +50–100% for "
                    "H. pylori/erosive esophagitis. Voriconazole (adults): choose an alternative (isavuconazole, "
                    "posaconazole, liposomal amphotericin B). SSRIs: standard start; if inadequate response, titrate or "
                    "switch to a non-CYP2C19 agent. Tertiary-amine TCAs: consider nortriptyline/desipramine.")
        if ph.startswith("Intermediate"):
            return ("Reduced clopidogrel activation. CPIC: for ACS/PCI avoid standard-dose clopidogrel if possible — "
                    "use prasugrel or ticagrelor if no contraindication (prasugrel contraindicated after stroke/TIA); "
                    "for neurovascular indications consider ticagrelor. Other listed CYP2C19 drugs: standard starting "
                    "doses for intermediate metabolizers (PPIs: consider 50% reduction for chronic >12-week use).")
    if gene == "CYP2C9" and ph.startswith("Intermediate") and "activity score 1.5" in ph:
        slco = other("SLCO1B1")
        fl = fluva(slco, ph)
        statin = ((f" Fluvastatin (with SLCO1B1 {short_pheno(slco).lower()})" if slco else " Fluvastatin")
                  + f" — {fl}." if fl else "")
        siponimod = (f" Siponimod (FDA label; no CPIC guideline) — standard 2 mg maintenance dose ({dip})."
                     if "*3" not in dip else "")
        return ("Activity score 1.5. CPIC: NSAIDs (celecoxib, ibuprofen…) and phenytoin — standard starting doses "
                "(adjust phenytoin by TDM). Warfarin — include in genotype-guided dosing algorithm (small reduction)."
                + statin + siponimod)
    if gene == "SLCO1B1":
        c9_reduced = any(other("CYP2C9").startswith(x) for x in ("Intermediate", "Poor"))
        abcg2 = other("ABCG2")
        rosu_lim = rosu_start(ph, abcg2)
        if is_poor(ph):
            rosu = (" (ABCG2 normal)" if abcg2.startswith("Normal") else
                    " (ABCG2 decreased function does not change this)" if abcg2.startswith("Decreased") else
                    " (ABCG2 poor function also present)" if is_poor(abcg2) else "")
            return ("Higher statin exposure and myopathy risk. CPIC (SLCO1B1 poor function): avoid simvastatin and "
                    "lovastatin" + (" (and fluvastatin, given CYP2C9 intermediate/poor)" if c9_reduced else "") +
                    f". Starting-dose limits: atorvastatin ≤20 mg, rosuvastatin {rosu_lim}{rosu}, pravastatin ≤40 mg, "
                    "pitavastatin ≤1 mg" + ("" if c9_reduced else ", fluvastatin ≤40 mg") +
                    ". If higher intensity is needed, consider rosuvastatin or add non-statin therapy (e.g. ezetimibe) "
                    "rather than withholding treatment. Based on rs4149056 only — confirm with a clinical PGx test.")
        if ph.startswith("Decreased"):
            fl = fluva(ph, other("CYP2C9"))
            rosu = (f"rosuvastatin {rosu_lim} starting dose (ABCG2 poor function also present)" if rosu_lim else
                    "rosuvastatin at the usual starting dose, with extra myopathy caution above 20 mg")
            return ("Increased simvastatin/lovastatin exposure and myopathy risk. CPIC: use an alternative statin, or "
                    "simvastatin <20 mg/day; lovastatin alternative or ≤20 mg/day. Atorvastatin ≤40 mg and "
                    f"pitavastatin ≤2 mg starting doses; {rosu}"
                    + (f"; with CYP2C9 {short_pheno(other('CYP2C9')).lower()}, fluvastatin {fl}" if fl else "")
                    + ". Based on rs4149056 only.")
    if gene == "ABCG2":
        slco = other("SLCO1B1")
        if ph.startswith("Decreased"):
            limit = (" (rosuvastatin limit for this patient comes from SLCO1B1, above)" if is_poor(slco) else "")
            return (f"Modestly higher rosuvastatin exposure; CPIC: no dose change for ABCG2 decreased function "
                    f"alone{limit} (CPIC lowers the rosuvastatin start to ≤20 mg for ABCG2 poor function, or ≤10 mg "
                    f"when that is combined with reduced SLCO1B1 function). Q141K is an established risk allele "
                    f"for raised urate/gout; no CPIC guideline for allopurinol.")
        if is_poor(ph):
            why = " (SLCO1B1 function also reduced)" if rosu_start(slco, ph) == "≤10 mg" else ""
            return (f"Higher rosuvastatin exposure. CPIC (optional): rosuvastatin {rosu_start(slco, ph)} starting "
                    f"dose or an alternative statin{why}. Q141K is an established risk allele for raised urate/gout.")
    if gene == "TPMT" and ph.startswith("Intermediate"):
        text = ("TPMT intermediate: CPIC — start azathioprine/mercaptopurine at 30–80% and thioguanine at 50–80% of "
                "normal dose (thioguanine: if the standard start is ≥40 mg/m²/day) and adjust to myelosuppression.")
        n = other("NUDT15")
        if n.startswith("Intermediate"):
            text = ("TPMT + NUDT15 intermediate (compound intermediate): CPIC — start azathioprine/mercaptopurine/"
                    "thioguanine at 20–50% of the standard starting dose (instead of the single-gene 30–80%) and "
                    "adjust to myelosuppression.")
        elif is_poor(n):
            text += (" Also NUDT15 poor metabolizer: CPIC — drastically reduced doses for malignancy, an alternative "
                     "agent for nonmalignant conditions.")
        elif n:
            text += " NUDT15 normal (*3 only tested)."
        return text + " Confirm with clinical TPMT/NUDT15 testing before starting."
    if gene == "NUDT15" and ph.startswith("Intermediate"):
        how = ("imputed call (variant not on this chip), so less certain than a direct genotype"
               if "rs116855232" in r.get("imputed_snps", "") else
               "array calls of rare variants are often false positives")
        return (f"Rare variant in Europeans (NUDT15*3 ~0.2–0.3%); {how} — confirm with clinical TPMT/NUDT15 "
                "genotyping. If confirmed, CPIC: start "
                "azathioprine/mercaptopurine at 30–80% and thioguanine at 50–80% of normal dose (thioguanine: if the "
                "standard start is ≥40 mg/m²/day)" + ("; 20–50% with TPMT also intermediate (see TPMT row)"
                                       if other("TPMT").startswith("Intermediate") else "")
                + ", adjust by myelosuppression.")
    if gene == "CYP3A5" and (ph.startswith("Normal") or ph.startswith("Intermediate")):
        return ("CYP3A5 expressor: CPIC — tacrolimus starting dose 1.5–2× standard (max 0.3 mg/kg/day), "
                "then therapeutic drug monitoring.")
    if gene == "G6PD" and ph.startswith("Variable"):
        return ("Heterozygous female: with random X-inactivation, red-cell G6PD activity can range from normal to "
                "deficient, and a screening (qualitative) test can miss intermediate activity. Use quantitative G6PD "
                "activity testing before rasburicase, primaquine/tafenoquine or dapsone; tafenoquine needs ≥70% of "
                "normal activity.")
    if gene == "HLA-B*58:01" and "Possible" in ph:
        why = ("Tag SNP (rs9263726) is validated in East Asians only; in Europeans the tag allele frequency is "
               "~15–23% (1000 Genomes EUR/TSI/IBS) vs ~1–2% for HLA-B*58:01, so this is most likely a false "
               "positive." if mostly_european(r.get("person")) else
               "Tag SNP only (rs9263726, validated in East Asians; its accuracy varies by ancestry) — confirm by "
               "HLA typing.")
        return (why + " If allopurinol is being considered, order HLA-B*58:01 typing; CPIC: allopurinol "
                "contraindicated only if typing confirms the allele.")
    if gene.startswith("HLA") and "Possible" in ph:
        action = HLA_ACTION.get(gene, f"{HLA_DRUG.get(gene, DRUGS[gene])}: see CPIC guideline")
        return f"Tag SNP only — confirm by HLA typing. If confirmed, CPIC: {action}."
    return note(gene, ph)


# What "nothing found" means for genes where the chip tests only a few variants.
LIMITED = {
    "DPYD": ("DPYD: none of the 4 CPIC core variants (*2A, *13, c.2846A>T, HapB3) detected; other variants not "
             "tested — does not exclude DPD deficiency; follow local DPD testing policy before fluoropyrimidines."),
    "TPMT": "TPMT: *2, *3A/B/C not detected",
    "NUDT15": "NUDT15: *3 not detected",
    "G6PD": ("G6PD: none of 3 tested variants (A- 202A/376G, Mediterranean) detected; other deficiency alleles not "
             "tested — does not exclude deficiency; use G6PD enzyme activity testing before rasburicase, "
             "primaquine/tafenoquine or dapsone."),
}


def limited_text(gene, r):
    """LIMITED note, adjusted when a core variant wasn't on this person's chip (DPYD HapB3 on the v4 chip)."""
    if gene == "DPYD" and "rs75017182" in r.get("missing_snps", ""):
        return ("DPYD: *2A, *13 and c.2846A>T not detected; HapB3 (the most common CPIC variant in Europeans) is not on "
                "this chip, so DPD status is less complete than usual — follow local DPD testing policy before "
                "fluoropyrimidines.")
    return LIMITED[gene]


def routine_label(gene, r):
    dip, ph = r["diplotype"], r["phenotype"]
    if gene == "CYP3A5" and ph.startswith("Poor metabolizer (non-expressor"):
        return ("CYP3A5 *3/*3 (poor metabolizer/non-expressor, the usual European genotype; CPIC: standard "
                "tacrolimus starting dose)") if dip == "*3/*3" else f"CYP3A5 {dip} (non-expressor; standard tacrolimus dosing)"
    if gene == "UGT1A1" and dip == "*1/*1":
        return "UGT1A1: *28 not detected (inferred from tag rs887829; TA repeat not directly assayed), *6 not detected"
    if gene == "UGT1A1":
        return f"UGT1A1 {dip} ({short_pheno(ph).lower()}; *28 inferred from tag rs887829)"
    if gene == "CYP2B6" and dip == "516G/G":
        return "CYP2B6: 516G>T not detected (other alleles, e.g. *4, *18, not tested)"
    if gene == "CYP2C9" and dip == "*1/*1":
        return "CYP2C9 *1/*1 (none of *2, *3, *5, *8, *11 detected)"
    if gene == "CYP2C19" and dip == "*1/*1":
        return "CYP2C19 *1/*1 (*2, *3, *17 not detected)"
    if gene in ("SLCO1B1", "ABCG2") and ph.startswith("Normal"):
        return f"{gene} {dip} (normal function)"
    if gene == "VKORC1":
        return f"VKORC1 {dip} (normal warfarin sensitivity)"
    if gene == "CYP4F2":
        return f"CYP4F2 {dip}" + (" (minor warfarin dose factor)" if "higher" in ph else "")
    return f"{gene} {dip} ({short_pheno(ph).lower()})"


def cyp2d6_line(r):
    dip, miss = r["diplotype"], r["missing_snps"]
    stars = dip.split("/")
    if len(stars) != 2 or any(s not in D6_VALUE for s in stars):
        return None, None
    score = sum(D6_VALUE[s] for s in stars)
    alt = ""
    if "rs1065852" in miss and "*1" in stars:
        alt_score = score - 1 + D6_VALUE["*10"]
        alt = f"; rs1065852 (*10) not called (if present, AS {alt_score:g} = {d6_pheno(alt_score)})"
    elif miss:
        alt = f"; {miss} not called"
    return score, (f"<b>CYP2D6</b> {esc(dip)}: activity score {score:g} = {d6_pheno(score)} if the "
                   f"{'*1' if '*1' in stars else 'diplotype'} call is correct (current CPIC activity values, *41 = 0.25){esc(alt)}; "
                   "deletions/duplications not detectable — get a clinical CYP2D6 test before codeine, tramadol, "
                   "tamoxifen, ondansetron or CYP2D6-dependent antidepressants.")


def pgx_section(rows):
    by_gene = {r["gene"]: r for r in rows}
    flagged, routine, limited, extra, untested = [], [], [], [], []
    warfarin_flag = "VKORC1" in by_gene and bool(note("VKORC1", by_gene["VKORC1"]["phenotype"]))
    for gene, r in by_gene.items():
        dip, ph = r["diplotype"], r["phenotype"]
        if gene == "CYP2D6":
            if "not covered" in dip or "Not determinable" in ph:
                extra.append("<b>CYP2D6: not tested</b> — key SNPs (*4, *10, *41) absent from this chip; "
                             "no inference about metabolizer status is possible.")
                continue
            score, line = cyp2d6_line(r)
            if score is not None and d6_pheno(score) == "normal metabolizer":
                extra.append(line)
            else:
                action = note(gene, ph) or "Non-normal CYP2D6 result; confirm with a clinical CYP2D6 test."
                shown = (f"{d6_pheno(score).capitalize()} (activity score {score:g})" if score is not None
                         else show_pheno(ph))
                flagged.append((gene, dip, shown, action + " Deletions/duplications not detectable."))
            continue
        if "Not determinable" in ph:  # none of the gene's SNPs were genotyped: say so, never "typical"
            untested.append(gene)
            continue
        action = pgx_action(gene, r, by_gene)
        if gene.startswith("HLA"):
            if action:  # the phenotype's "— typing needed…" aside is covered by the action column
                flagged.append((gene, dip, short_pheno(ph), action))
            else:
                routine.append(f"{gene} tag negative")
            continue
        if gene == "CYP4F2" and not action and "higher" in ph and warfarin_flag:
            continue  # folded into the warfarin row below
        if not action:
            if gene in LIMITED:
                limited.append(gene)
            else:
                caveat = ph.rsplit(" (", 1)[1].rstrip(")") if ph.endswith("not tested)") else ""
                routine.append(f"{gene} {dip} ({short_pheno(ph).lower()}; {caveat})" if caveat
                               else routine_label(gene, r))
            continue
        if gene == "VKORC1" and "CYP4F2" in by_gene:
            c9, f2 = by_gene.get("CYP2C9"), by_gene["CYP4F2"]
            if c9:
                c9_txt = (f"CYP2C9 {c9['diplotype']}" + (" (= none of *2, *3, *5, *8, *11 detected; other alleles not "
                                                        "tested)" if c9["diplotype"] == "*1/*1" else ""))
                action += (f" Enter with {c9_txt} and CYP4F2 {f2['diplotype']} into the calculator. Genotype guides "
                           "the initial dose only; titrate by INR.")
        flagged.append((gene, dip, show_pheno(ph), action))

    body = "".join(
        f"<tr><td><b>{esc(g)}</b></td><td>{esc(d)}</td><td>{esc(p)}</td>"
        f"<td>{esc(drugs(g))}</td><td>{esc(a)}</td></tr>" for g, d, p, a in flagged)
    table = (f'<table><thead><tr><th>Gene</th><th>Diplotype</th><th>Phenotype</th><th>Affected drugs</th>'
             f'<th>CPIC-based action</th></tr></thead><tbody>{body}</tbody></table>' if flagged
             else "<p>No non-routine pharmacogenomic result among the genes tested.</p>")
    hla_note = (" HLA results are tag SNP only — confirm by HLA typing before relying on a negative. "
                "HLA-B*15:02 (carbamazepine/phenytoin) not tested."
                if any(x.startswith("HLA") for x in routine) else "")
    lines = ([f"<p class=small><b>Typical / not detected:</b> {esc('; '.join(routine))}.{esc(hla_note)}</p>"]
             if routine else [])
    if limited:
        parts = [limited_text(g, by_gene[g]) for g in ("DPYD", "G6PD") if g in limited]
        thio = [limited_text(g, by_gene[g]) for g in ("TPMT", "NUDT15") if g in limited]
        if thio:
            parts.append("; ".join(thio) + " (other alleles not tested) — does not replace clinical TPMT/NUDT15 "
                         "testing before thiopurines.")
        lines.append("<p class=small><b>Limited panels (only a few variants tested):</b> "
                     + " ".join(esc(p) for p in parts) + "</p>")
    if untested:
        lines.append(f"<p class=small><b>Not tested on this chip</b> (no result, not a normal result): "
                     f"{esc(', '.join(untested))}.</p>")
    lines += [f"<p class=small>{x}</p>" for x in extra]
    return table + "".join(lines), flagged, by_gene, limited


def split_genotypes(geno):
    """'rs6025=CC rs1799963=GG*' -> [('rs6025', 'CC', False), ('rs1799963', 'GG', True)]"""
    out = []
    for tok in geno.split():
        rs, _, g = tok.partition("=")
        out.append((rs, g.rstrip("*"), g.endswith("*")))
    return out


def show_genotypes(geno):
    return " ".join(f"{rs}={g}" + ("†" if imp else "") for rs, g, imp in split_genotypes(geno))


SNP_NAME = {"rs6025": "Factor V Leiden (F5 rs6025)", "rs1799963": "prothrombin G20210A (F2 rs1799963)"}


def conditions_section(well, apoe_row, person=None):
    W = {r["report"]: r for r in well}
    rows = []

    def add(name, key, interp):
        r = W.get(key)
        if not r:
            return
        rows.append((name, r["gene"], show_genotypes(r["genotypes"]), interp(r)))

    def hfe(r):
        res = r["result"]
        if "No HFE" in res:
            return "Negative for C282Y and H63D."
        if "H63D carrier" in res:
            return "H63D heterozygote (carrier). No iron-overload risk from this alone."
        return res

    def clot(r):
        res = r["result"]
        if res.startswith("Not determinable"):
            return "Not determinable (Factor V Leiden and prothrombin G20210A both no-calls) — clinical test if indicated."
        calls = split_genotypes(r["genotypes"])
        untested = [rs for rs, g, _ in calls if not g.strip("-")]
        name = {"rs6025": "Factor V Leiden", "rs1799963": "Prothrombin G20210A"}
        gap = (" " + " and ".join(name.get(rs, rs) for rs in untested) + " not tested (chip no-call).") if untested else ""
        if not res.startswith("No "):
            found = res.split(":")[0]
            when = {"M": "surgery, immobilisation and long-haul travel",
                    "F": ("estrogen-containing contraception/HRT, pregnancy and the postpartum period, surgery, "
                          "immobilisation and long-haul travel")}.get(
                sex(person), "surgery, immobilisation, long-haul travel and, in women, estrogen-containing "
                             "contraception/HRT, pregnancy and the postpartum period")
            return f"{found}: raised venous thromboembolism risk — relevant to {when}.{gap}"
        tested = [(rs, imp) for rs, g, imp in calls if rs not in untested]
        if not any(imp for _, imp in tested) and not untested:
            text = "Factor V Leiden and prothrombin G20210A not detected (direct chip calls)."
        else:
            text = " ".join(
                f"{name.get(rs, rs)} "
                + ("not detected (imputed, chip no-call — not a tested negative)." if imp else
                   "not detected (direct chip call).") for rs, imp in tested) + gap
            text += " Confirm by clinical test if thrombophilia evaluation is indicated (e.g. personal/family VTE history)."
        return text + " Protein C, protein S and antithrombin deficiency not assessed."

    def a1at(r):
        res = r["result"]
        if "PiMM" in res:
            return ("No S or Z allele detected; other deficiency/null alleles (e.g. Mmalton) not tested — measure "
                    "serum AAT level and phenotype if clinically indicated (COPD, bronchiectasis, unexplained liver disease).")
        if "PiMS" in res:
            return ("One S allele, no Z allele (other alleles not tested). One S allele: AAT levels normal or mildly "
                    "reduced; no disease from this alone. Avoid smoking.")
        return res

    def celiac(r):
        res = r["result"]
        if res.startswith("No DQ2/DQ8"):
            return ("DQ2.5, DQ2.2 and DQ8 not detected by 4 tag SNPs (rs2187668, rs7454108, rs2395182/rs7775228; "
                    "not HLA typing; DQ7.5/half-DQ2 not assessed). "
                    "If confirmed by clinical HLA-DQ typing, celiac disease would be very unlikely. Do not use this "
                    "result in place of tTG-IgA serology (on a gluten-containing diet) when symptoms suggest celiac disease.")
        if "Carries" in res:
            dq = res.split(":")[0].replace("Carries ", "").replace("×1", " (one copy)").replace("×2", " (two copies)")
            return (f"Tag-SNP inference of {dq}; permissive genotype in ~30–40% of people and does not indicate "
                    "celiac disease. Not equivalent to clinical HLA-DQ typing — if HLA status is needed (equivocal "
                    "serology or already gluten-free), order clinical typing.")
        if res.startswith("Possible DQ2.2"):
            dq4 = " or DQ4 (rs4713586)" if "(or DQ4)" in res else ""
            return (f"DQ2.5 and DQ8 not detected; unphased tags (rs2395182/rs7775228) compatible with DQ2.2{dq4}, "
                    "not resolvable by tag SNPs. Not equivalent to clinical HLA-DQ typing — if HLA status is "
                    "needed (equivocal serology or already gluten-free), order clinical typing.")
        return res

    add("Hereditary hemochromatosis", "Hemochromatosis (iron)", hfe)
    add("Thrombophilia", "Inherited clotting", clot)
    add("Alpha-1 antitrypsin deficiency", "Alpha-1 antitrypsin", a1at)
    def lpa(r):
        if r["result"].startswith("No LPA"):
            return ("rs10455872 G and rs3798220 C not detected. These two variants identify only a minority of people "
                    "with high Lp(a) (most variation is KIV-2 copy number), so this does not rule it out.")
        return ("LPA risk allele(s) present (rs10455872 G / rs3798220 C; Clarke 2009: markedly higher Lp(a), "
                "~1.5–1.7× CHD risk per copy). Measure Lp(a) to confirm.")

    def b27(r):
        res = r["result"]
        if res.startswith("HLA-B27 not detected"):
            return ("rs4349859 A not detected (HLA-B*27 tag; 98% sensitive, 99% specific in Europeans, Evans 2011); "
                    "may miss non-European subtypes (B*27:03/04/06/07). Not HLA typing.")
        if res.startswith("HLA-B27 likely"):
            return ("rs4349859 A present (HLA-B*27 tag, confirmed by rs116488202). Present in ~8% of Europeans; "
                    "relevant only with inflammatory back pain or anterior uveitis. Not HLA typing.")
        return res

    def carrier(r):
        res = r["result"]
        if res.startswith("Not a carrier"):
            return res.replace("(only the variants listed; others not tested)",
                               "(listed alleles only; negative result does not exclude other variants)") + "."
        return res + "."

    add("Celiac disease susceptibility (HLA-DQ2/DQ8)", "Celiac genes (HLA-DQ2/DQ8)", celiac)
    add("Lipoprotein(a)", "Lipoprotein(a)", lpa)
    add("HLA-B27 (axial spondyloarthritis)", "HLA-B27 (spine arthritis)", b27)
    add("β-thalassemia / sickle cell / HbC carrier", "Beta-thalassemia & sickle cell", carrier)
    add("Cystic fibrosis carrier", "Cystic fibrosis", carrier)
    add("GJB2 hearing loss carrier", "Inherited hearing loss", carrier)
    add("Familial Mediterranean Fever (MEFV)", "Familial Mediterranean Fever", carrier)
    for name, key in (("Phenylketonuria carrier (PAH)", "Phenylketonuria (PKU)"), ("Gaucher carrier (GBA N370S)", "Gaucher disease"),
                      ("Hereditary fructose intolerance carrier (ALDOB)", "Hereditary fructose intolerance"),
                      ("MCAD deficiency carrier (ACADM)", "MCAD deficiency (fat metabolism)"), ("Pompe disease carrier (GAA)", "Pompe disease"),
                      ("PMM2-CDG carrier", "PMM2-CDG (metabolic disorder)"), ("Pendred syndrome carrier (SLC26A4)", "Pendred syndrome (hearing)"),
                      ("Pyruvate kinase deficiency carrier (PKLR)", "Pyruvate kinase deficiency")):
        add(name, key, carrier)

    def passthrough(note, always=False):
        return lambda r: (r["result"].rstrip(".") + ". " + (note if always or r["result"].startswith(
            ("Variants not", "G84E not")) else "")).rstrip()

    add("Parkinson disease risk variants", "Parkinson's genes",
        passthrough("Tested: LRRK2 G2019S, GBA N370S only (23andMe's panel)."))
    add("Hereditary transthyretin amyloidosis", "Hereditary amyloidosis",
        passthrough("Tested: TTR V30M (p.Val50Met), V122I (p.Val142Ile), T60A (p.Thr80Ala) only."))
    add("APOL1 kidney risk", "Kidney disease (APOL1)",
        passthrough("G1 called from rs73885319 alone (rs60910145 not on the chip); risk applies to two risk alleles.",
                    always=True))
    hoxb13_note = ("Single variant (population OR ~3–5, Karlsson 2014; ~20× in familial/early-onset cases, Ewing "
                   "2012); relevant to male relatives.")
    if sex(person) == "F":  # her own result carries no prostate risk; it matters for the men she shares it with
        def hoxb13(r):
            if r["result"].startswith("G84E not"):
                return passthrough(hoxb13_note)(r)
            return (r["result"].split(":")[0] + ": no prostate cancer risk for her; relevant to male relatives "
                    "(father, brothers, sons), who may carry it and could be offered testing and an earlier PSA "
                    "discussion (population OR ~3–5, Karlsson 2014; ~20× in familial/early-onset cases, Ewing 2012). "
                    "Confirm with a clinical test.")
        add("HOXB13 G84E (hereditary prostate cancer)", "Hereditary prostate cancer", hoxb13)
    else:
        add("HOXB13 G84E (hereditary prostate cancer)", "Hereditary prostate cancer", passthrough(hoxb13_note))
    # Two variants in a recessive gene is a possible affected genotype (phase unknown), not carrier status.
    rows = [(name.replace(" carrier", "") if " — two variants" in interp else name, gene, geno, interp)
            for name, gene, geno, interp in rows]
    body = "".join(f"<tr><td>{esc(a)}</td><td>{esc(b)}</td><td class=mono>{esc(c)}</td><td>{esc(d)}</td></tr>"
                   for a, b, c, d in rows)
    out = (f'<table><thead><tr><th>Condition</th><th>Gene</th><th>Genotype</th><th>Result</th></tr></thead>'
           f'<tbody>{body}</tbody></table><p class=small>Heterozygous carriers are usually asymptomatic (exceptions: '
           f'β-thalassemia trait causes microcytosis that can mimic iron deficiency; sickle trait carries rare '
           f'exertional/dehydration risks; MEFV heterozygotes, especially M694V, occasionally have mild FMF-like '
           f'symptoms). '
           f'Negative results cover the common variants listed only; rarer pathogenic variants are not tested. '
           f'† = that SNP is from 23andMe phased/imputed files (chip no-call), not a direct chip call.</p>')
    out += apoe_line(apoe_row)
    return out, rows


def apoe_line(apoe_row):
    if not apoe_row:
        return ""
    g, src = apoe_row["apoe"], apoe_row["source"]
    if not g:  # apoe.py writes an empty call when neither the chip nor imputation resolves it
        return ("<p><b>APOE (risk modifier, not an inherited condition):</b> not determined "
                f"({esc(src or 'no call')}). Clinical testing is needed for an APOE result.</p>")
    how = "imputed, chip no-call" if "imputed" in src else "direct consumer-chip call"
    n4 = g.count("ε4")
    if n4:
        odds = "~2.5×" if g == "ε2/ε4" else "~3×" if n4 == 1 else "~12–15×"
        interp = (("One ε4 allele: higher average risk of late-onset Alzheimer disease" if n4 == 1 else
                   "Two ε4 alleles: substantially higher average risk of late-onset Alzheimer disease")
                  + f" (odds {odds} vs ε3/ε3; Farrer 1997 clinic/autopsy meta-analysis, population-based estimates "
                  "lower)" + (" and modestly higher LDL." if g == "ε3/ε4" else ".") + " Not diagnostic and not recommended for predicting dementia (ACMG/NSGC). If "
                  "anti-amyloid therapy (lecanemab/donanemab) is ever considered, confirm APOE by clinical test: ε4 "
                  "carriers have a higher ARIA risk, and labels recommend testing and counseling before treatment.")
    elif g == "ε2/ε2":
        interp = "Relevant to type III hyperlipoproteinemia only if lipids are abnormal; not used for Alzheimer prediction."
    elif "ε2" in g:
        interp = ("ε2 carrier: slightly lower average Alzheimer risk (odds ~0.6× vs ε3/ε3) and LDL. Not recommended for Alzheimer risk "
                  "prediction (ACMG/NSGC); no clinical action. Confirm clinically if relevant to lipid evaluation.")
    else:
        interp = "Most common genotype, average risk; not recommended for Alzheimer risk prediction (ACMG/NSGC)."
    return (f"<p><b>APOE (risk modifier, not an inherited condition):</b> {esc(g)}"
            f"{'†' if 'imputed' in src else ''} ({esc(how)}). {esc(interp)}</p>")


PGS_ADVICE = {
    "glaucoma": "Periodic eye exam with IOP and optic-nerve assessment.",
    "urate": ("Polygenic tendency to higher urate (imprecise percentile, small reference). Check serum urate if gout "
              "symptoms; no guideline supports testing or treating asymptomatic hyperuricemia."),
    "breast cancer": "Discuss formal risk assessment (e.g. CanRisk) alongside family history to guide screening.",
    "prostate cancer": "Shared decision-making about PSA screening; consider alongside family history.",
    "coronary artery disease": "Attention to standard CV risk factors; lipid panel and overall risk estimate.",
    "atrial fibrillation": "Pulse/rhythm check at visits; evaluate palpitations promptly.",
    "type 2 diabetes": "Periodic HbA1c or fasting glucose; weight and activity counseling.",
    "ldl cholesterol": "Lipid panel per routine schedule.",
    "systolic blood pressure": "Regular blood-pressure measurement.",
    "colorectal cancer": ("Modest polygenic tendency; family history of colorectal cancer or advanced adenomas outweighs it "
                          "and sets screening age and interval. Not a reason to change screening on its own."),
    "benign neoplasm of colon": ("Modest polygenic tendency to colon polyps (adenomas). Family history of adenomas, "
                                 "their histology and age at removal should drive colonoscopy timing; the score does not."),
    "melanoma": ("Sun protection counseling. This score already includes the main MC1R red-hair variants, so don't "
                 "count an MC1R result as a separate risk factor on top of it. Consider a baseline full-skin exam with "
                 "fair/freckled skin, easy burning, many or atypical moles, or a family history of melanoma "
                 "(see skin section)."),
    "hip osteoarthritis": ("Weak polygenic tendency (OR ~1.2–1.3 per SD). Weight, joint injury and activity matter "
                           "more; no screening indicated."),
    "knee osteoarthritis": ("Weak polygenic tendency (OR ~1.2 per SD). Weight management and quadriceps "
                            "strengthening; no screening indicated."),
    "rheumatoid arthritis": ("Modest polygenic tendency (HR ~1.2 per SD; mostly HLA). Low absolute risk; evaluate "
                             "promptly if persistent small-joint swelling or morning stiffness (RF/anti-CCP). "
                             "Not smoking is the main modifiable factor."),
    "hemorrhoid": ("Weak polygenic tendency to hemorrhoids. Fiber, fluids and avoiding straining; rectal bleeding still "
                   "warrants evaluation rather than attribution to hemorrhoids."),
}
NOT_RISK = ("height", "body mass index")
# 1000 Genomes European sub-populations that pgs.py may also rank against (shown next to the EUR percentile).
EUR_SUBPOP = {"IBS": "Iberian", "TSI": "Tuscan", "GBR": "British", "CEU": "Utah CEPH", "FIN": "Finnish"}


def ordinal(n):
    n = round(n)
    return f"{n}{'th' if 10 <= n % 100 <= 20 else {1: 'st', 2: 'nd', 3: 'rd'}.get(n % 10, 'th')}"


def pgs_section(pgs, amd_row, personal=(), who=None):
    items = []
    for pid, r in pgs.items():
        trait = r["trait"]
        t = trait.lower().replace(" (t2d)", "").replace(" (bmi)", "")
        pct = r["groups"]["EUR"]["percentile"]
        if any(t.startswith(x) for x in NOT_RISK):
            continue
        if "speed of sound" in t:  # heel ultrasound bone score: higher = better
            if pct <= 20:
                items.append(("Bone density (heel ultrasound SOS) — low genetic score", pid, pct, r,
                              "Higher = better bone; this score is low (validated against incident fractures, OR ~0.7 "
                              "per SD higher). Consider DXA if other risk factors are present (e.g. early or surgical "
                              "menopause, aromatase inhibitors, steroids, low weight, prior fracture)."))
            continue
        if "glomerular" in t:
            if pct <= 20:
                items.append((f"Kidney function (eGFR, creatinine-based) — low genetic eGFR", pid, pct, r,
                              "Higher = better filtration; this score is low. Routine creatinine/eGFR and urine ACR."))
            continue
        if pct >= 80:
            advice = PGS_ADVICE.get(t, "Consider in the context of family history and routine screening.")
            if t == "breast cancer" and "breast_cancer" in personal:
                which = text(f"doctor.personal_cancer.{who}", "breast cancer")
                advice = (f"Consistent with {pron(who, 'poss')} personal history, but a polygenic score cannot explain "
                          f"{which}; germline testing (see family history) is the relevant step.")
            elif t == "breast cancer":
                pub = max(r.get("variants_published") or 0, r["variants_total"])  # as in pgs.py
                miss = pub - r["variants_used"]
                advice += (f" This research-grade score ({miss} of {pub} published variants missing, imputed "
                           "data) should not be entered into CanRisk as a PRS; use it only as a prompt for a "
                           "family-history-based risk assessment.")
            items.append((trait.replace(" (T2D)", "").replace("Cancer", "cancer"), pid, pct, r, advice))
    items.sort(key=lambda x: -x[2])
    if items:
        def eur(r):  # a 1000 Genomes European sub-population, if pgs.py ranked against one
            code = next((c for c in EUR_SUBPOP if c in r["groups"]), None)
            return f" · {EUR_SUBPOP[code]} ref. {ordinal(r['groups'][code]['percentile'])}" if code else ""
        body = "".join(
            f"<tr><td>{esc(name)}</td><td>{ordinal(pct)}</td><td class=small>{esc(pid)} · "
            f"{r['variants_used']}/{max(r.get('variants_published') or 0, r['variants_total'])} published variants used"
            f"{eur(r)}</td><td>{esc(advice)}</td></tr>"
            for name, pid, pct, r, advice in items)
        out = (f'<table><thead><tr><th>Score</th><th>Percentile (European ref.)</th><th>Source</th>'
               f'<th>Suggested consideration</th></tr></thead><tbody>{body}</tbody></table>')
    elif not pgs:
        out = ("<p>Polygenic scores not computed (no imputed data, or the scoring step was skipped). This is not a "
               "negative result.</p>")
    else:
        out = "<p>No disease-risk score at or above the 80th percentile, and no low eGFR or bone score.</p>"
    out += "" if not pgs else ("<p class=small>Polygenic scores (PGS Catalog) on imputed data, ranked against the ~500 European (1000 Genomes EUR) "
            "samples. Shown: disease/risk scores ≥80th percentile, plus eGFR and heel-ultrasound bone (SOS) scores "
            "≤20th, where higher is better; the 80th-percentile cut-off is a display highlight, not a clinical "
            "threshold. These are modest relative tendencies, not diagnoses: typically ~1.3–2.5× average at the top "
            "decile (higher end for the breast/prostate/colorectal/melanoma scores), less just above the 80th "
            "percentile and for weak scores; height and BMI omitted. Imputed dosages pull percentiles slightly toward "
            "the median, so extremes are if anything understated.</p>")
    if amd_row and "raised" in amd_row["result"].lower():
        risk = amd_row["result"].split(" — ")[0]
        common = "Common genotype; " if risk.startswith("Somewhat") else ""
        out += (f"<p class=small><b>Single-variant marker:</b> age-related macular degeneration, CFH/ARMS2 "
                f"({esc(show_genotypes(amd_row['genotypes']))}): {esc(risk.lower())}. {common}"
                f"routine eye exams and not smoking.</p>")
    return out, items


COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
# plus-strand chip letters -> MTHFR coding-strand alleles (gene is on the minus strand)
MTHFR = {"rs1801133": ("677", "C", "T"), "rs1801131": ("1298", "A", "C")}


def mthfr_genotype(geno):
    parts = []
    for rs, g, _ in split_genotypes(geno):
        if rs not in MTHFR or len(g) != 2:
            parts.append(f"{rs}={g}")
            continue
        pos, ref, alt = MTHFR[rs]
        coding = sorted((COMP.get(a, a) for a in g), key=lambda a: a != ref)
        n = sum(a == alt for a in coding)
        parts.append(f"{pos}{''.join(coding)} ({['no variant', 'heterozygous', 'homozygous'][n]})")
    return ", ".join(parts)


def nutrition_section(well):
    W = {r["report"]: r for r in well}
    rows, flags, typical = [], [], []
    m = W.get("Folate / homocysteine")
    if m:
        status = m["result"].split(":")[0]
        rows.append(("MTHFR", mthfr_genotype(m["genotypes"]),
                     f"{status}. ACMG (2013/2020) and ACOG: MTHFR genotype has no established clinical utility — "
                     f"not an indication for thrombophilia work-up or homocysteine testing; no management change "
                     f"beyond a normal folate-adequate diet."))
    v = W.get("Vitamin D")
    if v and "lower" in v["result"].lower():
        rows.append(("GC (vitamin D binding protein)", show_genotypes(v["genotypes"]),
                     "Genetic tendency to lower 25-OH vitamin D. Consider measuring levels if clinically indicated."))
        flags.append("vitamin D")
    elif v:
        typical.append("vitamin D (GC)")
    b = W.get("Vitamin B12")
    if b and "lower" in b["result"].lower():
        rows.append(("FUT2 (vitamin B12)", show_genotypes(b["genotypes"]),
                     "Genetic tendency to lower serum B12. Consider measuring levels if symptoms, low "
                     "animal-food intake, or metformin/PPI use."))
        flags.append("vitamin B12")
    elif b:
        typical.append("vitamin B12 (FUT2)")
    a = W.get("Vitamin A (from β-carotene)")
    if a and "reduced" in a["result"].lower():
        rows.append(("BCMO1 (vitamin A from β-carotene)", show_genotypes(a["genotypes"]),
                     "Somewhat reduced β-carotene → vitamin A conversion; not clinically actionable with a normal "
                     "diet that includes preformed vitamin A."))
    elif a:
        typical.append("β-carotene conversion (BCMO1)")
    body = "".join(f"<tr><td>{esc(x)}</td><td class=mono>{esc(y)}</td><td>{esc(z)}</td></tr>" for x, y, z in rows)
    out = (f'<table><thead><tr><th>Marker</th><th>Genotype</th><th>Interpretation</th></tr></thead>'
           f'<tbody>{body}</tbody></table>' if rows else "<p>No results.</p>")
    if typical:
        out += f"<p class=small>Typical: {esc(', '.join(typical))}.</p>"
    return out, rows, flags


def discussion_points(flagged, by_gene, cond_rows, pgs_items, nut_flags, limited, apoe_row, person=None):
    pts = []
    genes = {g: ph for g, _, ph, _ in flagged}
    other = lambda g: (by_gene.get(g) or {}).get("phenotype", "")
    if "SLCO1B1" in genes:
        rosu_lim = rosu_start(genes["SLCO1B1"], other("ABCG2"))
        if is_poor(genes["SLCO1B1"]):
            c9_reduced = any(other("CYP2C9").startswith(x) for x in ("Intermediate", "Poor"))
            pts.append("If a statin is indicated (SLCO1B1 poor function, CPIC 2022): avoid simvastatin and lovastatin"
                       + (" (and fluvastatin, given CYP2C9 intermediate/poor)" if c9_reduced else "")
                       + f"; starting-dose limits atorvastatin ≤20 mg, rosuvastatin {rosu_lim}"
                       + (" (ABCG2 poor function also present)" if rosu_lim == "≤10 mg" else "")
                       + ", pravastatin ≤40 mg, pitavastatin ≤1 mg" + ("" if c9_reduced else ", fluvastatin ≤40 mg")
                       + "; if higher intensity is needed, consider rosuvastatin or add non-statin therapy "
                       "(e.g. ezetimibe) rather than withholding treatment.")
        else:
            fl = fluva(genes["SLCO1B1"], other("CYP2C9"))
            pts.append("If a statin is indicated (SLCO1B1 decreased function, CPIC 2022): prefer a statin other than "
                       "simvastatin; if simvastatin is used, keep the dose below 20 mg/day. Starting-dose limits: "
                       "atorvastatin ≤40 mg, pitavastatin ≤2 mg; "
                       + (f"rosuvastatin {rosu_lim} (ABCG2 poor function also present)" if rosu_lim else
                          "rosuvastatin at the usual starting dose (myopathy caution above 20 mg)")
                       + (f"; with CYP2C9 {short_pheno(other('CYP2C9')).lower()}, fluvastatin {fl}" if fl else "")
                       + ". Avoid lovastatin or keep it ≤20 mg/day. Watch for myopathy.")
    elif "ABCG2" in genes:
        pts.append("If rosuvastatin is used: ≤20 mg starting dose (ABCG2 poor function; CPIC optional; ≤10 mg if "
                   "SLCO1B1 function is also reduced)." if is_poor(genes["ABCG2"])
                   else "ABCG2 decreased function: no rosuvastatin dose change (CPIC); associated with higher urate/gout risk.")
    if "CYP2C19" in genes:
        ph = genes["CYP2C19"]
        dip = by_gene["CYP2C19"]["diplotype"]
        if "Poor" in ph:
            pts.append("Clopidogrel: little or no activation (CYP2C19 poor metabolizer) — CPIC: avoid clopidogrel if "
                       "possible (for ACS/PCI use prasugrel or ticagrelor if no contraindication; prasugrel "
                       "contraindicated after stroke/TIA; for neurovascular indications consider ticagrelor).")
        elif "Intermediate" in ph:
            pts.append("Clopidogrel: reduced activation (CYP2C19) — for ACS/PCI avoid standard-dose clopidogrel if "
                       "possible (prasugrel or ticagrelor if no contraindication; prasugrel contraindicated after "
                       "stroke/TIA); for neurovascular indications consider ticagrelor.")
        elif ph.startswith("Ultrarapid"):
            pts.append("CYP2C19 ultrarapid: check citalopram/escitalopram, tertiary-amine TCA and voriconazole choice; "
                       "PPIs may need a higher dose; clopidogrel and sertraline at standard dosing.")
        else:
            pts.append(f"CYP2C19 rapid metabolizer ({dip}): avoid voriconazole as primary therapy; increase PPI dose "
                       "for H. pylori/erosive esophagitis; consider a non-CYP2C19 tricyclic; clopidogrel and SSRIs at "
                       "standard dosing.")
    if "TPMT" in genes or "NUDT15" in genes:
        which = [g for g in ("TPMT", "NUDT15") if g in genes]
        poor = [g for g in which if is_poor(genes[g])]
        if poor:
            pts.append(f"Thiopurines: CPIC — {' and '.join(poor)} poor metabolizer: drastically reduced doses for "
                       "malignancy, an alternative agent for nonmalignant conditions; confirm with clinical "
                       "TPMT/NUDT15 testing first.")
        elif len(which) == 2:
            pts.append("Thiopurines: CPIC — compound TPMT + NUDT15 intermediate: start azathioprine/mercaptopurine/"
                       "thioguanine at 20–50% of the standard starting dose; confirm with clinical TPMT/NUDT15 "
                       "testing first.")
        else:
            pts.append(f"Thiopurines: CPIC — start azathioprine/mercaptopurine at 30–80% and thioguanine at 50–80% "
                       f"of normal dose (thioguanine: if the standard start is ≥40 mg/m²/day) ({which[0]} intermediate); confirm "
                       f"with clinical TPMT/NUDT15 testing first.")
    if "G6PD" in genes:
        status = genes["G6PD"].lower().replace(" (", ", ").rstrip(")")
        pts.append(f"G6PD variant detected on a consumer chip ({status}): quantitative G6PD "
                   "activity test before rasburicase, primaquine/tafenoquine, dapsone or other oxidant drugs.")
    if any(g in genes for g in WARFARIN):
        pts.append("Warfarin: use genotype-guided dosing (VKORC1/CYP2C9/CYP4F2) for the initial dose; titrate by INR.")
    for g in genes:
        if g == "HLA-B*58:01":
            pts.append("HLA-B*58:01 tag SNP present" + (" (low predictive value in European ancestry)"
                                                         if mostly_european(person) else "")
                       + ": order HLA-B*58:01 typing before allopurinol.")
        elif g.startswith("HLA"):
            pts.append(f"{g} tag SNP present: confirm by HLA typing before any {HLA_DRUG.get(g, DRUGS[g])} prescription.")
    d6 = by_gene.get("CYP2D6") or {}
    if d6.get("diplotype", "").startswith("not covered"):
        pts.append("CYP2D6 not tested: a clinical PGx panel (with copy-number analysis) is needed before relying on "
                   "genotype for codeine, tramadol, tamoxifen or many antidepressants.")
    elif d6:
        pts.append("CYP2D6 incompletely tested (no copy-number analysis" +
                   ("; *10 not called" if "rs1065852" in d6.get("missing_snps", "") else "") +
                   "): clinical CYP2D6 test before codeine, tramadol, tamoxifen or CYP2D6-dependent antidepressants.")
    if limited:
        pre = []
        if "DPYD" in limited:
            pre.append("DPYD testing before fluoropyrimidines")
        if "TPMT" in limited or "NUDT15" in limited:
            pre.append("TPMT/NUDT15 testing before thiopurines")
        if "G6PD" in limited:
            pre.append("G6PD activity before rasburicase/primaquine/dapsone")
        if pre:
            pts.append("Only a few variants were tested in DPYD/TPMT/NUDT15/G6PD: use clinical "
                       + ", ".join(pre) + ".")
    for name, gene, geno, interp in cond_rows:
        if interp.startswith("One S allele") or "PiMZ" in interp:
            pts.append("Alpha-1 antitrypsin carrier: reinforce not smoking; serum AAT level only if lung/liver findings.")
        if gene == "LPA":
            pts.append("Lipoprotein(a): measure once if never done (EAS 2022 consensus, ESC guidelines)"
                       + (" — genetic risk allele present." if interp.startswith("LPA risk") else
                          "; the absence of the two common risk alleles does not rule out high Lp(a).")
                       + " Elevated Lp(a) is a risk enhancer (ACC/AHA 2018) that favors earlier, more intensive LDL lowering.")
        if gene == "F5 / F2" and "×" in interp:
            fvl = re.search(r"Factor V Leiden ×(\d)", interp)
            pt = re.search(r"Prothrombin G20210A ×(\d)", interp)
            n_fvl, n_pt = (int(m.group(1)) if m else 0 for m in (fvl, pt))
            high = n_fvl == 2 or n_pt == 2 or (n_fvl and n_pt)
            zyg = {1: "heterozygous", 2: "homozygous"}
            kind = " + ".join(f"{zyg[n]} {what}" for n, what in ((n_fvl, "Factor V Leiden"),
                                                                (n_pt, "prothrombin G20210A")) if n)
            estrogen = {"M": "", "F": " Weigh it before"}.get(sex(person), " In women, weigh it before")
            estrogen += ("" if not estrogen else " estrogen-containing contraception or HRT (progestin-only or "
                         "non-hormonal options carry less VTE risk) and in pregnancy/postpartum planning.")
            pts.append(f"Thrombophilia: {kind} on a consumer chip — confirm with a clinical test (F5/F2 genotyping). "
                       "No anticoagulation for the genotype alone; consider VTE prophylaxis around surgery, "
                       "immobilisation and long-haul travel." + estrogen
                       + (" This is a higher-risk genotype: hematology advice on prophylaxis"
                          + {"M": "", "F": " and on estrogen use"}.get(sex(person), " (and, in women, on estrogen use)")
                          + "." if high else ""))
        gba_only = gene == "LRRK2 / GBA" and "LRRK2" not in interp  # already covered by the Gaucher carrier point
        if gene == "HOXB13" and sex(person) == "F" and interp.startswith("G84E detected"):
            pts.append(f"{name}: {interp.split(':')[0]} on a consumer chip — no prostate cancer risk for "
                       f"{pron(person, 'obj')}; confirm by clinical genetic testing, then {pron(person, 'poss')} male "
                       "relatives (father, brothers, sons) can be offered targeted testing and, if positive, an earlier "
                       "PSA discussion.")
        elif (gene in ("LRRK2 / GBA", "TTR", "HOXB13") and "detected" in interp and "not detected" not in interp
                and not gba_only):
            pts.append(f"{name}: {interp.split(':')[0]} on a consumer chip — confirm by clinical genetic testing "
                       "and genetic counseling before acting.")
        if gene == "APOL1" and interp.startswith("Two APOL1"):
            pts.append("APOL1 high-risk genotype (two risk alleles): monitor blood pressure, urine albumin/creatinine "
                       "and eGFR; confirm by clinical testing.")
        if gene == "HLA-B" and interp.startswith("rs4349859 A present"):
            pts.append("HLA-B27 likely (tag SNP): no action without symptoms. With back pain lasting ≥3 months that "
                       "began before 45 and has inflammatory features (e.g. improves with activity), consider axial "
                       "spondyloarthritis (ASAS classification entry criteria); with acute anterior uveitis, consider "
                       "rheumatology referral (SpA feature). Confirm B27 by clinical typing.")
        if " — two variants" in interp:
            found = interp.split(" — ")[0]
            phase = "" if "×2" in found else " (phase unknown: they may sit on the same copy)"
            pts.append(f"{name}: {found} — two variants, possible affected status{phase}. Clinical confirmation and "
                       "evaluation; chip calls of rare variants can be false.")
        elif "Carrier:" in interp:
            pts.append(f"{name}: {interp.split('; matters for children')[0].rstrip('.')}. Confirm with clinical testing; partner carrier screening if planning "
                       "pregnancy." + (" A β-thalassemia trait causes microcytosis that can be mistaken for iron "
                                      "deficiency." if gene == "HBB" and "Carrier: β-thal" in interp
                                      and "iron deficiency" not in interp else ""))
    if apoe_row and "ε4" in apoe_row["apoe"]:
        pts.append("APOE ε4 carrier (consumer-chip result): not used to predict dementia; manage vascular risk factors. "
                   "If anti-amyloid therapy (lecanemab/donanemab) is ever considered, confirm APOE by clinical test: "
                   "ε4 carriers have a higher ARIA risk, and labels recommend testing and counseling before treatment.")
    for name, pid, pct, r, advice in pgs_items:
        extra = ""
        if name.lower() == "urate":
            if "ABCG2" in genes:
                extra += " ABCG2 Q141K (above) also raises urate."
            if "HLA-B*58:01" not in genes:
                extra += (" If allopurinol is considered, HLA-B*58:01 status here is a tag-SNP inference of limited "
                          "accuracy" + (" in European ancestry" if mostly_european(person) else "")
                          + " — use clinical HLA-B typing where indicated.")
        pts.append(f"{name} PGS {ordinal(pct)} percentile: {advice[0].lower() + advice[1:]}{extra}")
    if nut_flags:
        pts.append(f"Consider measuring {' and '.join(nut_flags)} levels if clinically indicated.")
    if flagged:
        pts.append("Confirm any medication-relevant finding with a CLIA-certified pharmacogenomic test before changing care.")
    return pts


# ---------------------------------------------------------------- page
CSS = """
:root{color-scheme:light;--fg:#1a1a1a;--dim:#555;--line:#bdbdbd;--head:#f2f2f2;--warn:#fff4dc;--warnline:#c98a00}
*{box-sizing:border-box}
html,body{background:#fff;color:var(--fg)}
body{font:13px/1.45 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;margin:0;padding:24px 16px}
main{max-width:860px;margin:0 auto}
h1{font-size:20px;margin:0 0 2px} h2{font-size:14.5px;margin:18px 0 6px;padding-bottom:2px;border-bottom:2px solid #333;
  break-after:avoid;page-break-after:avoid}
p{margin:4px 0} .small{font-size:11.5px;color:var(--dim)} .mono{font-family:ui-monospace,Menlo,Consolas,monospace;font-size:11.5px}
.meta{color:var(--dim);font-size:12px;margin-bottom:8px}
.disclaimer{border:1.5px solid var(--warnline);background:var(--warn);padding:8px 10px;margin:8px 0 4px;font-size:12.5px}
.tablewrap{overflow-x:auto}
table{border-collapse:collapse;width:100%;margin:4px 0;font-size:12px}
th,td{border:1px solid var(--line);padding:4px 6px;text-align:left;vertical-align:top}
thead th{background:var(--head);font-weight:600}
tr{break-inside:avoid;page-break-inside:avoid}
ul{margin:4px 0;padding-left:20px} li{margin:2px 0}
footer{margin-top:16px;border-top:1px solid var(--line);padding-top:6px}
a{color:#1f4f8a}
@page{size:auto;margin:14mm 12mm}
@media print{body{padding:0;font-size:11px} table{font-size:10.5px} .small,.mono{font-size:10px}
  h2{font-size:13px;margin-top:12px} .disclaimer{-webkit-print-color-adjust:exact;print-color-adjust:exact}
  thead th{-webkit-print-color-adjust:exact;print-color-adjust:exact} a{color:inherit;text-decoration:none}}
@media (max-width:600px){table{display:block;overflow-x:auto}}
"""


def skin_section(person):
    """MC1R red-hair/sun-sensitivity variants (results/pigmentation.tsv) plus family-reported skin type
    (raw/family_notes.json). Returns (html, discussion point or None)."""
    pig = next((r for r in read_tsv("pigmentation.tsv") if r["person"] == person), None)
    notes_path = os.path.join(RAW, "family_notes.json")
    observed = {}
    if os.path.exists(notes_path):
        with open(notes_path) as fh:
            observed = json.load(fh).get("observed", {}).get(person, {})
    skin = observed.get("skin")
    if not pig and not skin:
        return "", None
    strong = int(round(float(pig["mc1r_R"]))) if pig else 0
    weak = int(round(float(pig.get("mc1r_r") or 0))) if pig else 0
    geno = (f"MC1R: {strong} strong red-hair/sun-sensitivity (R) variant{'s' if strong != 1 else ''}, "
            f"{weak} weaker (r) variant{'s' if weak != 1 else ''} (R variants include R151C, R160W, D294H)"
            ) if pig else "MC1R: not assessed"
    reported = f"Family-reported (not clinically assessed): {skin}." if skin else ""
    point = None
    if strong >= 1:
        per = "per variant (higher, ~2–3×, with two)" if strong == 1 else "(two R variants: ~2–3× or more)"
        interp = (f"MC1R R-variant carriers have roughly 1.5–2× relative melanoma risk {per}, seen even in people with "
                  "darker skin or non-red hair, and higher keratinocyte-cancer (BCC/SCC) risk. Easy burning, freckling "
                  "and high mole count are additional risk factors. Absolute risk remains modest and depends mainly on "
                  "phenotype and sun exposure.")
        lower = (skin or "").lower()
        traits = [t for t, k in (("sun-sensitive", "sensitive"), ("freckled", "freckle")) if k in lower]
        point = ("Skin: MC1R R-variant carrier" + (f" with {' and '.join(traits)} skin (family-reported)" if traits else "")
                 + " — counsel sun protection; MC1R genotype alone is not a guideline indication for screening "
                   "(USPSTF: insufficient evidence for routine skin exams in average-risk adults), but consider a "
                   "baseline full-skin examination, with follow-up skin checks if the examination or other risk "
                   "factors (many/atypical moles, personal or family history of melanoma, red hair or easy burning) "
                   "suggest increased risk.")
    elif weak:
        interp = (f"No strong MC1R R variant detected; {weak} weaker r variant{'s' if weak != 1 else ''} (at most a "
                  "small effect: meta-analyses report per-variant melanoma ORs of roughly 1–1.6 for V60L, V92M and "
                  "R163Q, with inconsistent results across studies). The array covers the common European "
                  "R/r variants only.")
    else:
        interp = "No strong MC1R R variant detected (array covers the common European R/r variants only)."
    html_text = (f'<table><tbody><tr><th scope="row">Sun sensitivity / skin cancer risk factor</th>'
                 f'<td>{esc(geno)}. {esc(reported)}</td><td>{esc(interp)}</td></tr></tbody></table>')
    return html_text, point


def relatives(person):
    """Whose reported history counts for this person: [(person id, relation as seen from this person)].
    Children see their parents and siblings, parents see themselves and their children, the maternal grandmother
    her daughter."""
    r = role(person)
    out = [(person, "Self")] if r in ("mother", "father") else []
    if r in ("self", "sibling"):
        out += [(p, rel) for p, rel in ((role_person("mother"), "Mother"), (role_person("father"), "Father")) if p]
        out += [(c, kin("sibling", sex(c)).capitalize()) for c in children() if c != person]
    elif r in ("mother", "father"):
        out += [(c, kin("child", sex(c)).capitalize()) for c in [*siblings(), role_person("self")] if c]
    elif r == "maternal_grandmother" and role_person("mother"):
        out.append((role_person("mother"), "Daughter"))
    return out


# Display labels for common raw/family_notes.json "medical_history" keys; any other key is shown humanised.
HISTORY_LABELS = {
    "heart_disease": "Heart disease", "heart_attack": "Heart attack", "stroke": "Stroke",
    "high_blood_pressure": "High blood pressure", "high_cholesterol": "High cholesterol",
    "atrial_fibrillation": "Atrial fibrillation", "blood_clots": "Blood clots (VTE)", "diabetes": "Diabetes",
    "type_1_diabetes": "Type 1 diabetes", "type_2_diabetes": "Type 2 diabetes", "kidney_disease": "Kidney disease",
    "thyroid_disease": "Thyroid disease", "asthma": "Asthma", "copd": "COPD", "osteoporosis": "Osteoporosis",
    "arthritis": "Arthritis", "glaucoma": "Glaucoma", "macular_degeneration": "Macular degeneration",
    "dementia": "Dementia", "alzheimers": "Alzheimer disease", "parkinsons": "Parkinson disease",
    "depression": "Depression", "celiac_disease": "Celiac disease", "colon_polyps": "Colon polyps",
    "hemorrhoids": "Hemorrhoids", "colorectal_cancer": "Colorectal cancer", "breast_cancer": "Breast cancer",
    "ovarian_cancer": "Ovarian cancer", "prostate_cancer": "Prostate cancer", "pancreatic_cancer": "Pancreatic cancer",
    "endometrial_cancer": "Endometrial cancer", "lung_cancer": "Lung cancer", "skin_cancer": "Skin cancer",
    "melanoma": "Melanoma",
}
# History items that feed the hereditary breast–ovarian cancer (HBOC) discussion points.
HBOC_CANCERS = ("breast_cancer", "ovarian_cancer")


def history_label(item):
    return HISTORY_LABELS.get(item) or item.replace("_", " ").capitalize()
PANEL = ("e.g. BRCA1, BRCA2, PALB2, CDH1, PTEN, STK11, TP53, ATM, CHEK2, RAD51C, RAD51D, BRIP1 and the Lynch "
         "genes (MLH1, MSH2, MSH6, PMS2, EPCAM); the genetic counselor selects the final panel")
PANCREAS = ("pancreatic cancer screening (MRI/MRCP or EUS at an experienced center), from 50 or 10 years before the "
            "earliest family diagnosis: consider for BRCA2 carriers even without a family history, and for BRCA1 "
            "carriers only with a family history of exocrine pancreatic cancer (current NCCN); skin checks (BRCA2 and "
            "melanoma)")
# NCCN surveillance for a BRCA1/2 carrier, by sex, before the shared pancreatic/skin items.
CARRIER_SURVEILLANCE = {
    "M": ("breast self-exam training and annual clinical breast exam from 35; prostate cancer screening from 40 "
          "(recommended for BRCA2, consider for BRCA1); "),
    "F": ("breast MRI from 25 and annual mammography plus MRI from 30; discuss risk-reducing salpingo-oophorectomy "
          "(35–40 for BRCA1, 40–45 for BRCA2) and risk-reducing mastectomy; "),
    None: "sex-specific breast, ovarian or prostate surveillance; ",
}


def load_history():
    path = os.path.join(RAW, "family_notes.json")
    return json.load(open(path)).get("medical_history", {}) if os.path.exists(path) else {}


def chip_markers(person, panel):
    path = os.path.join(RESULTS, "hereditary_markers.tsv")
    rows = [r for r in read_tsv("hereditary_markers.tsv") if r["person"] == person and r["panel"] == panel] \
        if os.path.exists(path) else []
    if not rows:
        return ""
    tested = ", ".join(r["marker"] for r in rows if r["alt_copies"] != "")
    found = [f"{r['marker']} ({r['alt_copies']} cop{'y' if r['alt_copies'] == '1' else 'ies'})"
             for r in rows if r["alt_copies"] not in ("", "0")]
    return f"{'; '.join(found) if found else 'none detected'} (tested: {tested})"


def plural(word, n):
    return word if n == 1 else {"child": "children"}.get(word, word + "s")


def cascade_relatives(person):
    """'her sons', 'his children and mother'...: the relatives in the data who would be offered cascade testing."""
    kind = {"son": "child", "daughter": "child", "brother": "sibling", "sister": "sibling", "mother": "parent",
            "father": "parent"}
    groups = {}  # child/sibling/parent -> relation words; mixed sexes collapse to the neutral word
    rels = [rel for _, rel in relatives(person) if rel != "Self"]
    if role(person) == "mother" and role_person("maternal_grandmother"):
        rels.append("Mother")
    for rel in rels:
        groups.setdefault(kind.get(rel.lower(), rel.lower()), []).append(rel.lower())
    names = [plural(ws[0] if len(set(ws)) == 1 else k, len(ws)) for k, ws in groups.items()] or ["first-degree relatives"]
    joined = names[0] if len(names) == 1 else ", ".join(names[:-1]) + " and " + names[-1]
    return text(f"doctor.hboc_cascade.{person}", f"{pron(person, 'poss')} {joined}")


def hboc_points(person, cancers_self, cancers_mother, cancers_daughter, pgs):
    """Discussion points for a breast/ovarian cancer history (hereditary breast–ovarian cancer pattern).
    Each cancers_* list holds (history item, family-reported text)."""
    she, her = pron(person), pron(person, "poss")
    obj = pron(person, "obj")
    if cancers_self:
        kinds = {i for i, _ in cancers_self}
        both = {"breast_cancer", "ovarian_cancer"} <= kinds
        qualifies = text(f"doctor.hboc_breast_qualifies.{person}", f"breast cancer at ≤50 would qualify {obj}")
        ovary = ""
        if "ovarian_cancer" in kinds:
            ovary = (("If the ovarian tumor was invasive epithelial, breast plus ovarian cancer in one woman strongly "
                      "suggests a hereditary breast–ovarian cancer syndrome (most often BRCA1/2; also PALB2, RAD51C/D, "
                      "BRIP1, Lynch genes). " if both else "")
                     + "Confirm the ovarian tumor type from records: BRCA-associated ovarian cancer is uncommon before "
                     "40, and ovarian tumors at young ages are more often germ-cell, stromal or borderline, which are "
                     "not BRCA-associated" + (f"; {qualifies} for testing either way. " if both else ". "))
        return [f"Personal history (family-reported): {'; '.join(t for _, t in cancers_self)}. Breast cancer at ≤50 "
                "meets NCCN criteria for germline testing, and so does epithelial ovarian (including fallopian "
                f"tube/peritoneal) cancer at any age. {ovary}Ask whether testing was done at diagnosis; if not, or if "
                f"only BRCA1/2 were tested, refer for genetic counseling and a multigene panel ({PANEL}). A pathogenic "
                f"variant would guide {her} own surveillance (contralateral breast; risk-reducing salpingo-oophorectomy "
                f"if ovarian tissue remains) and cascade testing of {cascade_relatives(person)}. Bone health: if "
                + ("ovarian cancer treatment caused early or surgical menopause (or with aromatase-inhibitor use)"
                   if "ovarian_cancer" in kinds else
                   "cancer treatment caused early menopause, or with aromatase-inhibitor use")
                + ", consider DXA per osteoporosis guidelines."]
    if cancers_mother:
        kinds = {i for i, _ in cancers_mother}
        pattern = (", a pattern strongly suggestive of a hereditary breast–ovarian cancer syndrome"
                   if {"breast_cancer", "ovarian_cancer"} <= kinds else
                   "; check her history against NCCN germline-testing criteria")
        prostate = next((r for r in pgs.values() if r["trait"].lower() == "prostate cancer"), None)
        pro = (f" The polygenic prostate score ({ordinal(prostate['groups']['EUR']['percentile'])} percentile) is not "
               "reassuring here: a BRCA2 variant would outweigh it." if prostate and sex(person) == "M" else "")
        is_ = "are" if she == "they" else "is"
        self_ = {"he": "himself", "she": "herself"}.get(she, "themselves")
        return [f"Mother: {'; '.join(t for _, t in cancers_mother)} (family-reported){pattern}. Testing should start "
                f"with her; if she carries a pathogenic variant, offer {obj} targeted testing (50% chance of inheriting "
                f"it). If {she} {is_} a BRCA1/2 carrier (NCCN): {CARRIER_SURVEILLANCE.get(sex(person), CARRIER_SURVEILLANCE[None])}"
                f"{PANCREAS}. If a Lynch gene is found, NCCN colonoscopy timing depends on the gene: every 1–2 years "
                "from 20–25 for MLH1/MSH2/EPCAM, every 1–3 years from 30–35 for MSH6/PMS2 (earlier if a relative had "
                "colorectal cancer young), plus gene-specific surveillance for other cancers. If she cannot be tested, "
                f"{she} can have a multigene panel {self_}." + pro]
    if cancers_daughter:
        hist = load_history()
        sibs = [rel.lower() for k, rel in relatives(person) if rel in ("Son", "Daughter", "Child")
                and not (rel == "Daughter" and any(i in HBOC_CANCERS for i in hist.get(k, {})))]
        sib_words = (plural(sibs[0] if len(set(sibs)) == 1 else "child", len(sibs)) if sibs else "")
        others = text(f"doctor.hboc_daughter_siblings.{person}",
                      f"her siblings (this patient's other {sib_words})" if sibs else "any other children of this patient")
        others_short = text(f"doctor.hboc_daughter_siblings_short.{person}", "them" if sibs else "any other children")
        age = text(f"doctor.daughter_breast_age.{person}", ", especially before 50,")
        breast = (" Continue age-appropriate breast screening; a first-degree relative with breast cancer"
                  f"{age} is a recognized risk factor."
                  if sex(person) == "F" and any(i == "breast_cancer" for i, _ in cancers_daughter) else "")
        return [f"Daughter: {'; '.join(t for _, t in cancers_daughter)} (family-reported). Her germline testing comes "
                f"first. If she carries a pathogenic variant, {others} should be offered targeted testing whatever the "
                "result here (50% chance). Testing this patient shows which side of the family the variant came from, "
                "which tells more distant relatives (this patient's siblings and their families, or the other parent's "
                f"side) whether they need testing; a negative result here does not clear {others_short}.{breast}"]
    return []


def family_history_section(person, pgs):
    """Family-reported history (raw/family_notes.json "medical_history") plus chip markers
    (results/hereditary_markers.tsv). Returns (html, discussion points)."""
    hist = load_history()
    rows, points, fdr_polyps = [], [], []
    cancers = {"Self": [], "Mother": [], "Daughter": []}
    polyps_self = "colon_polyps" in hist.get(person, {})
    for key, relation in relatives(person):
        for item, details in hist.get(key, {}).items():
            # family_notes.json is hand-written: accept true/null/numbers as well as text
            details = details if isinstance(details, str) else "" if details in (None, True, False) else str(details)
            rows.append((relation, history_label(item), details))
            if item == "colon_polyps" and key != person:
                fdr_polyps.append(relation.lower())
            if item in HBOC_CANCERS and relation in cancers:
                cancers[relation].append((item, details[:1].lower() + details[1:]))
    if not rows:
        return "", []
    body = "".join(f"<tr><td>{esc(a)}</td><td>{esc(b)}</td><td>{esc(c)}</td></tr>" for a, b, c in rows)
    out = (f"<table><thead><tr><th>Relative</th><th>Reported history</th><th>Details</th></tr></thead>"
           f"<tbody>{body}</tbody></table>"
           f"<p class=small>{esc(hist.get('_about', 'Family-reported, not from medical records.'))}</p>")
    has_cancer = any(cancers.values())
    for panel, title, show in (("breast-ovarian", "Chip BRCA founder variants", has_cancer),
                               ("colorectal", "Chip colorectal markers", bool(fdr_polyps or polyps_self))):
        found = chip_markers(person, panel) if show else ""
        if found:
            out += f"<p class=small><b>{title}:</b> {esc(found)}.</p>"
    if out.count("<b>Chip"):
        out += ("<p class=small>Arrays test only a few founder variants. They do not test Lynch syndrome or the many "
                "other BRCA1/2, PALB2, APC and MUTYH variants, so a negative chip result excludes almost nothing.</p>")
    points += hboc_points(person, cancers["Self"], cancers["Mother"], cancers["Daughter"], pgs)

    crc = next((r for r in pgs.values() if r["trait"].lower() == "colorectal cancer"), None)
    crc_txt = (f" Polygenic colorectal cancer score: {ordinal(crc['groups']['EUR']['percentile'])} percentile "
               "(European reference), a modest tendency that does not by itself change screening." if crc else "")
    genetics = ("If any relative had ≥10 cumulative adenomas, refer for genetic counseling and a clinical multigene "
                "panel (APC, MUTYH and others), testing the affected relative first. If colorectal or endometrial "
                "cancer occurs in the family, the affected relative should have tumor MMR/MSI and/or germline testing; "
                "refer this patient for genetic counseling if the family history meets Lynch criteria (e.g. a first- "
                "or second-degree relative with a Lynch-spectrum cancer before 50, several relatives with "
                "Lynch-spectrum cancers, or PREMM5 ≥2.5%).")
    if polyps_self:
        points.append("Colon polyps (personal history, family-reported): obtain colonoscopy and pathology reports; "
                      "surveillance interval per USMSTF 2020 according to number, size and histology of the removed "
                      f"polyps. {(pron(person, 'poss').capitalize() + ' ' + ' and '.join(fdr_polyps) + ' reportedly also had polyps removed. ') if fdr_polyps else ''}"
                      f"{genetics}{crc_txt}")
    elif fdr_polyps:
        who = " and ".join(fdr_polyps)
        two = len(fdr_polyps) > 1
        points.append(f"Colorectal screening — family history: {who} reportedly had colon polyps removed "
                      f"({'two' if two else 'one'} first-degree relative{'s' if two else ''}). Obtain the pathology "
                      "(adenoma, sessile serrated/traditional serrated or hyperplastic, size, villous features, dysplasia) "
                      "and ages at removal. Per USMSTF 2017: colorectal cancer, an advanced adenoma (≥10 mm, villous or "
                      "high-grade dysplasia) or an advanced serrated lesion (sessile serrated lesion/traditional serrated "
                      "adenoma ≥10 mm or with dysplasia) "
                      + ("in one first-degree relative before 60, or in two at any age" if two else "in a first-degree relative before 60")
                      + " → colonoscopy every 5 years starting at 40 or 10 years before the youngest affected relative's "
                      "diagnosis, whichever is earlier; "
                      + ("" if two else "the same findings in one first-degree relative at 60 or older → start at 40 "
                         "at average-risk intervals; ")
                      + "non-advanced adenomas or small hyperplastic polyps only → average-risk screening, which "
                      "starts at 45 under USPSTF 2021 / the 2022 USMSTF update. " + genetics + crc_txt)
    if any(item == HISTORY_LABELS["hemorrhoids"] for _, item, _ in rows):
        points.append("Hemorrhoids in the family (family-reported) do not change colorectal screening; rectal bleeding "
                      "should still be evaluated rather than attributed to hemorrhoids.")
    return out, points


def page(title, body):
    return (f'<!doctype html><html lang="en"><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width,initial-scale=1">'
            f'<meta name="color-scheme" content="light"><title>{esc(title)}</title><style>{CSS}</style></head>'
            f'<body><main>{body}</main></body></html>')


def person_page(person, date, pgx, apoe, well, pgs, n_people=2):
    name = patient(person)
    chips, has_imputed = data_sources(person)
    pgx_html, flagged, by_gene, limited = pgx_section(pgx)
    cond_html, cond_rows = conditions_section(well, apoe, person)
    amd = next((r for r in well if r["report"] == "Macular degeneration"), None)
    hist = load_history()
    pgs_html, pgs_items = pgs_section(pgs, amd, hist.get(person, {}), person)
    nut_html, nut_rows, nut_flags = nutrition_section(well)
    pts = discussion_points(flagged, by_gene, cond_rows, pgs_items, nut_flags, limited, apoe, person)
    skin_html, skin_point = skin_section(person)
    if skin_point:
        pts.append(skin_point)
    fam_html, fam_points = family_history_section(person, pgs or {})
    hem = next((i for i, p in enumerate(pts) if p.startswith("Hemorrhoid PGS")), None)
    if hem is not None:  # one hemorrhoid bullet: the score's, noting the family history
        fam_points = [p for p in fam_points if not p.startswith("Hemorrhoids in the family")]
        # the relatives shown in this page's family table who reported hemorrhoids
        reported = [hist.get(k, {})["hemorrhoids"] for k, rel in relatives(person)
                    if rel != "Self" and "hemorrhoids" in hist.get(k, {})]
        what = text(f"doctor.hemorrhoid_history.{person}", "hemorrhoids" if reported else "")
        if fam_html and what:
            pts[hem] += f" Family history of {what} (family-reported) is consistent with this tendency."
    pts.extend(fam_points)

    imputed = []
    if pgs:
        imputed.append("all polygenic scores")
    if apoe and "imputed" in apoe["source"]:
        imputed.append("APOE (chip no-call)" if "no-call" in apoe["source"] else f"APOE ({apoe['source']})")
    filled = [SNP_NAME.get(rs, f"{r['gene']} {rs}") for r in well
              if r["category"] in ("Inherited conditions", "Nutrition & health", "Taste & diet")
              for rs, _, imp in split_genotypes(r["genotypes"]) if imp]
    filled += [f"{r['gene']} {rs}" for r in pgx if r["person"] == person
               for rs in r.get("imputed_snps", "").split(",") if rs]
    if filled:
        imputed.append("phased/imputed fill-in for chip no-calls: " + ", ".join(dict.fromkeys(filled)))
    src = "; ".join(chips) or "consumer genotyping chip"
    if has_imputed:
        src += "; 23andMe imputed data"
    meta = f"Data: {esc(src)}."
    if imputed:
        meta += f" Imputed items: {esc('; '.join(imputed))}. All other results are direct chip calls."
    if date:
        meta += f" Generated {esc(date)}."

    body = f"""<h1>{esc(name)} — DNA summary for clinicians</h1>
<p class="meta">{meta}</p>
<div class="disclaimer"><b>For learning and research only — not medical advice and not a clinical test.</b> Results come from a
direct-to-consumer SNP array analysed locally{' by family members' if n_people > 1 else ''} with open-source research code that has
not been clinically validated. A clinician must review these findings and confirm any actionable one with a
CLIA-certified/clinical test before changing care.</div>
<h2>1. Medication-relevant findings (pharmacogenomics)</h2><div class="tablewrap">{pgx_html}</div>
<h2>2. Inherited conditions screened (and APOE)</h2><div class="tablewrap">{cond_html}</div>
{f'<div class="tablewrap">{skin_html}</div>' if skin_html else ''}
<h2>3. Polygenic risk tendencies</h2><div class="tablewrap">{pgs_html}</div>
<h2>4. Nutrition-related markers</h2><div class="tablewrap">{nut_html}</div>
{f'<h2>5. Family history (reported by the family)</h2><div class="tablewrap">{fam_html}</div>' if fam_html else ''}
<h2>{6 if fam_html else 5}. Suggested discussion points</h2>
<ul>{''.join(f'<li>{esc(p)}</li>' for p in pts) or '<li>No specific items.</li>'}</ul>
<footer class="small"><b>Limitations.</b> Genotyping arrays test selected common variants only: rare pathogenic variants,
copy-number changes (e.g. CYP2D6 deletions/duplications) and HLA alleles are not directly assayed. Pharmacogenomic
phase is assumed where it cannot be determined. Polygenic scores were derived mainly in European-ancestry cohorts (many
in UK Biobank) and use imputed genotypes; predictive accuracy is lower in other ancestries and declines modestly even
within Europe with distance from the training population (Privé 2022); percentiles are relative to the 1000 Genomes
European reference sample. Allele calls follow CPIC definitions; actions summarise CPIC guidance (cpicpgx.org; the FDA
label where CPIC has no guideline, e.g. siponimod). Full family reports: <a href="../index.html">results/index.html</a>
(pharmacogenomics, polygenic scores, wellness tables).</footer>"""
    page_html = page(f"{name} — DNA summary for clinicians", body)
    counts = {"pgx_flagged": len(flagged), "conditions": len(cond_rows), "pgs_flagged": len(pgs_items),
              "nutrition": len(nut_rows), "discussion": len(pts), "pgs_run": int(bool(pgs))}
    return page_html, counts


def main():
    date = None
    if "--date" in sys.argv:
        i = sys.argv.index("--date")
        date = sys.argv[i + 1] if i + 1 < len(sys.argv) else None
    pgx, apoe, wellness = read_tsv("pharmacogenomics.tsv"), read_tsv("apoe.tsv"), read_tsv("wellness.tsv")
    found = people_with_results(pgx, apoe, wellness)
    os.makedirs(OUT, exist_ok=True)
    links = []
    for p in found:
        html_text, counts = person_page(
            p, date, [r for r in pgx if r["person"] == p], next((r for r in apoe if r["person"] == p), None),
            [r for r in wellness if r["person"] == p], load_pgs(p), len(found))
        with open(os.path.join(OUT, f"{p}.html"), "w") as fh:
            fh.write(html_text)
        print(f"{p:17s} " + " ".join(f"{k}={v}" for k, v in counts.items() if k != "pgs_run"))
        pgs_txt = f'{counts["pgs_flagged"]} polygenic flag(s)' if counts["pgs_run"] else "polygenic scores not computed"
        links.append(f'<li><a href="{esc(p)}.html">{esc(patient(p))}</a> — {counts["pgx_flagged"]} medication '
                     f'finding(s), {pgs_txt}</li>')
    body = (f"<h1>DNA summaries for clinicians</h1><p class=meta>One printable page per person"
            f"{' · generated ' + esc(date) if date else ''}.</p>"
            f'<div class="disclaimer"><b>For learning and research only — not medical advice and not a clinical test.</b> '
            f"A clinician must review these findings and confirm actionable ones with a clinical test before changing care.</div>"
            f"<ul>{''.join(links)}</ul><p class=small>Full reports: <a href=\"../index.html\">results/index.html</a></p>")
    with open(os.path.join(OUT, "index.html"), "w") as fh:
        fh.write(page("DNA summaries for clinicians", body))
    print(f"Wrote {len(found)} page(s) → results/doctor/")


if __name__ == "__main__":
    main()
