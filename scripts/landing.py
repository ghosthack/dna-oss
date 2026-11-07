"""Landing page (results/index.html): TL;DR of all analyses + links to the full reports.

Usage: python3 -I scripts/landing.py   (run after the other scripts)
Also writes one page per sibling (results/<sibling>.html, or family.json "pages"): the same page from that sibling's
point of view ("Me" = them), once their data is in.

Family-specific wording (genealogy, medical history, speculative notes) is optional and comes from family.json:
"text" keys "landing.*" ({viewer id: string} for anything first-person; a plain string is shown on every page, so
write it in neutral third person; "mine" in it means self's) and "ancestry" (regions, founders, places). Without them
the sections stay generic or are hidden.
"""
import csv
import glob
import html
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dnalib import DISCLAIMER_HTML, FAMILY, RAW, RESULTS, WORK, children, label, people as family_people, person, role, sex, text
from pgs import LOCAL_LABELS
from pgx import DRUGS, note
from pigmentation import MC1R_STRONG

esc = html.escape
SELF = person("self")
VIEWER = SELF  # whose page this is: the "Me" column and first-person text (main() sets it per page)
PAGES = FAMILY.get("pages") or {c: "index.html" if c == SELF else f"{c}.html" for c in children()}


def order(people):
    """People in family order with the viewer first, so the viewer's column comes first."""
    ids = family_people()
    return sorted(people, key=lambda p: (p != VIEWER, ids.index(p) if p in ids else 99))


def tag(p):
    """" (Name)" after a child whose label another child shares (two brothers), else ""."""
    if role(p) not in ("self", "sibling") or sum(label(c, VIEWER) == label(p, VIEWER) for c in children()) < 2:
        return ""
    return f" ({FAMILY['people'][p].get('name') or p})"


def full(p):
    """Full label with the tag: "Grandma (maternal)", "Brother (Name)"."""
    return label(p, VIEWER) + tag(p)


def col(p):
    """Column header for a person: "Me", "Mom", "Dad", "Brother", "Sister", "Grandma"."""
    return label(p, VIEWER).split(" (")[0] + tag(p)


def rel(p):
    """How running text refers to a person: "me", "mom", "dad", "grandma", "my brother", "my sister"."""
    if p == VIEWER:
        return "me"
    if role(p) in ("self", "sibling"):
        return "my " + label(p, VIEWER).lower() + tag(p)
    return {"mother": "mom", "father": "dad", "maternal_grandmother": "grandma"}.get(role(p), p)


def poss(p):
    """Capitalized possessive: "My", "My brother's", "Dad's"."""
    return "My" if p == VIEWER else (rel(p)[:1].upper() + rel(p)[1:] + "'s")


def pron(p):
    """his / her / their, from family.json sex."""
    return {"M": "his", "F": "her"}.get(sex(p), "their")


def and_join(words):
    return " and ".join(words) if len(words) < 3 else ", ".join(words[:-1]) + " and " + words[-1]


def as_viewer(s):
    """Family notes are written from self's point of view ("like mine"); on a sibling's page "mine" is self's."""
    return s if VIEWER == SELF else re.sub(r"\bmine\b", rel(SELF) + "'s", s)


def named(p):
    """A person in titles and links: their family.json name, else a role phrase ("the main tester", "the brother")."""
    n = FAMILY["people"].get(p, {}).get("name")
    if n:
        return n
    if p == SELF:
        return "the main tester"
    same = sum(label(c, SELF) == label(p, SELF) for c in children())  # labels as self sees them, same on every page
    return f"the {label(p, SELF).lower()}" + (f" ({p})" if same > 1 else "")


def vtext(key):
    """family.json text for this page: {viewer id: text}, or one neutral string shown on every page."""
    v = text(key)
    if isinstance(v, dict):
        return v.get(VIEWER, "")
    return as_viewer(v) if v else ""


def cap(s):
    return s[:1].upper() + s[1:]


def is_male(p):
    """family.json sex, else the sex inferred from the chip (summary.json qc)."""
    s = sex(p) or (load_json("summary.json") or {}).get("qc", {}).get(p, {}).get("sex")
    return s == "M"


def read_tsv(name):
    path = os.path.join(RESULTS, name)
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_json(name):
    path = os.path.join(RESULTS, name)
    return json.load(open(path)) if os.path.exists(path) else None


def load_json_path(path):
    return json.load(open(path)) if os.path.exists(path) else None


def load_raw_json(name):
    """Family-reported facts kept next to the raw data (raw/family_notes.json)."""
    path = os.path.join(RAW, name)
    return json.load(open(path)) if os.path.exists(path) else None


VENDOR_FILES = [("23andMe", ["genome_*.zip", "genome_*.txt"]), ("AncestryDNA", ["dna-data-*.zip", "AncestryDNA*.txt"])]


def vendors():
    """Chip vendors whose raw files are in raw/<person>/ (file names only), in a fixed order."""
    return [v for v, pats in VENDOR_FILES
            if any(glob.glob(os.path.join(RAW, p, pat)) for p in family_people() for pat in pats)]


def pgs_people():
    """People with polygenic scores (results/pgs_<id>.json), in family order."""
    return order(p for p in family_people() if os.path.exists(os.path.join(RESULTS, f"pgs_{p}.json")))


def pgs_groups():
    """[(code, short label)] of the 1000 Genomes comparison groups pgs.py used, read from its output, minus "ALL"."""
    codes = []
    for p in pgs_people():
        for r in (load_json(f"pgs_{p}.json") or {}).values():
            codes += [c for c in r.get("groups", {}) if c not in codes and c != "ALL"]
    names = {**LOCAL_LABELS, "EUR": "European", "AMR": "Admixed American"}
    return [(c, names.get(c, c).split(" (")[0]) for c in codes]


def pgs_trait(r):
    """Display name of a pgs_<id>.json record's trait ("Body mass index", "Colorectal cancer")."""
    t = r["trait"].replace(" (T2D)", "").replace(" (BMI)", "")
    t = t[0].upper() + t[1:].replace("Cancer", "cancer")
    return TRAIT_LABELS.get(t, t)


BIGY_MAJOR = ["U106", "L21", "U152", "DF19", "DF99", "L238"]  # the big R-P312/P311 branches bigy_preview.py counts
P311_ISOGG = "R1b1a2a1a"  # ISOGG 2016 code of R-P311; the preview only applies to men at or below it


def y_calls():
    """{person: y_haplogroups.tsv row} (haplogroup, isogg_2016), men only."""
    return {r["person"]: r for r in read_tsv("y_haplogroups.tsv") if r.get("haplogroup") and is_male(r["person"])}


def bigy_section():
    """Preview of a FamilyTreeDNA Big Y-700: chip facts next to best guesses (every guess ends in "?").

    The chip column is computed from bigy_preview.json and y_haplogroups.tsv; guesses are generic unless
    family.json text "landing.bigy.<row>" gives the family's own."""
    d = load_json("bigy_preview.json")
    ys = y_calls()
    men = [p for p in order(d["branches"]) if p in ys and ys[p].get("isogg_2016", "").startswith(P311_ISOGG)] if d else []
    if not men:
        return ""
    subj = max(men, key=lambda p: len(ys[p]["isogg_2016"]))  # deepest call; ties keep the viewer / family order
    br, hg = d["branches"][subj], ys[subj]["haplogroup"]
    same = [subj] + [p for p in order(ys) if p != subj and ys[p]["haplogroup"] == hg]
    major = [b for b in BIGY_MAJOR if b in br]
    neg = sum(br[b]["ancestral"] for b in major)
    pos = sum(br[b]["derived"] for b in major + ["DF27"] if b in br)
    negative = [b for b in major if br[b]["ancestral"] and not br[b]["derived"]]
    positive = [b for b in major if br[b]["derived"]]

    def df27(p):  # "negative" / "positive" / "no-call", or None where the chip doesn't test it
        k = d["key_snps"].get(p, {}).get("DF27")
        if k == ".":
            return "no-call"
        c = d["branches"].get(p, {}).get("DF27", {})
        return None if k in (None, "not on chip") else "positive" if c.get("derived") else (
            "negative" if c.get("ancestral") else None)
    by = {}
    for p in same:
        if df27(p):
            by.setdefault(df27(p), []).append(rel(p))
    p312 = d["key_snps"].get(subj, {}).get("P312") == "not on chip"
    now_hg = (f"{hg} ({', '.join(rel(p) for p in same)})."
              + (f" Negative for the big branches {', '.join(negative)}" if negative else "")
              + (f"{';' if negative else ''} derived calls on {', '.join(positive)}" if positive else "")
              + (f" ({pos} derived of {neg} SNPs)." if negative or positive else ""))
    if by:
        now_hg += (" The single DF27 marker reads "
                   + ", ".join(f"{s} for {' and '.join(w)}" for s, w in by.items())
                   + ("; P312 itself isn't on the chip." if p312 else "."))
    elif p312:
        now_hg += " P312 itself isn't on the chip."
    if positive:
        hg_guess = f"R-{positive[0]} ? (the chip already points there); Big Y names the sub-branch ?"
    elif by.get("negative") and rel(subj) in by["negative"]:
        hg_guess = "A branch of R-P312 outside DF27 ? (the chip reads DF27 negative); Big Y names it ?"
    else:
        hg_guess = "A branch below R-P312 ?; Big Y names it ?"
    # Paternal-line men in the family (everyone with the same Y), plus named paternal ancestors from family.json.
    dad = person("father")
    line = (["Me"] if is_male(VIEWER) else []) + (["dad"] if dad else []) + [
        label(c, VIEWER).lower() + tag(c) for c in children() if c != VIEWER and is_male(c)]
    ancestors = text("landing.bigy.ancestors")
    covers = " = ".join(line + ([ancestors] if ancestors else []))
    covers = covers[:1].upper() + covers[1:]
    rows = [
        ("Deep haplogroup", now_hg, text("landing.bigy.haplogroup", hg_guess)),
        ("Age of your branch", "—", text("landing.bigy.age", "Not predictable from the chip ?; Big Y dates the "
                                                             "branch once the result arrives")),
        ("New (private) SNPs", "—", "~5–20 unnamed variants ?, which could define a new branch for your line ?"),
        ("Y-STRs (111 of 700+)", "—", text("landing.bigy.strs", "Not predictable from the chip ?")),
        ("Matches", "—", text("landing.bigy.matches", "Unknown until the result arrives ?; depends on who else "
                                                      "has tested")),
        ("Deep-ancestry story", "R1b-P311: descends from Bronze Age steppe-related men who spread into Western "
                                "Europe ~2500 BCE",
         text("landing.bigy.story", "Placed on FamilyTreeDNA's time tree once the branch is known ?")),
        ("Who it covers", covers or "—",
         f"Same result for all of them; only one of {'us' if is_male(VIEWER) else 'them'} needs to test ?"),
    ]
    body = "".join(f"<tr><th scope=row>{esc(a)}</th><td>{esc(b)}</td><td>{esc(c)}</td></tr>" for a, b, c in rows)
    # Men whose (older or smaller) chip stops higher up the same line.
    lower = [f'{poss(p)} "{ys[p]["haplogroup"]}" is the same line at lower resolution.' for p in order(ys)
             if p not in same and ys[subj]["isogg_2016"].startswith(ys[p].get("isogg_2016") or "?")]
    return (f'<section class="card"><h2>🔮 Big Y-700 preview (FamilyTreeDNA / Gene by Gene)</h2>'
            f'<div class="tablewrap"><table class="fam bigy"><thead><tr><th></th><th>What we know now (chip)</th>'
            f'<th>Best guess for Big Y-700</th></tr></thead><tbody>{body}</tbody></table></div>'
            f'<p class="dim">Big Y-700 sequences the ~15–23 million usable bases of the Y chromosome (the part passed down '
            f'intact): it reads the direct father-to-son line '
            f'only (no health, no mtDNA, no other ancestors). Everything marked "?" is a guess until the result '
            f'arrives; the chip column is measured.{"".join(" " + esc(s) for s in lower)}</p></section>')


def grandpa_cells():
    return load_json("grandpa_cells.json") or {}


GP_NAME = FAMILY.get("maternal_grandfather", {}).get("name")
GP_TH = (f'<th title="Virtual reconstruction of mom\'s father{esc(", " + GP_NAME, quote=True) if GP_NAME else ""}">'
         'Grandpa*</th>')
GP_NA = '<td class=dim title="{}">—</td>'


def gp_td(c):
    """Grandpa* cell from a grandpa_cells.json entry ({cell, detail}); '—' when he has no entry."""
    if not c:
        return GP_NA.format("Not reconstructable from one copy")
    return f'<td class=dim title="{esc(c["detail"], quote=True)}">{esc(c["cell"])}</td>'


# Short drug names for the medications row headers (full lists stay in the hover text and the PGx report).
DRUG_SHORT = {"CYP2C19": "clopidogrel, PPIs", "CYP2C9": "warfarin, NSAIDs", "VKORC1": "warfarin",
              "CYP4F2": "warfarin", "SLCO1B1": "statins", "ABCG2": "rosuvastatin",
              "TPMT": "thiopurines", "NUDT15": "thiopurines", "HLA-B*58:01": "allopurinol",
              "HLA-B*57:01": "abacavir", "HLA-A*31:01": "carbamazepine", "DPYD": "5-FU, capecitabine",
              "CYP2D6": "codeine, tramadol, tamoxifen", "UGT1A1": "irinotecan", "CYP3A5": "tacrolimus",
              "CYP2B6": "efavirenz", "G6PD": "rasburicase, primaquine"}


def meds_section(pgx):
    """One table: a row per gene that is actionable for anyone, a column per person (siblings join automatically)."""
    people = order(dict.fromkeys(r["person"] for r in pgx))
    cell = {(r["person"], r["gene"]): r for r in pgx}
    genes = list(dict.fromkeys(r["gene"] for r in pgx if note(r["gene"], r["phenotype"])))
    if not genes:
        return ""

    def short(gene, ph):
        ph = ph.split(" (")[0].split(" — ")[0]
        for cut in (" metabolizer", " function", " warfarin sensitivity", gene + " ", " carrier"):
            ph = ph.replace(cut, "")
        return ph[0].upper() + ph[1:]

    def td(p, g):
        r = cell.get((p, g))
        if not r:
            return "<td class=dim>—</td>"
        text = note(g, r["phenotype"], p)
        tip = esc(f'{r["diplotype"]}: {r["phenotype"]}' + (f" — {text}" if text else ""), quote=True)
        label = esc(short(g, r["phenotype"]))
        return (f'<td title="{tip}" class=act>{label}</td>' if text else f'<td title="{tip}" class=dim>{label}</td>')
    gp = grandpa_cells().get("_meds")
    head = "".join(f"<th>{esc(col(p))}</th>" for p in people) + (GP_TH if gp else "")
    body = "".join(f'<tr><th scope=row title="{esc(DRUGS[g], quote=True)}"><b>{esc(g)}</b> '
                   f'<span class=dim>{esc(DRUG_SHORT.get(g, DRUGS[g]))}</span></th>'
                   + "".join(td(p, g) for p in people) + (gp_td(gp.get(g)) if gp else "") + "</tr>" for g in genes)
    return (f'<div class="tablewrap"><table class="fam meds"><thead><tr><th></th>{head}</tr></thead>'
            f'<tbody>{body}</tbody></table></div>')


# Lay labels that state the direction where "higher" isn't simply "more risk".
TRAIT_LABELS = {
    "Estimated glomerular filtration rate": "Kidney filtration (eGFR) — higher = better",
    "Urate": "Uric acid (gout risk factor)",
    "Systolic blood pressure": "Systolic blood pressure",
    "Glaucoma": "Glaucoma",
    "Hemorrhoid": "Hemorrhoids",
    "Benign neoplasm of colon": "Colon polyps (benign growths)",
    "Heel quantitative ultrasound speed of sound (SOS)": "Bone density (heel ultrasound) — higher = better",
}
NOT_DISEASE = {"Height", "Body mass index", "Kidney filtration (eGFR) — higher = better",
               "Bone density (heel ultrasound) — higher = better"}
LOW_IS_RISK = {"Kidney filtration (eGFR) — higher = better", "Bone density (heel ultrasound) — higher = better"}
PGS_GROUPS = [  # (subgroup label, [trait labels]); unlisted traits land in "Other conditions"
    ("Heart, kidney & metabolism", ["Coronary artery disease", "Atrial fibrillation", "Systolic blood pressure",
                                    "LDL cholesterol", "Type 2 diabetes", "Uric acid (gout risk factor)",
                                    "Kidney filtration (eGFR) — higher = better"]),
    ("Cancer & growths", ["Breast cancer", "Prostate cancer", "Colorectal cancer", "Colon polyps (benign growths)",
                          "Melanoma"]),
    ("Bones & joints", ["Hip osteoarthritis", "Knee osteoarthritis", "Rheumatoid arthritis",
                        "Bone density (heel ultrasound) — higher = better"]),
    ("Other conditions", ["Glaucoma", "Hemorrhoids"]),
    ("Body traits (not diseases)", ["Height", "Body mass index"]),
]


def family_pgs_table():
    """Trait × person table of European (1000 Genomes EUR) percentiles; blank where a score doesn't apply (sex-specific)."""
    people = pgs_people()
    if not people:
        return ""
    data = {p: load_json(f"pgs_{p}.json") for p in people}
    traits = {}
    for p in people:
        for r in data[p].values():
            traits.setdefault(pgs_trait(r), {})[p] = r["groups"]["EUR"]["percentile"]
    gp = bool(grandpa_cells())
    head = "".join(f"<th>{esc(col(p))}</th>" for p in people) + (GP_TH if gp else "")
    listed = {t for _, ts in PGS_GROUPS for t in ts}
    groups = [(group, ts + sorted(t for t in traits if t not in listed) if group == "Other conditions" else ts)
              for group, ts in PGS_GROUPS]
    ncol = 1 + len(people) + gp
    rows = []
    for group, ts in groups:
        present = [t for t in ts if t in traits]
        if present:
            rows.append(f'<tr class="grp"><th colspan="{ncol}">{esc(group)}</th></tr>')
        for trait in present:
            vals = traits[trait]
            cells = []
            for p in people:
                v = vals.get(p)
                if v is None:
                    cells.append('<td class="dim">—</td>')
                else:
                    tag = (' <span class="hi">↓ low</span>' if trait in LOW_IS_RISK and v <= 20 else
                           ' <span class="hi">↑ high</span>' if v >= 80 and trait not in LOW_IS_RISK else "")
                    cells.append(f"<td>{v:.0f}{tag}</td>")
            if gp:
                cells.append(GP_NA.format("A polygenic score needs both copies of every variant"))
            rows.append(f"<tr><th scope=row>{esc(trait)}</th>{''.join(cells)}</tr>")
    # Disease scores in the top 20% (height and BMI are traits, not diseases).
    def ordinal(n):
        n = round(n)
        return f"{n}{'th' if 10 <= n % 100 <= 20 else {1: 'st', 2: 'nd', 3: 'rd'}.get(n % 10, 'th')}"
    def whose(p):
        return "My" if p == VIEWER else f"{col(p)}'s"
    flagged = [f"{whose(p)} {t[0].lower() + t[1:]} ({ordinal(v)})" for t, vals in sorted(traits.items())
               for p, v in vals.items() if (v >= 80 and t not in NOT_DISEASE) or (t in LOW_IS_RISK and v <= 20)]
    weak = " (the hemorrhoid score is especially weak)" if any("hemorrhoid" in f.lower() for f in flagged) else ""
    note_text = ("Highlighted (top 20% — a display cut-off, not a clinical threshold): " + "; ".join(flagged) +
                 f". Most are small shifts{weak}; mention to a doctor mainly where "
                 "family history points the same way — a score shifts risk, it doesn't diagnose."
                 if flagged else "No disease score in the top 20% for anyone.")
    return (f'<div class="tablewrap"><table class="fam"><thead><tr><th></th>{head}</tr></thead>'
            f'<tbody>{"".join(rows)}</tbody></table></div><p class="dim">{esc(note_text)}</p>')


WELLNESS_SHORT = [  # (substring of full result, short label)
    ("HLA-B27 not detected", "Not detected"), ("HLA-B27 likely", "Likely carrier"),
    ("tag SNPs disagree", "Inconclusive"),
    ("No LPA risk variant", "No risk variant"), ("LPA risk variant", "Likely high"),
    ("Not a carrier of", "Not a carrier"), ("Carrier: ", "Carrier"), ("two variants:", "Two variants"),
    ("Unlikely to flush", "No flush"), ("flush", "Flushes"),
    ("0 variants", "Drinks a bit less"), ("caffeine", "Average / more"),
    ("Not likely to be an especially deep", "Not especially deep"), ("especially deep", "Deep sleeper"),
    ("Not likely to be lactose", "Tolerant"), ("lactose intolerant", "Intolerant"),
    ("power athletes", "Power type"), ("endurance athletes", "Endurance type"),
    ("Likely to weigh more", "Sensitive"), ("Not likely to weigh more", "Not sensitive"),
    ("move more", "Moves more"), ("About average movement", "Average"), ("move less", "Moves less"),
    ("Not determinable", "—"),
    ("Strong taster", "Strong taster"), ("Taster (PAV/AVI)", "Taster"), ("Non-taster", "Non-taster"),
    ("Rare TAS2R38", "Intermediate"),
    ("More likely to find cilantro", "Soapy-prone"), ("Average odds of cilantro", "Average"),
    ("Less likely to find cilantro", "Unlikely soapy"),
    ("G/G: ~57%", "Likely (~57%)"), ("one A: ~71%", "Likely (~71%)"), ("A/A: ~74%", "Likely (~74%)"),
    ("Variants not detected", "Not detected"), ("G84E not detected", "Not detected"), ("Variant detected", "Detected"),
    ("G84E detected", "Detected"), ("One APOL1", "One variant"), ("Two APOL1", "Two variants"),
    ("Wet earwax", "Wet"), ("Dry earwax", "Dry"),
    ("Efficient omega-3", "Efficient"), ("Moderate omega-3", "Moderate"), ("Slow omega-3", "Slow"),
    ("No DQ2/DQ8", "Very unlikely"), ("Carries DQ", "DQ carrier"), ("Possible DQ2.2", "Inconclusive"),
    ("Lower than average —", "Lower"), ("About average (", "About average"), ("Somewhat raised", "Somewhat raised"),
    ("Raised (", "Raised"),
    ("677TT", "Reduced (677TT)"), ("677CT + 1298", "Mildly reduced"), ("1298CC", "Mildly reduced"),
    ("One reduced copy", "Typical"),
    ("Typical MTHFR", "Typical"),
    ("Typical vitamin D", "Typical"), ("Slightly lower vitamin D", "Slightly lower"), ("Lower vitamin D", "Lower"),
    ("Higher B12", "Higher"), ("Slightly lower B12", "Slightly lower"), ("Lower B12", "Lower"),
    ("Typical β-carotene", "Typical"), ("Somewhat reduced β-carotene", "Somewhat reduced"),
    ("Reduced β-carotene", "Reduced"),
    ("No HFE variants", "None"), ("C282Y/C282Y", "At risk"), ("C282Y/H63D", "Mild risk"),
    ("H63D/H63D", "H63D ×2"), ("C282Y carrier", "Carrier"), ("H63D carrier", "Carrier (H63D)"),
    ("No Factor V Leiden", "None"), ("Factor V Leiden", "Variant"), ("Prothrombin G20210A", "Variant"),
    ("PiMM", "Normal (PiMM)"), ("PiMS", "Carrier (PiMS)"), ("PiMZ", "Carrier (PiMZ)"), ("PiSS", "PiSS"),
    ("PiSZ", "At risk (PiSZ)"), ("PiZZ", "At risk (PiZZ)"),
]


def wellness_notes(have_23, have_bmi, have_sleep):
    """The "About these results" notes; a sentence about the 23andMe comparison, the BMI score or the sleep-paralysis
    row only when that is on the page. Population context comes from family.json "landing.*" text keys."""
    anc = text("landing.carrier_ancestry")  # e.g. "Northern European": the ancestry the carrier panel was chosen for
    ttr = text("landing.ttr_note")
    return {
        "23andMe wellness": ("🌿 Wellness traits (23andMe reports, recreated)",
                             "Same single markers 23andMe uses"
                             + ("; " + text("landing.wellness_check", "where a 23andMe report exists, the comparison below "
                                                                       "shows its label next to ours") if have_23 else "")
                             + "."
                             + (" Genetic Weight ≈ BMI score above." if have_bmi else "")
                             + (" Sleep Paralysis is a private 23andMe polygenic model (11,348 "
                                "variants, weights unpublished; genetics alone AUROC ~0.60, top 5% ≈ 2× average odds), so its "
                                "row shows 23andMe's own result where we have the report, not a recreation." if have_sleep
                                else "")),
        "Taste & diet": ("🍽️ Taste & diet",
                         "Allele directions checked against the original studies (TAS2R38 via Ensembl amino-acid changes; "
                         "cilantro Eriksson 2012; asparagus Eriksson 2010 (Table 9), replicated by Pelchat 2011; FADS1 omega-3 studies). Effects are real but "
                         "modest except bitter taste. Celiac uses Monsuur 2008 HLA tag SNPs: not carrying DQ2/DQ8 "
                         "makes celiac very unlikely (about 95–99% of people with celiac carry one; DQ7.5 isn't tested), while carrying them is common and only a small minority develop it. "
                         "No marker reliably predicts spicy-food tolerance (TRPV1 evidence is "
                         "weak) or seafood allergy, so neither is shown."),
        "Nutrition & health": ("🥗 Nutrition & health markers",
                               text("landing.nutrition_source", "Markers often found in third-party wellness reports that "
                                                                "hold up in the literature.")
                               + " All are tendencies; a blood test "
                               "(vitamin D, B12, homocysteine) tells you your actual level."),
        "Carrier status": ("🧬 Carrier status",
                           "Carrier panels cover a few variants the chip measures reliably, each verified against Ensembl "
                           f"and ClinVar: β-thalassemia ({text('landing.thal_alleles', '5 alleles')}), sickle cell (HbS) "
                           "and HbC; cystic fibrosis (F508del, the commonest, plus 5 rarer variants; the next most common "
                           + (f"{anc} ones" if anc else "ones") + ", G542X and N1303K, are not tested); GJB2 hearing loss "
                           "(35delG plus 4); Familial Mediterranean Fever (M694V, M680I, V726A); plus "
                           + (f"the 23andMe carrier conditions most relevant to {anc} ancestry" if anc else
                              "a selection of the 23andMe carrier conditions")
                           + ": PKU (13 PAH variants, incl. the Mediterranean IVS10-11G>A), "
                           "Gaucher (N370S), hereditary fructose intolerance (3 ALDOB), MCAD (K304E), Pompe (IVS1-13T>G), "
                           "PMM2-CDG (R141H, F119L), Pendred (4 SLC26A4) and pyruvate kinase deficiency (R486W). The other "
                           f"~30 23andMe carrier reports are founder variants of {'other' if anc else 'specific'} populations "
                           "(Ashkenazi, Finnish, French-Canadian, Amish…)"
                           + (" and are listed in the 23andMe comparison" if have_23 else "")
                           + ". Being a carrier mostly matters only if both "
                           "parents carry the same condition (1 in 4 chance per child). Exceptions: one MEFV copy occasionally "
                           "causes mild FMF symptoms, β-thalassemia carriers have small red cells (thalassemia minor), "
                           "which can be mistaken for iron deficiency, and sickle-cell trait carries rare risks with "
                           "extreme exertion or dehydration. 'Not a carrier' covers only these variants. "
                           "Lp(a): rs10455872 and rs3798220 (Clarke 2009) flag about 1 in 6 Europeans with likely high Lp(a), "
                           "but many people with high Lp(a) carry neither (it is mainly set by the kringle IV-2 copy number), "
                           "so a one-time blood test is recommended for everyone."),
        "Inherited conditions": ("🩺 Inherited conditions (common, well-measured variants)",
                                 "Allele directions verified with Ensembl. Carriers usually have no disease from one copy"
                                 + (f"; {text('landing.pis_note')}." if text("landing.pis_note") else ".")
                                 + " These are the main variants only — a negative "
                                 "result doesn't exclude rarer ones. HLA-B27 is inferred from tag SNP rs4349859 "
                                 "(98% sensitive, 99% specific in Europeans; Evans 2011), cross-checked with rs116488202; "
                                 "about 90% of people with ankylosing spondylitis carry B27, but about 95% of carriers "
                                 "never develop it. Parkinson's (LRRK2 G2019S, GBA N370S), hereditary amyloidosis (TTR V30M, "
                                 f"V122I, T60A{f'; {ttr}' if ttr else ''}), APOL1 kidney risk (G1 tag, G2) and "
                                 "HOXB13 G84E (prostate) are the exact variants 23andMe tests. Many carriers never develop "
                                 "the disease, though for HOXB13 G84E a substantial minority of men do (lifetime prostate "
                                 "cancer roughly 33–60%). Colon polyp genes = MUTYH Y179C and G396D (two copies cause "
                                 "MUTYH-associated polyposis; only 2 of >100 MUTYH variants, no Lynch/FAP); BRCA = the 3 "
                                 "Ashkenazi founder variants only, so a negative excludes almost nothing — "
                                 + (f"with {text('landing.brca_context')}, clinical testing (see the doctor pages) is what "
                                    "counts." if text("landing.brca_context") else
                                    "family history decides whether clinical testing (see the doctor pages) is warranted.")),
    }


SLEEP_PARALYSIS = "Sleep paralysis (23andMe score)"
MUTYH_ROW, BRCA_ROW = "Colon polyp genes (MUTYH)", "BRCA founder variants"


def hereditary_rows():
    """{row: {person: (short, full, tested)}} for the MUTYH and BRCA chip markers (hereditary_markers.py)."""
    by = {}
    for r in read_tsv("hereditary_markers.tsv"):
        if r["alt_copies"] != "":
            by.setdefault((r["person"], r["panel"]), []).append((r["marker"], int(r["alt_copies"])))
    out = {}
    for (p, panel), marks in by.items():
        if panel == "colorectal":
            mutyh = [(m, n) for m, n in marks if m.startswith("MUTYH")]
            n = sum(c for _, c in mutyh)
            short = "Not detected" if not n else "Carrier (one copy)" if n == 1 else "Two variants: see a doctor"
            full = ("MUTYH Y179C and G396D, the two variants 23andMe tests: "
                    + ("none detected" if not n else f"{n} variant cop{'y' if n == 1 else 'ies'}; two (one from each "
                       "parent) cause MUTYH-associated polyposis, one makes a carrier")
                    + ". Only 2 of >100 known MUTYH variants; Lynch syndrome and FAP (APC) are not tested")
            out.setdefault(MUTYH_ROW, {})[p] = (short, full, ", ".join(f"{m}={c}" for m, c in mutyh))
        else:
            n = sum(c for _, c in marks)
            full = ("BRCA1 185delAG, BRCA1 5382insC, BRCA2 6174delT (the 3 Ashkenazi founder variants 23andMe tests): "
                    + ("none detected" if not n else "variant detected — confirm with a clinical test")
                    + ". Thousands of other BRCA1/2 variants are not tested, so a negative excludes almost nothing; "
                      "family history decides whether clinical testing is warranted")
            out.setdefault(BRCA_ROW, {})[p] = ("Not detected" if not n else "Detected: confirm clinically", full,
                                              ", ".join(f"{m}={c}" for m, c in marks))
    return out




def archive():
    """23andMe's own displayed results (report_archive.py), keyed by person."""
    path = os.path.join(RESULTS, "23andme_reports.json")
    return json.load(open(path)) if os.path.exists(path) else {}


# Landing-page sections: (title, [(subgroup label, [report names])], short note, notes for the "About" fold).
# Reports are matched by name; anything in wellness.tsv not listed here lands in the last subgroup.
WELLNESS_LAYOUT = [
    ("🩺 Health & wellness", [
        ("Health", ["Alzheimer's (APOE)", "Macular degeneration", "Lipoprotein(a)", "Hemochromatosis (iron)",
                    "Inherited clotting", "Alpha-1 antitrypsin", "HLA-B27 (spine arthritis)", "Parkinson's genes",
                    "Hereditary amyloidosis", "Kidney disease (APOL1)", "Hereditary prostate cancer", MUTYH_ROW,
                    BRCA_ROW]),
        ("Carrier status (for family planning)", ["Beta-thalassemia & sickle cell", "Cystic fibrosis",
                                                  "Inherited hearing loss", "Familial Mediterranean Fever",
                                                  "Phenylketonuria (PKU)", "Gaucher disease",
                                                  "Hereditary fructose intolerance", "MCAD deficiency (fat metabolism)", "Pompe disease",
                                                  "PMM2-CDG (metabolic disorder)", "Pendred syndrome (hearing)",
                                                  "Pyruvate kinase deficiency"]),
        ("Body & sleep", ["Muscle composition", "Deep sleep", "Sleep movement", SLEEP_PARALYSIS, "Earwax type"]),
    ], "APOE is the largest single genetic factor for late-onset Alzheimer's, not a prediction. Carriers of one "
       "copy usually have no disease (β-thalassemia carriers have small red cells; one MEFV copy occasionally gives "
       "mild FMF). No LPA risk variant doesn't rule out high Lp(a): everyone should have it measured once. HLA-B27 "
       "matters mainly if there is chronic inflammatory back pain or eye inflammation (uveitis).",
     ["Inherited conditions", "Carrier status", "23andMe wellness"]),
    ("🥗 Nutrition & taste", [
        ("Digestion & diet", ["Lactose intolerance", "Celiac genes (HLA-DQ2/DQ8)", "Saturated fat & weight",
                              "Alcohol flush", "Caffeine consumption"]),
        ("Vitamins & nutrients", ["Folate / homocysteine", "Vitamin D", "Vitamin B12", "Vitamin A (from β-carotene)",
                                  "Omega-3 conversion"]),
        ("Taste & smell", ["Bitter taste", "Cilantro soapiness", "Asparagus smell"]),
    ], "All are tendencies; a blood test (vitamin D, B12, homocysteine) tells you your actual level. "
       "Not carrying DQ2/DQ8 makes celiac very unlikely (about 95–99% of people with celiac carry one).",
     ["Nutrition & health", "Taste & diet", "23andMe wellness"]),
]


def implied_from_parents(table):
    """Fill a child's missing call where mom and dad have the same homozygous genotypes (so the child must too)
    with their result, marked "(from parents)"."""
    mom, dad = person("mother"), person("father")
    for vals in table.values():
        if mom not in vals or dad not in vals or vals[mom][2] != vals[dad][2]:
            continue
        gts = [g.partition("=")[2] for g in vals[mom][2].split()]
        if not gts or not all(len(g) == 2 and g[0] == g[1] and g[0] in "ACGTDI" for g in gts):
            continue
        short, full_, geno = vals[mom]
        for c in children():
            if c in vals and vals[c][0] != "—":
                continue
            vals[c] = (short + " (from parents)", f"Implied: mom and dad are both {geno}, so their children are too. "
                       + full_, geno + " (both parents)")


def wellness_section():
    rows = read_tsv("wellness.tsv")
    if not rows:
        return ""
    people = order(dict.fromkeys(r["person"] for r in rows))
    gp = grandpa_cells()
    head = "".join(f"<th>{esc(col(p))}</th>" for p in people) + (GP_TH if gp else "")
    ncol = len(people) + 1 + bool(gp)

    def gp_cell(report):
        if not gp:
            return ""
        c = gp.get(report)
        if not c or c["cell"] == "?":
            return f'<td class=dim title="{esc(c["detail"] if c else "Not reconstructed", quote=True)}">?</td>'
        # "!" = fully known (e.g. "Not sensitive!"); never for carrier hits, where his other copy stays unknown
        known = c["cell"].endswith("!") and not c["cell"].startswith("Carrier")
        full = c["cell"].rstrip("!") if known else c["cell"].rstrip("!") + " / ?"
        return f'<td class=dim title="{esc(c["detail"], quote=True)}">{esc(full)}</td>'

    table = {}
    for r in rows:
        short = next((s for k, s in WELLNESS_SHORT if k in r["result"]), r["result"])
        table.setdefault(r["report"], {})[r["person"]] = (esc(short), r["result"], r["genotypes"])
    for p, rec in archive().items():  # 23andMe's own polygenic result: shown, not recreated
        sp = rec.get("sleep_paralysis")
        if sp and p in people:
            label = sp["category"].split(" likelihood")[0]
            when = f" by their {sp['by_age']}" if sp.get("by_age") else ""
            pct = f"{sp['percent']}%{when.replace(' their', '')}" if sp.get("percent") is not None else ""
            full = (f"23andMe: {sp['category']}" + (f" ({sp['percent']}% of people with genetics like this have had it"
                    f"{when})" if pct else "") + ". Private 23andMe Premium polygenic model; variants unpublished, so not recomputed")
            table.setdefault(SLEEP_PARALYSIS, {})[p] = (esc(label + (f" ({pct})" if pct else "")), full, "23andMe report")
    for rep, vals in hereditary_rows().items():
        for p, (short, full, tested) in vals.items():
            if p in people:
                table.setdefault(rep, {})[p] = (esc(short), full, tested)
    implied_from_parents(table)
    listed = {rep for _, groups, _, _ in WELLNESS_LAYOUT for _, reps in groups for rep in reps}
    WELLNESS_LAYOUT[-1][1][-1][1].extend(rep for rep in table if rep not in listed)

    def row(rep):
        if rep == "Alzheimer's (APOE)":
            return apoe_row(people, gp_cell(rep))
        if rep not in table:
            return ""
        vals = table[rep]
        return (f"<tr><th scope=row>{esc(rep)}</th>"
                + "".join(f'<td title="{esc(vals[p][1])} · {esc(vals[p][2])}">{vals[p][0]}</td>' if p in vals
                          else "<td>—</td>" for p in people) + gp_cell(rep) + "</tr>")

    have_23 = any(archive().get(p, {}).get("reports") for p in family_people())  # twentythree_section() renders
    have_bmi = any(pgs_trait(r) == "Body mass index" for p in pgs_people() for r in load_json(f"pgs_{p}.json").values())
    notes = wellness_notes(have_23, have_bmi, SLEEP_PARALYSIS in table)
    shown = set()  # a note shared by several sections goes under the first one only
    sections = []
    for title, groups, short_note, about in WELLNESS_LAYOUT:
        body = "".join(f'<tr class="grp"><th colspan="{ncol}">{esc(label)}</th></tr>' + "".join(row(r) for r in reps)
                       for label, reps in groups)
        sib = apoe_family_note(people) if "Health" in [g for g, _ in groups] else ""
        keys = [k for k in about if k in notes and k not in shown]
        shown.update(keys)
        fold = " ".join(notes[k][1] for k in keys)
        sections.append(
            f'<section class="card"><h2>{esc(title)}</h2><div class="tablewrap">'
            f'<table class="fam"><thead><tr><th></th>{head}</tr></thead><tbody>{body}</tbody></table></div>'
            f'<p class="dim">{esc(short_note + sib)} Hover a cell for the full wording and genotype.</p>'
            f'<details><summary>About these results</summary><p class="dim">{esc(fold)}</p></details></section>')
    return "".join(sections)


# 23andMe report → where our own call for the same person lives: ("w", wellness report), ("pgs", trait),
# ("pgx", gene), ("hered", marker prefix), ("apoe",), ("eye",), ("hair",), ("red",), or (None, why not recreated).
TWENTYTHREE_MAP = {
    "Age-Related Macular Degeneration": ("w", "Macular degeneration"),
    "Alpha-1 Antitrypsin Deficiency": ("w", "Alpha-1 antitrypsin"),
    "Type 2 Diabetes": ("pgs", "Type 2 diabetes"),
    "Celiac Disease": ("w", "Celiac genes (HLA-DQ2/DQ8)"),
    "Chronic Kidney Disease (APOL1-Related)": ("w", "Kidney disease (APOL1)"),
    "Colorectal Cancer": ("pgs", "Colorectal cancer"),
    "MUTYH-Associated Polyposis": ("hered", "MUTYH"),
    "Familial Hypercholesterolemia": (None, "23andMe's 24-variant list isn't public, so it can't be recreated"),
    "G6PD Deficiency": ("pgx", "G6PD"),
    "Hereditary Amyloidosis (TTR-Related)": ("w", "Hereditary amyloidosis"),
    "Hereditary Hemochromatosis (HFE-Related)": ("w", "Hemochromatosis (iron)"),
    "Hereditary Thrombophilia": ("w", "Inherited clotting"),
    "Late-Onset Alzheimer's Disease": ("apoe",),
    "Parkinson's Disease": ("w", "Parkinson's genes"),
    "Prostate Cancer": ("pgs", "Prostate cancer"),
    "BRCA1/BRCA2 (Selected Variants)": ("hered", "BRCA"),
    "Hereditary Prostate Cancer (HOXB13-Related)": ("w", "Hereditary prostate cancer"),
    "Alcohol Flush Reaction": ("w", "Alcohol flush"),
    "Caffeine Consumption": ("w", "Caffeine consumption"),
    "Deep Sleep": ("w", "Deep sleep"),
    "Genetic Weight": ("pgs", "Body mass index"),
    "Muscle Composition": ("w", "Muscle composition"),
    "Sleep Movement": ("w", "Sleep movement"),
    "Sleep Paralysis": (None, "Private 23andMe polygenic model (weights unpublished); shown in Body & sleep as 23andMe's result"),
    "Lactose Intolerance": ("w", "Lactose intolerance"),
    "Saturated Fat and Weight": ("w", "Saturated fat & weight"),
    "Asparagus Odor Detection": ("w", "Asparagus smell"),
    "Bitter Taste": ("w", "Bitter taste"),
    "Cilantro Taste Aversion": ("w", "Cilantro soapiness"),
    "Earwax Type": ("w", "Earwax type"),
    "Eye Color": ("eye",),
    "Light or Dark Hair": ("hair",),
    "Red Hair": ("red",),
    "Beta Thalassemia and Related Hemoglobinopathies": ("w", "Beta-thalassemia & sickle cell"),
    "Sickle Cell Anemia": ("w", "Beta-thalassemia & sickle cell"),
    "Cystic Fibrosis": ("w", "Cystic fibrosis"),
    "Familial Mediterranean Fever": ("w", "Familial Mediterranean Fever"),
    "Nonsyndromic Hearing Loss and Deafness, DFNB1 (GJB2-Related)": ("w", "Inherited hearing loss"),
    "Phenylketonuria and Related Disorders": ("w", "Phenylketonuria (PKU)"),
    "Gaucher Disease Type 1": ("w", "Gaucher disease"),
    "Hereditary Fructose Intolerance": ("w", "Hereditary fructose intolerance"),
    "MCAD Deficiency": ("w", "MCAD deficiency (fat metabolism)"),
    "Pompe Disease": ("w", "Pompe disease"),
    "Congenital Disorder of Glycosylation Type 1a (PMM2-CDG)": ("w", "PMM2-CDG (metabolic disorder)"),
    "Pendred Syndrome and DFNB4 Hearing Loss (SLC26A4-Related)": ("w", "Pendred syndrome (hearing)"),
    "Pyruvate Kinase Deficiency": ("w", "Pyruvate kinase deficiency"),
    "Maternal Haplogroup": ("mt",),
    "Paternal Haplogroup": ("y",),
}
ROW_NOTES = {  # why a recreated call can read differently from 23andMe's
    "Celiac Disease": "23andMe tests DQ2.5/DQ8 only; ours also checks DQ2.2",
    "Genetic Weight": "23andMe's private model vs our BMI score: different methods",
    "Type 2 Diabetes": "23andMe's private model vs our polygenic score: different methods",
    "Cilantro Taste Aversion": "23andMe's current model differs from the single OR6A2 marker we use (Eriksson 2012); "
                               "both are weak, so labels can disagree",
    "Light or Dark Hair": "23andMe's survey model vs HIrisPlex; compare with the family's actual hair colours in the "
                          "looks section",
    "Paternal Haplogroup": "23andMe names branches at different depths (e.g. R-M269 sits a few steps above R-L21)",
}
NOT_RECREATED = {  # section → why reports without a mapping aren't recreated
    "Traits": "23andMe's own survey-based models; no public, validated equivalent",
    "Carrier Status": "not recreated: mostly founder variants of other populations (Ashkenazi, Finnish, "
                      "French-Canadian, Amish…); 23andMe's result shown",
}


def our_call(person, spec):
    """Our own short result for one 23andMe report, or None."""
    kind = spec[0]
    if kind == "w":
        r = next((r for r in read_tsv("wellness.tsv") if r["person"] == person and r["report"] == spec[1]), None)
        return r and next((v for k, v in WELLNESS_SHORT if k in r["result"]), r["result"])
    if kind == "pgs":
        path = os.path.join(RESULTS, f"pgs_{person}.json")
        if os.path.exists(path):
            r = next((r for r in json.load(open(path)).values() if r["trait"].startswith(spec[1])), None)
            return r and f"{round(r['groups']['EUR']['percentile'])}th percentile (our score)"
        return None
    if kind == "pgx":
        r = next((r for r in read_tsv("pharmacogenomics.tsv") if r["person"] == person and r["gene"] == spec[1]), None)
        return r and r["phenotype"]
    if kind == "hered":
        rs = [r for r in read_tsv("hereditary_markers.tsv") if r["person"] == person and r["marker"].startswith(spec[1])]
        called = [r for r in rs if r["alt_copies"] != ""]  # "" = marker not on this person's chip
        if not called:
            return rs and "Not tested"
        if any(float(r["alt_copies"]) > 0 for r in called):
            return "Detected"
        return "Not detected" + (f" ({len(rs) - len(called)} of {len(rs)} markers not tested)" if len(called) < len(rs) else "")
    if kind == "apoe":
        r = next((r for r in read_tsv("apoe.tsv") if r["person"] == person), None)
        return r and r["apoe"] and (r["apoe"] + (" (no ε4)" if "ε4" not in r["apoe"] else ""))
    if kind in ("mt", "y"):
        r = next((r for r in read_tsv("mt_haplogroups.tsv" if kind == "mt" else "y_haplogroups.tsv")
                  if r.get("SampleID", r.get("person")) == person), None)
        return r and (r.get("Haplogroup") or r.get("haplogroup"))
    pig = (load_json("pigmentation.json") or {}).get("people", {}).get(person, {})
    if kind == "eye":
        return pig.get("eye_call") and f"{pig['eye_call'].capitalize()} (HIrisPlex)"
    if kind == "hair":
        return pig.get("hair_call") and f"{pig['hair_call'].capitalize()} (HIrisPlex)"
    if kind == "red":
        red = pig.get("hair", {}).get("red")
        if red is None:
            return None
        return f"Red {'likely' if red >= 0.5 else 'unlikely'} ({red:.0%}, HIrisPlex)"
    return None


TWENTYTHREE_SECTIONS = ["Health Predisposition", "Carrier Status", "Wellness", "Traits", "Ancestry"]


def twentythree_section():
    arch = archive()
    people = order(p for p in family_people() if arch.get(p, {}).get("reports"))
    if not people:
        return ""
    head = "".join(f"<th>{esc(col(p))}</th>" for p in people)
    parts = []
    for section in TWENTYTHREE_SECTIONS:
        names = list(dict.fromkeys(n for p in people for n in arch[p]["reports"].get(section, {})
                                   if not (section == "Ancestry" and arch[p]["reports"][section][n].endswith("%"))))
        if not names:
            continue
        rows = []
        for name in names:
            spec = TWENTYTHREE_MAP.get(name)
            cells = []
            for p in people:
                says = arch[p]["reports"].get(section, {}).get(name)
                ours = our_call(p, spec) if spec and spec[0] else None
                cells.append(f"<td>{esc(says) if says else '<span class=dim>—</span>'}"
                             + (f'<div class="dim">ours: {esc(ours)}</div>' if ours else "") + "</td>")
            note = ROW_NOTES.get(name) or (spec[1] if spec and not spec[0] else "") or (
                "" if spec else NOT_RECREATED.get(section, ""))
            rows.append(f"<tr><th scope=row>{esc(name)}</th>{''.join(cells)}<td class=dim>{esc(note)}</td></tr>")
        parts.append(f'<details><summary>{esc(section)} ({len(names)})</summary><div class="tablewrap"><table class="fam">'
                     f'<thead><tr><th></th>{head}<th>Note</th></tr></thead><tbody>{"".join(rows)}</tbody></table>'
                     f'</div></details>')
    comp = []
    for p in people:
        tree = arch[p].get("ancestry_tree") or []
        top = [(n, v, l) for n, v, l in tree if l == 1]
        subs = {}
        cur = None
        for n, v, l in tree:
            if l == 1:
                cur = n
            elif l == 2 and cur:
                subs.setdefault(cur, []).append(f"{n} {v}")
        line = " · ".join(f"{n} {v}" + (f" ({', '.join(subs[n])})" if subs.get(n) else "") for n, v, _ in top)
        if line:
            comp.append(f"<li><b>{esc(col(p))}:</b> {esc(line)}</li>")
    if comp:
        parts.append('<details><summary>Ancestry Composition (23andMe)</summary><ul class="dim">' + "".join(comp) +
                     "</ul><p class=\"dim\">23andMe's own reference panels and labels; our fine-scale ancestry "
                     "(above) uses a different panel and method, so the two are not expected to match exactly.</p></details>")
    raw = [p for p in people if p not in {r["person"] for r in read_tsv("wellness.tsv")}]
    note = (f" {', '.join(full(p) for p in raw)}: only 23andMe's labels so far (raw DNA file not here yet)." if raw else "")
    return ('<section class="card"><h2>📋 23andMe\'s own reports vs ours</h2>'
            f'<p class="dim">23andMe\'s label for each person, with our own call from the raw data underneath ("ours:") '
            f'where we recreate it.{esc(note)} Accounts without Health/Carrier reports show "—" (not on their plan).</p>'
            + "".join(parts) +
            '<p class="dim">23andMe sometimes tests fewer variants (e.g. celiac: DQ2.5/DQ8 only) or uses its own '
            'survey-trained models, so wording can differ while the underlying call agrees.</p></section>')


def gp_inputs():
    """(children whose painting was used, whether 23andMe's one-parent phasing of mom was used) from grandpa.json."""
    inputs = (load_json("grandpa.json") or {}).get("inputs", {})
    return order(inputs.get("children_used", [SELF])), bool(inputs.get("phase_file"))


def gp_note():
    kids, phased = gp_inputs()
    src = ["mom", "grandma"] + (["mom's children's inheritance painting"] if kids else []) + (
        ["23andMe phasing of mom"] if phased else [])
    return ('<p class="dim" style="margin:12px 0 0">* Grandpa (mom\'s father) is a virtual reconstruction from '
            f'{esc(and_join(src), quote=False)}. We only '
            'know the copy he passed to mom, so each cell shows that copy and "?" for his other, unknown copy; "—" where a result needs both copies.</p>')


# Family-tree origins of the children's 16 great-great-grandparents (Ahnentafel numbers, GEDCOM via scripts/pedigree.py),
# from family.json "ancestry" (optional; the tree part of the ancestry card is hidden without it):
#   "regions": [region, ...] in display order (also the bar and map colour order),
#   "founders": {"16": {"region", "label", "assumed"}, ...}; "assumed" = birthplace missing in the GEDCOM, region
#               inferred from the line (same as the ancestry simulations),
#   "places": [{"lat", "lon", "place", "surnames", "region", "dx", "dy", "anchor"}, ...] documented birthplaces with
#             label placement; "*" after a surname = birthplace inferred from the line, not recorded.
# text "landing.ancestry_groups": [[label, [regions]], ...] optional summary for the lede (e.g. a country made of
# several regions); without it each region is its own group.
ANCESTRY = FAMILY.get("ancestry", {})
REGIONS = ANCESTRY.get("regions", [])
FOUNDERS = {int(k): (v["region"], v["label"], bool(v.get("assumed"))) for k, v in ANCESTRY.get("founders", {}).items()}
PLACES = [(p["lat"], p["lon"], p["place"], p["surnames"], p["region"], p.get("dx", 9), p.get("dy", -4),
           p.get("anchor", "start")) for p in ANCESTRY.get("places", [])]
# Whose great-great-grandparent slots each role covers (dad's are 16-23, mom's 24-31, ...).
ROLE_SLOTS = {"self": range(16, 32), "sibling": range(16, 32), "mother": range(24, 32), "father": range(16, 24),
              "maternal_grandmother": range(28, 32), "grandpa": range(24, 28)}
GRANDPA = "grandpa-maternal (inferred)"  # how ancestry.json / ancestry_ho.json name the virtual grandfather


def slots(p):
    return [i for i in ROLE_SLOTS.get("grandpa" if p == "grandpa" else role(p), []) if i in FOUNDERS]


def sw(region):
    """Colour class for a region (six colours, reused beyond six regions)."""
    return f"s{REGIONS.index(region) % 6 + 1}" if region in REGIONS else "s6"


def region_shares(p):
    known = slots(p)
    return {r: sum(FOUNDERS[i][0] == r for i in known) / len(known) if known else 0 for r in REGIONS}


def origins_bars(people):
    rows = []
    for p in people:
        shares = region_shares(p)
        segs = "".join(
            f'<span class="seg {sw(r)}" style="flex-basis:{v * 100:.2f}%" title="{esc(r)}: {v * 100:.0f}%">'
            f'{f"{v * 100:.0f}%" if v >= 0.12 else ""}</span>'
            for r, v in shares.items() if v)
        who = "Grandpa*" if p == "grandpa" else col(p)
        rows.append(f'<div class="obar"><span class="who">{esc(who)}</span><span class="segs">{segs}</span></div>')
    legend = "".join(f'<span><i class="sw {sw(r)}"></i>{esc(r)}</span>' for r in REGIONS)
    return f'<div class="olegend">{legend}</div>{"".join(rows)}'


def origins_map():
    base = load_json_path(os.path.join(WORK, "map", "basemap.json"))
    if not base or not PLACES:
        return ""
    pr = base["project"]

    def xy(lat, lon):
        return (lon - pr["lon0"]) * pr["cos"] * pr["k"], (pr["lat1"] - lat) * pr["k"]
    dots, labels = [], []
    for lat, lon, place, names, region, dx, dy, anchor in PLACES:
        x, y = xy(lat, lon)
        n = len(names.split(", "))
        dots.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{4 + 1.5 * n:.1f}" class="dot {sw(region)}">'
                    f'<title>{esc(place)}: {esc(names)}</title></circle>')
        labels.append(f'<text x="{x + dx:.1f}" y="{y + dy:.1f}" text-anchor="{anchor}">{esc(place)}</text>')
    hi = "".join(f'<path class="land hi" d="{d}"/>' for k, d in base["paths"].items() if k != "other")
    return (f'<svg class="omap" viewBox="{base["viewbox"]}" role="img" aria-label="Map of documented ancestral '
            f'birthplaces (* = inferred)"><path class="land" d="{base["paths"]["other"]}"/>'
            f'{hi}{"".join(dots)}<g class="olabels">{"".join(labels)}</g></svg>')


def dna_verdicts(people):
    """Short per-person cells derived from the verified fine-scale statements (ancestry_ho.json 'display')."""
    disp = (load_json("ancestry_ho.json") or {}).get("display", {}).get("people", {})
    adm = (load_json("ancestry.json") or {}).get("people", {})

    def said(p):
        v = disp.get(p, {})
        return (v.get("headline", "") + " " + " ".join(s["text"] for s in v.get("statements", []))).lower()

    def axis(p):
        t = said(p)
        return ("Mostly Spanish" if "mostly spanish" in t else "In between" if "mix of italian and spanish" in t
                else "Italian" if "clearly italian" in t or "italian (clearly" in t else "—")

    def italy(p):
        t = said(p)
        return ("Central/Northern" if "central/northern italians" in t else
                "Southern-type" if "southern-italian-type" in t else "—")

    def sicily(p):
        return "Not typical" if "not a typical sicilian" in said(p) or "not as a typical sicilian" in said(p) else "—"

    def balkan(p):
        return "Yes" if "balkan-like pull" in said(p) else "—"

    def tree(p):
        return "Fits ✓" if "consistent with" in said(p) else "—"

    def noneu(p):
        b = adm.get(p, {}).get("broad", {}).get("European")
        return "—" if b is None else "None" if b >= 0.95 else f"{(1 - b) * 100:.0f}% (likely noise)"
    keys = [p if p != "grandpa" else GRANDPA for p in people]
    # Non-European components that stay under 1% for everyone shown (computed, for the row's hover text).
    broad = [adm[k]["broad"] for k in keys if adm.get(k, {}).get("broad")]
    low = [c for c in dict.fromkeys(c for b in broad for c in b) if c != "European" and all(b.get(c, 0) < 0.01 for b in broad)]
    rows = [("Italian ↔ Spanish", axis, "Tuscan-vs-Spanish comparison; reference people placed correctly 99–100% of the time"),
            ("Within Italy", italy, "Southern-vs-Northern Italian comparison (88% accurate)"),
            ("Typical Sicilian?", sicily, "Four Sicily-vs-elsewhere comparisons (~90% accurate); partial Sicilian ancestry is still possible"),
            ("Balkan / E-Med pull", balkan, text("landing.balkan_tip", "Shifted toward Greeks and Albanians compared "
                                                                       "with Spaniards; see the evidence below")),
            ("Fits the family tree", tree, "Simulated families built from the GEDCOM tree"),
            ("Non-European", noneu, "Admixture vs 1000 Genomes" + (f"; {'/'.join(low)} all under 1%" if low else ""))]
    rows = [r for r in rows if any(r[1](k) != "—" for k in keys)]  # rows that say nothing for anyone are dropped
    body = "".join(f'<tr><th scope=row title="{esc(tip, quote=True)}">{esc(label)}</th>'
                   + "".join(f'<td{" class=dim" if p == "grandpa" else ""}>{esc(fn(k))}</td>'
                             for p, k in zip(people, keys)) + "</tr>"
                   for label, fn, tip in rows)
    if not body:
        return ""
    head = "".join(f"<th>{esc(col(p))}</th>" if p != "grandpa" else GP_TH for p in people)
    return (f'<div class="tablewrap"><table class="fam"><thead><tr><th></th>{head}</tr></thead>'
            f'<tbody>{body}</tbody></table></div>')


def ancestry_section():
    adm = (load_json("ancestry.json") or {}).get("people", {})
    # The viewer is the only child shown: the other children's tree is the same and their DNA verdicts are on their pages.
    people = order(p for p in family_people() if p in adm and (p == VIEWER or role(p) not in ("self", "sibling")))
    if person("mother") and person("maternal_grandmother"):
        people.append("grandpa")
    if len(people) < 2:
        return ""
    verdicts, fine = dna_verdicts(people), fine_scale_block()
    tree = ""
    if FOUNDERS and REGIONS:
        me = region_shares(VIEWER)
        groups = text("landing.ancestry_groups") or [[r, [r]] for r in REGIONS]
        shares = [(g, sum(me.get(r, 0) for r in rs)) for g, rs in groups]
        lede = and_join([f"<b>~{v * 100:.0f}% {esc(g)}</b>" for g, v in shares if v])
        where = text("landing.ancestry_note")
        dna = text("landing.ancestry_dna") or tree_fit_sentence()
        names = [FOUNDERS[i][1].split(" (")[0].replace("'s wife", "") for i in sorted(FOUNDERS) if FOUNDERS[i][2]]
        assumed = [f"{n} ×{names.count(n)}" if names.count(n) > 1 else n for n in dict.fromkeys(names)]
        n_assumed = {1: "One", 2: "Two", 3: "Three", 4: "Four", 5: "Five", 6: "Six", 7: "Seven", 8: "Eight"}.get(
            len(names), str(len(names)))
        missing = 16 - len(FOUNDERS)
        how = (f'{n_assumed} slot{"s have" if len(names) != 1 else " has"} no recorded birthplace and '
               f'{"are" if len(names) != 1 else "is"} inferred from their line ({", ".join(assumed)}; those shown on '
               f'the map are marked *). ' if names else "")
        tree = (f'<p>By the family tree I am {lede}{f" ({esc(where)})" if where else ""}. {esc(dna)}</p>'
                f'<div class="anc-grid"><div>{origins_map()}</div>'
                f'<div><h3 style="margin-top:0">By family tree (great-great-grandparents)</h3>'
                f'{origins_bars([p for p in people if slots(p)])}'
                f'<p class="dim">From the GEDCOM. {how}'
                + (f'{missing} of the 16 slots are unknown and left out. ' if missing > 0 else "")
                + 'Hover for exact shares.</p></div></div>')
    if not tree and not verdicts and not fine:
        return ""
    return (f'<section class="card"><h2>🌍 Ancestry</h2>{tree}'
            + (f'<h3>What the DNA says</h3>{verdicts}'
               f'<p class="dim">Hover a row name for how it was tested.{esc(untested_note(), quote=False)}</p>' if verdicts else "")
            + (f'<details><summary>Family lines in detail</summary>{genealogy_note()}</details>' if genealogy_note() else "")
            + (f'<details><summary>Evidence and limits of the DNA comparison</summary>{fine}</details>' if fine else "")
            + '</section>')


def ho_pops(d, *keys):
    """Population codes used by ancestry_ho.json, from the first of its sections (pairwise pairs, calibration) that has
    any, in first-seen order, spelled for display ("Italian_North" -> "Italian North")."""
    for k in keys:
        sec = d.get(k, {})
        codes = list(dict.fromkeys(c for key in (sec.get("pairs", {}) if k == "pairwise" else sec) for c in key.split("|")))
        if codes:
            return [c.replace("_", " ") for c in codes]
    return []


def untested_note():
    """Why admixture percentages aren't shown: ancestry_ho.py tested them on reference people of known origin
    (its leave-one-out calibration); text "landing.ho_tested_on" names those people in the family's own words."""
    pops = text("landing.ho_tested_on") or and_join(ho_pops(load_json("ancestry_ho.json") or {}, "calibration"))
    return (f" Country-level percentages from DNA aren't shown: they were tested on real {pops} and came out unreliable."
            if pops else "")


def tree_fit_sentence():
    """Whose DNA fits the documented family tree, from the ancestry_ho.json headlines ("" without them)."""
    disp = (load_json("ancestry_ho.json") or {}).get("display", {}).get("people", {})

    def said(p):
        v = disp.get(p, {})
        return v.get("headline", "") + " " + " ".join(s["text"] for s in v.get("statements", []))
    ids = order(p for p in disp if p in family_people())
    fits = [rel(p) for p in ids if "consistent with the documented family tree" in said(p)]
    off = [rel(p) for p in ids if "does not match the documented family tree" in said(p)]
    out = f"The DNA is consistent with the family tree for {and_join(fits)}." if fits else ""
    if off:
        out += (" It does not" if fits else "The DNA does not") + f" match the tree well for {and_join(off)}."
    return out


def fine_scale_block():
    """Verified statements from ancestry_ho.py (AADR Human Origins panel). Admixture percentages are biased there,
    so only the leave-one-out-calibrated statements are shown; evidence goes in the hover text."""
    d = (load_json("ancestry_ho.json") or {}).get("display")
    if not d:
        return ""
    items = []
    for p, v in d["people"].items():
        stmts = "".join(f'<li title="{esc(s["evidence"])}">{esc(s["text"])}</li>' for s in v["statements"])
        items.append(f'<li><b>{esc(v["label"])}:</b> {esc(v["headline"])}<ul>{stmts}</ul></li>')
    summary = text("landing.fine_scale_summary")  # the family's own reading of the statements, if any
    pops = text("landing.ho_panel") or ", ".join(ho_pops(load_json("ancestry_ho.json"), "pairwise", "calibration"))
    return (f'<p class="dim">Reference panel: AADR Human Origins{f" ({esc(pops)})" if pops else ""}. '
            'Hover a statement for its evidence.</p>'
            f'<ul>{"".join(items)}</ul>'
            + (f'<p class="dim">{esc(summary)}</p>' if summary else ""))


def genealogy_note():
    notes = (load_raw_json("family_notes.json") or {}).get("genealogy")
    if not notes:
        return ""
    return ((f'<p><b>Mom\'s side:</b> {esc(notes["maternal"])}</p>' if "maternal" in notes else "")
            + (f'<p><b>Dad\'s side:</b> {esc(notes["paternal"])}</p><p>{esc(notes.get("expected_me", ""))}</p>'
               if "paternal" in notes else ""))


def neanderthal_section():
    """Categories only: per-person percentiles are too uncertain to show (few segments, IBS-ascertained calls)."""
    d = load_json("neanderthal.json")
    if not d:
        return ""

    def category(z):
        return ("More than average" if z >= 1.5 else "Less than average" if z <= -1.5 else "Typical")
    gp = bool(grandpa_cells())
    people = order(d["people"])
    # Both comparisons when available: IBS is in-sample (Sprime found the segments in those Spaniards, so we look low
    # against them) and EUR_noIBS may run high for Southern Europeans; the truth is likely in between.
    both = all("EUR_noIBS" in d["people"][p]["groups"] for p in people)
    refs = [("IBS", "Spaniards"), ("EUR_noIBS", "Europeans outside Spain")] if both else [("IBS", "Spaniards")]

    def cell(p):
        cats = {lab: category(d["people"][p]["groups"][ref]["z"]) for ref, lab in refs}
        uniq = list(dict.fromkeys(cats.values()))
        text = uniq[0] + "".join(" or " + c[0].lower() + c[1:] for c in uniq[1:])
        tip = "; ".join(f"vs {lab}: {c}" for lab, c in cats.items())
        return f'<td title="{esc(tip, quote=True)}">{esc(text)}</td>'
    head = "".join(f"<th>{esc(col(p))}</th>" for p in people) + (GP_TH if gp else "")
    cells = ("<tr><th scope=row>Neanderthal segments</th>"
             + "".join(cell(p) for p in people)
             + (GP_NA.format("Needs both copies") if gp else "") + "</tr>")
    return (f'<section class="card"><h2>🦴 Neanderthal DNA</h2><div class="tablewrap"><table class="fam">'
            f'<thead><tr><th></th>{head}</tr></thead><tbody>{cells}</tbody></table></div>'
            f'<p class="dim">Relative amount of the Neanderthal segments found in Spaniards (Sprime calls, '
            f'Browning et al. 2018), compared with '
            + ('two groups: 1000 Genomes Spaniards (IBS) and 1000 Genomes Europeans outside Spain (British, Utah '
               'CEPH, Tuscans, Finns). The segments were discovered in those same Spaniards, who therefore score '
               'slightly high; Europeans outside Spain avoid that bias, but because the segments were found in '
               'Spaniards, Southern Europeans may score a little high against them. The truth is likely in '
               'between, so a cell names both results when the two comparisons disagree (hover for each). '
               if both else
               '1000 Genomes Spaniards (IBS). The segments were discovered in those same Spaniards, so they score '
               'slightly high and a result against them may read a little low. ')
            + 'All Europeans carry roughly 2% Neanderthal DNA; this measures who has a bit more or less than the '
            'comparison group, not a total %. Exact ranks are too uncertain to show.</p></section>')


def pigmentation_section():
    rows = read_tsv("pigmentation.tsv")
    if not rows:
        return ""
    rows = [r for p in order(dict.fromkeys(r["person"] for r in rows)) for r in rows if r["person"] == p]
    gp = grandpa_cells().get("_pigmentation")
    head = "".join(f"<th>{esc(col(r['person']))}</th>" for r in rows) + (GP_TH if gp else "")

    def eye(r):
        p = {"blue": r["p_blue"], "intermediate": r["p_intermediate"], "brown": r["p_brown"]}.get(r["eye_call"])
        if not p:  # too many eye SNPs missing to predict
            return f'<span class=dim>{esc(r["eye_call"].capitalize())}</span>'
        return f"{r['eye_call'].capitalize()} <span class=dim>{float(p) * 100:.0f}%</span>"

    def red(r):
        strong, weak = round(float(r["mc1r_R"])), round(float(r["mc1r_r"] or 0))  # round(), as on the doctor pages
        if strong >= 2:
            return "Likely"
        return ("Red tints possible (R + r)" if weak else "Tints possible") if strong == 1 else "Unlikely"
    observed = (load_raw_json("family_notes.json") or {}).get("observed", {})
    pig = (load_json("pigmentation.json") or {}).get("people", {})
    imputed = ["my" if p == VIEWER else f"{col(p)}'s"
               for p, v in pig.items() if set(v.get("imputed_or_phased", [])) & set(v.get("hair_snps_used", []))]
    who = esc(" and ".join(imputed))
    imputed_txt = (f" {who[:1].upper() + who[1:]} hair call{'s use' if len(imputed) > 1 else ' uses'} some imputed or "
                   "phased SNPs (statistical, not measured)." if imputed else "")
    lines = [("Eye color (predicted)", eye)]
    if observed:
        lines.append(("Eye color (actual)", lambda r: esc(as_viewer(observed.get(r["person"], {}).get("eye", "—")))))
    lines.append(("Hair color (predicted)", lambda r: r["hair_call"].capitalize()))
    if any("hair" in o for o in observed.values()):
        lines.append(("Hair color (actual)", lambda r: esc(as_viewer(observed.get(r["person"], {}).get("hair", "—")))))
    lines.append(("Red hair (MC1R)", red))
    body = "".join(f"<tr><th scope=row>{label}</th>" + "".join(f"<td>{fn(r)}</td>" for r in rows)
                   + (gp_td(gp.get(label)) if gp else "") + "</tr>" for label, fn in lines)
    return (f'<section class="card"><h2>🎨 Eye &amp; hair color (HIrisPlex)</h2><div class="tablewrap">'
            f'<table class="fam"><thead><tr><th></th>{head}</tr></thead><tbody>{body}</tbody></table></div>'
            f'<p class="dim">Published IrisPlex/HIrisPlex models (Walsh 2011/2013, coefficients from the paper\'s own '
            f'tool), independently recomputed and orientation-checked. Blue vs brown eyes is reliable; hair AUC is '
            f'≈0.72–0.93 (red and black predicted best, brown and blond weakest; Walsh 2013/2014) and green/hazel '
            f'eyes can\'t be predicted. A predicted hair shade that misses the reported one (e.g. blond vs brown) '
            f'sits on the model\'s weakest boundary, and reported shades are themselves uncertain.{imputed_txt}</p>'
            f'{red_hair_note()}</section>')


def red_hair_note():
    """Explains a red-haired maternal aunt (raw/family_notes.json "relatives") from the traced MC1R alleles
    (pigmentation.json + grandpa_cells.json)."""
    aunt = (load_raw_json("family_notes.json") or {}).get("relatives", {}).get("aunt_maternal")
    if not aunt or "red" not in json.dumps(aunt).lower():
        return ""
    pig = (load_json("pigmentation.json") or {}).get("people", {})
    gp = grandpa_cells().get("_pigmentation", {}).get("Red hair (MC1R)", {})

    def strong(p):  # strong R variants a person carries at least one copy of (imputed dosages rounded)
        v = pig.get(p, {}).get("mc1r", {}).get("variants", {})
        return {n: round(v[n]) for n in MC1R_STRONG.values() if round(v.get(n, 0)) >= 1}
    gma, mom, me, dad = (strong(person("maternal_grandmother")), strong(person("mother")), strong(VIEWER),
                         strong(person("father")))
    # Mom's single R copy came from grandma when grandma carries it and grandpa's transmitted copy shows none.
    # Only chip-measured sites in mom: grandpa's cell covers just the sites reconstructable from mom's chip.
    rs_of = {n: rs for rs, n in MC1R_STRONG.items()}
    src = pig.get(person("mother"), {}).get("genotype_sources", {})
    traced = ([n for n in mom if mom[n] == 1 and n in gma and src.get(rs_of[n]) == "chip"]
              if gp.get("cell", "").startswith("None found") else [])
    general = ('<p><b>The red-haired aunt:</b> red hair usually needs an MC1R variant on both copies (two strong R '
               'variants, or an R plus a weaker r), and about 6% of redheads carry only one.')
    if not traced:
        return general + '</p>'
    v = traced[0]
    mine = " I carry grandma's " + v + " too." if v in me and v not in dad else ""
    return (general + f' Mom\'s {v} came from grandma, and no strong red variant was found on the copy grandpa passed '
            'to mom (only the sites reconstructable from mom\'s chip; hover his cell for coverage). Assuming the aunt '
            f'is a full sister, she most likely got grandma\'s {v} plus an MC1R variant (strong R or weaker r) on '
            'grandpa\'s <i>other</i> copy, so grandpa was probably an MC1R variant carrier, unless the aunt got a '
            f'variant on grandma\'s other copy or is a single-variant redhead.{mine}</p>')


def grandpa_section():
    data = load_json("grandpa.json")
    rows = read_tsv("grandpa.tsv")
    if not data or not rows:
        return ""
    cov = data["coverage"]
    kids, phase = gp_inputs()
    how = ([f"mom's children's inheritance painting ({', '.join(rel(k) for k in kids)})"] if kids else []) + (
        ["23andMe's statistical phasing of mom"] if phase else [])
    auto = cov["autosomes"]
    phased = auto.get("by_method", {}).get("phased", 0) / auto["mom_called_snps"] * 100 if auto["mom_called_snps"] else 0
    agree = data.get("validation", {}).get("phased_vs_children_hethet_autosomes", {}).get("agree_pct")
    err = f", est. error at most ~{100 - agree:.1f}% at those sites" if agree is not None else ""
    body = "".join(
        f"<tr><th scope=row>{esc(r['gene'])}</th><td class=nw><b>{esc(r['grandpa_transmitted_allele'])}</b></td>"
        f"<td>{esc(r['meaning'])}{' <span class=dim>(imputed)</span>' if r.get('source', '').startswith('imputed') else ''}"
        f"</td></tr>" for r in rows)
    return ('<section class="card"><h2>👴 Virtual grandpa (mom\'s father)</h2>'
            '<p class="dim">Reconstructed from mom + grandma'
            + (f', with ambiguous sites resolved by {esc(" or ".join(how), quote=False)}' if how else "")
            + f'; rows marked "imputed" come from imputed data. His transmitted copy is '
            f'called at {auto["known_pct"]:.1f}% of autosomal SNPs'
            + (f' ({phased:.1f}% rest on statistical phasing alone{err})' if "phased" in auto.get("by_method", {}) else "")
            + ' '
            f'and {cov["X_nonPAR"]["known_pct"]:.1f}% of his X. '
            f'Each row is one of his two copies — he carried at least this.</p>'
            f'<div class="tablewrap"><table class="fam wrap"><tbody>{body}</tbody></table></div></section>')


def apoe_family_note(people):
    """Grandma→mom ε4 transmission and possible genotypes of siblings without data, derived from apoe.tsv."""
    calls = {r["person"]: r for r in read_tsv("apoe.tsv") if r["apoe"]}
    gma, mom, dad = (calls.get(person("maternal_grandmother")), calls.get(person("mother")),
                     calls.get(person("father")))
    imputed = [p for p, r in (("mom", mom), ("dad", dad)) if r and "imputed" in r["source"]]
    why = "mom's and dad's APOE are imputed" if len(imputed) == 2 else f"{''.join(imputed)}'s APOE is imputed"
    out = ""
    if gma and mom and "ε4" in gma["apoe"] and "ε4" not in mom["apoe"]:
        out += " Grandma's ε4 was not passed to mom" + (" (mom's APOE is imputed)" if "mom" in imputed else "") + "."
    missing = [c for c in children() if c not in people]
    if missing and mom and dad:
        kids = sorted({"/".join(sorted((a, b))) for a in mom["apoe"].split("/") for b in dad["apoe"].split("/")})
        extra = [] if any("ε4" in k for k in kids) else [f"no ε4 {'expected' if imputed else 'possible'}"]
        extra += [why] if imputed else []
        who = and_join([full(c) if i == 0 else label(c, VIEWER).lower() + tag(c) for i, c in enumerate(missing)])
        out += (f" {who} {'would most likely be' if imputed else 'can only be'} {' or '.join(kids)}"
                + (f" ({'; '.join(extra)})" if extra else "") + ".")
    return out


def apoe_row(people, extra=""):
    """Alzheimer's (APOE) as the first row of the health table."""
    rows = {r["person"]: r for r in read_tsv("apoe.tsv")}
    if not rows:
        return ""

    def td(p):
        r = rows.get(p)
        if not r:
            return "<td>—</td>"
        tip = esc(f"{r['risk']} · {r['source']}", quote=True)
        # Odds vs ε3/ε3 (Farrer 1997), by genotype so ε2/ε4 and ε4/ε4 don't get the one-copy figure.
        odds = {"ε2/ε4": "~2.5×", "ε3/ε4": "~3×", "ε4/ε4": "~12–15×"}.get(r["apoe"])
        flag = f' <span class=hi>↑ {odds} odds vs ε3/ε3</span>' if odds else ""
        return f'<td title="{tip}">{esc(r["apoe"] or "—")}{flag}</td>'
    return "<tr><th scope=row>Alzheimer's (APOE)</th>" + "".join(td(p) for p in people) + extra + "</tr>"


# Short "most important" medication pointers, in this order, for whoever is actionable for the gene.
MEDS_PRIORITY = [("SLCO1B1", "SLCO1B1 if a statin is ever prescribed ({})"), ("CYP2C19", "CYP2C19 (clopidogrel) for {}")]
MEDS_FILTER = {"CYP2C19": lambda ph: "Intermediate" in ph or "Poor" in ph}  # clopidogrel needs reduced function


def meds_note(pgx):
    """The medications card's note, from who is actionable for which gene (pharmacogenomics.tsv)."""
    act = {}
    for r in pgx:
        if note(r["gene"], r["phenotype"]) and MEDS_FILTER.get(r["gene"], lambda ph: True)(r["phenotype"]):
            act.setdefault(r["gene"], []).append(r["person"])
    everyone = {r["person"] for r in pgx}

    def who(g):
        ps = order(dict.fromkeys(act[g]))
        return "everyone" if set(ps) == everyone and len(ps) > 1 else and_join([rel(p) for p in ps])
    top = [tpl.format(who(g)) for g, tpl in MEDS_PRIORITY if g in act]
    # The tag SNP is low-specificity in Europeans (ancestry.json broad European share >= 0.9) but tracks the allele
    # well in some other groups, so "false positive" is only said when everyone with the hit is European.
    adm = (load_json("ancestry.json") or {}).get("people", {})
    eur = all((adm.get(p, {}).get("broad", {}).get("European") or 0) >= 0.9 for p in act.get("HLA-B*58:01", []))
    hla = ("An HLA-B*58:01 tag hit is a low-specificity\nmarker in Europeans (most likely a false positive), so "
           "HLA-B typing is needed before allopurinol.\n" if eur else
           "An HLA-B*58:01 tag hit needs confirming by HLA-B typing before allopurinol.\n")
    return (f"Most important: {', then '.join(top)}. " if top else "") + (hla if "HLA-B*58:01" in act else "")


def pair_word(a, b):
    """brothers / sisters / siblings"""
    return ({"M": "brothers", "F": "sisters"}.get(sex(a)) if sex(a) == sex(b) else None) or "siblings"


def find_trio(summary, viewer):
    """Mendel trio stats for viewer/mother/father (summary.json "mendel" key "trio <child>/<mother>/<father>")."""
    m = summary.get("mendel", {})
    return m.get(f"trio {viewer}/{person('mother')}/{person('father')}") or next(
        (v for k, v in m.items() if k.startswith(f"trio {viewer}/")), {})


def sibling_shares(kids, viewer):
    """[(other child, % of mom's half where both got the same grandparent's copy)] for the viewer.

    summary.json "children" gives "sibling_maternal_same_pct" as one number (two children) or {"<a>/<b>": pct};
    without it, each sibling's own "maternal_same_pct" (measured against self)."""
    v = kids.get("sibling_maternal_same_pct")
    if v is None:
        def own(c):
            return (kids.get(c) or {}).get("maternal_same_pct")
        if viewer == SELF:
            return [(s, own(s)) for s in children() if s != SELF and own(s) is not None]
        return [(SELF, own(viewer))] if own(viewer) is not None else []
    if not isinstance(v, dict):
        others = [c for c in children() if c != viewer and c in kids]
        return [(others[0], v)] if len(others) == 1 and viewer in kids else []
    out = []
    for k, pct in v.items():
        pair = k.split("/")
        if viewer in pair and pct is not None:
            out.append((next(c for c in pair if c != viewer) if len(pair) == 2 else k, pct))
    return out


def relatedness(summary, a, b):
    rel_ = summary.get("relatedness", {})
    return (rel_.get(f"{a} – {b}") or rel_.get(f"{b} – {a}") or {}).get("inferred")


def lineage_stats(summary, viewer):
    """The mtDNA and Y stat rows, computed from the haplogroup tables and summary.json."""
    mt = {r["SampleID"]: r["Haplogroup"] for r in read_tsv("mt_haplogroups.tsv")}
    ys = y_calls()
    mom, dad = person("mother"), person("father")
    out = []
    # mtDNA: the viewer's maternal line (siblings, mom, grandma) plus anyone else who happens to share it.
    hg = mt.get(viewer)

    def same_mt(p):  # chips resolve to different depths/private marks ("H1+@152"), so compare the base call
        a, b = (re.split(r"[+@ ]", x or "")[0] for x in (hg, mt.get(p)))
        return bool(a and b) and (a.startswith(b) or b.startswith(a))
    line = [p for p in order(mt) if same_mt(p) and role(p) in ("self", "sibling", "mother", "maternal_grandmother")]
    dad_same = bool(dad and same_mt(dad))
    who = ", ".join(rel(p) for p in line) + (", and dad too" if dad_same else "")
    unrelated = relatedness(summary, mom, dad) == "unrelated"
    out.append(("Maternal line (mtDNA)", hg or "—", (who[:1].upper() + who[1:]) or "—",
                "Passed unchanged from mother to child, so it traces mom's mother's mother… back thousands of years."
                + (f" {text('landing.mt_note')}" if text("landing.mt_note") else "")
                + (" Dad carries the same haplogroup." if dad_same else "")
                + (" Mom and dad show no close relationship (nothing closer than about 3rd cousins is detectable); "
                   "their shared maternal-line ancestor could lie anywhere from centuries to thousands of years "
                   "back, and the chip cannot say." if dad_same and unrelated else "")))
    # Y: the viewer's own call, or (for a daughter) the father's / a brother's, i.e. the line she belongs to.
    owner = viewer if viewer in ys else next((p for p in [dad] + children() if p in ys), None)
    if not owner:
        return out
    hg, isogg = ys[owner]["haplogroup"], ys[owner].get("isogg_2016", "")
    men = [p for p in order(ys) if role(p) in ("self", "sibling", "father")]
    names = [rel(p) for p in men if p != dad] + (["dad"] if dad in ys else [])
    who = and_join(names)
    r1b = isogg.startswith("R1b")
    surname = text("landing.surname")
    dad_hg = ys.get(dad, {}).get("haplogroup")
    finer = (viewer == owner and dad_hg and hg != dad_hg
             and ys[dad].get("isogg_2016", "").startswith(isogg or "?"))
    steps = len(ys[dad].get("isogg_2016", "")) - len(isogg) if finer else 0
    match = (summary.get("y", {}).get(f"{viewer} vs {dad}") or summary.get("y", {}).get(f"{dad} vs {viewer}") or {}
             ).get("concordance_pct")
    out.append(("Paternal line (Y)", hg, (who[:1].upper() + who[1:]) + (" · R1b, Western Europe" if r1b else ""),
                ("" if viewer == owner else f"My chip has no usable Y call; this is {rel(owner)}'s line, so mine too. "
                 if is_male(viewer) else "I don't carry a Y chromosome (only men do); this is "
                 + ("dad's" if owner == dad else f"{rel(owner)}'s (and so dad's)") + " line, my paternal line. ")
                + f"Passed from father to son, so it traces the {f'{surname} ' if surname else 'paternal '}surname line."
                + (" R1b is the most common Western European lineage" + ("; the Big Y-700 preview below guesses the "
                   "finer branch." if bigy_section() else ".") if r1b else
                   (" The Big Y-700 preview below guesses the finer branch." if bigy_section() else ""))
                + (f" My chip stops at {hg}, {'one step' if steps == 1 else 'a few steps'} above dad's {dad_hg}"
                   + (f"; the Y markers we share all match, so mine is {dad_hg} too." if match == 100 else ".")
                   if finer else "")))
    return out


def conflicts(m, confirmed):
    """The middle of the inheritance-conflicts text for a summary.json "mendel" record (trio or duo)."""
    if not m["errors"]:
        return (f"None in {m['snps']:,} SNPs: the data is clean, and {confirmed} (a non-parent gives thousands).")
    n = 100 / m["error_rate_pct"]
    return (f"Only {m['errors']} such conflicts, about 1 in {round(n, -3) if n >= 1000 else round(n):,.0f}: the data is "
            f"clean, and {confirmed} (a non-parent gives thousands).")


def crossover_note(n):
    """", a typical number for a mother (~40–45)" or why the count sits outside the usual range."""
    if not isinstance(n, (int, float)):
        return ""
    if 30 <= n <= 60:
        return ", a typical number for a mother (~40–45)"
    if n < 30:
        return (", fewer than the typical ~40–45 for a mother: the chip misses crossovers in short or marker-poor "
                "stretches")
    return ", more than the typical ~40–45 for a mother: some apparent switches may be chip misreads"


def main(viewer=None):
    global VIEWER
    VIEWER = viewer = viewer or SELF
    summary = load_json("summary.json") or {}
    pgx = read_tsv("pharmacogenomics.tsv")
    people = order(summary.get("qc", {}))
    mom, dad, gma = person("mother"), person("father"), person("maternal_grandmother")
    trio = find_trio(summary, viewer)
    meds = meds_section(pgx)

    kids = summary.get("children", {})
    mine = kids.get(viewer, {"crossovers": summary.get("crossovers", {}).get("total", "—"),
                             "grandma_share_pct": summary.get("grandma_share_pct", "—")})
    xo, share = {"total": mine["crossovers"]}, mine["grandma_share_pct"]
    # A child whose chip shares fewer informative sites with mom's (an older chip) misses short segments.
    sites = {c: kids[c]["informative_sites"] for c in children() if isinstance(kids.get(c), dict)
             and kids[c].get("informative_sites")}
    best = max(sites, key=sites.get) if sites else None
    sparse = viewer in sites and best != viewer and sites[viewer] < 0.75 * sites[best]  # clearly sparser, not chance
    me_qc = summary.get("qc", {}).get(viewer, {})
    chips = (f" Comparing my two chips gives a rough per-chip error of ~{me_qc['discordant'] / me_qc['shared'] * 50:.2f}% "
             "(very roughly, as it also includes differences in how the two platforms code some alleles)."
             if me_qc.get("shared") and me_qc.get("discordant") is not None else "")
    stats = []
    mendel = summary.get("mendel", {})
    duos = [(p, mendel[f"duo {viewer}/{p}"]) for p in (mom, dad) if p and f"duo {viewer}/{p}" in mendel]
    if "error_rate_pct" in trio:
        stats.append(
            ("Inheritance conflicts", f"{trio['error_rate_pct']:.3f}%", f"{trio['errors']} conflicts in {trio['snps']:,} SNPs",
             "SNPs where my result is impossible given mom's and dad's (e.g. I'm A/A but dad has no A). "
             + conflicts(trio, "mom and dad are confirmed as my biological parents") + " Most chip misreads don't create a "
             "conflict like this, so the true error rate is higher, though still very low." + chips))
    elif len(duos) == 1:  # one parent tested: the parent-child check (opposite homozygotes) analyze.py runs instead
        par, duo = duos[0]
        stats.append(
            ("Inheritance conflicts", f"{duo['error_rate_pct']:.3f}%", f"{duo['errors']} conflicts in {duo['snps']:,} SNPs",
             f"SNPs where {rel(par)} and I share no letter at all (e.g. I'm A/A and {rel(par)} is G/G), impossible for a "
             "parent and child. " + conflicts(duo, f"{rel(par)} is confirmed as my biological "
                                             + ("father" if par == dad else "mother"))
             + f" With only {rel(par)} tested, this check catches fewer chip misreads than a check with both parents, "
             "so the true error rate is higher, though still very low." + chips))
    elif mom and dad:
        stats.append(("Inheritance conflicts", "—", "not computed", "Not computed yet."))
    painted = share not in (None, "—") and xo.get("total") not in (None, "—")
    if mom and gma and not painted:  # grandma listed but her data (or the painting) isn't there yet
        stats.append(("Maternal crossovers", "—", "not computed", "Needs mom's and grandma's data; not computed yet."))
    elif mom and gma:
        stats.append(
            ("Maternal crossovers", str(xo.get("total", "—")), f"{share}% of mom's half from grandma",
             "When mom made the egg I came from, her two copies of each chromosome (one from grandma, one from grandpa) "
             f"swapped over {xo.get('total', '—')} times{crossover_note(xo.get('total'))}. The result is a patchwork: "
             f"{share}% of the DNA I got from mom is grandma's and the rest grandpa's (50% on average; roughly 40–60% "
             "is normal). This patchwork is what lets us rebuild grandpa's half."
             + (" My chip shares fewer markers with mom's, so short segments are missed: the true crossover count "
                f"is likely higher and the share is rougher than {rel(best)}'s." if sparse else "")))
    stats += lineage_stats(summary, viewer)
    for other, pct in sibling_shares(kids, viewer):
        stats.append(("Sibling sharing (mom's side)" + (f" · {col(other)}" if len(children()) > 2 else ""),
                      f"{pct:.0f}%", "of mom's half: same grandparent's copy",
                      f"Where {rel(other)} and I inherited the same grandparent's copy from mom (grandma's or "
                      "grandpa's). Full siblings average 50% with a wide spread. Dad's side can't be painted this way, "
                      "because dad's parents aren't tested."))
    stat_html = "".join(f'<div class="stat"><div><div class="label">{esc(a)}</div><div class="num">{esc(b)}</div>'
                        f'<div class="sub">{esc(c)}</div></div><p class="dim">{esc(d)}</p></div>' for a, b, c, d in stats)

    painting = mom in people and gma in people  # sibling sharing compares maternal paintings (mom + grandma)
    pending = [f"{full(s)}'s file — unlocks {'sibling sharing and ' if painting else ''}{pron(s)} own page"
               for s in children() if s not in people]
    missing_pgs = [full(p) for p in people if not os.path.exists(os.path.join(RESULTS, f"pgs_{p}.json"))]
    if missing_pgs:
        pending.append(f"Imputed downloads for {', '.join(missing_pgs)} — needed for polygenic scores")
    others = [c for c in children() if c != viewer and c in PAGES and c in people]
    switch = (f'This is {esc(named(viewer))}\'s page ("Me" = {esc(named(viewer))}). '
              + " · ".join(f'<a href="{esc(PAGES[c])}">{esc(cap(named(c)))}\'s version →</a>' for c in others)
              if others and viewer in PAGES else "")
    height = vtext("landing.height_check")  # e.g. measured height/build vs the height and BMI scores
    colon = vtext("landing.colon_history")  # e.g. who had polyps removed
    bone = vtext("landing.bone_density_context")
    qc = summary.get("qc", {})
    painted = bool(summary.get("crossovers")) or any(isinstance(v, dict) and v.get("crossovers") is not None
                                                     for v in kids.values())
    lines = [n for n, k in (("Y", "y"), ("mtDNA", "mt")) if summary.get(k)]
    report_sub = ", ".join(["QC"] + (["relatedness"] if len(qc) > 1 else [])
                           + (["Mendelian checks"] if summary.get("mendel") else [])
                           + (["maternal grandparent painting"] if painted else [])
                           + ([f"{'/'.join(lines)} matches"] if lines else []))
    whole = " × whole family" if len(qc) > 1 else ""
    pgs_links = " · ".join(f'<a href="pgs_{esc(p)}.html">{esc(col(p))}</a>' for p in pgs_people())
    pgs_card = (f'<div class="card report"><b>Polygenic scores →</b><div class="dim">Percentiles vs '
                f'{esc(" / ".join(n for _, n in pgs_groups()))}:\n{pgs_links}</div></div>' if pgs_links else "")
    local = [n for c, n in pgs_groups() if c not in ("EUR", "AMR")]
    local_note = f" ({esc(local[0])}-only values are in the full reports)" if local else ""
    pgs_table = family_pgs_table()
    chips_from = " + ".join(vendors()) or "Consumer DNA"
    diffs = "".join(
        f'<a class="card report" href="diffs_{esc(s)}.html"><b>{esc(cap(named(SELF)))} vs {esc(named(s))} →</b>'
        f'<div class="dim">Only the results where the two {pair_word(SELF, s)} differ</div></a>'
        for s in children() if s != SELF and os.path.exists(os.path.join(RESULTS, f"diffs_{s}.html")))

    doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{f"Family DNA ({esc(named(viewer))})" if viewer != SELF or FAMILY["people"].get(SELF, {}).get("name") else "Family DNA"}</title>
<style>
:root{{--bg:#fbfaf7;--fg:#1d1d1f;--dim:#6b6b70;--line:#e3e1db;--card:#fff;--accent:#2f6fb0;--track:#ecebe6;--warn:#b26b00;
--c1:#d6404e;--c2:#2a78d6;--c3:#1baf7a;--c4:#dca500;--c5:#8a63d2;--c6:#d94c1c;--sea:#eef2f5;--land:#e2dfd8;--land-hi:#cfcbc1}}
@media (prefers-color-scheme:dark){{:root{{--bg:#151517;--fg:#ececec;--dim:#9a9aa0;--line:#2c2c30;--card:#1d1d20;--accent:#6aa6e6;--track:#2c2c31;--warn:#e0a040;
--c1:#e05561;--c2:#3987e5;--c3:#199e70;--c4:#b58500;--c5:#9b74e0;--c6:#c9421a;--sea:#18191c;--land:#2e2f34;--land-hi:#45464d}}}}
*{{box-sizing:border-box}}
body{{background:var(--bg);color:var(--fg);font:15px/1.5 -apple-system,system-ui,sans-serif;margin:0;padding:28px 16px 48px}}
main{{max-width:980px;margin:0 auto}} h1{{font-size:30px;margin:0;letter-spacing:-.01em}} h2{{font-size:13px;text-transform:uppercase;letter-spacing:.06em;color:var(--dim);margin:0 0 12px}}
.dim{{color:var(--dim);font-size:13px}} .lede{{margin:4px 0 24px}}
.stack{{display:grid;gap:16px}}
.row2{{display:grid;gap:16px;grid-template-columns:repeat(auto-fit,minmax(min(360px,100%),1fr))}}
.row3{{display:grid;gap:16px;grid-template-columns:repeat(auto-fit,minmax(min(240px,100%),1fr))}}
.card{{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:18px;min-width:0;overflow-wrap:anywhere}}
.stat{{display:grid;grid-template-columns:minmax(0,220px) 1fr;gap:6px 24px;padding:14px 0;border-bottom:1px solid var(--line)}}
.stat:first-child{{padding-top:0}}
.stat .sub{{font-size:13px;font-weight:600}} .stat p{{margin:0;align-self:center}}
@media (max-width:600px){{.stat{{grid-template-columns:1fr}}}}
.person{{margin-bottom:12px}} .person:last-child{{margin-bottom:0}} .who{{font-weight:650;margin-bottom:4px}}
.person li{{margin-bottom:2px;font-size:14px}}
.stat .label{{font-size:12px;color:var(--dim)}} .stat .num{{font-size:24px;font-weight:650;margin:2px 0}}
ul{{margin:0;padding-left:18px}} li{{margin-bottom:10px}}
.flagged li::marker{{color:var(--warn)}}
h3{{font-size:14px;margin:22px 0 8px}} .tablewrap{{overflow-x:auto}}
table.fam{{border-collapse:collapse;width:100%;min-width:760px;table-layout:fixed;font-size:14px;font-variant-numeric:tabular-nums;overflow-wrap:normal}}
table.fam th,table.fam td{{text-align:left;padding:6px 8px;border-bottom:1px solid var(--line);vertical-align:top}}
table.fam th:first-child{{width:24%}} td.act{{color:var(--warn);font-weight:600}}
table.fam thead th{{color:var(--dim);font-weight:600;font-size:12px}} table.fam tbody th{{font-weight:400}}
.hi{{color:var(--warn);font-size:12px;font-weight:600}}
table.fam td.nw{{white-space:nowrap}} table.bigy th:first-child{{width:18%}}
table.fam tr.grp th{{font-size:11px;text-transform:uppercase;letter-spacing:.06em;color:var(--dim);font-weight:600;padding-top:14px}}
.anc-grid{{display:grid;gap:20px;grid-template-columns:repeat(auto-fit,minmax(min(380px,100%),1fr));align-items:start;margin:8px 0 4px}}
.omap{{width:100%;height:auto;display:block}} .omap{{background:var(--sea);border-radius:8px}} .omap .land{{fill:var(--land);stroke:var(--sea);stroke-width:.6}}
.omap .land.hi{{fill:var(--land-hi)}}
.omap .dot{{stroke:var(--card);stroke-width:2}} .olabels text{{font-size:11px;fill:var(--fg);paint-order:stroke;stroke:var(--card);stroke-width:3px;stroke-linejoin:round}}
.s1{{--c:var(--c1)}} .s2{{--c:var(--c2)}} .s3{{--c:var(--c3)}} .s4{{--c:var(--c4)}} .s5{{--c:var(--c5)}} .s6{{--c:var(--c6)}}
.dot{{fill:var(--c)}} .seg,.sw{{background:var(--c)}}
.olegend{{display:flex;flex-wrap:wrap;gap:4px 14px;font-size:12px;color:var(--dim);margin-bottom:10px}}
.olegend span{{display:inline-flex;align-items:center;gap:6px}} .sw{{display:inline-block;width:10px;height:10px;border-radius:3px}}
.obar{{display:grid;grid-template-columns:80px 1fr;align-items:center;gap:10px;margin:7px 0}} .obar .who{{font-size:13px}}
.segs{{display:flex;gap:2px;height:22px}} .seg{{display:flex;align-items:center;justify-content:center;font-size:11px;font-weight:600;color:#111;min-width:0;overflow:hidden}}
.seg:first-child{{border-radius:4px 0 0 4px}} .seg:last-child{{border-radius:0 4px 4px 0}} .seg:only-child{{border-radius:4px}}
details{{margin-top:10px;border-top:1px solid var(--line);padding-top:8px}} summary{{cursor:pointer;font-weight:600;font-size:14px}}
details[open] summary{{margin-bottom:8px}}
a.report{{display:block;text-decoration:none;color:inherit}} a.report:hover{{border-color:var(--accent)}}
.report b,.report a{{color:var(--accent)}}
.note{{font-size:13px;color:var(--dim);margin-top:24px;border-top:1px solid var(--line);padding-top:12px}}
</style></head><body><main>
{DISCLAIMER_HTML}
<h1>Family DNA</h1>
<p class="lede dim">{chips_from} raw data, analysed locally. The short version is below; full reports are linked at the bottom.
{switch}</p>

<div class="stack">
<section class="card"><div class="stats">{stat_html}</div>{gp_note() if grandpa_cells() else ""}</section>

<section class="card"><h2>💊 Medications</h2>{meds or "<p>Nothing non-routine</p>"}
<p class="dim">Highlighted = worth telling a doctor; hover a cell for the genotype and the CPIC guidance. {esc(meds_note(pgx))}
Confirm with a clinical PGx test before acting.</p></section>

{f'''<section class="card"><h2>📊 Genetic tendency scores (rank among ~500 Europeans)</h2>
{pgs_table}
<p class="dim">Each number is a percentile among the ~500 Europeans of 1000 Genomes (Tuscans, Spaniards, British, Finns,
Utah CEPH): 50 = average. {esc(height)}
Percentiles move ±~10 points depending on the reference group{local_note}. New scores (glaucoma, eGFR, blood pressure, uric acid)
were independently recomputed and direction-checked at known genes.
Colorectal cancer (PGS003850, 205 published variants), colon polyps (PGS001811, 2,231) and hemorrhoids (PGS001024, 786)
{f"were added because {esc(colon)}. For colon screening, that family history counts for more than these scores." if colon else
 "are included; for colon screening, family history counts for more than these scores."}
Bone density (PGS000657, 21,716 published variants; validated against real fractures){f" matters for {esc(bone)}" if bone else ""};
melanoma (PGS002247, 68 variants) already includes the MC1R red-hair variants (and other pigmentation genes), so it is
not extra evidence beyond them; MC1R status and skin type are the main drivers.
Hip and knee osteoarthritis (PGS002728, 44 variants; PGS002729, 24 variants) are weak predictors: weight, injuries and
physical work matter more. Rheumatoid arthritis (PGS004326, 154 variants, mostly immune HLA genes) is about as weak per SD
(HR ~1.2), and RA is rare (~0.5–1%), so even a high score means low absolute risk. Variant counts are the published
score sizes; the per-person pages show how many were usable.</p></section>''' if pgs_table else
'<section class="card"><h2>📊 Genetic tendency scores</h2><p class="dim">Not computed yet: polygenic scores need imputed '
'data and a run without --offline, --quick or --no-pgs.</p></section>'}

{wellness_section()}

{twentythree_section()}

{ancestry_section()}

{bigy_section()}

{neanderthal_section()}

{pigmentation_section()}

{grandpa_section()}

{f'<section class="card"><h2>⏳ Still to do</h2><ul>{"".join(f"<li>{esc(p)}</li>" for p in pending)}</ul></section>' if pending else ""}

<div class="row3">
<a class="card report" href="report.html"><b>Family report →</b><div class="dim">{report_sub}</div></a>
<a class="card report" href="pharmacogenomics.html"><b>Pharmacogenomics →</b><div class="dim">14 genes + 3 HLA alleles{whole}, CPIC-based calls and drug notes</div></a>
{diffs}
<a class="card report" href="doctor/index.html"><b>Doctor summaries →</b><div class="dim">{"One printable page per person" if len(qc) > 1 else "A printable page"}
with only the clinically relevant findings and suggested discussion points</div></a>
{pgs_card}
</div>
</div>

<p class="note">For learning and research only — not medical advice or a diagnosis. Consumer chips miss rare variants and gene
copy-number changes; a doctor or genetic counselor must review anything health-relevant and confirm it with a clinical test.</p>
</main></body></html>"""
    path = os.path.join(RESULTS, PAGES[viewer])
    with open(path, "w") as fh:
        fh.write(doc)
    print(path)


if __name__ == "__main__":
    have = (load_json("summary.json") or {}).get("qc", {})
    for v in PAGES:
        if v == SELF or (v in children() and v in have):
            main(v)
