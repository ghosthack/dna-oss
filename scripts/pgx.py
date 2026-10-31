"""Pharmacogenomic calls from measured chip genotypes, following CPIC allele definitions.

Usage: python3 -I scripts/pgx.py
Writes results/pharmacogenomics.tsv and results/pharmacogenomics.html.

Chips only see a handful of SNPs per gene and no copy-number changes, so every call
here is a screening result: confirm with a clinical PGx test before acting on it.
"""
import functools
import html
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dnalib import RESULTS, imputed_genotypes, label, load_person, name, people as family, person

# rsid -> (plus-strand variant allele, star allele / label)
SNPS = {
    "CYP2C19": {"rs4244285": ("A", "*2"), "rs4986893": ("A", "*3"), "rs12248560": ("T", "*17")},
    "CYP2C9": {"rs1799853": ("T", "*2"), "rs1057910": ("C", "*3"), "rs28371686": ("G", "*5"),
               "rs7900194": ("A", "*8"), "rs28371685": ("T", "*11")},
    "VKORC1": {"rs9923231": ("T", "-1639A")},
    "CYP4F2": {"rs2108622": ("T", "*3")},
    "SLCO1B1": {"rs4149056": ("C", "521C (*5/*15)")},
    "ABCG2": {"rs2231142": ("T", "Q141K")},
    "CYP3A5": {"rs776746": ("C", "*3"), "rs10264272": ("T", "*6")},
    "DPYD": {"rs3918290": ("T", "*2A"), "rs55886062": ("C", "*13"),
             "rs67376798": ("A", "c.2846A>T"), "rs75017182": ("C", "HapB3")},
    "TPMT": {"rs1800462": ("G", "*2"), "rs1800460": ("T", "*3B"), "rs1142345": ("C", "*3C")},
    "NUDT15": {"rs116855232": ("T", "*3")},
    "CYP2B6": {"rs3745274": ("T", "*6/*9")},
    "UGT1A1": {"rs887829": ("T", "*80 (tags *28)"), "rs4148323": ("A", "*6")},
    "G6PD": {"rs1050828": ("T", "A- (202A)"), "rs1050829": ("C", "376G"), "rs5030868": ("A", "Mediterranean")},
    "CYP2D6": {"rs3892097": ("T", "*4"), "rs5030655": ("D", "*6"), "rs1065852": ("A", "*10"),
               "rs28371725": ("T", "*41")},
    "HLA-B*57:01": {"rs2395029": ("G", "tag")},
    "HLA-A*31:01": {"rs1061235": ("T", "tag")},
    "HLA-B*58:01": {"rs9263726": ("A", "tag")},
}

DRUGS = {
    "CYP2C19": "clopidogrel, PPIs (omeprazole…), SSRIs (citalopram, escitalopram, sertraline), voriconazole, tricyclics",
    "CYP2C9": "warfarin, phenytoin, NSAIDs (celecoxib, ibuprofen…), siponimod",
    "VKORC1": "warfarin",
    "CYP4F2": "warfarin",
    "SLCO1B1": "statins (simvastatin, atorvastatin, rosuvastatin…)",
    "ABCG2": "rosuvastatin (CPIC); allopurinol response: research association only, no CPIC guidance",
    "CYP3A5": "tacrolimus",
    "DPYD": "fluorouracil, capecitabine (chemotherapy)",
    "TPMT": "azathioprine, mercaptopurine, thioguanine",
    "NUDT15": "azathioprine, mercaptopurine, thioguanine",
    "CYP2B6": "efavirenz, sertraline (CPIC); methadone: pharmacokinetic effect, no CPIC recommendation",
    "UGT1A1": "irinotecan, atazanavir; Gilbert syndrome",
    "G6PD": "rasburicase, primaquine, dapsone, high-dose vitamin C",
    "CYP2D6": "codeine, tramadol, tamoxifen, ondansetron, many antidepressants",
    "HLA-B*57:01": "abacavir",
    "HLA-A*31:01": "carbamazepine",
    "HLA-B*58:01": "allopurinol",
}


def counts(data, gene):
    """{rsid: variant allele count or None if not genotyped}."""
    out = {}
    for rs, (alt, _) in SNPS[gene].items():
        rec = data.get(rs)
        out[rs] = None if rec is None else rec[2].count(alt)
    return out


def star_list(gene, c):
    stars = []
    for rs, n in c.items():
        if n:
            stars += [SNPS[gene][rs][1]] * n
    return stars


def diplotype(stars, ref="*1"):
    stars = stars + [ref] * max(0, 2 - len(stars))
    def order(s):
        num = "".join(ch for ch in s if ch.isdigit())
        return (s != ref, int(num) if num else 0, s)
    stars = sorted(stars, key=order)
    if len(stars) == 2:
        return "/".join(stars)
    # More than two variant alleles on a diploid gene: phase (which copy carries what) is unknown, so no diplotype.
    seen = list(dict.fromkeys(stars))
    return "ambiguous (variants at " + ", ".join(s + (f" ×{stars.count(s)}" if stars.count(s) > 1 else "")
                                                for s in seen) + ")"


def missing(c):
    return [rs for rs, n in c.items() if n is None]


def call(gene, data, male):
    c = counts(data, gene)
    g = lambda rs: c.get(rs) or 0
    if gene == "CYP2C19":
        nf, inc = g("rs4244285") + g("rs4986893"), g("rs12248560")
        pheno = ("Poor metabolizer" if nf >= 2 else "Intermediate metabolizer" if nf == 1
                 else "Ultrarapid metabolizer" if inc == 2 else "Rapid metabolizer" if inc == 1
                 else "Normal metabolizer")
        return diplotype(star_list(gene, c)), pheno, c
    if gene == "CYP2C9":
        value = {"*2": 0.5, "*3": 0, "*5": 0.5, "*8": 0.5, "*11": 0.5}  # CPIC: *5 is decreased function (0.5)
        stars = star_list(gene, c)
        score = sum(value[s] for s in stars) + (2 - len(stars))
        pheno = "Normal metabolizer" if score == 2 else "Intermediate metabolizer" if score >= 1 else "Poor metabolizer"
        return diplotype(stars), f"{pheno} (activity score {score:g})", c
    if gene == "VKORC1":
        n = g("rs9923231")
        return {0: "G/G", 1: "G/A", 2: "A/A"}[n], ["Normal warfarin sensitivity", "Increased warfarin sensitivity",
                                                     "High warfarin sensitivity"][n], c
    if gene == "CYP4F2":
        n = g("rs2108622")
        return diplotype(star_list(gene, c)), "Typical" if n == 0 else "Slightly higher warfarin dose need", c
    if gene == "SLCO1B1":
        n = g("rs4149056")
        return ["T/T", "T/C", "C/C"][n], ["Normal function", "Decreased function", "Poor function"][n], c
    if gene == "ABCG2":
        n = g("rs2231142")
        return ["G/G", "G/T", "T/T"][n], ["Normal function", "Decreased function", "Poor function"][n], c
    if gene == "CYP3A5":
        nf = g("rs776746") + g("rs10264272")
        return diplotype(star_list(gene, c)), ["Normal metabolizer (expressor)", "Intermediate metabolizer",
                                               "Poor metabolizer (non-expressor, most common in Europeans)"][min(nf, 2)], c
    if gene == "DPYD":
        lost = {"rs3918290": 1, "rs55886062": 1, "rs67376798": 0.5, "rs75017182": 0.5}
        score = max(0, 2 - sum(g(rs) * v for rs, v in lost.items()))
        pheno = "Normal metabolizer" if score == 2 else "Intermediate metabolizer" if score >= 1 else "Poor metabolizer"
        untested = [SNPS[gene][rs][1] for rs in lost if c[rs] is None]
        caveat = f"; {', '.join(untested)} not tested" if untested else ""
        return diplotype(star_list(gene, c)), f"{pheno} (activity score {score:g}{caveat})", c
    if gene == "TPMT":
        nf = g("rs1800462") + max(g("rs1800460"), g("rs1142345"))
        stars = ["*2"] * g("rs1800462")
        both = min(g("rs1800460"), g("rs1142345"))
        stars += ["*3A"] * both + ["*3B"] * (g("rs1800460") - both) + ["*3C"] * (g("rs1142345") - both)
        return diplotype(stars), ["Normal metabolizer", "Intermediate metabolizer", "Poor metabolizer"][min(nf, 2)], c
    if gene == "NUDT15":
        if c["rs116855232"] is None:
            return "not covered", "Not determinable (*3 not on this chip)", c
        n = g("rs116855232")
        return diplotype(star_list(gene, c)), ["Normal metabolizer", "Intermediate metabolizer", "Poor metabolizer"][n], c
    if gene == "CYP2B6":
        n = g("rs3745274")
        return ["516G/G", "516G/T", "516T/T"][n], ["Normal metabolizer", "Intermediate metabolizer", "Poor metabolizer"][n], c
    if gene == "UGT1A1":
        n = g("rs887829") + g("rs4148323")
        return diplotype(star_list(gene, c)), ["Normal function", "Intermediate (Gilbert carrier)",
                                               "Poor function (likely Gilbert syndrome)"][min(n, 2)], c
    if gene == "G6PD":
        a_minus = g("rs1050828") if g("rs1050829") else 0  # A- copies = 202A count, on a 376G background
        med = g("rs5030868")
        deficient = a_minus or med
        if male:
            pheno = "Deficient" if deficient else "Normal"
        else:
            pheno = "Normal" if not deficient else "Variable (heterozygous female)" if max(a_minus, med) == 1 else "Deficient"
        return "variant" if deficient else "B (reference)", pheno, c
    if gene == "CYP2D6":
        if c["rs3892097"] is None:  # *4 is the most common non-functional allele; without it there's no call
            return "not covered", "Not determinable from this chip", c
        value = {"*4": 0, "*6": 0, "*10": 0.25, "*41": 0.25}  # current CPIC allele functionality table: *10 and *41 = 0.25
        # The *4 core haplotype includes 100C>T (rs1065852); count *10 only for 100T not explained by *4.
        n4 = g("rs3892097")
        stars = (["*4"] * n4 + ["*6"] * g("rs5030655") + ["*10"] * max(0, g("rs1065852") - n4)
                 + ["*41"] * g("rs28371725"))
        score = max(0, sum(value[s] for s in stars) + (2 - len(stars)))
        pheno = ("Normal metabolizer" if score >= 1.25 else "Intermediate metabolizer" if score > 0
                 else "Poor metabolizer")
        return diplotype(stars), f"Probably {pheno.lower()} (activity score ~{score:g}; deletions/duplications not detectable)", c
    if gene == "HLA-B*58:01":
        # rs9263726 is a good tag only in East Asians; in Europeans it is far more common than B*58:01 itself.
        n = g("rs9263726")
        return ("tag present" if n else "tag absent"), ("Possible HLA-B*58:01 carrier (tag SNP, low predictive value "
                                                        "in Europeans) — HLA-B*58:01 typing needed before allopurinol"
                                                        if n else "HLA-B*58:01 unlikely (tag SNP only)"), c
    if gene.startswith("HLA"):
        n = next(iter(c.values())) or 0
        return ("tag present" if n else "tag absent"), (f"Possible {gene} carrier — confirm by HLA typing"
                                                         if n else f"{gene} unlikely (tag SNP only)"), c
    raise KeyError(gene)


# Plain-language notes for non-routine results. Keep in sync with CPIC guidelines (cpicpgx.org).
NOTES = {
    ("CYP2C19", "Ultrarapid"): "Faster clearance of several CYP2C19 drugs. CPIC: PPIs (omeprazole, lansoprazole…) may need a higher dose; for citalopram/escitalopram and tricyclics, consider an alternative not mainly cleared by CYP2C19; voriconazole: choose an alternative. Clopidogrel works normally.",
    ("CYP2C19", "Rapid"): "Somewhat faster clearance of several CYP2C19 drugs. CPIC: clopidogrel works normally; PPIs at standard dose, but may need a higher dose for H. pylori/erosive esophagitis; voriconazole (adults): choose an alternative; tertiary-amine tricyclics (amitriptyline, clomipramine…): consider an alternative not mainly cleared by CYP2C19.",
    ("CYP2C19", "Intermediate"): "Reduced clopidogrel activation; CPIC suggests an alternative antiplatelet for cardiac indications.",
    ("CYP2C19", "Poor"): "Clopidogrel largely ineffective; several SSRIs/PPIs reach higher levels.",
    ("SLCO1B1", "Poor"): "Higher statin blood levels and myopathy risk. CPIC 2022: avoid simvastatin and lovastatin; atorvastatin or rosuvastatin ≤20 mg starting dose (rosuvastatin ≤10 mg if ABCG2 is also poor function); pravastatin ≤40 mg, pitavastatin ≤1 mg, fluvastatin ≤40 mg (an alternative statin instead if CYP2C9 activity is also reduced). The single most actionable finding if a statin is ever prescribed (rs4149056 only — confirm clinically).",
    ("SLCO1B1", "Decreased"): "Moderately higher statin exposure and myopathy risk. CPIC 2022: prefer an alternative to simvastatin (if used, keep it below 20 mg/day); lovastatin: alternative or ≤20 mg/day; starting-dose limits atorvastatin ≤40 mg, pitavastatin ≤2 mg; rosuvastatin, pravastatin and fluvastatin: usual starting dose, with myopathy awareness at higher doses (fluvastatin ≤20 mg if CYP2C9 activity is also reduced, an alternative statin if CYP2C9 is poor; rosuvastatin ≤10 mg if ABCG2 is poor function). Based on rs4149056 only — confirm clinically.",
    ("ABCG2", "Decreased"): "Modestly higher rosuvastatin exposure; CPIC: no dose change for ABCG2 decreased function alone (CPIC lowers the rosuvastatin start to ≤20 mg for ABCG2 poor function, or ≤10 mg when that is combined with reduced SLCO1B1 function). Slightly higher uric acid/gout tendency.",
    ("ABCG2", "Poor"): "Higher rosuvastatin exposure. CPIC: rosuvastatin ≤20 mg starting dose (≤10 mg if SLCO1B1 function is also reduced). Higher uric acid/gout tendency.",
    ("TPMT", "Intermediate"): "Thiopurines (azathioprine, mercaptopurine) build up: CPIC recommends reduced starting doses.",
    ("TPMT", "Poor"): "Thiopurines can be life-threatening at standard doses.",
    ("NUDT15", "Intermediate"): "Thiopurines need reduced starting doses. Combined with a TPMT intermediate result, CPIC recommends starting at 20–50% of the standard dose.",
    ("NUDT15", "Poor"): "Thiopurines can be life-threatening at standard doses.",
    ("VKORC1", "Increased"): "Needs somewhat less warfarin; a genotype-guided dosing calculator (warfarindosing.org) handles this together with CYP2C9 and CYP4F2.",
    ("VKORC1", "High"): "Needs considerably less warfarin.",
    ("CYP2C9", "Intermediate"): "Slower clearance of warfarin, phenytoin and some NSAIDs.",
    ("CYP2C9", "Poor"): "Much slower clearance of warfarin, phenytoin and some NSAIDs.",
    ("DPYD", "Intermediate"): "Fluoropyrimidine chemo needs a reduced dose — important to tell an oncologist.",
    ("DPYD", "Poor"): "Fluoropyrimidine chemo should be avoided.",
    ("UGT1A1", "Poor"): "Likely Gilbert syndrome; irinotecan dose reduction.",
    ("G6PD", "Deficient"): "Avoid rasburicase, primaquine, dapsone and other oxidant drugs.",
    ("G6PD", "Variable"): "One variant copy on the X chromosome: red-cell G6PD activity can range from normal to deficient "
                          "(X-inactivation), so the genotype can't predict it. CPIC: measure G6PD enzyme activity before "
                          "rasburicase, primaquine, tafenoquine, dapsone or other oxidant drugs.",
    ("CYP2D6", "intermediate"): "Codeine/tramadol less effective; tamoxifen less active.",
    ("CYP2D6", "poor"): "Avoid codeine and tramadol; tamoxifen less active.",
}


@functools.lru_cache(maxsize=None)
def ancestry_json():
    """results/ancestry.json from scripts/ancestry.py ({} until it has run)."""
    path = os.path.join(RESULTS, "ancestry.json")
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        return json.load(fh)


def mostly_european(entry):
    """True when an ancestry.json entry (a person, or "grandpa_inferred") is >=90% broad European."""
    return ((entry or {}).get("broad") or {}).get("European", 0) >= 0.9


def note(gene, pheno, person=None):
    """Plain-language note for a non-routine result ('' = routine). `person` picks ancestry-dependent wording."""
    for (g, key), text in NOTES.items():
        if g == gene and key in pheno:
            return text
    if gene == "HLA-B*58:01" and "Possible" in pheno:
        # The tag is validated in East Asians; "most likely a false positive" only holds for European ancestry.
        why = ("rs9263726 tags HLA-B*58:01 reliably only in East Asians; in Europeans the tag is common (~15–23%) while "
               "B*58:01 itself is ~1–2%, so this is most likely a false positive."
               if person and mostly_european(ancestry_json().get("people", {}).get(person)) else
               "rs9263726 tags HLA-B*58:01 reliably in East Asians; its accuracy varies by ancestry, so this is only a "
               "hint.")
        return why + (" If allopurinol is ever considered, ask for HLA-B*58:01 typing first (CPIC: contraindicated only "
                      "if typing confirms the allele).")
    if gene.startswith("HLA") and "Possible" in pheno:
        return f"Tag SNP suggests {gene}; risk of severe reaction to {DRUGS[gene]}. Confirm by HLA typing first."
    return ""


def main():
    people = {}
    for p in family():
        data, _ = load_person(p)
        if data:
            male = sum(1 for r in data.values() if r[0] == "Y") > 500
            people[p] = (data, male)
    imputed = {}
    for p, (d, _) in people.items():  # fill chip gaps (e.g. NUDT15 on the v4 chip) from confident imputed calls
        imputed[p] = imputed_genotypes(p, {rs for snps in SNPS.values() for rs in snps if rs not in d})
        d.update({rs: (None, None, gt) for rs, gt in imputed[p].items()})
    calls = {p: {gene: call(gene, d, male) for gene in SNPS} for p, (d, male) in people.items()}

    with open(os.path.join(RESULTS, "pharmacogenomics.tsv"), "w") as fh:
        fh.write("person\tgene\tdiplotype\tphenotype\tmissing_snps\timputed_snps\n")
        for p, genes in calls.items():
            for gene, (dip, pheno, c) in genes.items():
                imp = ",".join(rs for rs in c if rs in imputed[p])
                fh.write(f"{p}\t{gene}\t{dip}\t{pheno}\t{','.join(missing(c))}\t{imp}\n")
                print(f"{p:17s} {gene:12s} {dip:22s} {pheno}" + (f"  [missing {missing(c)}]" if missing(c) else ""))

    esc = html.escape
    def header(p):  # "Me", "Mom", "Sister"… as on the other pages, with the family.json name (or id) underneath
        same = sum(label(q) == label(p) for q in calls) > 1
        sub = name(p) if name(p) != label(p) else (p if same else "")
        return f"<th>{esc(label(p))}" + (f'<div class="muted">{esc(sub)}</div>' if sub else "") + "</th>"
    head = "".join(header(p) for p in calls)
    rows = []
    for gene in SNPS:
        cells = []
        for p in calls:
            dip, pheno, c = calls[p][gene]
            flag = " class=\"flag\"" if note(gene, pheno) else ""
            miss = f'<div class="muted">missing {len(missing(c))}/{len(c)} SNPs</div>' if missing(c) else ""
            miss += ('<div class="muted">imputed: ' + esc(", ".join(rs for rs in c if rs in imputed[p])) + "</div>"
                     if any(rs in imputed[p] for rs in c) else "")
            cells.append(f"<td{flag}><b>{esc(dip)}</b><br>{esc(pheno)}{miss}</td>")
        rows.append(f"<tr><th scope=\"row\">{esc(gene)}<div class=\"muted\">{esc(DRUGS[gene])}</div></th>{''.join(cells)}</tr>")
    me = person("self")
    mine = [(gene, *calls[me][gene][:2]) for gene in SNPS if note(gene, calls[me][gene][1])] if me in calls else []
    findings = "".join(f"<li><b>{esc(g)} — {esc(ph)}</b> ({esc(dip)}). {esc(note(g, ph, me))}</li>" for g, dip, ph in mine)

    mine_html = (f'<h2>What stands out for me</h2><div class="card"><ul>{findings or "<li>Nothing non-routine.</li>"}</ul></div>\n'
                 if me in calls else "")
    doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Family Pharmacogenomics</title>
<style>
:root{{--bg:#fbfaf7;--fg:#1d1d1f;--muted:#6b6b70;--line:#e3e1db;--card:#fff;--flag:#fff3d6}}
@media (prefers-color-scheme:dark){{:root{{--bg:#151517;--fg:#ececec;--muted:#9a9aa0;--line:#2c2c30;--card:#1d1d20;--flag:#3a3020}}}}
body{{background:var(--bg);color:var(--fg);font:15px/1.5 -apple-system,system-ui,sans-serif;margin:0;padding:24px 16px}}
main{{max-width:1100px;margin:0 auto}} h1{{font-size:26px;margin:0 0 4px}} h2{{font-size:18px;margin:28px 0 8px}}
.muted{{color:var(--muted);font-size:12px;font-weight:400}} .card{{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:16px;overflow-x:auto}}
table{{border-collapse:collapse;width:100%;font-size:13px}} th,td{{text-align:left;vertical-align:top;padding:8px;border-bottom:1px solid var(--line)}}
thead th{{color:var(--muted)}} td.flag{{background:var(--flag)}} li{{margin-bottom:8px}}
.warn{{border-left:4px solid #d08a00;padding:8px 12px;background:var(--card)}}
</style></head><body><main>
<h1>Family pharmacogenomics</h1>
<p class="muted">Measured chip genotypes, CPIC allele definitions. Generated locally.</p>
<p class="warn"><b>Screening only.</b> Consumer chips test a few SNPs per gene and miss rare variants, gene deletions
and duplications. Don't start, stop or change any medication based on this; bring it to a doctor or pharmacist,
who can order a clinical pharmacogenomic panel to confirm.</p>
{mine_html}<h2>Whole family</h2><div class="card"><table><thead><tr><th>Gene · affected drugs</th>{head}</tr></thead>
<tbody>{''.join(rows)}</tbody></table></div>
<p class="muted">Highlighted cells have a non-routine result. Where phase can't be determined, variants are assumed to be
on separate copies, except TPMT: c.460A and c.719G together are assumed to be on one copy (*3A, the common allele). If they
were on separate copies (*3B/*3C), the result would be Poor rather than Intermediate metabolizer, which a clinical test can resolve.</p>
</main></body></html>"""
    with open(os.path.join(RESULTS, "pharmacogenomics.html"), "w") as fh:
        fh.write(doc)


if __name__ == "__main__":
    main()
