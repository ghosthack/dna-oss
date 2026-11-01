"""Self vs each sibling: only the results where the two differ (results/diffs_<sibling>.html, one page per sibling).

Usage: python3 -I scripts/sibling_diffs.py   (run after the other scripts; a sibling is skipped until they have results)

Everything is read from results/ (pharmacogenomics.tsv, wellness.tsv, pgs_*.json, pigmentation.tsv/.json, summary.json,
23andme_reports.json) plus raw/family_notes.json. Results are compared by outcome (the landing page's short label), so
a genotype-only difference is listed on its own line. Differences that only reflect a missing call (a marker an older
23andMe chip lacks, or a failed read) are re-scored with the genotype the parents imply where both carry two identical
copies; the ones the parents can't settle are listed apart.
Also re-checks each sibling's imputed data: at polygenic-score positions where both parents are confidently
homozygous, a child's genotype is fixed, so their imputed calls must match it.
"""
import csv
import glob
import html
import json
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dnalib import DISCLAIMER_HTML, BCFTOOLS, FAMILY, RAW, RESULTS, WORK, find, load_person, name, person, sex, siblings, text
from landing import TRAIT_LABELS, WELLNESS_SHORT
from pgs import load_scores, match
from pgx import DRUGS, note
from pigmentation import MC1R_STRONG
import wellness as W

esc = html.escape
A = person("self")
MOM, DAD, GRANDMA = person("mother"), person("father"), person("maternal_grandmother")
CHECK_CHROMS = ("1", "6", "11", "19")  # imputation check: a sample of chromosomes is enough (~11k score positions)
IMPUTE_OK_PCT = 2.0  # a mismatch rate below this at parent-fixed positions counts as accurate imputation
NUM = {1: "one", 2: "two", 3: "three", 4: "four"}
MHC = ("chr6", 28_510_120, 33_480_577)  # GRCh38 extended MHC: the HLA immune-gene region
HLA_SHARE = 0.5  # a score gap counts as HLA-driven when at least this share of it comes from the MHC
NEG = ("Not a carrier", "Variants not detected", "No ")  # results that report no counted variant


def tsv(name):
    path = os.path.join(RESULTS, name)
    return list(csv.DictReader(open(path), delimiter="\t")) if os.path.exists(path) else []


def load(name):
    path = os.path.join(RESULTS, name)
    return json.load(open(path)) if os.path.exists(path) else None


def page(p):
    """A child's own page: family.json "pages", else index.html for self and <id>.html for a sibling."""
    return (FAMILY.get("pages") or {}).get(p) or ("index.html" if p == A else f"{p}.html")


def chip_version(p):
    """23andMe chip version from the raw file name (..._v4_... / ..._v5_...), or None."""
    f = find(p, "genome_*")
    m = re.search(r"_v(\d+)_", os.path.basename(f)) if f else None
    return int(m.group(1)) if m else None


def older(p, q):
    """True when p's 23andMe chip is an older version than q's."""
    vp, vq = chip_version(p), chip_version(q)
    return bool(vp and vq and vp < vq)


def pair_word(a, b):
    """"brothers", "sisters" or "siblings" for two children."""
    return {"M": "brothers", "F": "sisters"}.get(sex(a), "siblings") if sex(a) == sex(b) else "siblings"


def poss(p):
    return {"M": "his", "F": "her"}.get(sex(p), "their")


def chip_only(a, b):
    """True when two wellness results differ only because a marker is missing from one chip."""
    neg = ("Not a carrier", "Variants not detected")
    return ("Not determinable" in a + b or "not tested on this chip" in a + b
            or (a.startswith(neg) and b.startswith(neg)))


def outcome(result):
    """A result's outcome: the landing page's short label, plus any copy counts (×1/×2) it states."""
    return next((s for k, s in WELLNESS_SHORT if k in result), result), tuple(re.findall(r"×\d", result))


def genos(r):
    """{rsid: genotype or None} from a wellness.tsv row ("rs1=AG rs2=-- rs3=CT*"; * = filled from phased/imputed)."""
    return {rs: (g.rstrip("*") if not g.startswith("--") else None)
            for rs, g in (x.split("=", 1) for x in r["genotypes"].split())}


def parents_fix(rsids):
    """{rsid: child genotype} where mom and dad both carry two identical copies, so a child's genotype is certain."""
    if not (MOM and DAD):
        return {}
    mom, dad = load_person(MOM)[0], load_person(DAD)[0]
    if not (mom and dad):
        return {}
    return {rs: "".join(sorted(mom[rs][2][0] + dad[rs][2][0])) for rs in rsids if rs in mom and rs in dad
            and len(set(mom[rs][2])) == 1 and len(set(dad[rs][2])) == 1}


REPORT_FN = {rep: (cat, snps, fn) for cat, reports in (("23andMe wellness", W.REPORTS), ("Taste & diet", W.TASTE),
                                                       ("Nutrition & health", W.NUTRITION),
                                                       ("Inherited conditions", W.CONDITIONS),
                                                       ("Carrier status", W.CARRIERS))
             for rep, _, snps, fn in reports}


def rescore(rep, p, g):
    """The wellness result for p from genotypes g, scored the way wellness.py does, or None for an unknown report."""
    if rep not in REPORT_FN:
        return None
    cat, snps, fn = REPORT_FN[rep]
    g = {rs: g.get(rs) for rs in snps}
    ok = (all(v for rs, v in g.items() if rs not in W.OPTIONAL)
          or ((cat == "Carrier status" or rep in W.ANY_CALLED) and any(g.values())))
    return (fn(g, sex(p)) if fn in W.BY_SEX else fn(g)) if ok else "Not determinable (no-call)"


_HDS = {}  # person -> {CHROM:POS:REF:ALT: imputed dosage}, or None without imputed data; loaded once per run


def load_hds(p):
    if p not in _HDS:
        bcfdir = glob.glob(os.path.join(WORK, "imputed", p, "*"))
        _HDS[p] = None
        if bcfdir:
            _HDS[p] = {}
            for c in CHECK_CHROMS:
                reg = os.path.join(WORK, "pgs", f"regions.chr{c}.tsv")
                out = subprocess.run([BCFTOOLS, "query", "-R", reg, "-f",
                                      "%CHROM:%POS:%REF:%ALT\t[%HDS]\n", os.path.join(bcfdir[0], f"chr{c}.bcf")],
                                     capture_output=True, text=True).stdout
                for line in out.splitlines():
                    k, h = line.split("\t")
                    _HDS[p][k] = sum(float(x) for x in h.split(","))
    return _HDS[p]


def imputation_check(kids):
    """{child: (positions, mismatch %)} for score positions fixed by two confidently homozygous parents."""
    if not (MOM and DAD):
        return None
    hds = {p: load_hds(p) for p in (MOM, DAD, *kids)}
    if any(v is None for v in hds.values()):
        return None
    res = {}
    for child in kids:
        n = bad = 0
        for k, m in hds[MOM].items():
            d, c = hds[DAD].get(k), hds[child].get(k)
            if d is None or c is None:
                continue
            rm, rd = round(m), round(d)
            if rm in (0, 2) and rd in (0, 2) and abs(m - rm) < .1 and abs(d - rd) < .1:
                n += 1
                bad += abs(round(c) - (rm + rd) / 2) >= 1
        res[child] = (n, bad / n * 100 if n else 0)
    return res


def far_note(pair):
    return (" ⚠ marks a score more than 2.5 of its expected spread away from mom and dad's average, which inheritance "
            f"rarely produces: possible, but treat that number as less certain. A large gap where both {pair} sit "
            "equally far from the average in opposite directions instead points to a "
            "few large-effect genes, such as HLA, inherited differently.")


_SCORES = {}  # PGS id -> [(chrom, pos, effect, other, weight)] inside the MHC; loaded once per run


def hla_part(pid, a, b):
    """Score difference b - a (raw score units) that comes from the MHC, from the imputed dosages; None without data."""
    ha, hb = load_hds(a), load_hds(b)
    if not (ha and hb) or MHC[0][3:] not in CHECK_CHROMS:
        return None
    if not _SCORES:
        _SCORES.update({k: [v for v in s["variants"] if v[0] == MHC[0] and MHC[1] <= v[1] <= MHC[2]]
                        for k, s in load_scores().items()})
    if pid not in _SCORES:
        return None
    sites = {}
    for k in ha.keys() & hb.keys():
        c, pos, ref, alt = k.split(":")
        c = c if c.startswith("chr") else "chr" + c
        sites.setdefault((c, int(pos)), []).append((ref, alt, ha[k], hb[k]))
    diff = 0.0
    for chrom, pos, effect, other, w in _SCORES[pid]:
        rec, is_alt = match(sites.get((chrom, pos), []), effect, other)
        if rec:
            da, db = (d if is_alt else 2 - d for d in rec[2:])
            diff += w * (db - da)
    return diff


def table(head, rows):
    return ('<div class="tablewrap"><table><thead><tr>' + "".join(f"<th>{h}</th>" for h in head) + "</tr></thead><tbody>"
            + "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows) + "</tbody></table></div>")


def meds(B):
    """(rows that differ, genes that are the same, same genes that carry medication guidance,
    genes with the same phenotype from a different diplotype as [gene, diplotype A, diplotype B, phenotype])."""
    pg = {(r["person"], r["gene"]): r for r in tsv("pharmacogenomics.tsv")}
    rows, same, flagged, geno_only = [], [], [], []
    for gene in dict.fromkeys(g for _, g in pg):
        a, b = pg.get((A, gene)), pg.get((B, gene))
        if not a or not b:
            continue
        if a["phenotype"].split(" (")[0] == b["phenotype"].split(" (")[0]:  # same outcome, as the landing table shows
            if a["diplotype"] != b["diplotype"]:
                geno_only.append([gene, a["diplotype"], b["diplotype"], a["phenotype"].split(" (")[0]])
                continue
            same.append(gene)
            if note(gene, a["phenotype"]) and not gene.startswith("HLA"):  # HLA tag SNPs are unconfirmed guesses
                flagged.append(gene)
            continue
        cell = lambda r: (f"<b>{esc(r['diplotype'])}</b><br>{esc(r['phenotype'].split(' (')[0])}"
                          + (f'<div class="dim">{esc(note(gene, r["phenotype"], r["person"]))}</div>' if note(gene, r["phenotype"]) else ""))
        rows.append([f"<b>{esc(gene)}</b>", cell(a), cell(b)])
    return rows, same, flagged, geno_only


def wellness(B):
    """(rows that differ, missing-call pairs, count of same results, (report, outcome, [(rsid, genotype A, genotype B)])
    for the same outcome from a different genotype). Wording that differs only by sex, on the same genotypes, counts as same."""
    w = {(r["person"], r["report"]): r for r in tsv("wellness.tsv")}
    real, chip, same, geno_only = [], [], 0, []
    for rep in dict.fromkeys(r for _, r in w):
        a, b = w.get((A, rep)), w.get((B, rep))
        if not a or not b:
            continue
        if a["result"] == b["result"]:
            same += 1
        elif chip_only(a["result"], b["result"]):
            chip.append((rep, a, b))
        elif outcome(a["result"]) == outcome(b["result"]):
            if genos(a) == genos(b):
                same += 1
            else:
                ga, gb = genos(a), genos(b)
                geno_only.append((rep, outcome(a["result"])[0],
                                  [(rs, ga[rs] or "--", gb.get(rs) or "--") for rs in ga if ga[rs] != gb.get(rs)]))
        else:
            real.append([esc(rep), esc(a["result"]), esc(b["result"])])
    return real, chip, same, geno_only


MIDPARENT_SD = 0.71  # a child's score varies around the parents' average by ~sqrt(1/2) SD (many small effects)


def pgs(B):
    """(big, small) score gaps; each item is (label, pct A, pct B, z gap, children far from the parents' average,
    PGS id, raw score gap B - A)."""
    da, db = load(f"pgs_{A}.json"), load(f"pgs_{B}.json")
    if not da or not db:
        return [], []
    par = {p: load(f"pgs_{p}.json") or {} for p in (MOM, DAD) if p}
    out = []
    for k in da:
        if k not in db:
            continue
        t = da[k]["trait"].replace(" (T2D)", "").replace(" (BMI)", "")
        t = TRAIT_LABELS.get(t[0].upper() + t[1:].replace("Cancer", "cancer"), t[0].upper() + t[1:])
        ga, gb = da[k]["groups"]["EUR"], db[k]["groups"]["EUR"]
        far = []
        if len(par) == 2 and all(k in par[p] for p in par):  # sex-specific scores (prostate) have no parent pair
            mid = sum(par[p][k]["groups"]["EUR"]["z"] for p in par) / 2
            far = [c for c, g in ((A, ga), (B, gb)) if abs(g["z"] - mid) / MIDPARENT_SD > 2.5]
        out.append((t, ga["percentile"], gb["percentile"], gb["z"] - ga["z"], far, k, db[k]["score"] - da[k]["score"]))
    out.sort(key=lambda x: -abs(x[3]))
    big = [x for x in out if abs(x[2] - x[1]) >= 20]
    return big, [x for x in out if x not in big]


def ancestry_23(B):
    """23andMe ancestry rows. family.json "ancestry": {"sibling_diff_rows": [label or [label, name prefix], ...]}
    picks the populations; without it, every population whose percentage differs between the two is shown."""
    d = load("23andme_reports.json") or {}
    trees = {p: d.get(p, {}).get("ancestry_tree", []) for p in (A, B)}

    def pct(p, prefix):
        hit = next((x for x in trees[p] if x[0].startswith(prefix)), None)
        return hit[1] if hit else "—"
    picks = (FAMILY.get("ancestry") or {}).get("sibling_diff_rows")
    if picks:
        return [[esc(lab), pct(A, pre), pct(B, pre)]
                for lab, pre in ((x, x) if isinstance(x, str) else x for x in picks)]
    if not (trees[A] and trees[B]):
        return []
    val = {p: {x[0]: x[1] for x in trees[p]} for p in (A, B)}
    return [[esc(lab), val[A].get(lab, "—"), val[B].get(lab, "—")]
            for lab in dict.fromkeys(x[0] for p in (A, B) for x in trees[p]) if val[A].get(lab) != val[B].get(lab)]


def same_looks(B, hair=""):
    """Appearance markers the two share, as a phrase ("predicted brown eyes and one copy of ... each"), or "".
    hair: a hair-color prediction the two share ("black"), named alongside the eye color."""
    pj = (load("pigmentation.json") or {}).get("people", {})
    a, b = pj.get(A), pj.get(B)
    if not (a and b):
        return ""
    predicted = [f"{hair} hair"] if hair else []
    if a.get("eye_call") == b.get("eye_call") and a.get("eye_call") not in (None, "not predicted"):
        predicted.append(f"{a['eye_call']} eyes")
    parts = [f"predicted {' and '.join(predicted)}"] if predicted else []
    carried = lambda r: {v for v, dose in r.get("mc1r", {}).get("variants", {}).items() if round(dose) >= 1}
    n_a, n_b = (round(r.get("mc1r", {}).get("strong_R_allele_count") or 0) for r in (a, b))
    w_a, w_b = (round(r.get("mc1r", {}).get("weak_r_allele_count") or 0) for r in (a, b))
    if n_a == n_b and n_a:
        shared = [v for v in MC1R_STRONG.values() if v in carried(a) & carried(b)]
        src = ""
        if len(shared) == 1:  # say whose copy it is when the tested elders show one route only
            has = {p: shared[0] in carried(pj[p]) for p in (GRANDMA, MOM, DAD) if p in pj}
            if has.get(DAD) is False and has.get(GRANDMA):
                src = "grandma's "
            elif has.get(DAD) is False and has.get(MOM):
                src = "mom's "
            elif has.get(MOM) is False and has.get(DAD):
                src = "dad's "
        variant = (f"MC1R {' and '.join(shared)} red-hair variant{'s' if len(shared) > 1 else ''}" if shared
                   else "strong MC1R red-hair variants")
        phrase = f"{NUM.get(n_a, n_a)} cop{'y' if n_a == 1 else 'ies'} of {src}{variant} each"
        if w_a == w_b and w_a:
            phrase += f" (plus {NUM.get(w_a, w_a)} weaker variant{'s' if w_a > 1 else ''})"
        parts.append(phrase)
    # "predicted black hair and brown eyes, and one copy of ...": a comma keeps the two "and"s apart.
    return (", and " if len(predicted) > 1 and len(parts) > 1 else " and ").join(parts)


def write_page(B):
    NAME = {A: name(A), B: name(B)}
    pair = pair_word(A, B)
    med_rows, med_same, med_flagged, med_geno = meds(B)
    real, chip, same_n, well_geno = wellness(B)
    big, small = pgs(B)
    kids = (load("summary.json") or {}).get("children", {})
    check = imputation_check((A, B))

    # Missing-call gaps: when both parents carry two identical copies, a child can only inherit those, so the missing
    # call is filled with that genotype and the report re-scored.
    snps = {rs for _, a, b in chip for r in (a, b) for rs in genos(r)}
    fixed = parents_fix(snps)
    chip_items, open_items, lacking = [], [], []
    for rep, a, b in chip:
        g = {A: genos(a), B: genos(b)}
        # Who is missing a call the other has (the sibling's row is checked first); else anyone missing one.
        who, miss = B, []
        for p, q in ((B, A), (A, B), (B, None), (A, None)):
            m = [rs for rs, v in g[p].items() if v is None and (q is None or g[q].get(rs))]
            if m:
                who, miss = p, m
                break
        other = A if who == B else B
        res = {A: a["result"], B: b["result"]}
        filled = {**g[who], **{rs: fixed[rs] for rs in miss if rs in fixed}}
        implied = rescore(rep, who, filled) if miss and all(rs in fixed for rs in miss) else None
        if implied is None or implied.startswith("Not determinable") or res[other].startswith("Not determinable"):
            have = (f"{NAME[other]}: {res[other].split(' (')[0].split(' — ')[0]}" if not res[other].startswith("Not determinable")
                    else f"{NAME[other]} has no result either")
            open_items.append((who, f"<li><b>{esc(rep)}</b>: {esc(NAME[who])} has no call at "
                                    f"{esc(', '.join(miss) or 'the marker')}, and mom and dad's genotypes don't settle "
                                    f"it ({esc(have)}).</li>"))
            continue
        gts = " and ".join(f"{fixed[rs]} at {rs}" for rs in miss)
        if outcome(implied) != outcome(res[other]):  # the parents show the two really differ
            res[who] = f"{implied} (from mom and dad: {gts}; {poss(who)} own file has no call)"
            real.append([esc(rep), esc(res[A]), esc(res[B])])
            continue
        lacking.append(who)
        how = (f"settled by mom and dad: both carry two identical, non-variant copies, so {NAME[who]} does too"
               if implied.startswith(NEG) else
               f"settled by mom and dad: each carries two identical copies, so {NAME[who]} must be {gts}: "
               f"{implied[0].lower() + implied[1:]}")
        chip_items.append(f"<li><b>{esc(rep)}</b>: {esc(how)}.</li>")

    pgs_rows = [[esc(t) + ("".join(f' <span class="warn">⚠ {esc(NAME[c])}\'s is unusually far from the parents\' '
                               f'average</span>' for c in far)), f"{a:.0f}", f"{b:.0f}", f"{dz:+.1f} SD"]
                for t, a, b, dz, far, _, _ in big]
    small_txt = ", ".join(f"{esc(t)} ({a:.0f} vs {b:.0f})" for t, a, b, *_ in small)
    far_any = any(x[4] for x in big + small)
    # Gaps mostly explained by the HLA region, measured from the two children's imputed dosages there.
    hla = []
    for t, _, _, _, _, pid, gap in big:
        part = hla_part(pid, A, B) if gap else None
        if part is not None and part / gap >= HLA_SHARE:
            hla.append((t, min(part / gap, 1), B if part > 0 else A))
    hla_note = ""
    if hla:
        each = "; ".join(f"{t[0].lower() + t[1:]}: {'nearly all' if share >= .95 else f'about {share * 100:.0f}%'} of "
                         f"the gap, {NAME[hi]}'s HLA copies pushing {poss(hi)} score up" for t, share, hi in hla)
        hla_note = (f" {'These gaps come' if len(hla) > 1 else 'This gap comes'} mostly from the HLA immune-gene region "
                    f"on chromosome 6 ({esc(each)}). Siblings often inherit different HLA blocks (a 1 in 4 chance of "
                    "sharing neither of a parent pair's), which is why HLA-heavy scores can differ a lot.")
    chk = ""
    if check:
        chip_s = "older chip's " if older(B, A) else ""
        verdict = (f"so {poss(B)} {chip_s}imputation is accurate overall and does not explain these gaps."
                   if check[B][1] < IMPUTE_OK_PCT else "so imputation errors may explain part of these gaps.")
        chk = (f" Checked: at {check[B][0]:,} score positions where mom and dad fix the genotype, {esc(NAME[B])}'s "
               f"imputed calls disagree {check[B][1]:.2f}% of the time ({esc(NAME[A])}'s: {check[A][1]:.2f}%), {verdict}")

    notes = (json.load(open(os.path.join(RAW, "family_notes.json"))) if os.path.exists(os.path.join(RAW, "family_notes.json"))
             else {}).get("observed", {})
    pig = {r["person"]: r for r in tsv("pigmentation.tsv")}
    # Family notes are written by self ("mine", "than me"); name them on this two-person page.
    fam = lambda p, k: esc(notes.get(p, {}).get(k, "—").replace("mine", f"{NAME[A]}'s").replace("than me", f"than {NAME[A]}"))
    hair = {p: (pig[p].get("hair_call") or "—").capitalize() for p in (A, B) if p in pig}
    looks = [["Hair color (predicted)", esc(hair[A]), esc(hair[B])],
             ["Hair color (reported)", fam(A, "hair"), fam(B, "hair")],
             ["Hair texture (reported)", fam(A, "hair_texture"), fam(B, "hair_texture")],
             ["Eye color (reported)", fam(A, "eye"), fam(B, "eye")]] if len(hair) == 2 else []
    # Only rows where the two differ and at least one value is known; a shared prediction joins the "Same:" line.
    looks = [r for r in looks if r[1] != r[2] and (r[1], r[2]) != ("—", "—")]
    same_hair = hair[A].lower() if len(hair) == 2 and hair[A] == hair[B] and hair[A] not in ("—", "Not predicted") else ""
    looks_same = same_looks(B, same_hair)
    # How the shared markers fit what the family sees: family.json text, else any reported reddish hair.
    fit = text(f"sibling_diffs.appearance_fit.{B}")
    reds = [p for p in (A, B) if re.search(r"red|auburn|ginger|copper", notes.get(p, {}).get("hair", ""), re.I)]
    if not fit and reds and "MC1R" in looks_same:
        fit = ", which fits " + " and ".join(f"{NAME[p]}'s reported reddish tint" for p in reds)
    looks_txt = f'<p class="dim">Same: {esc(looks_same + fit, quote=False)}.</p>' if looks_same else ""

    stats = []
    same_pct = kids.get("sibling_maternal_same_pct")
    if isinstance(same_pct, dict):  # one value per sibling (vs self)
        same_pct = same_pct.get(f"{A}/{B}", same_pct.get(f"{B}/{A}"))
    elif len(siblings()) > 1:  # a single number is ambiguous with several siblings
        same_pct = None
    if same_pct is None:
        same_pct = (kids.get(B) or {}).get("maternal_same_pct")
    if same_pct is not None:
        stats.append((f"{same_pct:.0f}%", "Same grandparent's copy from mom",
                      "Where you both got grandma's copy, or both grandpa's. Siblings average 50%."))
    if A in kids and B in kids:
        gp = (FAMILY.get("maternal_grandfather") or {}).get("name", "")
        gp = gp.rsplit(" ", 1)[0] if " " in gp else gp  # first name(s), without the surname
        rough = next((f"; {NAME[x]}'s is rougher (older chip)" for x, y in ((B, A), (A, B)) if older(x, y)), "")
        stats.append((f"{kids[A]['grandma_share_pct']:.0f}% · {kids[B]['grandma_share_pct']:.0f}%",
                      f"Mom's half that is grandma's ({NAME[A]} · {NAME[B]})",
                      f"The rest is grandpa{' ' + gp if gp else ''}'s. 40–60% is normal{rough}."))
    stat_html = "".join(f'<div class="stat"><div class="num">{esc(a)}</div><div class="label">{esc(b)}</div>'
                        f'<div class="dim">{esc(c)}</div></div>' for a, b, c in stats)

    # Same outcome, different genotype: one line under the matching table ("Report (outcome; rs1 CT vs CC)").
    med_geno_txt = [f"{g} ({ph}; {da} vs {db})" for g, da, db, ph in med_geno]
    well_geno_txt = [f"{r} ({lab}; " + ", ".join(f"{rs} {x} vs {y}" for rs, x, y in d) + ")" for r, lab, d in well_geno]
    geno_line = lambda items: (f'<p class="dim">Same outcome, different genotype ({esc(NAME[A])} vs {esc(NAME[B])}): '
                               f'{esc("; ".join(items))}.</p>' if items else "")
    sections = []
    if med_rows:
        # The shared results worth pointing out: family.json text, else every shared one with medication guidance.
        flagged = text(f"sibling_diffs.meds_shared.{B}") or (
            f"That includes the shared result{'s' if len(med_flagged) > 1 else ''} with medication guidance, "
            + ", ".join(f"{g} ({DRUGS[g].split(' (')[0]})" if g in DRUGS else g for g in med_flagged)
            + ", which you both inherited." if med_flagged else "")
        sections.append(f'<section class="card"><h2>💊 Medications</h2>{table(["Gene", NAME[A], NAME[B]], med_rows)}'
                        + (f'<p class="dim">Same for both: {esc(", ".join(med_same))}.'
                           f'{" " + esc(flagged, quote=False) if flagged else ""}</p>' if med_same else "")
                        + geno_line(med_geno_txt) + "</section>")
    if real:
        sections.append(f'<section class="card"><h2>🩺 Health, diet &amp; traits</h2>{table(["Report", NAME[A], NAME[B]], real)}'
                        f'<p class="dim">{same_n} other reports give the same result for both of you.</p>'
                        + geno_line(well_geno_txt) + "</section>")
    if pgs_rows:
        sections.append(f'<section class="card"><h2>📊 Genetic tendency scores (rank among ~500 Europeans)</h2>'
                        f'{table(["Score", NAME[A], NAME[B], "Gap"], pgs_rows)}'
                        f'<p class="dim">Shown: scores at least 20 percentile points apart, largest gap first. The gap '
                        f'is in standard deviations; full siblings typically differ by about 1 SD and gaps over 2 SD are '
                        f'uncommon, because each of you got a different random half of mom\'s and dad\'s variants. '
                        f'These are tendencies, not diagnoses.{hla_note}{chk}{far_note(pair) if far_any else ""}</p>'
                        + (f'<p class="dim">Close (under 20 points apart): {small_txt}.</p>' if small_txt else "")
                        + "</section>")
    if looks or looks_txt:
        sections.append(f'<section class="card"><h2>🎨 Appearance</h2>'
                        f'{table(["", NAME[A], NAME[B]], looks) if looks else ""}{looks_txt}</section>')
    va, vb = chip_version(A), chip_version(B)
    anc = ancestry_23(B)
    if anc:
        how = ("; its two reports also use differently structured breakdowns, so small differences mean little"
               if va and vb and va != vb else ", so small differences mean little")
        sections.append(f'<section class="card"><h2>🌍 23andMe ancestry estimate</h2>{table(["", NAME[A], NAME[B]], anc)}'
                        f'<p class="dim">Your family trees are identical, so this gap is partly real (each of you inherited '
                        f'a different mix of the same ancestors) and partly how 23andMe estimates{how}.</p></section>')
    if chip_items:
        who = list(dict.fromkeys(lacking))
        if va and vb and va != vb:  # different chip versions: the marker is missing from one chip
            whose = (f"{NAME[who[0]]}'s {'older ' if older(who[0], A if who[0] == B else B) else ''}chip lacks the marker"
                     if len(who) == 1 else "one of the two chips lacks the marker")
        else:
            whose = (f"{NAME[who[0]]}'s file" if len(who) == 1 else "one of the two files") + \
                " has no call at the marker (a failed read)"
        sections.append(f'<section class="card"><h2>🔍 Looks different, but isn\'t</h2><ul>{"".join(chip_items)}</ul>'
                        f'<p class="dim">These differ on paper only because {esc(whose, quote=False)}.</p>'
                        f'</section>')
    if open_items:
        sections.append(f'<section class="card"><h2>❔ Not comparable yet</h2><ul>{"".join(i for _, i in open_items)}</ul>'
                        f'<p class="dim">A missing call leaves these open: they may or may not differ. '
                        + ("A clinical test, or a newer chip, would settle them." if va and vb and va != vb
                           else "Re-testing, or a clinical test, would settle them.") + '</p></section>')

    doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{pair.capitalize()} Compared</title>
<style>
:root{{--bg:#fbfaf7;--fg:#1d1d1f;--dim:#6b6b70;--line:#e3e1db;--card:#fff;--accent:#2f6fb0}}
@media (prefers-color-scheme:dark){{:root:not([data-theme="light"]){{--bg:#151517;--fg:#ececec;--dim:#9a9aa0;--line:#2c2c30;--card:#1d1d20;--accent:#6aa6e6}}}}
*{{box-sizing:border-box}}
body{{background:var(--bg);color:var(--fg);font:15px/1.5 -apple-system,system-ui,sans-serif;margin:0;padding:28px 16px 48px}}
main{{max-width:980px;margin:0 auto}} h1{{font-size:30px;margin:0;letter-spacing:-.01em}}
h2{{font-size:13px;text-transform:uppercase;letter-spacing:.06em;color:var(--dim);margin:0 0 12px}}
.dim{{color:var(--dim);font-size:13px}} .lede{{margin:4px 0 24px}} a{{color:var(--accent)}}
.stack{{display:grid;gap:16px}}
.card{{background:var(--card);border:1px solid var(--line);border-radius:12px;padding:18px;min-width:0;overflow-wrap:anywhere}}
.stats{{display:grid;gap:16px;grid-template-columns:repeat(auto-fit,minmax(min(260px,100%),1fr))}}
.stat .num{{font-size:24px;font-weight:650}} .stat .label{{font-weight:600;font-size:14px}}
.tablewrap{{overflow-x:auto}}
table{{border-collapse:collapse;width:100%;font-size:14px;font-variant-numeric:tabular-nums}}
th,td{{text-align:left;padding:6px 8px;border-bottom:1px solid var(--line);vertical-align:top}}
thead th{{color:var(--dim);font-weight:600;font-size:12px}} td:first-child{{width:26%}} .warn{{color:#b26b00;font-size:12px;font-weight:600;display:block}}
@media (prefers-color-scheme:dark){{.warn{{color:#e0a040}}}}
ul{{margin:0;padding-left:18px}} li{{margin-bottom:8px}}
</style></head><body><main>
{DISCLAIMER_HTML}
<h1>{esc(NAME[A])} vs {esc(NAME[B])}</h1>
<p class="lede dim">Only the results where the two of you differ; everything not listed is the same. As full {pair} you
share about half your DNA, with each of you getting a different random half of mom's and dad's.
<a href="{esc(page(A))}">{esc(NAME[A])}'s page</a> · <a href="{esc(page(B))}">{esc(NAME[B])}'s page</a></p>
<div class="stack">
{f'<section class="card"><div class="stats">{stat_html}</div></section>' if stat_html else ""}
{"".join(sections)}
</div>
<p class="dim" style="margin-top:24px">For learning and research only — not medical advice or a diagnosis. A clinician
must review anything health-relevant and confirm it with a clinical test.</p>
</main></body></html>"""
    path = os.path.join(RESULTS, f"diffs_{B}.html")
    with open(path, "w") as fh:
        fh.write(doc)
    return path


def main():
    sibs = siblings()
    if not A or not sibs:
        print("sibling_diffs: no siblings configured, skipped")
        return
    for B in sibs:
        if not any(r["person"] == B for r in tsv("wellness.tsv")):  # PGS are optional (need imputed data)
            print(f"sibling_diffs: {B} has no results yet, skipped")
            continue
        path = write_page(B)
        print(path)


if __name__ == "__main__":
    main()
