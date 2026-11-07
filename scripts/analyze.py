"""Family QC, relatedness, trio phasing and maternal grandparent painting.

Usage: python3 -I scripts/analyze.py
Writes to work/ and results/. Uses only measured (not imputed) genotypes.
"""
import gzip
import html
import itertools
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dnalib import (DISCLAIMER_HTML, CHROM_LEN, RESULTS, WORK, children, find, is_autosome, is_het, load_person, people as family,
                    person, role, sex, siblings, text)

EXPECTED = {  # kinship coefficient, by role pair (self and siblings are all "child")
    frozenset(["child", "mother"]): ("parent-child", 0.25),
    frozenset(["child", "father"]): ("parent-child", 0.25),
    frozenset(["mother", "maternal_grandmother"]): ("parent-child", 0.25),
    frozenset(["child"]): ("full siblings", 0.25),
    frozenset(["child", "maternal_grandmother"]): ("grandparent", 0.125),
    frozenset(["mother", "father"]): ("unrelated", 0.0),
    frozenset(["father", "maternal_grandmother"]): ("unrelated", 0.0),
}
GENOTYPE_ERROR = 0.005
SWITCH_PER_BP = 1.5e-8  # ~1.5 cM/Mb, female autosomal average (Kong 2002: ~4,460 cM / ~2,880 Mb, ~42 crossovers/egg)


def hemi(gt):
    return gt[0] if len(gt) == 2 and gt[0] == gt[1] else gt


def expected(a, b):
    """(relationship, kinship) expected for two people from their family.json roles."""
    kind = lambda p: "child" if role(p) in ("self", "sibling") else role(p)
    return EXPECTED.get(frozenset([kind(a), kind(b)]), ("?", None))


def x_y_calls(data):
    x = [r[2] for r in data.values() if r[0] == "X"]
    y = sum(1 for r in data.values() if r[0] == "Y")
    return sum(map(is_het, x)) / max(len(x), 1), y


def is_male(data):
    """Genetic sex from the chip: Y calls and an (almost) homozygous X."""
    x_het, y = x_y_calls(data)
    return y > 500 and x_het < 0.05


def qc(data, stats):
    auto = [r for r in data.values() if is_autosome(r[0])]
    x_het, y = x_y_calls(data)
    return {**stats, "autosomal": len(auto), "x_het_rate": round(x_het, 4),
            "y_calls": y, "sex": "M" if is_male(data) else "F"}


def kinship(a, b):
    """KING-robust kinship over shared autosomal SNPs."""
    het_a = het_b = hethet = ibs0 = n = 0
    for rs in a.keys() & b.keys():
        ca, _, ga = a[rs]
        if not is_autosome(ca):
            continue
        gb = b[rs][2]
        if len(ga) != 2 or len(gb) != 2:
            continue
        n += 1
        ha, hb = is_het(ga), is_het(gb)
        het_a += ha
        het_b += hb
        hethet += ha and hb
        ibs0 += not ha and not hb and ga != gb
    phi = (hethet - 2 * ibs0) / max(1, het_a + het_b)  # max(): no shared heterozygous sites → kinship 0
    return {"snps": n, "kinship": round(phi, 4), "ibs0_pct": round(ibs0 / max(1, n) * 100, 3)}


def king_degree(phi, ibs0_pct):
    """Relationship class from KING degree bounds (Manichaikul 2010).
    Parent-child vs full siblings is split by IBS0, not kinship; the 0.1% cut-off
    is set from this chip data (parent-child pairs sit at genotyping-error level)."""
    if phi > 0.354:
        return "duplicate"
    if phi > 0.177:
        return "parent-child" if ibs0_pct < 0.1 else "full siblings"
    if phi > 0.0884:
        return "2nd degree"
    if phi > 0.0442:
        return "3rd degree"
    return "unrelated"


def mendel_trio(child, mom, dad):
    errors = n = 0
    for rs in child.keys() & mom.keys() & dad.keys():
        ch, _, c = child[rs]
        if not is_autosome(ch) or len(c) != 2:
            continue
        m, d = mom[rs][2], dad[rs][2]
        n += 1
        ok = (c[0] in m and c[1] in d) or (c[1] in m and c[0] in d)
        errors += not ok
    return {"snps": n, "errors": errors, "error_rate_pct": round(errors / max(1, n) * 100, 4)}


def mendel_duo(child, parent):
    errors = n = 0
    for rs in child.keys() & parent.keys():
        ch, _, c = child[rs]
        if not is_autosome(ch):
            continue
        n += 1
        errors += not set(c) & set(parent[rs][2])
    return {"snps": n, "errors": errors, "error_rate_pct": round(errors / max(1, n) * 100, 4)}


def lineage_match(people, chrom):
    """Pairwise concordance on Y or MT SNPs (haploid)."""
    out = {}
    for a, b in itertools.combinations(people, 2):
        da, db = people[a], people[b]
        shared = [rs for rs in da.keys() & db.keys() if da[rs][0] == chrom]
        if not shared:
            continue
        same = sum(hemi(da[rs][2]) == hemi(db[rs][2]) for rs in shared)
        out[f"{a} vs {b}"] = {"snps": len(shared), "concordance_pct": round(same / len(shared) * 100, 2)}
    return out


def maternal_allele(rs, child, mom, dad, male):
    """Allele the child got from mom, or None if not determinable. dad may be {} (untested)."""
    ch, _, c = child[rs]
    m = mom[rs][2]
    if ch == "X" and male:  # a son's X is entirely maternal (a daughter's X is handled like an autosome)
        c = hemi(c)
        return c if len(c) == 1 and c in m else None
    if len(c) != 2:
        return None
    if c[0] == c[1]:
        return c[0] if c[0] in m else None
    if not is_het(m):
        return m[0] if m[0] in c else None
    d = dad.get(rs)
    if d and not is_het(d[2]):
        other = c.replace(d[2][0], "", 1)
        return other if len(other) == 1 and other in m else None
    return None


def grandparent_sites(child, mom, dad, gma):
    """Per chromosome: sorted [(pos, 0=from grandma | 1=from maternal grandpa)]."""
    sites, male = {}, is_male(child)
    for rs in child.keys() & mom.keys() & gma.keys():
        ch, pos, _ = child[rs]
        if ch not in CHROM_LEN:
            continue
        m, g = mom[rs][2], gma[rs][2]
        if not is_het(m) or is_het(g) or len(g) != 2 or g[0] not in m:
            continue
        a = maternal_allele(rs, child, mom, dad, male)
        if a is None:
            continue
        sites.setdefault(ch, []).append((pos, 0 if a == g[0] else 1))
    for v in sites.values():
        v.sort()
    return sites


def viterbi(obs):
    """Two-state HMM smoothing of noisy grandparent-of-origin calls."""
    le, lk = math.log(GENOTYPE_ERROR), math.log(1 - GENOTYPE_ERROR)
    score = [lk if obs[0][1] == s else le for s in (0, 1)]
    back = []
    for (p0, _), (p1, o) in zip(obs, obs[1:]):
        sw = min(0.5, max(1e-12, (p1 - p0) * SWITCH_PER_BP))
        ls, lst = math.log(sw), math.log(1 - sw)
        new, bp = [], []
        for s in (0, 1):
            stay, move = score[s] + lst, score[1 - s] + ls
            bp.append(s if stay >= move else 1 - s)
            new.append(max(stay, move) + (lk if o == s else le))
        score, back = new, back + [bp]
    state = 0 if score[0] >= score[1] else 1
    path = [state]
    for bp in reversed(back):
        state = bp[state]
        path.append(state)
    return path[::-1]


def segments(obs, path):
    segs, start = [], 0
    for i in range(1, len(obs) + 1):
        if i == len(obs) or path[i] != path[start]:
            calls = obs[start:i]
            agree = sum(o == path[start] for _, o in calls)
            segs.append({"start": calls[0][0], "end": calls[-1][0], "source": path[start],
                         "sites": len(calls), "agree_pct": round(agree / len(calls) * 100, 1)})
            start = i
    return segs


def write_merged(data, path):
    with gzip.open(path, "wt") as fh:
        fh.write("rsid\tchromosome\tposition\tgenotype\n")
        for rs, (ch, pos, gt) in sorted(data.items(), key=lambda kv: (kv[1][0], kv[1][1])):
            fh.write(f"{rs}\t{ch}\t{pos}\t{gt}\n")


def chips_note(me):
    """The chips merged for self, e.g. "23andMe + AncestryDNA for me" (family.json text can name the versions)."""
    chips = [n for n, pats in (("23andMe", ("genome_*.zip", "genome_*.txt")),
                               ("AncestryDNA", ("dna-data-*.zip", "AncestryDNA*.txt")))
             if any(find(me, pat) for pat in pats)]
    return text("analyze.chips", " + ".join(chips) + " for me")


def crossover_note(n):
    """Report wording for the maternal crossover count, judged against the ~40–45 expected per egg."""
    if n < 30:
        return (f"{n} maternal crossovers detected, fewer than the typical ~40–45 autosomal per egg "
                "(at chip density short segments, and crossovers near chromosome ends, are missed)")
    if n > 60:
        return (f"{n} maternal crossovers detected, more than the typical ~40–45 autosomal per egg "
                "(genotype errors can split a segment into spurious short ones)")
    return f"{n} maternal crossovers detected (typical: ~40–45 autosomal per egg)"


def shares_line(a, b, chrom, sexes):
    """True if a and b are expected to carry the same Y (father, sons) or mtDNA line (maternal grandmother, mother,
    the children). sexes: id -> "M"/"F" (chip QC, else family.json)."""
    kind = lambda p: "child" if role(p) in ("self", "sibling") else role(p)
    if chrom == "Y":
        return all(sexes.get(p) == "M" and kind(p) in ("child", "father") for p in (a, b))
    return all(kind(p) in ("child", "mother", "maternal_grandmother") for p in (a, b))


def render_report(summary, painting, path):
    me = person("self")
    labels = {0: ("Grandma (maternal)", "#c2417a"), 1: ("Grandpa (maternal, inferred)", "#2f7fc1")}
    rows, bar_w, row_h = [], 760, 34
    scale = bar_w / CHROM_LEN["1"]
    chroms = [c for c in CHROM_LEN if c in painting]
    for i, ch in enumerate(chroms):
        y = 10 + i * row_h
        w = CHROM_LEN[ch] * scale
        rows.append(f'<text x="34" y="{y + 15}" class="lbl" text-anchor="end">{ch}</text>')
        rows.append(f'<rect x="44" y="{y}" width="{w:.1f}" height="10" rx="3" class="pat"/>')
        rows.append(f'<rect x="44" y="{y + 12}" width="{w:.1f}" height="10" rx="3" class="none"/>')
        for s in painting[ch]["segments"]:
            x0, x1 = 44 + s["start"] * scale, 44 + s["end"] * scale
            rows.append(f'<rect x="{x0:.1f}" y="{y + 12}" width="{max(x1 - x0, 1):.1f}" height="10" '
                        f'fill="{labels[s["source"]][1]}"><title>chr{ch}:{s["start"]:,}-{s["end"]:,} '
                        f'{labels[s["source"]][0]} ({s["sites"]} sites)</title></rect>')
    svg_h = 20 + len(chroms) * row_h
    svg = f'<svg viewBox="0 0 820 {svg_h}" role="img" aria-label="Chromosome painting">{"".join(rows)}</svg>'

    def table(d, cols):
        head = "".join(f"<th>{html.escape(c)}</th>" for c in ["", *cols])
        body = "".join("<tr><td>" + html.escape(k) + "</td>" +
                       "".join(f"<td>{html.escape(str(v.get(c, '')))}</td>" for c in cols) + "</tr>"
                       for k, v in d.items())
        return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"

    me_qc = summary["qc"].get(me, {})
    chip_err = ""
    if me_qc.get("shared") and me_qc.get("discordant") is not None:
        disc = me_qc["discordant"] / me_qc["shared"] * 100
        chip_err = (f" My two chips disagree at {disc:.3f}% of shared SNPs (very roughly {disc / 2:.3f}% per chip,"
                    " a rough indication that also includes platform allele-coding differences).")
    painted = ""
    if painting:  # needs mom and her mother
        x_note = ("" if "X" not in painting else
                  "\nX chromosome included: as a son, my single X is entirely from mom." if me_qc.get("sex") == "M" else
                  "\nX chromosome included: as a daughter, the X I got from mom is split like an autosome.")
        painted = f"""
<h2>Maternal grandparent painting</h2>
<p>Each chromosome shows my two copies: <b>top</b> from dad (his parents untested, so unsplit), <b>bottom</b> from mom,
split by which of <i>her</i> parents it came from. {crossover_note(summary["crossovers"]["total"])}.</p>
<div class="card"><div class="key muted"><span style="background:var(--pat)"></span>From dad
<span style="background:#c2417a"></span>Grandma (maternal)<span style="background:#2f7fc1"></span>Grandpa (maternal, inferred)
<span style="background:var(--none)"></span>Not determinable</div>{svg}</div>
<p class="muted">Share of my maternal genome from grandma: {summary["grandma_share_pct"]}% (expected ~50%).{x_note}</p>
"""
    sexes = {p: q.get("sex") or sex(p) for p, q in summary["qc"].items()}
    lines = {}  # Y/mtDNA rows only for pairs expected to share that line (e.g. not dad's mtDNA vs mom's)
    for chrom, key, tag in (("Y", "y", "Y "), ("MT", "mt", "mtDNA ")):
        for pair, v in summary[key].items():
            if shares_line(*pair.split(" vs "), chrom, sexes):
                lines[tag + pair] = v
    has_y = any(k.startswith("Y ") for k in lines)
    has_mt = any(k.startswith("mtDNA ") for k in lines)
    compare = ""
    if len(summary["qc"]) < 2:
        compare = f'\n<p class="muted">Only one person tested: no relatedness, inheritance or lineage comparisons.{chip_err}</p>\n'
    if summary["relatedness"]:
        compare += f"""
<h2>Relatedness</h2><div class="card">{table(summary["relatedness"], ["expected", "expected_kinship", "kinship", "ibs0_pct", "snps", "inferred", "verdict"])}</div>
<p class="muted">KING-robust kinship: 0.25 parent/child or full siblings, 0.125 grandparent, 0 unrelated.
Verdict uses KING degree bounds (Manichaikul 2010: 1st degree 0.177–0.354, 2nd 0.088–0.177, 3rd 0.044–0.088,
unrelated &lt;0.044). Parent/child and full siblings have the same kinship and are told apart by IBS0
(parent/child ≈0; cut-off 0.1%, set from this data).</p>
"""
    if summary["mendel"]:
        compare += f"""
<h2>Inheritance checks</h2><div class="card">{table(summary["mendel"], ["snps", "errors", "error_rate_pct"])}</div>
<p class="muted">Errors are SNPs where a child's genotype can't come from the parent(s); mostly chip errors
(occasionally a small deletion). Most genotyping errors don't produce a conflict, so error_rate_pct is a
Mendelian conflict rate and a lower bound on the true genotype error rate.{chip_err}</p>
"""
    elif len(summary["qc"]) > 1 and chip_err:
        compare += f'\n<p class="muted">{chip_err.strip()}</p>\n'
    if lines:
        title = " and ".join(t for t, on in (("Y", has_y), ("mitochondrial DNA", has_mt)) if on)
        how = "; ".join(t for t, on in (("Y passes father→son", has_y), ("mtDNA passes mother→children", has_mt)) if on)
        compare += f"""
<h2>Lineages ({title})</h2><div class="card">{table(lines, ["snps", "concordance_pct"])}</div>
<p class="muted">{how}. Only pairs expected to share a line are listed; they should match at ~100%.</p>
"""
    doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Family DNA Report</title>
<style>
:root{{--bg:#fbfaf7;--fg:#1d1d1f;--muted:#6b6b70;--line:#e3e1db;--card:#fff;--pat:#b9b7b0;--none:#ecebe6}}
@media (prefers-color-scheme:dark){{:root{{--bg:#151517;--fg:#ececec;--muted:#9a9aa0;--line:#2c2c30;--card:#1d1d20;--pat:#55545a;--none:#2a2a2e}}}}
body{{background:var(--bg);color:var(--fg);font:15px/1.5 -apple-system,system-ui,sans-serif;margin:0;padding:24px 16px}}
main{{max-width:900px;margin:0 auto}} h1{{font-size:26px;margin:0 0 4px}} h2{{font-size:18px;margin:32px 0 8px}}
.muted{{color:var(--muted)}} .card{{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:16px;overflow-x:auto}}
table{{border-collapse:collapse;width:100%;font-size:13px}} th,td{{text-align:left;padding:6px 8px;border-bottom:1px solid var(--line)}}
th{{color:var(--muted);font-weight:600}} svg{{width:100%;height:auto}} .lbl{{font-size:11px;fill:var(--muted)}}
.pat{{fill:var(--pat)}} .none{{fill:var(--none)}} .key span{{display:inline-block;width:12px;height:12px;border-radius:3px;vertical-align:-1px;margin:0 6px 0 14px}}
</style></head><body><main>
{DISCLAIMER_HTML}
<h1>Family DNA report</h1>
<p class="muted">Measured genotypes only ({html.escape(chips_note(me))}). Generated locally.</p>
{painted}
<h2>Sample QC</h2><div class="card">{table(summary["qc"], ["sex", "autosomal", "total", "x_het_rate", "y_calls", "shared", "discordant"])}</div>
{compare}</main></body></html>"""
    with open(path, "w") as fh:
        fh.write(doc)


def main():
    os.makedirs(WORK, exist_ok=True)
    os.makedirs(RESULTS, exist_ok=True)
    people, qcs = {}, {}
    for p in family():
        data, stats = load_person(p)
        if data is None:
            print(f"- {p}: no files, skipped")
            continue
        people[p] = data
        qcs[p] = qc(data, stats)
        print(f"- {p}: {qcs[p]}")
    me, mom, dad, gma = person("self"), person("mother"), person("father"), person("maternal_grandmother")
    if me not in people:
        sys.exit(f"No genotype file for self in raw/{me or '<self>'}/ (see family.json)")
    kids = [c for c in children() if c in people]  # self first
    write_merged(people[me], os.path.join(WORK, f"{me}_merged.tsv.gz"))

    rel = {}
    for a, b in itertools.combinations(people, 2):
        exp, exp_k = expected(a, b)
        k = kinship(people[a], people[b])
        inferred = king_degree(k["kinship"], k["ibs0_pct"])
        ok = exp_k is not None and inferred == {"grandparent": "2nd degree"}.get(exp, exp)
        rel[f"{a} – {b}"] = {"expected": exp, "expected_kinship": exp_k, **k, "inferred": inferred,
                             "verdict": "OK" if ok else "CHECK"}
        print(f"  {a:17s} {b:17s} {k}  {exp} -> {'OK' if ok else 'CHECK'}")

    # Trio checks for every child when both parents are tested; parent duos for self (and for siblings without a trio).
    trio = mom in people and dad in people
    mendel = {}
    if trio:
        mendel[f"trio {me}/{mom}/{dad}"] = mendel_trio(people[me], people[mom], people[dad])
    if mom in people and gma in people:
        mendel[f"duo {mom}/{gma}"] = mendel_duo(people[mom], people[gma])
    for par in (mom, dad):
        for c in kids:
            if par in people and (c == me or not trio):
                mendel[f"duo {c}/{par}"] = mendel_duo(people[c], people[par])
    for c in kids[1:]:
        if trio:
            mendel[f"trio {c}/{mom}/{dad}"] = mendel_trio(people[c], people[mom], people[dad])
    print("  mendel:", mendel)

    def maternal_painting(child, tsv=None):
        """Paint a child's maternal copy as grandma/grandpa segments; optional TSV of the segments."""
        sites = grandparent_sites(people[child], people[mom], people.get(dad, {}), people[gma])
        painting, total_xo, bp_g, bp_all = {}, 0, 0, 0
        fh = open(tsv, "w") if tsv else None
        if fh:
            fh.write("chromosome\tstart\tend\tsource\tsites\tagree_pct\n")
        for ch in CHROM_LEN:
            obs = sites.get(ch)
            if not obs:
                continue
            segs = segments(obs, viterbi(obs))
            painting[ch] = {"segments": segs, "crossovers": len(segs) - 1}
            if ch != "X":
                total_xo += len(segs) - 1
            for s in segs:
                span = s["end"] - s["start"]
                bp_all += span
                bp_g += span if s["source"] == 0 else 0
                if fh:
                    fh.write(f"{ch}\t{s['start']}\t{s['end']}\t{'grandma' if s['source'] == 0 else 'grandpa'}"
                             f"\t{s['sites']}\t{s['agree_pct']}\n")
        if fh:
            fh.close()
        return sites, painting, total_xo, bp_g, bp_all

    def same_pct(pa, pb):
        """Siblings: share of the maternal genome where both got the same grandparent's copy."""
        same = both = 0
        for ch in set(pa) & set(pb):
            if ch == "X":
                continue
            for a in pa[ch]["segments"]:
                for b in pb[ch]["segments"]:
                    ov = min(a["end"], b["end"]) - max(a["start"], b["start"])
                    if ov > 0:
                        both += ov
                        same += ov if a["source"] == b["source"] else 0
        return round(same / both * 100, 1) if both else None

    # Painting needs mom and her mother; dad is optional (he only helps resolve sites where mom and child are both het).
    # A child on an older chip has fewer informative sites (less overlap with the others' chips), so short segments
    # can be missed.
    paints, kids_out = {}, {}
    if mom in people and gma in people:
        for c in kids:
            sites, paints[c], xo, bp_g, bp_all = maternal_painting(c, os.path.join(RESULTS, f"{c}_maternal_segments.tsv"))
            kids_out[c] = {"crossovers": xo, "grandma_share_pct": round(bp_g / bp_all * 100, 1) if bp_all else None,
                           "informative_sites": sum(len(v) for v in sites.values())}
        for c in kids[1:]:
            kids_out[c]["maternal_same_pct"] = same_pct(paints[me], paints[c])  # vs self
        if len(kids) == 2:  # one sibling: a single number
            kids_out["sibling_maternal_same_pct"] = kids_out[kids[1]]["maternal_same_pct"]
        elif len(kids) > 2:  # several siblings: every pair, keyed "<a>/<b>"
            kids_out["sibling_maternal_same_pct"] = {f"{a}/{b}": same_pct(paints[a], paints[b])
                                                     for a, b in itertools.combinations(kids, 2)}
        if len(kids) > 1:
            print("  maternal painting per child:", kids_out)
    else:
        print("  maternal painting skipped: needs mother and maternal grandmother")

    summary = {"qc": qcs, "relatedness": rel, "mendel": mendel,
               "y": lineage_match(people, "Y"), "mt": lineage_match(people, "MT")}
    painting = paints.get(me, {})
    if painting:
        summary.update({"crossovers": {"total": kids_out[me]["crossovers"],
                                       **{ch: v["crossovers"] for ch, v in painting.items()}},
                        "informative_sites": kids_out[me]["informative_sites"],
                        "grandma_share_pct": kids_out[me]["grandma_share_pct"]})
    summary["children"] = kids_out
    with open(os.path.join(RESULTS, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    render_report(summary, painting, os.path.join(RESULTS, "report.html"))
    print("  y:", summary["y"])
    print("  mt:", summary["mt"])
    if painting:
        print(f"  informative sites={summary['informative_sites']} maternal crossovers={summary['crossovers']['total']} "
              f"grandma share={summary['grandma_share_pct']}%")
        print("  per-chrom crossovers:", {ch: v["crossovers"] for ch, v in painting.items()})


if __name__ == "__main__":
    main()
