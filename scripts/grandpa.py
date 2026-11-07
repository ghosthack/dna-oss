"""Virtual maternal grandfather: reconstruct the haplotype grandpa transmitted to mom.

Usage: python3 -I scripts/grandpa.py
Writes work/grandpa/grandpa_haplotype.tsv.gz (full per-SNP haplotype, GRCh37),
results/grandpa.json (coverage, validation, markers) and results/grandpa.tsv (markers).

Grandpa never tested. Mom = grandma's haplotype + grandpa's haplotype, so at each SNP:
  - mom homozygous X                    -> grandpa gave X            (method mom_hom)
  - mom X/Y, grandma homozygous X        -> grandpa gave Y            (method duo)
  - mom X/Y, grandma X/Y (or untyped)    -> ambiguous; resolved by
      * mom's tested children (self, siblings): within a maternal-painting segment, the child's
        maternal allele is grandpa's (grandpa segment) or grandma's (grandma segment)
                                                                       (method children)
      * 23andMe "one_parent" phased file for mom (phased with grandma; allele2 =
        non-grandma haplotype)                                          (method phased)
X: grandpa's single X went entirely to mom, so the same logic recovers his whole X.
Only one of grandpa's two autosomal copies is ever recoverable this way.
APOE (and markers missing from mom's chip) use the imputed GRCh38 haplotype dosages, with
the grandpa haplotype identified from nearby grandma-informative sites.
Needs mom and grandma (maternal); dad and the children are optional but resolve more sites.
"""
import bisect
import glob
import gzip
import json
import os
import subprocess
import sys
import zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dnalib import BCFTOOLS, RAW, RESULTS, WORK, children, is_het, load_person, person
import analyze  # read-only reuse of maternal_allele / grandparent painting (no side effects on import)
from apoe import call as apoe_call

OUT = os.path.join(WORK, "grandpa")
# GRCh37 pseudo-autosomal regions on X: excluded from the "grandpa's X" statistics.
PAR = [(60001, 2699520), (154931044, 155260560)]

# rsid: (gene, label, GRCh38 (chrom, pos, ref, alt) for imputed fallback/cross-check, {allele: meaning})
MARKERS = {
    "rs4988235": ("LCT/MCM6", "Lactase persistence", ("chr2", 135851076, "G", "A"),
                  {"A": "lactase persistence (-13910T)", "G": "non-persistence (ancestral)"}),
    "rs12913832": ("HERC2/OCA2", "Eye colour", ("chr15", 28120472, "A", "G"),
                   {"G": "blue-eye associated", "A": "brown-eye associated"}),
    "rs1815739": ("ACTN3", "R577X muscle fibre", ("chr11", 66560624, "C", "T"),
                  {"T": "577X (null; endurance-associated)", "C": "R577 (functional; sprint/power)"}),
    "rs1801133": ("MTHFR", "C677T", ("chr1", 11796321, "G", "A"),
                  {"A": "677T (reduced activity)", "G": "677C (normal)"}),
    "rs4149056": ("SLCO1B1", "*5/*15 (statin myopathy)", ("chr12", 21178615, "T", "C"),
                  {"C": "521C (reduced function)", "T": "521T (normal)"}),
    "rs4244285": ("CYP2C19", "*2", ("chr10", 94781859, "G", "A"),
                  {"A": "*2 (no function)", "G": "not *2"}),
    "rs12248560": ("CYP2C19", "*17", ("chr10", 94761900, "C", "T"),
                   {"T": "*17 (increased function)", "C": "not *17"}),
    "rs1800562": ("HFE", "C282Y", ("chr6", 26092913, "G", "A"),
                  {"A": "C282Y (haemochromatosis risk)", "G": "normal"}),
    "rs1799945": ("HFE", "H63D", ("chr6", 26090951, "C", "G"),
                  {"G": "H63D (low-penetrance HFE variant)", "C": "normal"}),
    "rs1061170": ("CFH", "Y402H (macular degeneration)", ("chr1", 196690107, "T", "C"),
                  {"C": "402H (AMD risk)", "T": "Y402 (normal)"}),
    "rs713598": ("TAS2R38", "Bitter taste (A49P)", ("chr7", 141973545, "C", "G"),
                 {"G": "taster (PAV-associated)", "C": "non-taster (AVI-associated)"}),
}
APOE38 = {"rs429358": ("chr19", 44908684, "T", "C"), "rs7412": ("chr19", 44908822, "C", "T")}
ISOFORM = {("T", "T"): "ε2", ("T", "C"): "ε3", ("C", "C"): "ε4", ("C", "T"): "ε1"}


def in_par(pos):
    return any(a <= pos <= b for a, b in PAR)


def load_one_parent_phase(p):
    """23andMe phased-with-one-parent file: {rsid: (allele1, allele2)} or {}."""
    zips = glob.glob(os.path.join(RAW, p, "phased_genotype_data_*.zip"))
    for zp in zips:
        z = zipfile.ZipFile(zp)
        names = [n for n in z.namelist() if "one_parent" in n]
        if not names:
            continue
        out = {}
        with z.open(names[0]) as fh:
            for line in fh:
                if line[:1] == b"#":
                    continue
                f = line.decode().split()
                if len(f) == 5 and f[3] in "ACGTDI" and f[4] in "ACGTDI" and f[3] != "-":
                    out[f[0]] = (f[3], f[4])
        return out, names[0]
    return {}, None


def painting(child, mom, dad, gma):
    """{chrom: (starts, ends, states)} for one child, using analyze.py's HMM."""
    out = {}
    for ch, obs in analyze.grandparent_sites(child, mom, dad, gma).items():
        segs = analyze.segments(obs, analyze.viterbi(obs))
        out[ch] = ([s["start"] for s in segs], [s["end"] for s in segs], [s["source"] for s in segs])
    return out


def paint_state(paint, ch, pos):
    """0 grandma / 1 grandpa if pos lies inside a segment, else None (gap/crossover zone)."""
    if ch not in paint:
        return None
    starts, ends, states = paint[ch]
    i = bisect.bisect_right(starts, pos) - 1
    return states[i] if i >= 0 and pos <= ends[i] else None


def other(gt, a):
    rest = gt.replace(a, "", 1)
    return rest if len(rest) == 1 else None


def reconstruct(mom, gma, kids, phase):
    """Per-SNP grandpa allele. Returns rows and counters."""
    rows = []
    for rs, (ch, pos, m) in mom.items():
        if not (ch.isdigit() or ch == "X") or len(m) != 2:
            continue
        g = gma[rs][2] if rs in gma else None
        if g is not None and len(g) != 2:
            g = None
        ph = phase.get(rs)
        ph_allele = ph[1] if ph and sorted(ph) == sorted(m) else None
        kid_calls = {}
        for name, (child, dad, paint, male) in kids.items():
            if rs not in child:
                continue
            st = paint_state(paint, ch, pos)
            if st is None:
                continue
            a = analyze.maternal_allele(rs, child, mom, dad, male)
            if a is None or a not in m:
                continue
            kid_calls[name] = a if st == 1 else other(m, a)
        kid_vals = {v for v in kid_calls.values() if v}
        kid_allele = kid_vals.pop() if len(kid_vals) == 1 else None
        kid_conflict = len({v for v in kid_calls.values() if v}) > 1

        allele, method = None, None
        if m[0] == m[1]:
            if g is not None and m[0] not in g:
                method = "mendel_error"
            else:
                allele, method = m[0], "mom_hom"
        else:
            shared = set(m) & set(g) if g else set(m)
            if g is not None and not shared:
                method = "mendel_error"
            elif g is not None and len(shared) == 1 and g[0] == g[1]:
                allele, method = other(m, shared.pop()), "duo"
            elif g is not None and len(shared) == 1:  # tri-allelic: grandma carries only one of mom's alleles
                allele, method = other(m, shared.pop()), "duo"
            elif kid_allele:
                allele, method = kid_allele, "children"
            elif ph_allele:
                allele, method = ph_allele, "phased"
            else:
                method = "kid_conflict" if kid_conflict else "unresolved"
        rows.append((ch, pos, rs, m, g or "", allele or "", method, ph_allele or "",
                     ",".join(f"{k}:{v}" for k, v in sorted(kid_calls.items()) if v),
                     bool(g is not None and is_het(g) and is_het(m)) or (g is None and is_het(m))))
    rows.sort(key=lambda r: (r[0] == "X", int(r[0]) if r[0].isdigit() else 0, r[1]))
    return rows


def pct(a, b):
    return round(100 * a / b, 3) if b else None


def coverage(rows, label):
    n = len(rows)
    known = sum(1 for r in rows if r[5])
    by = {}
    for r in rows:
        by[r[6]] = by.get(r[6], 0) + 1
    het = [r for r in rows if is_het(r[3])]
    return {"region": label, "mom_called_snps": n, "grandpa_allele_known": known, "known_pct": pct(known, n),
            "by_method": dict(sorted(by.items(), key=lambda kv: -kv[1])),
            "mom_het_snps": len(het), "mom_het_known_pct": pct(sum(1 for r in het if r[5]), len(het))}


def validate(rows, kids, mom, gma):
    v = {}
    auto = [r for r in rows if r[0].isdigit()]
    # 1. Mendelian consistency mom vs grandma
    typed = [r for r in rows if r[4]]
    err = sum(1 for r in typed if r[6] == "mendel_error")
    v["mendel_mom_grandma"] = {"snps_both_typed": len(typed), "errors": err, "error_rate_pct": pct(err, len(typed))}
    # 2. Independent check: 23andMe one-parent phasing vs children/painting at ambiguous het/het sites
    for label, sel in (("autosomes", auto), ("X_nonPAR", [r for r in rows if r[0] == "X" and not in_par(r[1])])):
        both = [r for r in sel if r[9] and r[7] and r[8] and len({x.split(":")[1] for x in r[8].split(",")}) == 1]
        agree = sum(1 for r in both if r[7] == r[8].split(",")[0].split(":")[1])  # children all agree (above)
        v[f"phased_vs_children_hethet_{label}"] = {
            "sites": len(both), "agree": agree, "agree_pct": pct(agree, len(both)),
            "note": "independent: 23andMe statistical+grandma phasing vs Mendelian inheritance through child + painting "
                    "(crossover-zone gaps excluded by the painting); children calls are preferred because they come from "
                    "direct inheritance. 100 - agree_pct is a rough upper bound on the error rate of 'phased' calls, "
                    "assuming the children calls are correct; where disagreements fall (switch-error runs vs painting "
                    "errors) is not measured here"}
    # 3. Phased file respects duo constraints (sanity; not independent)
    duo = [r for r in rows if r[6] == "duo" and r[7]]
    v["phased_vs_duo_sites"] = {"sites": len(duo), "agree_pct": pct(sum(r[7] == r[5] for r in duo), len(duo)),
                                "note": "sanity only: 23andMe used grandma's genotypes"}
    # 4. Child maternal allele in grandpa-painted segments vs reconstructed grandpa allele
    for name, (child, dad, paint, male) in kids.items():
        res = {}
        for method_set, lab in (({"duo"}, "duo_sites_circular"), ({"phased"}, "phased_only_sites"),
                                ({"mom_hom"}, "mom_hom_sites")):
            n = agree = 0
            for r in rows:
                if r[6] not in method_set or r[2] not in child or not is_het(r[3]) and lab != "mom_hom_sites":
                    continue
                if paint_state(paint, r[0], r[1]) != 1:
                    continue
                a = analyze.maternal_allele(r[2], child, mom, dad, male)
                if a is None:
                    continue
                n += 1
                agree += a == r[5]
            res[lab] = {"sites": n, "agree_pct": pct(agree, n)}
        res["note"] = ("in segments where this child inherited grandpa's haplotype, the child's maternal allele "
                       "should equal grandpa's; duo sites built the painting (circular), phased_only sites are "
                       "those the child could not resolve itself (expected ~0)")
        v[f"{name}_grandpa_segments"] = res
    # 5. Sibling consistency (sites where two or more children made a call)
    if len(kids) > 1:
        both = [r for r in rows if r[8].count(":") >= 2]
        agree = sum(1 for r in both if len({x.split(":")[1] for x in r[8].split(",")}) == 1)
        v["sibling_concordance"] = {"sites": len(both), "agree_pct": pct(agree, len(both))}
    return v


# ---------- imputed (GRCh38) fallback ----------

def _imputed(person, chrom, start, end):
    files = glob.glob(os.path.join(WORK, "imputed", person, "*", f"{chrom}.bcf"))
    if not files:
        return {}
    out = subprocess.run([BCFTOOLS, "query", "-r", f"{chrom}:{max(1, start)}-{end}",
                          "-f", "%POS\t%REF\t%ALT\t[%HDS]\n", files[0]],
                         capture_output=True, text=True).stdout
    res = {}
    for line in out.splitlines():
        pos, ref, alt, h = line.split("\t")
        hs = h.split(",")
        if len(hs) != 2 or "." in hs:
            continue
        hs = [float(x) for x in hs]
        res[(int(pos), ref, alt)] = None if any(0.1 < x < 0.9 for x in hs) else [int(x >= 0.5) for x in hs]
    return res


def imputed_grandpa(chrom, pos, ref, alt, window=3_000_000, flank=30):
    """Grandpa's transmitted allele at a GRCh38 site from mom's + grandma's imputed haplotypes."""
    mom = _imputed(person("mother"), chrom, pos - window, pos + window)
    gma = _imputed(person("maternal_grandmother"), chrom, pos - window, pos + window)
    key = (pos, ref, alt)
    if key not in mom and (pos, alt, ref) in mom:  # catalogue REF/ALT order may differ from the BCF's
        key, ref, alt = (pos, alt, ref), alt, ref
    if key not in mom:
        return {"status": "site not in imputed data (or REF/ALT mismatch)"}
    mh = mom[key]
    if mh is None:
        return {"status": "mom's imputed haplotypes uncertain"}
    alleles = [alt if x else ref for x in mh]
    if mh[0] == mh[1]:
        return {"status": "ok", "allele": alleles[0], "how": "mom homozygous (imputed)", "mom_haps": alleles}
    gh = gma.get(key)
    if gh and gh[0] == gh[1]:
        g = alt if gh[0] else ref
        return {"status": "ok", "allele": alleles[1 - alleles.index(g)] if g in alleles else None,
                "how": "grandma homozygous (imputed)", "mom_haps": alleles}
    # orientation: which of mom's imputed haplotypes is grandma's, from nearby informative sites
    inf = []
    for k, h in mom.items():
        g = gma.get(k)
        if h is None or g is None or h[0] == h[1] or g[0] != g[1]:
            continue
        inf.append((k[0], 0 if h[0] == g[0] else 1))
    inf.sort()
    left = [s for p, s in inf if p < pos][-flank:]
    right = [s for p, s in inf if p > pos][:flank]
    if len(left) < flank // 2 or len(right) < flank // 2:
        return {"status": "too few informative imputed sites nearby"}
    fl, fr = sum(left) / len(left), sum(right) / len(right)
    if fl <= 0.1 and fr <= 0.1:
        gm_hap = 0
    elif fl >= 0.9 and fr >= 0.9:
        gm_hap = 1
    else:
        return {"status": "phase orientation inconsistent around site", "left_frac_hap1": fl, "right_frac_hap1": fr}
    return {"status": "ok", "allele": alleles[1 - gm_hap], "mom_haps": alleles,
            "how": f"imputed phase; grandma = hap{gm_hap} at {len(left)}+{len(right)} flanking informative sites "
                   f"({round(100 * max(1 - fl, fl))}%/{round(100 * max(1 - fr, fr))}% consistent)"}


def apoe(kids, mom_chip):
    res = {}
    calls = {rs: imputed_grandpa(*spec) for rs, spec in APOE38.items()}
    out = {"per_snp": calls, "mom_chip_missing": [rs for rs in APOE38 if rs not in mom_chip]}
    if all(c.get("status") == "ok" and c.get("allele") for c in calls.values()):
        hap = (calls["rs429358"]["allele"], calls["rs7412"]["allele"])
        out["grandpa_transmitted"] = ISOFORM.get(hap)
        moms = list(zip(calls["rs429358"]["mom_haps"], calls["rs7412"]["mom_haps"]))
        out["mom_imputed"] = "/".join(sorted(ISOFORM.get(h, "?") for h in moms))
    # cross-check via children: child's maternal APOE isoform + painting state at chr19:45.41 Mb (GRCh37)
    for name, (child, dad, paint, _) in kids.items():
        st = paint_state(paint, "19", 45411941)
        g1, g2 = (child.get(rs, (0, 0, ""))[2] for rs in ("rs429358", "rs7412"))
        d1, d2 = (dad.get(rs, (0, 0, ""))[2] for rs in ("rs429358", "rs7412"))
        mat = None
        if g1 and g2 and d1 and d2 and not is_het(d1) and not is_het(d2):
            m1, m2 = other(g1, d1[0]), other(g2, d2[0])
            if m1 and m2:
                mat = ISOFORM.get((m1, m2))
        res[name] = {"painting_at_APOE": {0: "grandma", 1: "grandpa", None: "unknown"}[st],
                     "child_maternal_isoform": mat}
        isos = (out.get("mom_imputed") or "").split("/")
        # The two APOE SNPs are resolved separately, so grandpa's combination may not be one of mom's isoforms.
        if mat and st is not None and out.get("grandpa_transmitted") in isos:
            isos.remove(out["grandpa_transmitted"])
            res[name]["consistent"] = mat == (out["grandpa_transmitted"] if st == 1 else isos[0])
    out["children_check"] = res
    return out


def describe(allele, meanings):
    return meanings.get(allele, "?")


def apoe_notes(ap):
    """(source note, meaning note) for the APOE row: whether mom's chip lacks the SNPs, and which SNP(s) decide
    grandpa's isoform (those where mom's imputed haplotypes differ)."""
    no_call = bool(ap.get("mom_chip_missing"))
    drivers = [rs for rs, c in ap["per_snp"].items()
               if c.get("status") == "ok" and not c["how"].startswith("mom homozygous")]
    why = ["mom's chip no-call"] * no_call + [f"driven by {' + '.join(drivers)}"] * bool(
        drivers and ap.get("grandpa_transmitted"))
    return ("imputed (mom chip no-call)" if no_call else "imputed",
            " — imputed only" + (f" ({'; '.join(why)})" if why else ""))


def main():
    mom_id, gma_id = person("mother"), person("maternal_grandmother")
    mom, _ = load_person(mom_id) if mom_id else (None, None)
    gma, _ = load_person(gma_id) if gma_id else (None, None)
    if mom is None or gma is None:
        print("- virtual grandpa skipped: needs mom's and grandma's (maternal) genotypes")
        return
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(RESULTS, exist_ok=True)
    dad = (load_person(person("father"))[0] if person("father") else None) or {}
    kids = {}
    for name in children():
        if not glob.glob(os.path.join(RAW, name, "*")):
            print(f"- {name}: no files, skipped")
            continue
        child, _ = load_person(name)
        if child is None:
            continue
        kids[name] = (child, dad, painting(child, mom, dad, gma), analyze.is_male(child))
        print(f"- {name}: painted {sum(len(p[0]) for p in kids[name][2].values())} maternal segments")
    phase, phase_file = load_one_parent_phase(mom_id)
    print(f"- mom one-parent phase file: {phase_file} ({len(phase)} called sites)")

    rows = reconstruct(mom, gma, kids, phase)
    hap_path = os.path.join(OUT, "grandpa_haplotype.tsv.gz")
    with gzip.open(hap_path, "wt") as fh:
        fh.write("chromosome\tposition_grch37\trsid\tmom\tgrandma\tgrandpa_transmitted\tmethod\t"
                 "phased_file_allele\tchildren_allele\tambiguous_hethet\n")
        for r in rows:
            fh.write("\t".join(map(str, r[:9])) + f"\t{int(r[9])}\n")

    auto = [r for r in rows if r[0].isdigit()]
    xnp = [r for r in rows if r[0] == "X" and not in_par(r[1])]
    cov = {"autosomes": coverage(auto, "autosomes"), "X_nonPAR": coverage(xnp, "X (non-PAR)")}
    val = validate(rows, kids, mom, gma)

    # notable markers
    by_rs = {r[2]: r for r in rows}
    table = []
    for rs, (gene, label, g38, meanings) in MARKERS.items():
        r = by_rs.get(rs)
        imp = imputed_grandpa(*g38)
        entry = {"rsid": rs, "gene": gene, "marker": label}
        if r and r[5]:
            entry.update(allele=r[5], source=f"chip ({r[6]})", mom=r[3], grandma=r[4] or None)
        elif imp.get("status") == "ok" and imp.get("allele"):
            entry.update(allele=imp["allele"], source=f"imputed ({imp['how']})", mom=None, grandma=None)
        else:
            entry.update(allele=None, source="not determinable", mom=r[3] if r else None)
        entry["imputed_crosscheck"] = imp.get("allele") if imp.get("status") == "ok" else imp.get("status")
        entry["meaning"] = describe(entry["allele"], meanings) if entry["allele"] else ""
        table.append(entry)

    ap = apoe(kids, mom)
    gma_apoe, gma_how = apoe_call(gma_id)
    src_note, meaning_note = apoe_notes(ap)
    table.insert(0, {"rsid": "rs429358+rs7412", "gene": "APOE", "marker": "ε2/ε3/ε4",
                     "allele": ap.get("grandpa_transmitted"), "source": f"{src_note}; "
                     + "; ".join(f"{k}: {v.get('how', v.get('status'))}" for k, v in ap["per_snp"].items()),
                     "mom": ap.get("mom_imputed"), "grandma": f"{gma_apoe} ({gma_how})" if gma_apoe else None,
                     "meaning": {"ε2": "ε2 (lower Alzheimer's risk)", "ε3": "ε3 (neutral)",
                                 "ε4": "ε4 (higher Alzheimer's risk)"}.get(ap.get("grandpa_transmitted"), "")
                                + meaning_note,
                     "imputed_crosscheck": ap.get("grandpa_transmitted")})

    # ABO: rs8176719 (23andMe D = 261delG = O allele; I = functional A/B), rs8176746 (GRCh37 plus strand T = B)
    r19, r46 = by_rs.get("rs8176719"), by_rs.get("rs8176746")
    abo = {"rsid": "rs8176719+rs8176746", "gene": "ABO", "marker": "blood-group allele", "allele": None,
           "source": "chip", "mom": None, "grandma": None, "meaning": "not determinable", "imputed_crosscheck": None}
    if r19 and r19[5] and r46 and r46[5]:
        o = r19[5] == "D"
        b = r46[5] == "T"
        abo.update(allele=f"{r19[5]}/{r46[5]}", source=f"chip ({r19[6]}, {r46[6]})", mom=f"{r19[3]} {r46[3]}",
                   grandma=f"{r19[4]} {r46[4]}",
                   meaning="O allele" if o else ("B allele" if b else
                                                 "A allele (A1 vs A2 and rare O2/O03 variants not resolvable)"))
    table.append(abo)

    tsv = os.path.join(RESULTS, "grandpa.tsv")
    with open(tsv, "w") as fh:
        fh.write("gene\tmarker\trsid\tgrandpa_transmitted_allele\tmeaning\tsource\timputed_crosscheck\n")
        for e in table:
            fh.write(f"{e['gene']}\t{e['marker']}\t{e['rsid']}\t{e['allele'] or ''}\t{e['meaning']}\t"
                     f"{e['source']}\t{e['imputed_crosscheck'] or ''}\n")
    summary = {
        "description": "Virtual maternal grandfather: the haplotype grandpa transmitted to mom, reconstructed from "
                       "mom + grandma (Mendelian), mom's children via the maternal grandparent painting, and 23andMe's "
                       "one-parent phasing of mom. Each allele is ONE of grandpa's two copies: he carried at least one "
                       "copy; his other (untransmitted) copy is unknown. X: his only X, fully transmitted to mom.",
        "inputs": {"children_used": list(kids), "phase_file": phase_file},
        "haplotype_file": os.path.relpath(hap_path, os.path.dirname(WORK)),
        "coverage": cov, "validation": val, "apoe_detail": ap, "markers": table,
        "caveats": [
            "Only grandpa's transmitted haplotype (to mom) is recoverable; his second copy is unknown, so "
            "'grandpa carried allele X' means at least one copy, not his genotype.",
            "Ambiguous het/het sites rely on children painting (crossover-zone gaps excluded) or on 23andMe "
            "statistical phasing; both can err near crossovers / switch errors.",
            "Imputed calls (APOE) are statistical; orientation depends on flanking grandma-informative sites.",
            "Chip genotyping errors propagate directly into single-site calls (rare-variant calls are least reliable). "
            "Mom-grandma Mendelian conflicts are only a lower bound on that error, since most misreads don't "
            "create a conflict.",
        ],
    }
    with open(os.path.join(RESULTS, "grandpa.json"), "w") as fh:
        json.dump(summary, fh, indent=2, ensure_ascii=False)

    for k, c in cov.items():
        print(f"  coverage {k}: {c['known_pct']}% of {c['mom_called_snps']} ({c['by_method']})")
    print("  validation:", json.dumps(val, ensure_ascii=False))
    print("  APOE:", json.dumps({k: v for k, v in ap.items() if k != 'per_snp'}, ensure_ascii=False))
    for e in table:
        print(f"  {e['gene']:10s} {e['rsid']:20s} {str(e['allele']):6s} {e['meaning']:45s} {e['source'][:60]} | imp={e['imputed_crosscheck']}")


if __name__ == "__main__":
    main()
