"""Polygenic scores from 23andMe imputed data (GRCh38), with 1000 Genomes percentiles.

Usage: python3 -I scripts/pgs.py [person ...]   (default: everyone in work/imputed)
Needs: tools/pgs/*_hmPOS_GRCh38.txt.gz (PGS Catalog harmonized), bcftools (Homebrew),
       work/imputed/<person>/*/chr*.bcf (indexed).
Reference genotypes are fetched by remote range query from the 1000 Genomes
high-coverage GRCh38 panel and cached in work/pgs/ref/.
"""
import bisect
import concurrent.futures
import glob
import gzip
import html
import json
import os
import subprocess
import sys
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dnalib  # noqa: E402
from dnalib import BCFTOOLS, RESULTS, TOOLS, WORK, load_person

KG = ("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/"
      "working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.{chrom}.filtered."
      "SNV_INDEL_SV_phased_panel.vcf.gz")
KG_PED = ("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/"
          "20130606_g1k_3202_samples_ped_population.txt")
CHROMS = [f"chr{c}" for c in range(1, 23)]
PGS_DIR = os.path.join(TOOLS, "pgs")
OUT = os.path.join(WORK, "pgs")
SEX_SPECIFIC = {"prostate cancer": "M", "breast cancer": "F"}  # reported only for that sex (lowercase)
# Optional extra 1000G comparison group, family.json "ancestry" -> "pgs_local_pop" (e.g. "GBR"); none by default.
LOCAL_LABELS = {"IBS": "Iberian (Spain)", "TSI": "Tuscan (Italy)", "GBR": "British", "CEU": "Utah CEPH (NW European)",
                "FIN": "Finnish"}
LOCAL = dnalib.FAMILY.get("ancestry", {}).get("pgs_local_pop")
GROUPS = (([(LOCAL, LOCAL_LABELS.get(LOCAL, LOCAL))] if LOCAL else [])
          + [("EUR", "European"), ("AMR", "Admixed American"), ("ALL", "All 1000 Genomes")])


def load_scores():
    scores = {}
    for path in sorted(glob.glob(os.path.join(PGS_DIR, "*_hmPOS_GRCh38.txt.gz"))):
        meta, rows, cols = {}, [], None
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    k, _, v = line[1:].strip().partition("=")
                    meta[k] = v
                    continue
                f = line.rstrip("\n").split("\t")
                if cols is None:
                    cols = {c: i for i, c in enumerate(f)}
                    continue
                get = lambda c: f[cols[c]] if c in cols and cols[c] < len(f) else ""
                chrom, pos = get("hm_chr"), get("hm_pos")
                if not chrom or not pos.isdigit() or chrom not in [c[3:] for c in CHROMS]:
                    continue
                other = get("other_allele") or get("hm_inferOtherAllele")
                if "/" in other:
                    other = ""
                rows.append((f"chr{chrom}", int(pos), get("effect_allele").upper(), other.upper(),
                             float(get("effect_weight"))))
        # published size from the header; rows on chrX or without a GRCh38 harmonized position are dropped above
        published = int(meta.get("variants_number", "") or 0) or len(rows)
        scores[meta["pgs_id"]] = {"trait": meta.get("trait_reported", ""), "name": meta.get("pgs_name", ""),
                                  "variants": rows, "published": published}
    return scores


def write_regions(scores):
    sites = sorted({(v[0], v[1]) for s in scores.values() for v in s["variants"]},
                   key=lambda x: (int(x[0][3:]), x[1]))
    os.makedirs(OUT, exist_ok=True)
    per_chrom = {}
    for chrom, pos in sites:
        per_chrom.setdefault(chrom, []).append(pos)
    for chrom, positions in per_chrom.items():
        with open(os.path.join(OUT, f"regions.{chrom}.tsv"), "w") as fh:
            fh.writelines(f"{chrom}\t{p}\n" for p in positions)
    return per_chrom


def bcf_query(source, regions, fmt, cwd):
    out = subprocess.run([BCFTOOLS, "query", "-R", regions, "-f", fmt, source],
                         capture_output=True, text=True, cwd=cwd)
    if out.returncode != 0:
        raise RuntimeError(out.stderr[-500:])
    return out.stdout


def person_dosages(person, chroms):
    """{(chrom, pos): [(ref, alt, alt_dosage)]} from 23andMe imputed HDS."""
    out = {}
    for chrom in chroms:
        bcf = glob.glob(os.path.join(WORK, "imputed", person, "*", f"{chrom}.bcf"))[0]
        text = bcf_query(bcf, os.path.join(OUT, f"regions.{chrom}.tsv"),
                         "%CHROM\t%POS\t%REF\t%ALT\t[%HDS]\n", OUT)
        for line in text.splitlines():
            ch, pos, ref, alt, hds = line.split("\t")
            dose = sum(float(x) for x in hds.split(",") if x not in (".", ""))
            out.setdefault((ch, int(pos)), []).append((ref, alt, dose))
    return out


def fetch_reference(chroms):
    ref_dir = os.path.join(OUT, "ref")
    os.makedirs(ref_dir, exist_ok=True)

    def one(chrom):
        """Fetch only positions not queried before; append them to the cache."""
        path = os.path.join(ref_dir, f"{chrom}.tsv.gz")
        queried_path = os.path.join(ref_dir, f"{chrom}.queried.txt")
        wanted = [l.split("\t")[1].strip() for l in open(os.path.join(OUT, f"regions.{chrom}.tsv"))]
        queried = set()
        if os.path.exists(queried_path):
            queried = set(open(queried_path).read().split())
        elif os.path.exists(path):  # cache from before queried-tracking: trust the positions it holds
            with gzip.open(path, "rt") as fh:
                queried = {line.split("\t", 2)[1] for line in fh}
        new = [p for p in wanted if p not in queried]
        if not new:
            return chrom, "cached"
        regions = os.path.join(OUT, f"regions.{chrom}.new.tsv")
        with open(regions, "w") as fh:
            fh.writelines(f"{chrom}\t{p}\n" for p in new)
        text = bcf_query(KG.format(chrom=chrom), regions, "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", ref_dir)
        with gzip.open(path, "at") as fh:  # appends a gzip member; gzip.open reads them all
            fh.write(text)
        with open(queried_path, "w") as fh:
            fh.write("\n".join(sorted(queried | set(new), key=int)) + "\n")
        return chrom, f"+{text.count(chr(10))} records for {len(new)} new positions"

    samples_path = os.path.join(ref_dir, "samples.txt")
    if not os.path.exists(samples_path):
        names = subprocess.run([BCFTOOLS, "query", "-l", KG.format(chrom="chr22")], capture_output=True,
                               text=True, cwd=ref_dir, check=True).stdout
        with open(samples_path, "w") as fh:
            fh.write(names)
    ped_path = os.path.join(ref_dir, "ped.txt")
    if not os.path.exists(ped_path):
        urllib.request.urlretrieve(KG_PED, ped_path)
    with concurrent.futures.ThreadPoolExecutor(4) as pool:
        for chrom, status in pool.map(one, chroms):
            print(f"  1000G {chrom}: {status}", flush=True)


def reference_samples():
    ref_dir = os.path.join(OUT, "ref")
    samples = open(os.path.join(ref_dir, "samples.txt")).read().split()
    pop = {}
    with open(os.path.join(ref_dir, "ped.txt")) as fh:
        next(fh)
        for line in fh:
            fam, sid, father, mother, sex, p, superpop = line.split()
            if father == "0" and mother == "0":  # unrelated founders only
                pop[sid] = (p, superpop)
    return samples, pop


def match(cands, effect, other):
    """Pick the record whose alleles match the score variant; return (record, effect_is_alt)."""
    for rec in cands:
        ref, alt = rec[0], rec[1]
        if other and {ref, alt} != {effect, other}:
            continue
        if effect in (ref, alt):
            return rec, effect == alt
    return None, None


def score_all(scores, people):
    """people: {person: dosages}. Uses only variants present for every person and in the reference,
    so everyone (and the reference) is scored on the same set."""
    samples, pop = reference_samples()
    keep = [i for i, s in enumerate(samples) if s in pop]
    ref_sites = {}
    for path in sorted(glob.glob(os.path.join(OUT, "ref", "chr*.tsv.gz"))):
        with gzip.open(path, "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                ref_sites.setdefault((f[0], int(f[1])), []).append((f[2], f[3], f[4:]))
    results = {p: {} for p in people}
    for pid, s in scores.items():
        totals, ref_scores, used = {p: 0.0 for p in people}, [0.0] * len(keep), 0
        for chrom, pos, effect, other, w in s["variants"]:
            r, r_alt = match(ref_sites.get((chrom, pos), []), effect, other)
            hits = {p: match(d.get((chrom, pos), []), effect, other) for p, d in people.items()}
            if r is None or any(m is None for m, _ in hits.values()):
                continue
            used += 1
            for p, (m, m_alt) in hits.items():
                totals[p] += w * (m[2] if m_alt else 2 - m[2])
            gts = r[2]
            for j, i in enumerate(keep):
                g = gts[i]
                n = (g[0] == "1") + (len(g) > 2 and g[2] == "1")
                ref_scores[j] += w * (n if r_alt else 2 - n)
        dist = {}
        for code, _ in GROUPS:
            vals = sorted(v for j, v in enumerate(ref_scores)
                          if code == "ALL" or code in pop[samples[keep[j]]])
            mean = sum(vals) / len(vals)
            sd = (sum((v - mean) ** 2 for v in vals) / (len(vals) - 1)) ** 0.5
            dist[code] = (vals, mean, sd)
        for p, total in totals.items():
            groups = {code: {"n": len(vals), "percentile": round(bisect.bisect_left(vals, total) / len(vals) * 100, 1),
                             "z": round((total - mean) / sd, 2)}
                      for code, (vals, mean, sd) in dist.items()}
            results[p][pid] = {"trait": s["trait"], "variants_total": len(s["variants"]), "variants_used": used,
                               "variants_published": s["published"], "score": total, "groups": groups}
    return results


def render(person, results, path):
    esc = html.escape
    order = sorted(results.items(), key=lambda kv: kv[1]["trait"])
    rows = []
    for pid, r in order:
        published = max(r.get("variants_published") or 0, r["variants_total"])
        cov = r["variants_used"] / published * 100
        dropped = published - r["variants_total"]
        note = f" · {dropped:,} on chrX or not mapped to GRCh38, not scored" if dropped > 0 else ""
        cells = "".join(f"<td><b>{r['groups'][c]['percentile']:.0f}</b><span class=muted> (z {r['groups'][c]['z']:+.2f})</span></td>"
                        for c, _ in GROUPS)
        rows.append(f"<tr><th scope=row>{esc(r['trait'])}<div class=muted>{pid} · {r['variants_used']:,} of "
                    f"{published:,} published variants used ({cov:.0f}%){note}</div></th>{cells}</tr>")
    head = "".join(f"<th>{esc(label)}<div class=muted>n={next(iter(results.values()))['groups'][c]['n']}</div></th>"
                   for c, label in GROUPS)
    local_note = (f" The {esc(LOCAL_LABELS.get(LOCAL, LOCAL).split(' (')[0])} ({esc(LOCAL)}) column changes the comparison "
                  "group, not the score's accuracy." if LOCAL else "")
    doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Polygenic Scores</title>
<style>
:root{{--bg:#fbfaf7;--fg:#1d1d1f;--muted:#6b6b70;--line:#e3e1db;--card:#fff}}
@media (prefers-color-scheme:dark){{:root{{--bg:#151517;--fg:#ececec;--muted:#9a9aa0;--line:#2c2c30;--card:#1d1d20}}}}
body{{background:var(--bg);color:var(--fg);font:15px/1.5 -apple-system,system-ui,sans-serif;margin:0;padding:24px 16px}}
main{{max-width:960px;margin:0 auto}} h1{{font-size:26px;margin:0 0 4px}} h2{{font-size:18px;margin:28px 0 8px}}
.muted{{color:var(--muted);font-size:12px;font-weight:400}} .card{{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:16px;overflow-x:auto}}
table{{border-collapse:collapse;width:100%;font-size:14px}} th,td{{text-align:left;vertical-align:top;padding:8px;border-bottom:1px solid var(--line)}}
thead th{{color:var(--muted)}} .warn{{border-left:4px solid #d08a00;padding:8px 12px;background:var(--card)}}
</style></head><body><main>
<h1>Polygenic scores — {esc(person)}</h1>
<p class="muted">23andMe imputed genotypes (GRCh38) scored with PGS Catalog harmonized weights; percentiles against
unrelated 1000 Genomes high-coverage samples scored on the same variants. Your scores use imputed dosages, which are
slightly less certain than the directly sequenced reference genotypes; this tends to pull percentiles a little toward
the middle, so a high or low percentile, if anything, understates how far from average the score is. Generated locally.</p>
<p class="warn"><b>How to read this.</b> A percentile ranks your genetic predisposition against a reference group;
it is not a diagnosis or a lifetime risk. Most of these scores shift risk modestly (for these scores, people in the top 10%
typically have ~1.3–2.5× the average risk, depending on the score; the cancer scores sit at the higher end),
and lifestyle, family history and clinical measurements matter as much or more.
Scores were built mainly from European-ancestry data (many from UK Biobank). They are less accurate for other
ancestries, and somewhat less accurate even for Southern Europeans, because accuracy falls with genetic distance from
the training group (Privé et al. 2022).{local_note}</p>
<h2>Percentile (z-score) by reference group</h2>
<div class="card"><table><thead><tr><th>Trait</th>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table></div>
</main></body></html>"""
    with open(path, "w") as fh:
        fh.write(doc)


def sex_of(person):
    """From the raw chip (Y calls); family.json "sex" if there is no raw file; None if unknown."""
    data, _ = load_person(person)
    if not data:
        return dnalib.sex(person)
    return "M" if sum(1 for r in data.values() if r[0] == "Y") > 500 else "F"


def main():
    imp = os.path.join(WORK, "imputed")
    people = sys.argv[1:] or [p for p in sorted(os.listdir(imp) if os.path.isdir(imp) else [])
                              if glob.glob(os.path.join(imp, p, "*", "chr*.bcf"))]
    complete = [p for p in people if all(glob.glob(os.path.join(imp, p, "*", f"{c}.bcf")) for c in CHROMS)]
    for p in people:
        if p not in complete:
            print(f"Skipping {p}: incomplete imputed data in work/imputed.")
    people = complete
    if not people:
        print("No one has imputed data in work/imputed; nothing to score.")
        return
    scores = load_scores()
    if not scores:
        print(f"No PGS Catalog score files in {PGS_DIR}; nothing to score.")
        return
    for pid, s in scores.items():
        print(f"{pid}: {s['trait']} ({len(s['variants'])} variants)")
    chroms = list(write_regions(scores))
    print("Reading imputed dosages…", flush=True)
    dosages = {}
    for p in people:
        dosages[p] = person_dosages(p, chroms)
        print(f"  {p}: {len(dosages[p])} score positions", flush=True)
    print("Fetching 1000 Genomes reference genotypes…", flush=True)
    fetch_reference(chroms)
    print("Scoring…", flush=True)
    results = score_all(scores, dosages)
    for p, res in results.items():
        sex = sex_of(p)
        res = {pid: r for pid, r in res.items() if SEX_SPECIFIC.get(r["trait"].lower(), sex) == sex}
        results[p] = res
        with open(os.path.join(RESULTS, f"pgs_{p}.json"), "w") as fh:
            json.dump(res, fh, indent=2)
        render(p, res, os.path.join(RESULTS, f"pgs_{p}.html"))
        print(f"== {p} ({sex or '?'})")
        for pid, r in sorted(res.items(), key=lambda kv: kv[1]["trait"]):
            g = r["groups"]
            print(f"  {r['trait'][:28]:28s} {pid} used {r['variants_used']}/{r['variants_total']}  "
                  + "  ".join(f"{c}:{g[c]['percentile']:5.1f}%" for c, _ in GROUPS))


if __name__ == "__main__":
    main()
