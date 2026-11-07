"""Relative Neanderthal ancestry from 23andMe imputed data (GRCh38), vs 1000 Genomes.

Method
  Browning et al. 2018 (Cell 173:53) ran Sprime on 1000 Genomes phase 3 and published, per
  population, the introgressed segments it detected, with the putatively archaic allele at each
  site and whether the Altai Neanderthal genome carries it (NMATCH) and the Altai Denisovan (DMATCH).
  We use the Iberian (IBS) calls:
    1. keep segments that look Neanderthal (Neanderthal match rate >= MIN_NMATCH over >= MIN_COMP
       comparable sites and higher than the Denisovan match rate);
    2. inside each segment, take up to TAGS_PER_SEG SNVs where the archaic allele matches the Altai
       Neanderthal, lifted from GRCh37 to GRCh38 (UCSC chain) and present in the imputed data;
    3. per person and per 1000 Genomes sample, estimate the number of archaic copies of each segment
       (mean archaic-allele dosage over its tags, 0-2), weight by segment length and sum;
    4. report percentile / z / ratio vs unrelated 1000 Genomes IBS and EUR samples scored the same way,
       and vs Europeans outside Spain (EUR_noIBS = CEU/GBR/TSI/FIN). The IBS reference samples are the
       same people the segments were discovered in, so their own low-frequency segments are counted
       while an outside person's private segments are missed; scores vs IBS (and slightly vs EUR) are
       therefore biased low for the family. EUR_noIBS is out-of-sample like the family, but because
       the segments were found in Iberians, Southern European ancestry may score somewhat higher
       against it; the truth likely lies between the IBS and EUR_noIBS comparisons.
  The result is "Sprime-detectable Neanderthal sequence relative to Europeans", not an absolute
  Neanderthal percentage.

Usage: python3 -I scripts/neanderthal.py [person ...]   (default: everyone in work/imputed)
Needs: work/neanderthal/sprime/mendeley_data/IBS.chr*.ND_match (Mendeley y7hyt83vxr),
       work/neanderthal/hg19ToHg38.over.chain.gz (UCSC), bcftools, imputed BCFs.
Reference genotypes are cached incrementally in work/neanderthal/ref/.
"""
import bisect
import concurrent.futures
import glob
import gzip
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dnalib import BCFTOOLS, KG_PANEL as KG, RESULTS, WORK, append_gzip, kg_sample_files, write_atomic

OUT = os.path.join(WORK, "neanderthal")
SPRIME_GLOB = os.path.join(OUT, "sprime", "mendeley_data", "IBS.chr*.ND_match")
CHAIN = os.path.join(OUT, "hg19ToHg38.over.chain.gz")
REF_DIR = os.path.join(OUT, "ref")
POP = "IBS"
MIN_NMATCH, MIN_COMP, TAGS_PER_SEG = 0.6, 5, 5
GENOME_AUTOSOMAL_BP = 2_881_033_286  # GRCh37 autosomes (sum of chr1-22)
COMP = str.maketrans("ACGT", "TGCA")
GROUPS = [("IBS", "Iberian (Spain)"), ("EUR", "European"),
          ("EUR_noIBS", "European outside Spain (CEU/GBR/TSI/FIN)"), ("AMR", "Admixed American"),
          ("EAS", "East Asian"), ("SAS", "South Asian"), ("AFR", "African")]


# ---------- liftover (UCSC chain) ----------
def load_chain(path):
    """{chrom: sorted [(tStart, tEnd, qChrom, qStart, qStrand, qSize)]} for ungapped blocks."""
    blocks = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.split()
            if not f:
                continue
            if f[0] == "chain":
                t_chrom, t_pos = f[2], int(f[5])
                q_chrom, q_size, q_strand, q_pos = f[7], int(f[8]), f[9], int(f[10])
                continue
            size = int(f[0])
            blocks.setdefault(t_chrom, []).append((t_pos, t_pos + size, q_chrom, q_pos, q_strand, q_size))
            if len(f) == 3:
                t_pos += size + int(f[1])
                q_pos += size + int(f[2])
    for v in blocks.values():
        v.sort()
    starts = {c: [b[0] for b in v] for c, v in blocks.items()}
    return blocks, starts


def lift(chain, chrom, pos1):
    """GRCh37 1-based position -> (chrom, pos1, strand) on GRCh38, or None."""
    blocks, starts = chain
    v = blocks.get(chrom)
    if not v:
        return None
    p0 = pos1 - 1
    i = bisect.bisect_right(starts[chrom], p0) - 1
    for j in (i, i - 1):  # tolerate rare overlaps between chains
        if j >= 0 and v[j][0] <= p0 < v[j][1]:
            t0, _, qc, q0, strand, qsize = v[j]
            q = q0 + (p0 - t0)
            if strand == "-":
                q = qsize - q - 1
            return qc, q + 1, strand
    return None


# ---------- Sprime segments and tags ----------
def load_sprime():
    segs = {}
    for path in glob.glob(SPRIME_GLOB):
        with open(path) as fh:
            for line in fh:
                f = line.split()
                if f[0] == "CHROM":
                    continue
                chrom, pos, rsid, ref, alt, seg, allele, score, nm, dm = f[:10]
                s = segs.setdefault((chrom, seg), {"chrom": chrom, "sites": [], "n": [0, 0], "d": [0, 0]})
                s["sites"].append((int(pos), ref, alt, int(allele), nm))
                if nm in ("match", "mismatch"):
                    s["n"][nm == "match"] += 1
                if dm in ("match", "mismatch"):
                    s["d"][dm == "match"] += 1
    return segs


def neanderthal_segments(segs):
    keep = {}
    for key, s in segs.items():
        n_comp, d_comp = sum(s["n"]), sum(s["d"])
        n_rate = s["n"][1] / n_comp if n_comp else 0
        d_rate = s["d"][1] / d_comp if d_comp else 0
        if n_comp >= MIN_COMP and n_rate >= MIN_NMATCH and n_rate > d_rate:
            pos = [x[0] for x in s["sites"]]
            s.update(start=min(pos), end=max(pos), length=max(pos) - min(pos) + 1,
                     n_rate=round(n_rate, 3), d_rate=round(d_rate, 3))
            keep[key] = s
    return keep


def candidate_tags(segs, chain):
    """{seg_key: [(chr38, pos38, archaic_allele_38, other_allele_38)]} for SNVs matching Altai."""
    out = {}
    for key, s in segs.items():
        for pos, ref, alt, allele, nm in s["sites"]:
            if nm != "match" or len(ref) != 1 or len(alt) != 1:
                continue
            if {ref, alt} in ({"A", "T"}, {"C", "G"}):  # strand-ambiguous; skip for safety
                continue
            hit = lift(chain, "chr" + s["chrom"], pos)
            if not hit or hit[0] != "chr" + s["chrom"]:
                continue
            arch, other = (alt, ref) if allele == 1 else (ref, alt)
            if hit[2] == "-":
                arch, other = arch.translate(COMP), other.translate(COMP)
            out.setdefault(key, []).append((hit[0], hit[1], arch, other))
    return out


# ---------- genotype access ----------
def bcf_query(source, regions, fmt, cwd):
    out = subprocess.run([BCFTOOLS, "query", "-R", regions, "-f", fmt, source],
                         capture_output=True, text=True, cwd=cwd)
    if out.returncode != 0:
        raise RuntimeError(out.stderr[-500:])
    return out.stdout


def write_regions(sites, name):
    per = {}
    for chrom, pos in sorted(set(sites), key=lambda x: (int(x[0][3:]), x[1])):
        per.setdefault(chrom, []).append(pos)
    os.makedirs(OUT, exist_ok=True)
    for chrom, positions in per.items():
        with open(os.path.join(OUT, f"{name}.{chrom}.tsv"), "w") as fh:
            fh.writelines(f"{chrom}\t{p}\n" for p in positions)
    return per


def person_dosages(person, chroms, name):
    """{(chrom, pos): [(ref, alt, alt_dosage)]} from imputed HDS."""
    out = {}
    for chrom in chroms:
        bcf = glob.glob(os.path.join(WORK, "imputed", person, "*", f"{chrom}.bcf"))[0]
        text = bcf_query(bcf, os.path.join(OUT, f"{name}.{chrom}.tsv"), "%CHROM\t%POS\t%REF\t%ALT\t[%HDS]\n", OUT)
        for line in text.splitlines():
            ch, pos, ref, alt, hds = line.split("\t")
            dose = sum(float(x) for x in hds.split(",") if x not in (".", ""))
            out.setdefault((ch, int(pos)), []).append((ref, alt, dose))
    return out


def ref_dosages(gts, keep, count_alt):
    """Archaic-allele dosages of the reference samples `keep` from 1000G GTs ("0|1"), or None if any is missing."""
    out = []
    for i in keep:
        g = gts[i].replace("|", "").replace("/", "")
        if len(g) != 2 or not set(g) <= {"0", "1"}:
            return None
        n = g.count("1")
        out.append(n if count_alt else 2 - n)
    return out


def fetch_reference(chroms, name):
    os.makedirs(REF_DIR, exist_ok=True)

    def one(chrom):
        path = os.path.join(REF_DIR, f"{chrom}.tsv.gz")
        queried_path = os.path.join(REF_DIR, f"{chrom}.queried.txt")
        wanted = [l.split("\t")[1].strip() for l in open(os.path.join(OUT, f"{name}.{chrom}.tsv"))]
        queried = set(open(queried_path).read().split()) if os.path.exists(queried_path) else set()
        new = [p for p in wanted if p not in queried]
        if not new:
            return chrom, "cached"
        regions = os.path.join(REF_DIR, f"new.{chrom}.tsv")
        with open(regions, "w") as fh:
            fh.writelines(f"{chrom}\t{p}\n" for p in new)
        for attempt in range(3):  # remote reads occasionally fail mid-stream; retry the chromosome
            try:
                text = bcf_query(KG.format(chrom=chrom), regions, "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", REF_DIR)
                break
            except RuntimeError:
                if attempt == 2:
                    raise
        append_gzip(path, text)
        write_atomic(queried_path, "\n".join(sorted(queried | set(new), key=int)) + "\n")
        return chrom, f"+{text.count(chr(10))} records for {len(new)} new positions"

    kg_sample_files(REF_DIR)
    with concurrent.futures.ThreadPoolExecutor(4) as pool:
        for chrom, status in pool.map(one, chroms):
            print(f"  1000G {chrom}: {status}", flush=True)


def load_reference():
    samples = open(os.path.join(REF_DIR, "samples.txt")).read().split()
    pop = {}
    with open(os.path.join(REF_DIR, "ped.txt")) as fh:
        next(fh)
        for line in fh:
            fam, sid, father, mother, sex, p, superpop = line.split()
            if father == "0" and mother == "0":
                pop[sid] = (p, superpop)
    sites = {}
    for path in glob.glob(os.path.join(REF_DIR, "chr*.tsv.gz")):
        with gzip.open(path, "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                sites.setdefault((f[0], int(f[1])), []).append((f[2], f[3], f[4:]))
    keep = [i for i, s in enumerate(samples) if s in pop]
    return samples, pop, keep, sites


def match(cands, arch, other):
    for rec in cands:
        if {rec[0], rec[1]} == {arch, other}:
            return rec, arch == rec[1]
    return None, None


def pick_even(items, k):
    if len(items) <= k:
        return items
    step = (len(items) - 1) / (k - 1)
    return [items[round(i * step)] for i in range(k)]


def stats(vals, x):
    vals = sorted(vals)
    mean = sum(vals) / len(vals)
    sd = (sum((v - mean) ** 2 for v in vals) / (len(vals) - 1)) ** 0.5
    return {"n": len(vals), "mean": mean, "sd": sd,
            "percentile": round(bisect.bisect_left(vals, x) / len(vals) * 100, 1),
            "z": round((x - mean) / sd, 2), "ratio_to_mean": round(x / mean, 3)}


def imputed_people():
    """Everyone in work/imputed with an imputed BCF for every autosome."""
    imp = os.path.join(WORK, "imputed")
    return [p for p in sorted(os.listdir(imp) if os.path.isdir(imp) else [])
            if all(glob.glob(os.path.join(imp, p, "*", f"chr{c}.bcf")) for c in range(1, 23))]


def main():
    people = sys.argv[1:] or imputed_people()
    if not people:
        print("No one has imputed data in work/imputed; nothing to score.")
        return
    chain = load_chain(CHAIN)
    all_segs = load_sprime()
    segs = neanderthal_segments(all_segs)
    print(f"Sprime {POP}: {len(all_segs)} segments, {len(segs)} Neanderthal-like", flush=True)
    cands = candidate_tags(segs, chain)
    # keep candidates the imputation covers with matching alleles (same site set for everyone)
    cand_sites = [(c, p) for v in cands.values() for c, p, _, _ in v]
    chroms = list(write_regions(cand_sites, "cand"))
    first = person_dosages(people[0], chroms, "cand")
    tags = {}
    for key, v in cands.items():
        ok = [t for t in sorted(v, key=lambda t: t[1]) if match(first.get((t[0], t[1]), []), t[2], t[3])[0]]
        if ok:
            tags[key] = pick_even(ok, TAGS_PER_SEG)
    tag_sites = [(c, p) for v in tags.values() for c, p, _, _ in v]
    print(f"Candidate Altai-matching SNVs: {len(cand_sites)}; segments with imputed tags: {len(tags)}; "
          f"tags: {len(tag_sites)}", flush=True)
    chroms = list(write_regions(tag_sites, "tags"))
    dosages = {p: person_dosages(p, chroms, "tags") for p in people}
    print("Fetching 1000 Genomes reference genotypes…", flush=True)
    fetch_reference(chroms, "tags")
    samples, pop, keep, ref_sites = load_reference()

    # score: sum over segments of length * mean archaic dosage over tags (same tags for everyone)
    person_score = {p: 0.0 for p in people}
    ref_score = [0.0] * len(keep)
    person_alleles = {p: 0.0 for p in people}
    ref_alleles = [0.0] * len(keep)
    used_segs = used_tags = 0
    bp_used = 0
    for key, tlist in tags.items():
        rows = []
        for chrom, pos, arch, other in tlist:
            r, r_alt = match(ref_sites.get((chrom, pos), []), arch, other)
            hits = {p: match(d.get((chrom, pos), []), arch, other) for p, d in dosages.items()}
            if r is None or any(m is None for m, _ in hits.values()):
                continue
            pd = {p: (m[2] if a else 2 - m[2]) for p, (m, a) in hits.items()}
            rd = ref_dosages(r[2], keep, r_alt)
            if rd is not None:  # None: a reference genotype is missing, so the tag is skipped for everyone
                rows.append((pd, rd))
        if not rows:
            continue
        used_segs += 1
        used_tags += len(rows)
        length = segs[key]["length"]
        bp_used += length
        for p in people:
            m = sum(r[0][p] for r in rows) / len(rows)
            person_score[p] += length * m
            person_alleles[p] += sum(r[0][p] for r in rows)
        for j in range(len(keep)):
            m = sum(r[1][j] for r in rows) / len(rows)
            ref_score[j] += length * m
            ref_alleles[j] += sum(r[1][j] for r in rows)
    denom = 2 * GENOME_AUTOSOMAL_BP / 100  # -> percent of diploid autosomal genome

    def in_group(pp, code):
        if code == "EUR_noIBS":           # out-of-sample Europeans: not part of Sprime discovery
            return pp[1] == "EUR" and pp[0] != "IBS"
        return code in pp

    def group_vals(vals, code):
        return [v for j, v in enumerate(vals) if in_group(pop[samples[keep[j]]], code)]

    ref_means = {code: round(sum(group_vals(ref_score, code)) / len(group_vals(ref_score, code)) / denom, 3)
                 for code, _ in GROUPS}
    results = {
        "method": ("Sprime (Browning et al. 2018) Iberian introgressed segments classified as Neanderthal "
                   "via Altai Neanderthal/Denisovan allele matching; per-segment archaic copy number "
                   "estimated from imputed dosages at Altai-matching tag SNVs; length-weighted sum compared "
                   "with unrelated 1000 Genomes high-coverage samples scored on the same tags."),
        "sources": {
            "sprime": "Browning SR et al. 2018, Cell 173:53-61, doi:10.1016/j.cell.2018.02.031; "
                      "data: Mendeley Data doi:10.17632/y7hyt83vxr.1 (IBS_sprime_results)",
            "liftover": "UCSC hg19ToHg38.over.chain.gz",
            "reference": "1000 Genomes high-coverage phased panel (Byrska-Bishop et al. 2022), unrelated samples",
        },
        "parameters": {"population_calls": POP, "min_neanderthal_match_rate": MIN_NMATCH,
                       "min_comparable_sites": MIN_COMP, "tags_per_segment": TAGS_PER_SEG},
        "segments_sprime": len(all_segs), "segments_neanderthal_like": len(segs),
        "segments_used": used_segs, "tags_used": used_tags,
        "segment_bp_used": bp_used,
        "reference_group_mean_coverage_pct": ref_means,
        "note": ("coverage_pct = estimated % of the diploid autosomal genome lying in the Neanderthal "
                 "segments Sprime detected in Iberians. It captures only detectable, IBS-ascertained "
                 "segments, so it is NOT total Neanderthal ancestry; use the percentile/ratio vs Europeans. "
                 "The IBS reference samples are the discovery sample, so they score slightly high relative to "
                 "an out-of-sample person. EUR_noIBS (CEU/GBR/TSI/FIN) avoids this discovery-sample bias, but "
                 "because the segments were found in Iberians, Southern European ancestry may score somewhat "
                 "higher against it. The truth likely lies between the IBS and EUR_noIBS comparisons."),
        "people": {},
    }
    print(f"Segments used {used_segs}, tags used {used_tags}, span {bp_used/1e6:.1f} Mb")
    print("Reference mean coverage % by group:", ref_means)
    for p in people:
        x = person_score[p]
        groups = {}
        for code, label in GROUPS:
            st = stats(group_vals(ref_score, code), x)
            sa = stats(group_vals(ref_alleles, code), person_alleles[p])
            groups[code] = {"label": label, "n": st["n"], "percentile": st["percentile"], "z": st["z"],
                            "ratio_to_mean": st["ratio_to_mean"],
                            "tag_allele_count_percentile": sa["percentile"]}
        results["people"][p] = {"coverage_pct": round(x / denom, 3),
                                "archaic_tag_alleles": round(person_alleles[p], 1), "groups": groups}
        g = groups
        print(f"  {p:18s} cov {x/denom:.3f}%  IBS pct {g['IBS']['percentile']:5.1f} (z {g['IBS']['z']:+.2f}, "
              f"x{g['IBS']['ratio_to_mean']:.2f})  EUR pct {g['EUR']['percentile']:5.1f} "
              f"(x{g['EUR']['ratio_to_mean']:.2f})  EUR_noIBS z {g['EUR_noIBS']['z']:+.2f}  "
              f"tag-count IBS pct {g['IBS']['tag_allele_count_percentile']:.1f}")
    with open(os.path.join(RESULTS, "neanderthal.json"), "w") as fh:
        json.dump(results, fh, indent=2)


if __name__ == "__main__":
    main()
