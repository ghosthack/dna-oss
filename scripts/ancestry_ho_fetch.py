"""Fine-scale reference panel: AADR Human Origins (HO) present-day samples x family chip SNPs.

Usage: python3 -I scripts/ancestry_ho_fetch.py     (fetch what's missing, build the cache)
       from ancestry_ho_fetch import load_cache    (downstream loader)
Re-running after the cache exists re-downloads the selected records (parts/ is a resume
cache and may be deleted once ref_geno.bin.gz is written).

Source: Allen Ancient DNA Resource (AADR) v66.p1, Mallick et al. 2024, Harvard Dataverse
doi:10.7910/DVN/FFIDCW, "HO" release (v66.p1_HO.aadr.patch.PUB.*; 584,131 SNPs, 27,594
individuals). The .geno is TGENO (transposed PACKEDANCESTRYMAP): a 48-byte header
"TGENO nind nsnp hashind hashsnp", then one record per *individual* of ceil(nsnp/4) bytes,
2 bits per SNP (high bits first) = copies of the .snp column-5 allele, 3 = missing.
Because records are per individual, only the selected present-day individuals are fetched
(HTTP Range requests, ~150 KB each) instead of the 4 GB file. Their integrity is checked via
the header hashes (EIGENSOFT hashit/hasharr over the .ind/.snp IDs) plus md5 of the metadata
files; the full-file md5 can't be checked without downloading all of it.

Selection: date mean BP == 0 (present-day), Group ID in POPS (exact; so "-o", "-QCremove",
"-oRelative", "Ignore_" outliers are excluded), ASSESSMENT in Pass/PROVISIONAL_PASS/
MERGE_PASS, data type HumanOrigins or Shotgun.diploid (pseudo-haploid "Shotgun" dropped),
deduplicated by Individual ID (HO array copy preferred over the shotgun copy of the same person).
Close relatives: AADR's "Family relations" field is blank for the related pairs found in the
reference groups used here, so it cannot be relied on. A KING-robust kinship scan within every
population (scripts/ancestry_ho.py, step 0; results/ancestry_ho.json validation.reference_relatives)
found two first-degree pairs, both in Italian_South: PuG17/PuG20 (Soleto/Calimera, kinship ~0.25,
parent-child pattern) and ITS2/ITS5 (Naples/Salerno, kinship ~0.29, full siblings). One member of
each pair is dropped (EXCLUDE_RELATIVES: the one with fewer HO SNPs). load_cache() applies the
exclusion to an existing cache, so no re-download is needed; select_samples() applies it too.

SNPs: autosomal HO SNPs where every family member who has the rsID has it at the same
GRCh37 chrom:pos, alleles match ref/alt directly or (consistently, for all members) after
complementing; palindromic A/T, C/G SNPs dropped; anything else dropped.

Cache (work/ancestry_ho/):
  snps.tsv           rsid chrom pos a1 a2 flipped in_family (a1 = counted allele, GRCh37 +)
  samples.tsv        genetic_id iid pop region locality country lat lon datatype publication
                     sex hosnps aadr_index
  ref_geno.bin.gz    uint8 [n_samples x n_snps] row-major; value = copies of a1, 3 = missing
  family_geno.tsv.gz rsid + one column per person (copies of a1, '.' missing)
  pop_counts.tsv.gz  rsid then per pop "a1count/called_alleles"
  manifest.json      source, URLs, checksums, counts, validation
Nothing about the family is uploaded; only public reference files are downloaded.
"""
import array
import collections
import concurrent.futures
import csv
import gzip
import hashlib
import json
import operator
import os
import sys
import threading
import time
import urllib.error
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dnalib import RAW, WORK, load_person, people as family_ids  # noqa: E402

OUT = os.path.join(WORK, "ancestry_ho")
DL = os.path.join(OUT, "download")
PARTS = os.path.join(OUT, "parts")

DOI = "doi:10.7910/DVN/FFIDCW"
API = "https://dataverse.harvard.edu/api"
VERSION = "v66.p1"
FILES = {  # Dataverse file id, expected md5 (from dataset metadata / v66.p1__files.md5sum)
    "geno": ("v66.p1_HO.aadr.patch.PUB.geno", 13994808, "9f9f0629f803835985bc38856e989995"),
    "snp": ("v66.p1_HO.aadr.patch.PUB.snp", 13994527, "fbd2cfe09ec6bad794f7607768b57a3a"),
    "ind": ("v66.p1_compatibility_HO.aadr.patch.PUB.ind", 13994523,
            "a14985ffae6dd33a00f4fd5bfefe58c5"),  # identical to HO .ind (same md5)
    "anno": ("v66.p1_compatibility_HO.aadr.PUB.anno", 13994525,
             "548b011b4adbcd46f0ff66b957c59c3f"),  # identical to HO .anno (same md5)
    "md5sum": ("v66.p1__files.md5sum", 13995422, None),
}

DV_SIZES = {"geno": 4029634650, "snp": 36800253, "ind": 1176253, "anno": 15465906,
            "md5sum": 1769}
BYTES_FETCHED_GENO = "~326 MB (2,231 records x 146,033 B + header)"

# label -> region (for grouping downstream). Exact AADR Group IDs.
POPS = {
    # Italy & central Mediterranean
    "Sicilian": "Italy", "Italian_South": "Italy", "Italian_Central": "Italy",
    "Italian_North": "Italy", "TSI": "Italy", "Sardinian": "Italy", "Maltese": "CentralMed",
    # Balkans / Greece
    "Greek": "Balkans", "Greek_WGA": "Balkans", "Greek_Crete": "Balkans", "Albanian": "Balkans",
    "Cypriot": "EastMed", "Bulgarian": "Balkans", "Romanian": "Balkans", "Croatian": "Balkans",
    "Serbian_Serb": "Balkans", "Montenegrin_Serb": "Balkans", "Bosnian_Serb": "Balkans",
    "Croatian_Serb": "Balkans",
    # Iberia / West & North Europe
    "Spanish": "Iberia", "Spanish_North": "Iberia", "IBS": "Iberia",
    "IBS_CanaryIslands": "Iberia", "Basque": "Iberia",
    "French": "WestEurope", "English": "WestEurope", "GBR": "WestEurope",
    "Scottish": "WestEurope", "Orcadian": "WestEurope", "CEU": "WestEurope",
    "Norwegian": "NorthEurope", "Icelandic": "NorthEurope", "Finnish": "NorthEurope",
    "Hungarian": "EastEurope", "Czech": "EastEurope", "Ukrainian": "EastEurope",
    "Russian": "EastEurope",
    # Anatolia / Caucasus / Iran
    "Turkish": "Anatolia", "Turkish_Trabzon": "Anatolia", "Armenian": "Caucasus",
    "Georgian": "Caucasus", "Adygei": "Caucasus", "Circassian": "Caucasus",
    "Chechen": "Caucasus", "Assyrian": "NearEast", "Kurd": "NearEast", "Iranian": "Iran",
    # Levant / Arabia
    "Lebanese": "Levant", "Lebanese_Christian": "Levant", "Lebanese_Muslim": "Levant",
    "Druze": "Levant", "Palestinian": "Levant", "Jordanian": "Levant", "Syrian": "Levant",
    "Saudi": "Arabia", "BedouinA": "Arabia", "BedouinB": "Arabia",
    # Jewish groups
    "Jew_Ashkenazi": "Jewish", "Jew_Moroccan": "Jewish", "Jew_Libyan": "Jewish",
    "Jew_Tunisian": "Jewish", "Jew_Turkish": "Jewish", "Jew_Iraqi": "Jewish",
    "Jew_Iranian": "Jewish", "Jew_Georgian": "Jewish", "Jew_Yemenite": "Jewish",
    # North Africa
    "Moroccan": "NorthAfrica", "Tunisian": "NorthAfrica", "Algerian": "NorthAfrica",
    "Libyan": "NorthAfrica", "Mozabite": "NorthAfrica", "Egyptian": "NorthAfrica",
    # Out-groups
    "Yoruba": "SubSaharan", "YRI": "SubSaharan", "Mbuti": "SubSaharan",
    "Somali": "SubSaharan", "LWK": "SubSaharan", "GWD": "SubSaharan",
    "Han": "EastAsia", "CHB": "EastAsia",
    "Pima": "America", "Karitiana": "America", "Mexico_Zapotec": "America", "MXL": "America",
}
REQUESTED_MISSING_CHECK = ["Sicilian", "Italian_South", "Italian_North", "TSI", "Sardinian",
                           "Maltese", "Greek", "Greek_Crete", "Albanian", "Spanish",
                           "Spanish_North", "Basque", "French", "English", "Cypriot",
                           "Turkish", "Lebanese", "Druze", "Palestinian", "Moroccan",
                           "Tunisian", "Algerian", "Mozabite", "Yoruba", "Han", "Pima",
                           "Karitiana", "MXL"]
# One member of each first-degree pair inside a reference population (kept: PuG17 535,350 vs
# PuG20 534,964 HO SNPs; ITS5 573,269 vs ITS2 572,866). See the module docstring.
EXCLUDE_RELATIVES = {"PuG20": "parent-child of PuG17 (Italian_South)",
                     "ITS2": "full sibling of ITS5 (Italian_South)"}
OK_ASSESS = {"Pass", "PROVISIONAL_PASS", "MERGE_PASS"}
OK_TYPES = {"HumanOrigins": 0, "Shotgun.diploid": 1}   # value = preference (lower wins)
VALIDATION_SNPS = ["rs1042602", "rs260690"]   # the classic LCT/HERC2/SLC24A5 SNPs are not on HO
N_DUP_CHECK = 24     # HO-array vs shotgun copies of the same person, for a concordance test
COMP = str.maketrans("ACGT", "TGCA")
WORKERS = 6
CHUNK = 16          # individuals per range request (contiguous runs only)
UA = {"User-Agent": "family-dna-ref-fetch/1.0 (python urllib)"}  # default urllib UA gets 403


# ---------------------------------------------------------------- downloads
def path(key):
    return os.path.join(DL, FILES[key][0])


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def fetch_meta():
    os.makedirs(DL, exist_ok=True)
    for key in ("snp", "ind", "anno", "md5sum"):
        name, fid, want = FILES[key]
        p = path(key)
        if not os.path.exists(p) or (want and md5(p) != want):
            print(f"downloading {name}", flush=True)
            req = urllib.request.Request(f"{API}/access/datafile/{fid}", headers=UA)
            with urllib.request.urlopen(req, timeout=600) as r, open(p + ".part", "wb") as f:
                for b in iter(lambda: r.read(1 << 20), b""):
                    f.write(b)
            os.replace(p + ".part", p)
        if want and md5(p) != want:
            sys.exit(f"md5 mismatch for {name}")


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, *a, **k):
        return None


_SIGNED = {"url": None, "t": 0.0}
_LOCK = threading.Lock()


def signed_url(force=False):
    """Resolve Dataverse's redirect to the presigned S3 URL once and reuse it (Dataverse
    rate-limits the access API with 403s when every Range request goes through it)."""
    with _LOCK:
        if force or not _SIGNED["url"] or time.time() - _SIGNED["t"] > 1800:
            opener = urllib.request.build_opener(_NoRedirect)
            req = urllib.request.Request(f"{API}/access/datafile/{FILES['geno'][1]}", headers=UA)
            for t in range(8):
                try:
                    opener.open(req, timeout=60)
                    raise IOError("expected a redirect")
                except urllib.error.HTTPError as e:
                    if e.code in (301, 302, 303, 307) and e.headers.get("Location"):
                        _SIGNED.update(url=e.headers["Location"], t=time.time())
                        break
                    print(f"  resolve: HTTP {e.code}, waiting", flush=True)
                    time.sleep(60 * (t + 1))
            else:
                raise IOError("could not resolve signed URL")
        return _SIGNED["url"]


def get_range(start, length, tries=6):
    for t in range(tries):
        url = signed_url(force=t > 0)
        try:
            req = urllib.request.Request(url, headers={**UA, "Range": f"bytes={start}-{start + length - 1}"})
            with urllib.request.urlopen(req, timeout=300) as r:
                if r.status != 206:  # a server that ignores Range would stream the whole multi-GB file
                    raise IOError(f"expected HTTP 206 Partial Content, got {r.status}")
                data = r.read(length + 1)
            if len(data) == length:
                return data
            raise IOError(f"short read {len(data)}/{length}")
        except Exception as e:  # noqa: BLE001
            if t == tries - 1:
                raise
            print(f"  retry {start}: {e}", flush=True)
            time.sleep(5 * (t + 1))


# ---------------------------------------------------------------- EIGENSOFT hashes
def _hashit(s):
    h = 0
    for c in s.encode():
        h = (h * 23 + c) & 0xFFFFFFFF
    return h


def hasharr(ids):
    h = 0
    for s in ids:
        h = ((h * 17) & 0xFFFFFFFF) ^ _hashit(s)
    return format(h, "x")


# ---------------------------------------------------------------- metadata
def read_snp():
    snps = []
    with open(path("snp")) as f:
        for line in f:
            rs, ch, _gp, pos, a1, a2 = line.split()
            snps.append((rs, ch, int(pos), a1, a2))
    return snps


def read_ind():
    with open(path("ind")) as f:
        return [line.split() for line in f]


def read_anno():
    with open(path("anno"), newline="") as f:
        rows = list(csv.reader(f, delimiter="\t"))
    head = rows[0]

    def col(prefix):
        return next(i for i, h in enumerate(head) if h.startswith(prefix))
    c = {"gid": 0, "iid": col("Individual ID"), "pub": col("Publication abbreviation"),
         "date": col("Date mean in BP"), "group": col("Group ID"), "loc": col("Locality"),
         "country": col("Political Entity"), "lat": col("Latitude"), "lon": col("Longitude"),
         "dtype": col("Data type"), "sex": col("Molecular Sex"),
         "rel": col("Family relations"), "assess": col("ASSESSMENT"),
         "hosnps": col("SNPs hit on autosomal targets (Computed using easystats on HO")}
    return [{k: r[i] for k, i in c.items()} for r in rows[1:]]


def select_samples(ind, anno):
    assert len(ind) == len(anno) and all(a[0] == b["gid"] for a, b in zip(ind, anno)), \
        ".ind and .anno disagree"
    cand = []
    for i, a in enumerate(anno):
        if a["date"] != "0" or a["group"] not in POPS:
            continue
        if a["assess"] not in OK_ASSESS or a["dtype"] not in OK_TYPES:
            continue
        if a["iid"] in EXCLUDE_RELATIVES:
            continue
        if a["rel"] not in ("..", ""):
            a["rel_flag"] = True      # informational only; relatives are found by kinship scan
        cand.append((OK_TYPES[a["dtype"]], -int(a["hosnps"] or 0) if a["hosnps"].isdigit() else 0, i))
    best, dups = {}, []
    for pref, negcov, i in sorted(cand):
        iid = anno[i]["iid"]
        if iid in best:                              # keep preferred copy per person
            if anno[best[iid]]["dtype"] != anno[i]["dtype"]:
                dups.append((best[iid], i))
        else:
            best[iid] = i
    return sorted(best.values()), dups


# ---------------------------------------------------------------- family overlap
def family_people():
    """(person id, genotypes) for every person configured in family.json who has data, in id order."""
    out = []
    for d in sorted(family_ids()):
        if os.path.isdir(os.path.join(RAW, d)):
            data, _ = load_person(d)
            if data:
                out.append((d, data))
    return out


def match_snps(snps, people):
    """Return kept [(idx, rs, ch, pos, a1, a2, flipped)], per-person counts, drop stats."""
    drop = collections.Counter()
    kept = []
    per = collections.Counter()
    need = set(VALIDATION_SNPS)
    for idx, (rs, ch, pos, a1, a2) in enumerate(snps):
        if not ch.isdigit() or int(ch) > 22:
            continue
        have = [(p, d[rs]) for p, d in people if rs in d]
        if not have:
            if rs in need:
                kept.append((idx, rs, ch, pos, a1, a2, False, []))
            continue
        if {a1, a2} in ({"A", "T"}, {"C", "G"}):
            drop["palindromic"] += 1
            continue
        if any(c != ch or p_ != pos for _p, (c, p_, _g) in have):
            drop["position_mismatch"] += 1
            continue
        modes = set()
        for _p, (_c, _q, g) in have:
            al = set(g)
            if len(g) != 2 or not al <= set("ACGT"):
                modes.add("bad")
            elif al <= {a1, a2}:
                modes.add("direct" if not al <= {x.translate(COMP) for x in (a1, a2)} else "either")
            elif {x.translate(COMP) for x in al} <= {a1, a2}:
                modes.add("flip")
            else:
                modes.add("bad")
        modes.discard("either")
        if "bad" in modes or len(modes) > 1:
            drop["allele_mismatch"] += 1
            continue
        flip = modes == {"flip"}
        drop["kept_flipped" if flip else "kept_direct"] += 1
        kept.append((idx, rs, ch, pos, a1, a2, flip, [p for p, _ in have]))
        for p, _ in have:
            per[p] += 1
    return kept, per, drop


def dosage(gt, a1, flip):
    if flip:
        gt = gt.translate(COMP)
    return gt.count(a1)


# ---------------------------------------------------------------- genotype extraction
_M3 = {}


def unpack(rec, nsnp):
    """TGENO record (bytes) -> bytearray of 2-bit values, one per SNP (stdlib, big-int trick)."""
    n = len(rec)
    m = _M3.get(n)
    if m is None:
        m = _M3[n] = int.from_bytes(b"\x03" * n, "big")
    v = int.from_bytes(rec, "big")
    out = bytearray(4 * n)
    for k, sh in enumerate((6, 4, 2, 0)):
        out[k::4] = ((v >> sh) & m).to_bytes(n, "big")
    return out[:nsnp]


def fetch_genotypes(sample_idx, snp_idx, nsnp, rlen, hdr=48):
    """Download records of the selected individuals, keep the selected SNP columns.

    Returns {sample_idx: bytes(len(snp_idx))}; resumable via parts/."""
    os.makedirs(PARTS, exist_ok=True)
    pick = operator.itemgetter(*snp_idx)
    tag = hashlib.md5(json.dumps(snp_idx).encode()).hexdigest()[:10]
    runs, cur = [], []
    for i in sample_idx:
        if cur and (i != cur[-1] + 1 or len(cur) >= CHUNK):
            runs.append(cur)
            cur = []
        cur.append(i)
    if cur:
        runs.append(cur)

    def job(run):
        pf = os.path.join(PARTS, f"{tag}_{run[0]}_{run[-1]}.bin")
        if os.path.exists(pf):
            return run, open(pf, "rb").read()
        data = get_range(hdr + run[0] * rlen, len(run) * rlen)
        buf = b"".join(bytes(pick(unpack(data[k * rlen:(k + 1) * rlen], nsnp)))
                       for k in range(len(run)))
        with open(pf + ".tmp", "wb") as f:
            f.write(buf)
        os.replace(pf + ".tmp", pf)
        return run, buf

    res, done, t0 = {}, 0, time.time()
    with concurrent.futures.ThreadPoolExecutor(WORKERS) as ex:
        for run, buf in ex.map(job, runs):
            k = len(snp_idx)
            for j, i in enumerate(run):
                res[i] = buf[j * k:(j + 1) * k]
            done += 1
            if done % 20 == 0 or done == len(runs):
                print(f"  {done}/{len(runs)} chunks, {len(res)} individuals, "
                      f"{time.time() - t0:.0f}s", flush=True)
    return res


# ---------------------------------------------------------------- loader for downstream use
def load_cache(exclude_relatives=True):
    """Return dict(snps, samples, geno, family, pops, excluded) from work/ancestry_ho/.

    geno: list of bytes per sample (value per SNP = copies of a1, 3 = missing).
    family: {person: list of int|None}.  pops: {pop: [sample positions]}.
    With exclude_relatives, samples in EXCLUDE_RELATIVES stay in samples/geno (positions are
    unchanged) but are left out of pops, so no analysis uses them."""
    with open(os.path.join(OUT, "snps.tsv")) as f:
        snps = [line.rstrip("\n").split("\t") for line in f][1:]
    with open(os.path.join(OUT, "samples.tsv")) as f:
        rd = csv.DictReader(f, delimiter="\t")
        samples = list(rd)
    k = len(snps)
    with gzip.open(os.path.join(OUT, "ref_geno.bin.gz"), "rb") as f:
        raw = f.read()
    geno = [raw[j * k:(j + 1) * k] for j in range(len(samples))]
    family = {}
    with gzip.open(os.path.join(OUT, "family_geno.tsv.gz"), "rt") as f:
        head = f.readline().rstrip("\n").split("\t")[1:]
        cols = [[] for _ in head]
        for line in f:
            for c, v in zip(cols, line.rstrip("\n").split("\t")[1:]):
                c.append(None if v == "." else int(v))
    family = dict(zip(head, cols))
    pops = collections.defaultdict(list)
    excluded = []
    for j, s in enumerate(samples):
        if exclude_relatives and s["iid"] in EXCLUDE_RELATIVES:
            excluded.append(j)
            continue
        pops[s["pop"]].append(j)
    return {"snps": snps, "samples": samples, "geno": geno, "family": family, "pops": dict(pops),
            "excluded": excluded}


# ---------------------------------------------------------------- main
def main():
    os.makedirs(OUT, exist_ok=True)
    fetch_meta()
    snps, ind, anno = read_snp(), read_ind(), read_anno()
    hdr = get_range(0, 48)
    magic, nind, nsnp, hind, hsnp = hdr.rstrip(b"\0").decode().split()
    nind, nsnp = int(nind), int(nsnp)
    assert magic == "TGENO" and nind == len(ind) and nsnp == len(snps), hdr
    my_hind, my_hsnp = hasharr([r[0] for r in ind]), hasharr([s[0] for s in snps])
    hash_ok = (my_hind == hind and my_hsnp == hsnp)
    print(f"geno header {magic} nind={nind} nsnp={nsnp} hashes {hind}/{hsnp} "
          f"(computed {my_hind}/{my_hsnp}: {'OK' if hash_ok else 'MISMATCH'})")
    if not hash_ok:  # the hashes tie each range-fetched record to its sample and SNP; never read past a mismatch
        sys.exit("The AADR genotype file does not match the downloaded .ind/.snp files (EIGENSOFT header hashes "
                 "differ): the release was probably updated. Delete work/ancestry_ho/ and rerun, or update the "
                 "pinned file IDs.")
    rlen = (nsnp + 3) // 4

    sel, dups = select_samples(ind, anno)
    print(f"selected {len(sel)} present-day individuals ({len(dups)} cross-platform duplicates)")
    step = max(1, len(dups) // N_DUP_CHECK)
    dups = dups[::step][:N_DUP_CHECK]

    people = family_people()
    print("family:", ", ".join(f"{p} ({len(d)})" for p, d in people))
    kept, per, drop = match_snps(snps, people)
    fam_kept = [k for k in kept if k[7]]
    print(f"overlap SNPs {len(fam_kept)}; per person {dict(per)}; {dict(drop)}")

    snp_idx = [k[0] for k in kept]
    geno = fetch_genotypes(sorted(set(sel) | {b for _a, b in dups}), snp_idx, nsnp, rlen)
    same = diff = 0
    for a, b in dups:
        for x, y in zip(geno[a], geno[b]):
            if x != 3 and y != 3:
                same += x == y
                diff += x != y
    dup_check = {"pairs": len(dups), "compared_genotypes": same + diff,
                 "concordance": round(same / max(1, same + diff), 5),
                 "note": "HO array vs Shotgun.diploid copy of the same individual"}
    print("duplicate concordance", dup_check)

    # flip so a1 is the allele counted on the plus strand of the family's chips
    out_snps = []
    for k in kept:
        idx, rs, ch, pos, a1, a2, flip, _ = k
        if flip:
            a1, a2 = a1.translate(COMP), a2.translate(COMP)
        out_snps.append((rs, ch, pos, a1, a2, flip))

    with open(os.path.join(OUT, "snps.tsv"), "w") as f:
        f.write("rsid\tchrom\tpos\ta1\ta2\tflipped\tin_family\n")
        for (rs, ch, pos, a1, a2, flip), k in zip(out_snps, kept):
            f.write(f"{rs}\t{ch}\t{pos}\t{a1}\t{a2}\t{int(flip)}\t{int(bool(k[7]))}\n")
    cols = ["genetic_id", "iid", "pop", "region", "locality", "country", "lat", "lon",
            "datatype", "publication", "sex", "hosnps", "aadr_index"]
    with open(os.path.join(OUT, "samples.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for i in sel:
            a = anno[i]
            f.write("\t".join([a["gid"], a["iid"], a["group"], POPS[a["group"]], a["loc"],
                               a["country"], a["lat"], a["lon"], a["dtype"], a["pub"],
                               a["sex"], a["hosnps"], str(i)]) + "\n")
    with gzip.open(os.path.join(OUT, "ref_geno.bin.gz"), "wb", compresslevel=6) as f:
        for i in sel:
            f.write(geno[i])

    with gzip.open(os.path.join(OUT, "family_geno.tsv.gz"), "wt") as f:
        f.write("rsid\t" + "\t".join(p for p, _ in people) + "\n")
        for (rs, ch, pos, a1, a2, flip) in out_snps:
            vals = []
            for _p, d in people:
                r = d.get(rs)
                vals.append(str(dosage(r[2], a1, False)) if r else ".")
            f.write(rs + "\t" + "\t".join(vals) + "\n")

    pops = collections.defaultdict(list)
    for i in sel:
        pops[anno[i]["group"]].append(i)
    pop_names = sorted(pops)
    sums = {p: pop_sums([geno[i] for i in pops[p]]) for p in pop_names}
    with gzip.open(os.path.join(OUT, "pop_counts.tsv.gz"), "wt") as f:
        f.write("rsid\t" + "\t".join(pop_names) + "\n")
        for j, s in enumerate(out_snps):
            f.write(s[0] + "\t" + "\t".join(f"{sums[p][0][j]}/{sums[p][1][j]}"
                                             for p in pop_names) + "\n")

    validation = validate(out_snps, geno, pops)
    validation["duplicate_concordance"] = dup_check
    validation["geno_header_hash_check"] = hash_ok
    pop_n = {p: len(v) for p, v in sorted(pops.items())}
    pop_detail = {p: dict(collections.Counter(anno[i]["dtype"] for i in v))
                  for p, v in sorted(pops.items())}
    miss = {p: pop_n.get(p, 0) for p in REQUESTED_MISSING_CHECK if pop_n.get(p, 0) < 5}
    manifest = {
        "source": "Allen Ancient DNA Resource (AADR), Mallick et al. 2024, Sci Data 11:182",
        "doi": DOI, "aadr_version": VERSION, "dataverse_dataset_version": dataverse_version(),
        "release": "HO (Human Origins) present-day + shotgun-diploid genomes on HO SNP set",
        "files": {k: {"name": v[0], "url": f"{API}/access/datafile/{v[1]}", "md5": v[2],
                      "size": os.path.getsize(path(k)) if os.path.exists(path(k)) else None,
                      "size_on_dataverse": DV_SIZES.get(k),
                      "downloaded": "range-requests for selected individual records only, "
                                    + BYTES_FETCHED_GENO
                      if k == "geno" else "full"} for k, v in FILES.items()},
        "geno_format": "TGENO, 48-byte header, record per individual of ceil(nsnp/4) bytes",
        "geno_header": {"nind": nind, "nsnp": nsnp, "hashind": hind, "hashsnp": hsnp,
                        "hash_check": hash_ok},
        "selection": {"present_day": "date mean BP == 0", "assessment": sorted(OK_ASSESS),
                      "data_types": sorted(OK_TYPES), "dedup": "by Individual ID, HO preferred",
                      "outlier_labels": "excluded by exact Group ID match"},
        "n_samples": len(sel), "populations": pop_n, "population_datatypes": pop_detail,
        "low_or_missing_requested": miss,
        "family_people": [p for p, _ in people],
        "snps": {"total_in_cache": len(out_snps), "family_overlap": len(fam_kept),
                 "validation_only": len(out_snps) - len(fam_kept),
                 "per_person": dict(per), "filters": dict(drop),
                 "matching": "rsID + GRCh37 chrom:pos; autosomes; no A/T or C/G; consistent "
                             "complement flips; others dropped"},
        "validation": validation,
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    with open(os.path.join(OUT, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1)
    print(json.dumps({k: manifest[k] for k in ("populations", "low_or_missing_requested",
                                               "snps", "validation")}, indent=1))


_DOS = bytes([0, 1, 2, 0] + [0] * 252)
_CALL = bytes([2, 2, 2, 0] + [0] * 252)


def pop_sums(rows):
    """Per-SNP (a1 count, called alleles) over rows of genotype bytes (16-bit big-int lanes)."""
    k = len(rows[0])
    tot_c = tot_n = 0
    lane = bytearray(2 * k)
    for g in rows:
        lane[1::2] = g.translate(_DOS)
        tot_c += int.from_bytes(lane, "big")
        lane[1::2] = g.translate(_CALL)
        tot_n += int.from_bytes(lane, "big")
    out = []
    for t in (tot_c, tot_n):
        a = array.array("H")
        a.frombytes(t.to_bytes(2 * k, "big"))
        if sys.byteorder == "little":
            a.byteswap()
        out.append(a)
    return out


def dataverse_version():
    try:
        req = urllib.request.Request(f"{API}/datasets/:persistentId/?persistentId={DOI}",
                                     headers=UA)
        with urllib.request.urlopen(req, timeout=60) as r:
            v = json.load(r)["data"]["latestVersion"]
        return f"{v['versionNumber']}.{v['versionMinorNumber']} ({v.get('releaseTime')})"
    except Exception:  # noqa: BLE001
        return None


def validate(out_snps, geno, pops):
    """Allele frequencies of well-known pigmentation/lactase SNPs per population."""
    want = {"rs1042602": "A", "rs260690": "A"}
    show = ["CEU", "GBR", "TSI", "Sicilian", "IBS", "Spanish", "Sardinian", "Greek",
            "Turkish", "Palestinian", "Mozabite", "Yoruba", "YRI", "Han", "CHB", "Pima"]
    res = {}
    for j, (rs, ch, pos, a1, a2, flip) in enumerate(out_snps):
        if rs not in want:
            continue
        allele = want[rs]
        if allele not in (a1, a2):
            res[rs] = f"allele {allele} not in {a1}/{a2}"
            continue
        row = {}
        for p in show:
            c = n = 0
            for i in pops.get(p, []):
                g = geno[i][j]
                if g != 3:
                    c += g if a1 == allele else 2 - g
                    n += 2
            row[p] = round(c / n, 3) if n else None
        res[f"{rs} {allele}"] = row
    return res


if __name__ == "__main__":
    main()
