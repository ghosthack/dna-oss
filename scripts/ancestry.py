"""Ancestry composition (supervised admixture) from 23andMe imputed data vs 1000 Genomes.

Usage: python3 -I scripts/ancestry.py            (fetch reference if needed, then analyse)
       python3 -I scripts/ancestry.py --fetch    (only download/cache the reference windows)

Method
  * Reference: 1000 Genomes high-coverage GRCh38 phased panel, queried remotely in ~2,000
    windows of 30 kb spread evenly over chr1-22 (biallelic SNVs, MAF >= 5%). Cached in
    work/ancestry/ref/win/ (one gzip per batch of windows), so reruns need no network.
  * SNPs: kept if present (same REF/ALT) in the family's imputed BCFs, thinned to >= 5 kb
    apart. Per person only confident imputed calls are used (each HDS < 0.1 or > 0.9).
  * Frequencies: per 1000G population from unrelated founders, excluding validation hold-outs.
  * Admixture: maximum likelihood of mixture proportions with fixed source allele
    frequencies (the ADMIXTURE likelihood with P fixed; concave, solved exactly by
    active-set Newton, same optimum as EM). Standard errors: leave-one-
    chromosome-out jackknife. Best single population = highest log-likelihood of 26.
  * Native American proxy: all PEL founders, each one's estimated European/African/East
    Asian share subtracted from the allele frequencies (shares estimated against an
    iteratively refined least-admixed PEL subset). ANC_NAT=subset uses that subset instead.
Nothing personal is sent anywhere; only public reference data is downloaded.
"""
import concurrent.futures
import glob
import gzip
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dnalib  # noqa: E402
from dnalib import BCFTOOLS, KG_PANEL, RESULTS, WORK, kg_sample_files  # noqa: E402
OUT = os.path.join(WORK, "ancestry")
REF = os.path.join(OUT, "ref")
WIN = os.path.join(REF, "win")
PGS_REF = os.path.join(WORK, "pgs", "ref")

KG = KG_PANEL
GT_TOKEN = re.compile(r"([0-9.])[|/]([0-9.])|\.")  # one diploid GT ("0|1") or a bare missing "."
CHROMS = [f"chr{c}" for c in range(1, 23)]
GRCH38_LEN = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
              145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
              101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468]

N_WINDOWS = 2000
WIN_SIZE = 30000
BATCH = 10            # windows per bcftools call
THREADS = 4
MIN_MAF = 0.05
THIN_BP = 5000
CONF = 0.1            # haplotype dosage must be < CONF or > 1-CONF

NAT_METHOD = os.environ.get("ANC_NAT", "deconv")   # "deconv" (all PEL) or "subset"
PEOPLE_ORDER = dnalib.people()   # self, siblings, mother, father, maternal grandmother

# Admixture sources: label -> 1000G populations (PEL is filtered, see native_american_subset)
SOURCES = [
    ("Iberian (IBS)", ["IBS"]),
    ("Northwest European (GBR+CEU)", ["GBR", "CEU"]),
    ("Southern European / Italian (TSI)", ["TSI"]),
    ("Finnish (FIN)", ["FIN"]),
    ("West African (YRI+ESN)", ["YRI", "ESN"]),
    ("East Asian (CHB+JPT)", ["CHB", "JPT"]),
    ("South Asian (GIH)", ["GIH"]),
    ("Native American proxy (PEL, admixture-corrected)", ["PEL*"]),
]
BROAD = {
    "European": ["Iberian (IBS)", "Northwest European (GBR+CEU)",
                 "Southern European / Italian (TSI)", "Finnish (FIN)"],
    "African": ["West African (YRI+ESN)"],
    "East Asian": ["East Asian (CHB+JPT)"],
    "South Asian": ["South Asian (GIH)"],
    "Native American": ["Native American proxy (PEL, admixture-corrected)"],
}
# Coarser European split, reported with its own SE (the 4-way European split is noisy).
EURO_REGIONS = {
    "Southern European (IBS+TSI)": ["Iberian (IBS)", "Southern European / Italian (TSI)"],
    "Northern European (GBR+CEU+FIN)": ["Northwest European (GBR+CEU)", "Finnish (FIN)"],
}
# Held out of every frequency estimate and used for validation.
HOLDOUT_POPS = {"IBS": 3, "GBR": 3, "TSI": 3, "FIN": 2, "YRI": 3, "CHB": 2, "GIH": 2,
                "MXL": 10, "PUR": 10, "CLM": 10, "ASW": 10}


def log(*a):
    print(*a, flush=True)


# ---------------------------------------------------------------- reference fetch

def windows():
    total = sum(GRCH38_LEN)
    step = total / N_WINDOWS
    out = []
    for chrom, length in zip(CHROMS, GRCH38_LEN):
        x = step / 2
        while x + WIN_SIZE < length - 10000:
            out.append((chrom, int(x), int(x) + WIN_SIZE))
            x += step
    return out


def batches():
    by_chrom = {}
    for w in windows():
        by_chrom.setdefault(w[0], []).append(w)
    out = []
    for chrom in CHROMS:
        ws = by_chrom.get(chrom, [])
        for i in range(0, len(ws), BATCH):
            out.append((f"{chrom}_{i // BATCH:03d}", ws[i:i + BATCH]))
    return out


def fetch_batch(name, wins):
    path = os.path.join(WIN, name + ".tsv.gz")
    if os.path.exists(path):
        return name, None
    chrom = wins[0][0]
    reg = os.path.join(WIN, name + ".regions")
    with open(reg, "w") as fh:
        for c, s, e in wins:
            fh.write(f"{c}\t{s}\t{e}\n")
    err = ""
    for attempt in range(4):
        t0 = time.time()
        view = subprocess.Popen([BCFTOOLS, "view", "-R", reg, "-v", "snps", "-m2", "-M2",
                                 "-q", f"{MIN_MAF}:minor", "-Ou", KG.format(chrom=chrom)],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=REF)
        query = subprocess.run([BCFTOOLS, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n"],
                               stdin=view.stdout, capture_output=True, cwd=REF)
        view.stdout.close()
        verr = view.stderr.read().decode()
        view.wait()
        if view.returncode == 0 and query.returncode == 0:
            tmp = path + ".tmp"
            n = 0
            with gzip.open(tmp, "wt", compresslevel=6) as fh:
                for line in query.stdout.decode().splitlines():
                    c, pos, ref, alt, gts = line.split("\t")
                    # "0|10|0..." -> two allele chars per sample; a bare "." (missing) becomes ".." to keep alignment
                    fh.write(f"{c}\t{pos}\t{ref}\t{alt}\t{GT_TOKEN.sub(lambda m: m.group(1) + m.group(2) if m.group(1) else '..', gts)}\n")
                    n += 1
            os.replace(tmp, path)
            os.remove(reg)
            return name, f"{n} SNVs in {time.time() - t0:.0f}s"
        err = (verr + query.stderr.decode())[-300:]
        time.sleep(5 * (attempt + 1))
    return name, "FAILED: " + err


def fetch_reference():
    os.makedirs(WIN, exist_ok=True)
    for chrom in CHROMS:  # reuse the tabix indexes already downloaded by pgs.py
        tbi = os.path.basename(KG.format(chrom=chrom)) + ".tbi"
        src, dst = os.path.join(PGS_REF, tbi), os.path.join(REF, tbi)
        if not os.path.exists(dst) and os.path.exists(src):
            with open(src, "rb") as a, open(dst, "wb") as b:
                b.write(a.read())
    todo = [b for b in batches() if not os.path.exists(os.path.join(WIN, b[0] + ".tsv.gz"))]
    if not todo:
        return
    log(f"Fetching {len(todo)} batches ({BATCH} x {WIN_SIZE // 1000} kb windows each) from 1000G…")
    t0 = time.time()
    failed = 0
    with concurrent.futures.ThreadPoolExecutor(THREADS) as pool:
        for i, (name, status) in enumerate(pool.map(lambda b: fetch_batch(*b), todo), 1):
            failed += status.startswith("FAILED")
            log(f"  [{i}/{len(todo)} {time.time() - t0:.0f}s] {name}: {status}")
    if failed:
        sys.exit(f"{failed} batches failed; rerun to retry (finished batches are cached)")


def load_samples():
    kg_sample_files(PGS_REF)  # shared with pgs.py; fetched here when pgs.py hasn't run (e.g. --no-pgs)
    samples = open(os.path.join(PGS_REF, "samples.txt")).read().split()
    pop, founder = {}, set()
    with open(os.path.join(PGS_REF, "ped.txt")) as fh:
        next(fh)
        for line in fh:
            f = line.split()
            pop[f[1]] = f[5]
            if f[2] == "0" and f[3] == "0":
                founder.add(f[1])
    return samples, pop, founder


def load_reference():
    """[(chrom, pos, ref, alt, allele-string)] for unique-position MAF>=5% biallelic SNVs."""
    rows, seen, dup = [], set(), set()
    for name, _ in batches():
        path = os.path.join(WIN, name + ".tsv.gz")
        with gzip.open(path, "rt") as fh:
            for line in fh:
                c, pos, ref, alt, al = line.rstrip("\n").split("\t")
                key = (c, int(pos))
                if key in seen:
                    dup.add(key)
                seen.add(key)
                rows.append((c, int(pos), ref, alt, al))
    return [r for r in rows if (r[0], r[1]) not in dup]


# ---------------------------------------------------------------- family genotypes

def people():
    out = []
    imp = os.path.join(WORK, "imputed")
    for p in sorted(os.listdir(imp) if os.path.isdir(imp) else [],
                    key=lambda p: (PEOPLE_ORDER.index(p) if p in PEOPLE_ORDER else 99, p)):
        if all(glob.glob(os.path.join(WORK, "imputed", p, "*", f"{c}.bcf")) for c in CHROMS):
            out.append(p)
    return out


def family_calls(person, ref_rows):
    """{(chrom,pos): (ref, alt, hds-string)} for reference positions, cached per candidate set."""
    os.makedirs(os.path.join(OUT, "family"), exist_ok=True)
    regfile = os.path.join(OUT, "family", "candidates.tsv")
    text = "".join(f"{r[0]}\t{r[1]}\n" for r in ref_rows)
    digest = hashlib.sha1(text.encode()).hexdigest()[:12]
    if not os.path.exists(regfile) or open(regfile).read() != text:
        with open(regfile, "w") as fh:
            fh.write(text)
    cache = os.path.join(OUT, "family", f"{person}.{digest}.tsv.gz")
    bcfs = [glob.glob(os.path.join(WORK, "imputed", person, "*", f"{c}.bcf"))[0] for c in CHROMS]
    newest = max(os.path.getmtime(b) for b in bcfs)
    if not os.path.exists(cache) or os.path.getmtime(cache) < newest:
        def q(bcf):
            r = subprocess.run([BCFTOOLS, "query", "-R", regfile, "-f",
                                "%CHROM\t%POS\t%REF\t%ALT\t[%HDS]\n", bcf],
                               capture_output=True, text=True)
            if r.returncode:
                raise RuntimeError(r.stderr[-300:])
            return r.stdout
        with concurrent.futures.ThreadPoolExecutor(THREADS) as pool:
            parts = list(pool.map(q, bcfs))
        with gzip.open(cache + ".tmp", "wt") as fh:
            fh.write("".join(parts))
        os.replace(cache + ".tmp", cache)
        for old in glob.glob(os.path.join(OUT, "family", f"{person}.*.tsv.gz")):
            if old != cache:
                os.remove(old)
    out = {}
    with gzip.open(cache, "rt") as fh:
        for line in fh:
            c, pos, ref, alt, hds = line.rstrip("\n").split("\t")
            out[(c, int(pos))] = (ref, alt, hds)
    return out


def call_genotype(ref, alt, hds, want_ref, want_alt):
    """ALT-allele count w.r.t. the 1000G alleles, or None if not confident / mismatched."""
    if (ref, alt) == (want_ref, want_alt):
        flip = False
    elif (ref, alt) == (want_alt, want_ref):
        flip = True
    else:
        return None
    try:
        h = [float(x) for x in hds.split(",")]
    except ValueError:
        return None
    if len(h) != 2 or any(CONF <= x <= 1 - CONF for x in h):
        return None
    g = sum(1 for x in h if x > 0.5)
    return 2 - g if flip else g


# ---------------------------------------------------------------- admixture likelihood fit
# Log-likelihood  L(pi) = sum_j g_j log q_j + (2-g_j) log(1-q_j),  q_j = sum_k pi_k p_jk
# (the ADMIXTURE model with allele frequencies held fixed). L is concave in pi, so the
# maximum over the simplex is found exactly by an active-set Newton method.

def _derivs(pi, rows, K, free):
    ll = 0.0
    grad = [0.0] * K
    H = [[0.0] * K for _ in range(K)]
    for g, p in rows:
        q = 0.0
        for k in range(K):
            q += pi[k] * p[k]
        a = g / q
        b = (2 - g) / (1 - q)
        ll += g * math.log(q) + (2 - g) * math.log(1 - q)
        w = a / q + b / (1 - q)
        for k in range(K):
            pk = p[k]
            grad[k] += (a - b) * pk
            if k in free:
                wk = w * pk
                Hk = H[k]
                for l in free:
                    if l >= k:
                        Hk[l] -= wk * p[l]
    for k in free:
        for l in free:
            if l < k:
                H[k][l] = H[l][k]
    return ll, grad, H


def _loglik(pi, rows):
    ll = 0.0
    K = len(pi)
    for g, p in rows:
        q = 0.0
        for k in range(K):
            q += pi[k] * p[k]
        ll += g * math.log(q) + (2 - g) * math.log(1 - q)
    return ll


def _solve(A, b):
    n = len(b)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        piv = max(range(c, n), key=lambda r: abs(M[r][c]))
        M[c], M[piv] = M[piv], M[c]
        if abs(M[c][c]) < 1e-300:
            raise ZeroDivisionError
        for r in range(n):
            if r != c:
                f = M[r][c] / M[c][c]
                if f:
                    for cc in range(c, n + 1):
                        M[r][cc] -= f * M[c][cc]
    return [M[i][n] / M[i][i] for i in range(n)]


def admix(rows, K, start=None, max_iter=200):
    """Maximum-likelihood mixture proportions on the simplex. rows: [(g, [p_k])]."""
    pi = list(start) if start else [1.0 / K] * K
    s = sum(pi)
    pi = [x / s for x in pi]
    free = [k for k in range(K) if pi[k] > 0]
    for _ in range(max_iter):
        ll, grad, H = _derivs(pi, rows, K, set(free))
        # Newton step on the free set with sum(d) = 0:  H d = mu*1 - grad
        F = free
        try:
            hg = _solve([[H[k][l] for l in F] for k in F], [grad[k] for k in F])
            h1 = _solve([[H[k][l] for l in F] for k in F], [1.0] * len(F))
            mu = sum(hg) / sum(h1)
            d = [mu * x - y for x, y in zip(h1, hg)]
        except ZeroDivisionError:
            d = [grad[k] - sum(grad[l] for l in F) / len(F) for k in F]
            mu = sum(grad[k] for k in F) / len(F)
        step = 1.0
        for x, k in zip(d, F):
            if x < 0 and pi[k] + x < 0:
                step = min(step, -pi[k] / x)
        new = pi[:]
        while True:
            new = pi[:]
            for x, k in zip(d, F):
                new[k] = max(0.0, pi[k] + step * x)
            ll_new = _loglik(new, rows)
            if ll_new >= ll - 1e-9 or step < 1e-8:
                break
            step /= 2
        moved = max((abs(new[k] - pi[k]) for k in F), default=0.0)
        pi = new
        free = [k for k in range(K) if pi[k] > 1e-12]
        for k in range(K):
            if k not in free:
                pi[k] = 0.0
        if moved < 1e-8:
            # KKT: a zero component should be re-activated if its gradient beats the free level.
            ll, grad, _ = _derivs(pi, rows, K, set())
            level = sum(grad[k] for k in free) / len(free)
            cand = [k for k in range(K) if k not in free and grad[k] > level + 1e-6 * abs(level) + 1e-6]
            if not cand:
                break
            k = max(cand, key=lambda k: grad[k])
            pi = [x * (1 - 1e-4) for x in pi]
            pi[k] = 1e-4
            free = sorted(free + [k])
    s = sum(pi)
    pi = [x / s for x in pi]
    return pi, _loglik(pi, rows)

# ---------------------------------------------------------------- main analysis

def main():
    fetch_only = "--fetch" in sys.argv
    t_start = time.time()
    fetch_reference()
    if fetch_only:
        return
    samples, pop_of, founders = load_samples()
    n_samp = len(samples)
    log("Loading cached reference…")
    ref_rows = load_reference()
    log(f"  {len(ref_rows)} common biallelic SNVs in {len(windows())} windows")
    fam = people()
    if not fam:
        log("No one has complete imputed data in work/imputed; nothing to analyse.")
        return
    log("Family members with imputed data: " + ", ".join(fam))
    calls = {p: family_calls(p, ref_rows) for p in fam}

    # SNPs present with matching alleles in every family member's imputed data, thinned.
    cand = []
    for r in ref_rows:
        ok = True
        for p in fam:
            c = calls[p].get((r[0], r[1]))
            if c is None or sorted(c[:2]) != sorted(r[2:4]):
                ok = False
                break
        if ok:
            cand.append(r)
    cand.sort(key=lambda r: (CHROMS.index(r[0]), r[1]))
    snps, last = [], (None, -10 ** 9)
    for r in cand:
        if r[0] != last[0] or r[1] - last[1] >= THIN_BP:
            snps.append(r)
            last = (r[0], r[1])
    log(f"  {len(cand)} present in imputed data; {len(snps)} after thinning to >= {THIN_BP // 1000} kb")

    # Reference dosages (0/1/2, or None if any allele missing), per sample index.
    def dosages(al):
        out = []
        for i in range(n_samp):
            a, b = al[2 * i], al[2 * i + 1]
            out.append(None if a == "." or b == "." else (a != "0") + (b != "0"))
        return out
    D = [dosages(r[4]) for r in snps]
    M = len(snps)
    idx = {s: i for i, s in enumerate(samples)}

    # Hold-outs: first N founders (sorted ID) of selected populations.
    holdout = []
    for popn, n in HOLDOUT_POPS.items():
        ids = sorted(s for s in samples if s in founders and pop_of.get(s) == popn)
        holdout += ids[:n]
    hold = set(holdout)
    pops = sorted({pop_of[s] for s in samples if s in founders})
    members = {pp: [idx[s] for s in samples if s in founders and s not in hold and pop_of[s] == pp]
               for pp in pops}

    def counts(ids):
        cs, ns = [], []
        for j in range(M):
            c = n = 0
            dj = D[j]
            for i in ids:
                g = dj[i]
                if g is not None:
                    c += g
                    n += 2
            cs.append(c)
            ns.append(n)
        return cs, ns

    def freqs(ids, exclude=None):
        cs, ns = counts(ids)
        if exclude is not None:  # leave one individual out
            for j in range(M):
                g = D[j][exclude]
                if g is not None:
                    cs[j] -= g
                    ns[j] -= 2
        return [(c + 0.5) / (n + 1.0) for c, n in zip(cs, ns)]

    log("Computing population allele frequencies…")
    popfreq = {pp: freqs(members[pp]) for pp in pops}

    def geno_rows(gvec, flist):
        return [(g, [f[j] for f in flist]) for j, g in enumerate(gvec) if g is not None]

    # ---- Native American proxy: iteratively keep the least European/African-admixed PEL.
    eur_ids = members["IBS"] + members["GBR"] + members["CEU"] + members["TSI"]
    afr_ids = members["YRI"] + members["ESN"]
    eas_ids = members["CHB"] + members["JPT"]
    f_eur, f_afr, f_eas = freqs(eur_ids), freqs(afr_ids), freqs(eas_ids)
    pel = list(members["PEL"])
    nat = list(pel)
    log(f"Refining Native American proxy from {len(pel)} PEL founders…")
    pel_est = {}
    for it in range(4):
        f_nat = freqs(nat)
        cs, ns = counts(nat)
        for i in pel:
            fn = f_nat
            if i in nat:
                fn = [(c - (D[j][i] or 0) + 0.5) / (n - (2 if D[j][i] is not None else 0) + 1.0)
                      for j, (c, n) in enumerate(zip(cs, ns))]
            pi, _ = admix(geno_rows([D[j][i] for j in range(M)], [f_eur, f_afr, f_eas, fn]), 4)
            pel_est[i] = pi
        ranked = sorted(pel, key=lambda i: pel_est[i][0] + pel_est[i][1])
        if it < 2:
            nat = ranked[:max(10, len(pel) // 2)]
        else:
            nat = [i for i in ranked if pel_est[i][0] + pel_est[i][1] < 0.03] or ranked[:10]
    log(f"  using {len(nat)} PEL individuals (non-Native estimate "
        f"{max(pel_est[i][0] + pel_est[i][1] for i in nat):.1%} max); all-PEL mean European "
        f"{sum(pel_est[i][0] for i in pel) / len(pel):.1%}")
    if NAT_METHOD == "subset":
        f_nat = freqs(nat)
    else:
        # Deconvolve all PEL founders: subtract each one's estimated European/African/East Asian
        # share, leaving the Native American allele frequency (uses all PEL, less sampling noise).
        f_nat = []
        for j in range(M):
            num = den = 0.0
            for i in pel:
                g = D[j][i]
                if g is None:
                    continue
                a, b, c, n = pel_est[i]
                num += g / 2 - a * f_eur[j] - b * f_afr[j] - c * f_eas[j]
                den += n
            lo = 0.5 / (2 * den + 1)
            f_nat.append(min(1 - lo, max(lo, (num + lo) / (den + 2 * lo))))

    names = [s[0] for s in SOURCES]
    src_f = []
    for label, plist in SOURCES:
        if plist == ["PEL*"]:
            src_f.append(f_nat)
        else:
            src_f.append(freqs([i for pp in plist for i in members[pp]]))
    K = len(names)

    def best_pop(gvec):
        lls = {}
        for pp in pops:
            f = popfreq[pp]
            ll = 0.0
            for j, g in enumerate(gvec):
                if g is not None:
                    ll += g * math.log(f[j]) + (2 - g) * math.log(1 - f[j])
            lls[pp] = ll
        order = sorted(lls, key=lls.get, reverse=True)
        return order[0], {pp: round(lls[pp] - lls[order[0]], 1) for pp in order[:5]}

    def analyse(gvec, jackknife=False):
        rows = geno_rows(gvec, src_f)
        pi, ll = admix(rows, K)
        def groups(v, table):
            return {b: sum(v[names.index(n)] for n in ns) for b, ns in table.items()}
        res = {"components": {n: round(x, 4) for n, x in zip(names, pi)},
               "broad": {b: round(x, 4) for b, x in groups(pi, BROAD).items()},
               "european_regions": {b: round(x, 4) for b, x in groups(pi, EURO_REGIONS).items()},
               "n_snps": len(rows)}
        if jackknife:
            # Weighted delete-one-chromosome jackknife (Busing et al. 1999): chromosomes differ several-fold in
            # SNP count, so each replicate is weighted by the SNPs it removes (as in ancestry_ho.py).
            ests, sizes = [], []
            for c in CHROMS:
                m = sum(1 for j, g in enumerate(gvec) if g is not None and snps[j][0] == c)
                if not m:
                    continue
                keep = [g if snps[j][0] != c else None for j, g in enumerate(gvec)]
                ests.append(admix(geno_rows(keep, src_f), K, start=pi)[0])
                sizes.append(m)
            total = sum(sizes)

            def se(vals, full):
                if len(vals) < 2:
                    return None
                h = [total / m for m in sizes]
                est = len(vals) * full - sum((1 - m / total) * v for m, v in zip(sizes, vals))
                var = sum((hj * full - (hj - 1) * v - est) ** 2 / (hj - 1)
                          for hj, v in zip(h, vals)) / len(vals)
                return round(math.sqrt(var), 4)
            res["se"] = {nm: se([e[k] for e in ests], pi[k]) for k, nm in enumerate(names)}
            full_broad, full_eu = groups(pi, BROAD), groups(pi, EURO_REGIONS)
            res["broad_se"] = {b: se([groups(e, BROAD)[b] for e in ests], full_broad[b]) for b in BROAD}
            res["european_regions_se"] = {b: se([groups(e, EURO_REGIONS)[b] for e in ests], full_eu[b])
                                          for b in EURO_REGIONS}
        bp, top = best_pop(gvec)
        res["best_population"] = bp
        res["best_population_loglik_vs_best"] = top
        return res

    # ---- validation on held-out 1000G individuals
    log(f"Validating on {len(holdout)} held-out 1000G individuals…")
    validation = []
    for s in holdout:
        i = idx[s]
        r = analyse([D[j][i] for j in range(M)])
        top = sorted(r["components"].items(), key=lambda kv: -kv[1])
        validation.append({"sample": s, "population": pop_of[s], "best_population": r["best_population"],
                           "broad": r["broad"],
                           "top_components": {k: v for k, v in top if v >= 0.02}})
        log(f"  {s} {pop_of[s]}: best={r['best_population']} "
            + ", ".join(f"{k.split(' (')[0]} {v:.0%}" for k, v in top if v >= 0.05))

    # ---- family
    out_people = {}
    for p in fam:
        gvec = [call_genotype(*calls[p][(r[0], r[1])], r[2], r[3]) for r in snps]
        log(f"Analysing {p} ({sum(g is not None for g in gvec)} confident SNPs)…")
        out_people[p] = analyse(gvec, jackknife=True)
        log("  " + ", ".join(f"{k} {v:.1%}" for k, v in out_people[p]["components"].items() if v >= 0.005)
            + f" | best={out_people[p]['best_population']}")

    def south_asian_caveat(res):
        # Worded for small GIH shares; optional follow-up text (e.g. what a finer panel
        # attributes them to) comes from family.json "text" -> "ancestry.south_asian_followup".
        sa = "South Asian (GIH)"
        small = [(r["components"][sa], r["se"][sa]) for r in res.values() if 0.005 <= r["components"][sa] < 0.1]
        if small and all(x <= 1.5 * se for x, se in small):
            head = ("Components < ~1% are within noise; family members' small South Asian (GIH) shares (a few %) "
                    "are at the ~1-SE level and do not indicate recent South Asian ancestry. They may")
        elif small and all(x <= 2 * se for x, se in small):
            head = ("Components < ~1% are within noise; family members' small South Asian (GIH) shares (a few %) "
                    "are within ~2 SE (not statistically distinguishable from zero). They may")
        elif small:
            head = ("Components < ~1% are within noise; family members' small South Asian (GIH) shares (a few %) "
                    "should be read against their SE before suggesting recent South Asian ancestry. They may")
        else:
            head = "Components < ~1% are within noise. Small South Asian (GIH) shares (a few %) may"
        follow = dnalib.text("ancestry.south_asian_followup")
        return (head + " partly reflect West Eurasian / Near-Eastern-related ancestry that has no source in 1000G "
                "(no Middle Eastern, North African or Ashkenazi populations), but held-out GIH themselves fit as "
                "only ~75% GIH (+~25% European), so the GIH/European boundary is poorly resolved in this panel and "
                "a few-% GIH share is not specific evidence for any ancestry" + (f"; {follow}." if follow else "."))

    def mix(a, b, wa, wb):
        return {n: wa * a[n] + wb * b[n] for n in names}

    mom, dad, gma = dnalib.person("mother"), dnalib.person("father"), dnalib.person("maternal_grandmother")
    consistency = {}
    for child in dnalib.children():
        if child in out_people and mom in out_people and dad in out_people:
            exp = mix(out_people[mom]["components"], out_people[dad]["components"], 0.5, 0.5)
            diff = {n: round(out_people[child]["components"][n] - exp[n], 4) for n in names}
            consistency[child] = {"expected_mid_parent": {n: round(v, 4) for n, v in exp.items()},
                                  "observed_minus_expected": diff,
                                  "max_abs_diff": max(abs(v) for v in diff.values())}
    grandpa = None
    if mom in out_people and gma in out_people:
        raw = mix(out_people[mom]["components"], out_people[gma]["components"], 2, -1)
        clipped = {n: max(0.0, v) for n, v in raw.items()}
        s = sum(clipped.values())
        comp = {n: round(v / s, 4) for n, v in clipped.items()}
        grandpa = {"label": "inferred maternal grandfather (2 x mom - grandma, clipped at 0, renormalised)",
                   "inferred": True,
                   "components": comp,
                   "broad": {b: round(sum(comp[n] for n in ns), 4) for b, ns in BROAD.items()},
                   "raw_before_clipping": {n: round(v, 4) for n, v in raw.items()},
                   "negative_mass_clipped": round(-sum(v for v in raw.values() if v < 0), 4),
                   "note": "Mom inherits a random ~50% of each parent's genome, so this is approximate; "
                           "uncertainty is roughly 2x mom's plus grandma's SE plus segregation noise "
                           "(a few % per component)."}

    expected_broad = {"IBS": "European", "GBR": "European", "TSI": "European", "FIN": "European",
                      "YRI": "African", "CHB": "East Asian", "GIH": "South Asian"}
    ok_best = sum(1 for v in validation if v["best_population"] == v["population"])
    ok_broad = sum(1 for v in validation if v["population"] in expected_broad
                   and v["broad"][expected_broad[v["population"]]] >= 0.9)
    n_unadm = sum(1 for v in validation if v["population"] in expected_broad)
    summary = {
        "best_population_correct": f"{ok_best}/{len(validation)}",
        "unadmixed_reference_holdouts_>=90%_expected_continent": f"{ok_broad}/{n_unadm}",
        "admixed_means": {pp: {b: round(sum(v["broad"][b] for v in validation if v["population"] == pp)
                                         / max(1, sum(1 for v in validation if v["population"] == pp)), 3)
                               for b in BROAD}
                          for pp in ("MXL", "PUR", "CLM", "ASW")},
        "admixed_n": {pp: sum(1 for v in validation if v["population"] == pp) for pp in ("MXL", "PUR", "CLM", "ASW")},
        "admixed_published_ballpark": "MXL ~45-50% Native American (wide individual range); PUR ~65-70% European, "
                                      "~15-20% African, ~12-15% Native; CLM ~65-70% European, ~25% Native, ~8% "
                                      "African; ASW ~75-80% African, ~20% European",
    }

    result = {
        "people": out_people,
        "grandpa_inferred": grandpa,
        "family_consistency": consistency,
        "validation": {"summary": summary, "individuals": validation},
        "n_snps_panel": M,
        "sources": {n: (pl if pl != ["PEL*"] else
                        [f"PEL subset ({len(nat)} of {len(pel)} founders)" if NAT_METHOD == "subset" else
                         f"all {len(pel)} PEL founders, admixture-corrected"])
                    for n, pl in SOURCES},
        "method": {
            "reference": "1000 Genomes high-coverage GRCh38 phased panel (20220422), unrelated founders, "
                         "validation hold-outs excluded from all frequencies",
            "windows": f"{len(windows())} x {WIN_SIZE // 1000} kb windows evenly spaced over chr1-22",
            "snp_filter": f"biallelic SNVs, 1000G MAF >= {MIN_MAF:.0%}, unique position, present with same "
                          f"alleles in all family imputed BCFs, thinned >= {THIN_BP // 1000} kb; per person "
                          f"only confident imputed calls (each HDS < {CONF} or > {1 - CONF})",
            "model": "supervised admixture: maximum-likelihood mixture proportions with fixed source allele "
                     "frequencies (pseudocount 0.5), concave likelihood maximised exactly by active-set Newton "
                     "(checked equal to plain EM); SE = leave-one-chromosome-out jackknife; best_population "
                     "= max binomial log-likelihood over 26 1000G populations",
            "native_american_proxy": (
                f"Each PEL founder's European/African/East Asian share was estimated (sources IBS+GBR+CEU+TSI, "
                f"YRI+ESN, CHB+JPT, and a Native proxy from the {len(nat)} least-admixed PEL, refined iteratively; "
                f"mean PEL European {sum(pel_est[i][0] for i in pel) / len(pel):.1%}). "
                + ("Native frequencies = those least-admixed PEL." if NAT_METHOD == "subset" else
                   f"Native frequencies = all {len(pel)} PEL with those shares subtracted.")),
            "caveats": [
                "Sub-European splits (Iberian vs Italian vs NW European vs Finnish) are fine-scale; with "
                "~10k SNPs from ~2,000 30-kb windows (SNPs within a window are linked, so the effective number "
                "of independent loci is far smaller) they are noisy (see 'se') and reflect similarity to these 1000G samples, "
                "not national origin.",
                "Native American proxy is Peruvian (PEL); Native ancestry from other regions (e.g. Mexico, "
                "Caribbean) is captured approximately; residual East Asian may partly reflect Native American "
                "ancestry.",
                "African is represented only by West Africans (YRI+ESN); South Asian only by GIH.",
                south_asian_caveat(out_people),
                "Jackknife SEs for the 4-way European split are large (~10-18 points); the Southern vs Northern "
                "European split ('european_regions') and best_population log-likelihoods are more robust.",
                "Imputed genotypes, not directly typed; only confident calls used.",
            ],
            "runtime_seconds": round(time.time() - t_start),
        },
    }
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "ancestry.json"), "w") as fh:
        json.dump(result, fh, indent=1)
    with open(os.path.join(RESULTS, "ancestry.tsv"), "w") as fh:
        fh.write("person\t" + "\t".join(names) + "\tbest_population\tn_snps\n")
        for p, r in out_people.items():
            fh.write(p + "\t" + "\t".join(f"{r['components'][n]:.4f}" for n in names)
                     + f"\t{r['best_population']}\t{r['n_snps']}\n")
        if grandpa:
            fh.write("grandpa-maternal (inferred)\t" + "\t".join(f"{grandpa['components'][n]:.4f}" for n in names)
                     + "\t\t\n")
    log(f"Wrote results/ancestry.json and results/ancestry.tsv ({time.time() - t_start:.0f}s)")


if __name__ == "__main__":
    main()
