"""Fine-scale ancestry of the family against the AADR Human Origins (HO) reference panel.

Usage: .venv/bin/python -I scripts/ancestry_ho.py
       (needs numpy from the project venv; reads the cache built by scripts/ancestry_ho_fetch.py
        in work/ancestry_ho/ and never downloads anything)
Writes results/ancestry_ho.json, results/ancestry_ho.tsv and work/ancestry_ho/loo_rows.tsv.

Methods
  0 QC     (a) KING-robust kinship within every reference population (and across populations of
           the same region) on all HO SNPs; first/second-degree relatives are dropped (the member
           with fewer SNPs; the known pairs are listed in ancestry_ho_fetch.EXCLUDE_RELATIVES).
           (b) Strand/allele calibration of each family member against pooled European
           frequencies (slope, frequency bins, opposite homozygotes at near-fixed SNPs).
           (c) Family relationships: opposite homozygotes per pair, trio Mendelian errors.
  SNPs     Main analysis on the autosomal HO SNPs called in *every* family member with data,
           so family members, leave-one-out reference individuals and simulated relatives are
           all measured on the same SNPs. Self is also run on all of their own SNPs as a
           sensitivity check.
  1 Admix  Supervised admixture: maximum likelihood of mixture proportions with fixed source
           allele frequencies (ADMIXTURE likelihood with P fixed; concave; active-set Newton on
           the simplex). Frequencies = pooled counts + LAMBDA pseudo-alleles of the panel mean.
           Several source configurations; leave-one-chromosome-out jackknife SEs.
           Known bias: a source estimated from few people has noisier frequencies and is
           under-assigned (real Sicilians get mostly 'NorthCentral_Italy'). Equal-n
           subsampling and empirical-Bayes shrinkage were tried and were worse (flat or more
           biased components), so the components are interpreted only through steps 3/3b.
  2 Pops   Nearest single population by log-likelihood (beta-binomial predictive) with every
           population subsampled to the same size (NEQ individuals, R_NEQ random replicates,
           per-SNP log-likelihoods averaged) so that large panels such as Spanish (n=173) do not
           win by sample size. Populations with fewer than NEQ people are skipped. Margins to
           the best population carry delete-one-block jackknife SEs over 5-Mb windows.
  3 LOO    Calibration: every reference individual of the key populations is removed from its
           own population/source counts and run through 1 and 2 like a family member. Gives the
           component distribution of real Sicilians, Greeks, Spaniards ..., how often each is
           assigned to its own population, and the yardstick for the family.
  3b Sims  Pedigree simulation: the family tree (16 great-great-grandparents of self) is
           rebuilt from reference individuals of proxy populations chosen from the documented
           birthplaces (family.json "ancestry": "founders" ref_pop, plus optional alternative
           "hypotheses"; founders removed from the source counts), passed through meiosis with
           Poisson crossovers at sex-specific rates (~0.95 cM/Mb paternal, ~1.5 cM/Mb maternal;
           Kong 2010 / Halldorsson 2019; random founder phase), and the simulated self/mother/
           father/maternal grandmother are analysed exactly like the real ones. This gives the
           expected component distribution under each pedigree hypothesis, including chance
           variation. Without founders in family.json the pedigree tests are skipped.
  2b Pair  Pairwise contrasts (the headline evidence): for a population pair A/B, per-SNP
           beta-binomial log-likelihood ratio log P(g|A) - log P(g|B) with both populations
           subsampled to the same size m = min(nA, nB) - 1 (R_PAIR replicates), reported per 10k
           SNPs. Calibrated by leave-one-out (the target is never in its own subsample): LOO
           accuracy and each family member's percentile among the LOO members of A and of B.
           Family values carry 5-Mb block-jackknife SEs. Also run on simulated relatives
           (step 3b) so observed values can be compared with each pedigree hypothesis.
           The 65-way nearest-population ranking (step 2) has 'attractor' populations (IBS,
           Bulgarian) and is kept only as a diagnostic, never shown as a person's population.
  4 Family self vs (mother+father)/2 per component (chromosome-jackknife SE); maternal grandfather
           inferred as 2*mother - grandmother; observed vs simulated pedigree hypotheses.
  4b Grandpa  The maternal grandfather's transmitted haplotype (work/grandpa/, reconstructed from
           the mother, grandmother and the children) is analysed as a haploid genome: admixture on 2h and
           haploid pairwise contrasts, against pseudo-haploid LOO reference individuals (one random
           allele per SNP) and against the simulated transmitted haplotype of each hypothesis.
  6 Display  results/ancestry_ho.json 'display': plain-language per-person statements built only
           from calibrated results (pairwise contrasts, clouds with >= MIN_CLOUD members,
           simulations). `--display-only` rebuilds it from an existing results file.
  5 'South Asian' test: the coarse 1000G analysis had GIH as its only eastern West-Eurasian
           source. Here the analogue is an Iranian source fitted with only W/NW European
           sources (E1) and again with Eastern-Mediterranean/Levant/Anatolian sources (E2).
Scope: the reference panel (source configurations, PAIRS/SIM_PAIRS, SUBPOPS, the leave-one-out and
typicality populations) is designed for Southern European ancestry, i.e. Italian and Iberian lines with
their Balkan/Eastern-Mediterranean neighbours. The plain-language Italy/Spain/Sicily statements are only
made for people whose overall profile places them in those reference clouds (or whose documented tree has
Italian/Iberian founders); anyone else is reported as not close to the groups this analysis tests.
Nothing is uploaded; genotypes are never printed or written for family members.
"""
import collections
import csv
import json
import math
import os
import re
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ancestry_ho_fetch as F  # noqa: E402

import dnalib as D  # noqa: E402
from dnalib import RESULTS, WORK  # noqa: E402
OUT_JSON = os.path.join(RESULTS, "ancestry_ho.json")
OUT_TSV = os.path.join(RESULTS, "ancestry_ho.tsv")
SELF, MOM, DAD, GMA = (D.person(r) for r in ("self", "mother", "father", "maternal_grandmother"))

LAMBDA = 2.0          # pseudo-alleles toward the panel mean (admixture source frequencies)
PRIOR = 2.0           # beta prior strength (alleles) for the predictive likelihood
NEQ = 8               # equal-n population size for the nearest-population likelihood
R_NEQ = 10            # subsample replicates
WIN_BP = 5_000_000    # jackknife window for log-likelihood margins
N_SIM = int(os.environ.get("HO_NSIM", 200))   # simulated families per pedigree hypothesis
R_PAIR = 10           # subsample replicates for pairwise contrasts (R_PAIR_SIM for simulations)
R_PAIR_SIM = 5
MIN_CLOUD = 20        # Mahalanobis cloud fits only against clouds with at least this many members
KIN_DROP = 0.0884     # KING-robust kinship above which (2nd degree or closer) one member is dropped
NO_RANK = 0.15        # hypothesis fits differing by less than this are not ranked
RIDGE = 0.05 ** 2     # ridge added to covariance matrices for Mahalanobis distances
EPS = 1e-4
SEED = int(D.FAMILY.get("ancestry", {}).get("seed", 12345))  # simulations; fixed per family for reproducibility
NON_WEST_EURASIAN = {"Yoruba", "YRI", "GWD", "LWK", "Mbuti", "Somali", "Han", "CHB", "Pima",
                     "Karitiana", "Mexico_Zapotec", "MXL"}
FAR = {"West_African", "East_Asian", "Native_American"}   # left out of Mahalanobis distances

# ---------------------------------------------------------------- sources
S_ITALY = ["Sicilian", "Italian_South"]
NC_ITALY = ["Italian_North", "Italian_Central", "TSI"]
IBERIA = ["Spanish", "IBS", "Spanish_North"]
BALKAN_AEGEAN = ["Greek", "Greek_WGA", "Greek_Crete", "Albanian"]
NW_EUROPE = ["English", "GBR", "CEU", "French", "Scottish"]
LEVANT = ["Lebanese", "Lebanese_Christian", "Lebanese_Muslim", "Palestinian", "Jordanian", "Syrian"]
N_AFRICA = ["Moroccan", "Tunisian", "Algerian", "Libyan", "Mozabite"]
W_AFRICA = ["Yoruba", "YRI"]
E_ASIA = ["Han", "CHB"]
NAT_AM = ["Pima", "Karitiana", "Mexico_Zapotec"]
E_EUROPE = ["Russian", "Ukrainian", "Czech", "Hungarian", "Croatian", "Serbian_Serb",
            "Bosnian_Serb", "Montenegrin_Serb", "Croatian_Serb"]
ANAT_CAUC_IRAN = ["Turkish", "Armenian", "Georgian", "Iranian", "Kurd"]
ANAT_CAUC = ["Turkish", "Armenian", "Georgian"]
# Sub-groups of Italian_South by AADR locality (they stay pooled in the admixture sources).
# Italian_South = 5 Campania + 6 Calabria + 10 Griko-speaking Salento samples + 1 other Apulian sample
# (in neither sub-group); minus excluded relatives.
SUBPOPS = {"Italian_South_CampCal": ("Italian_South", ("Campania", "Calabria", "Naples", "Salerno",
                                                       "Cosenza")),
           "Salento_Griko": ("Italian_South", ("Lecce",))}
PAIRS = [("Sicilian", "TSI"), ("Sicilian", "Italian_North"), ("Sicilian", "Italian_Central"),
         ("Sicilian", "Spanish"), ("Sicilian", "Greek"), ("Sicilian", "Italian_South"),
         ("Italian_South", "Italian_North"), ("Italian_South", "TSI"), ("Italian_South", "Spanish"),
         ("Italian_South_CampCal", "Italian_North"), ("Italian_South_CampCal", "TSI"),
         ("Italian_South_CampCal", "Sicilian"), ("Italian_South_CampCal", "Salento_Griko"),
         ("Italian_Central", "Italian_North"), ("TSI", "Spanish"), ("Italian_North", "Spanish"),
         ("TSI", "Greek"), ("Spanish", "Greek"), ("Albanian", "Greek"), ("Albanian", "Italian_South"),
         ("Croatian", "Italian_North"), ("IBS", "Spanish")]
SIM_PAIRS = [("Sicilian", "TSI"), ("Sicilian", "Italian_North"), ("Italian_South", "Italian_North"),
             ("Italian_South_CampCal", "Italian_North"), ("TSI", "Spanish"), ("Italian_North", "Spanish"),
             ("Spanish", "Greek"), ("Albanian", "Italian_South"), ("Croatian", "Italian_North")]
OUTGROUPS = [("Levant", LEVANT), ("North_African", N_AFRICA), ("West_African", W_AFRICA),
             ("East_Asian", E_ASIA), ("Native_American", NAT_AM)]

CONFIGS = {
    "A_core": {
        "about": "Main configuration: 5 European/Mediterranean sources + 5 out-groups",
        "sources": [("South_Italy_Sicily", S_ITALY), ("NorthCentral_Italy", NC_ITALY),
                    ("Iberia", IBERIA), ("Balkan_Aegean", BALKAN_AEGEAN),
                    ("NW_Europe", NW_EUROPE)] + OUTGROUPS},
    "B_extended": {
        "about": "A + Sardinian, Basque, East European/Slavic and Anatolia-Caucasus-Iran sources",
        "sources": [("South_Italy_Sicily", S_ITALY), ("NorthCentral_Italy", NC_ITALY),
                    ("Iberia", ["Spanish", "IBS"]), ("Basque", ["Basque", "Spanish_North"]),
                    ("Sardinian", ["Sardinian"]), ("Balkan_Aegean", BALKAN_AEGEAN),
                    ("NW_Europe", NW_EUROPE), ("East_Europe", E_EUROPE),
                    ("Anatolia_Caucasus_Iran", ANAT_CAUC_IRAN)] + OUTGROUPS},
    "C_coarse": {
        "about": "All Italy merged into one source; Iberia, Balkan-Aegean, NW Europe + out-groups",
        "sources": [("Italy", S_ITALY + NC_ITALY), ("Iberia", IBERIA),
                    ("Balkan_Aegean", BALKAN_AEGEAN), ("NW_Europe", NW_EUROPE)] + OUTGROUPS},
    "D_sicily_split": {
        "about": "A with Sicilian (n=11) and mainland Italian_South (n=20 after removing relatives) as separate sources",
        "sources": [("Sicilian", ["Sicilian"]), ("Italian_South_mainland", ["Italian_South"]),
                    ("NorthCentral_Italy", NC_ITALY), ("Iberia", IBERIA),
                    ("Balkan_Aegean", BALKAN_AEGEAN), ("NW_Europe", NW_EUROPE)] + OUTGROUPS},
    "E1_sas_like": {
        "about": "South-Asian test, analogue of the coarse 1000G panel: Iberia, Tuscan, NW Europe + "
                 "Iranian as the only eastern West-Eurasian source (GIH's role) + W Africa, E Asia, Nat Am",
        "sources": [("Iberia", ["IBS", "Spanish"]), ("Tuscan", ["TSI"]), ("NW_Europe", NW_EUROPE),
                    ("Iranian", ["Iranian"]), ("West_African", W_AFRICA), ("East_Asian", E_ASIA),
                    ("Native_American", NAT_AM)]},
    "E2_sas_eastmed": {
        "about": "South-Asian test: E1 + Southern Italy, Balkan-Aegean, Levant, Anatolia-Caucasus, N Africa",
        "sources": [("South_Italy_Sicily", S_ITALY), ("NorthCentral_Italy", NC_ITALY),
                    ("Iberia", IBERIA), ("Balkan_Aegean", BALKAN_AEGEAN), ("NW_Europe", NW_EUROPE),
                    ("Levant", LEVANT), ("Anatolia_Caucasus", ANAT_CAUC), ("Iranian", ["Iranian"]),
                    ("North_African", N_AFRICA), ("West_African", W_AFRICA),
                    ("East_Asian", E_ASIA), ("Native_American", NAT_AM)]},
}
MAIN_CONFIGS = ["A_core", "B_extended", "C_coarse", "D_sicily_split"]
LOO_POPS = ["Sicilian", "Italian_South", "TSI", "Italian_Central", "Italian_North", "Spanish",
            "Greek", "Greek_WGA", "Albanian", "Maltese"]
POP_REGION = {
    "Sicilian": "S_Italy", "Italian_South": "S_Italy",
    "TSI": "NC_Italy", "Italian_Central": "NC_Italy", "Italian_North": "NC_Italy",
    "Spanish": "Iberia", "IBS": "Iberia", "Spanish_North": "Iberia",
    "Greek": "Aegean", "Greek_WGA": "Aegean", "Greek_Crete": "Aegean", "Albanian": "Aegean",
    "Maltese": "Malta",
}
OWN_SOURCE = {
    "A_core": {"Sicilian": "South_Italy_Sicily", "Italian_South": "South_Italy_Sicily",
               "TSI": "NorthCentral_Italy", "Italian_Central": "NorthCentral_Italy",
               "Italian_North": "NorthCentral_Italy", "Spanish": "Iberia", "Greek": "Balkan_Aegean",
               "Greek_WGA": "Balkan_Aegean", "Albanian": "Balkan_Aegean"},
    "C_coarse": {"Sicilian": "Italy", "Italian_South": "Italy", "TSI": "Italy", "Italian_Central": "Italy",
                 "Italian_North": "Italy", "Spanish": "Iberia", "Greek": "Balkan_Aegean",
                 "Greek_WGA": "Balkan_Aegean", "Albanian": "Balkan_Aegean"},
}
OWN_SOURCE["B_extended"] = OWN_SOURCE["A_core"]
OWN_SOURCE["E2_sas_eastmed"] = OWN_SOURCE["A_core"]
OWN_SOURCE["D_sicily_split"] = dict(OWN_SOURCE["A_core"], Sicilian="Sicilian",
                                    Italian_South="Italian_South_mainland")
ITALIAN_SOURCES = {"South_Italy_Sicily", "NorthCentral_Italy", "Italy", "Sicilian", "Italian_South_mainland"}

# ---------------------------------------------------------------- pedigree hypotheses
# Ahnentafel numbers as in work/genealogy/pedigree.json: 1 self, 2 father, 3 mother, 7 maternal grandmother,
# 16..31 = self's 16 great-great-grandparents (founders), each mapped to a proxy reference population.
# All of it comes from family.json "ancestry" (optional; without it the pedigree tests are skipped):
#   "founders": {"16": {"ref_pop": "Spanish", ...}, ... "31": {...}}  -> the "gedcom" hypothesis
#   "hypotheses": {"gedcom": {"about": ...},                            -> optional description
#                  "<name>": {"about": ..., "founders": {"23": "<pop>"},  -> overrides on top of "gedcom"
#                             "control": {"role": "father", "tree": ..., "as_parent": ..., "contrast": "A|B",
#                                         "contrast_words": [A-side, B-side]},   -> optional display notes
#                             "brief": ...}}                                     -> "<brief> is rejected ..."
#   "proxy_notes": [...]                                                -> caveats on the proxy choice
FOUNDER_SLOTS = range(16, 32)
ANC = D.FAMILY.get("ancestry", {})


def _hypotheses():
    fd = {int(k): v.get("ref_pop") for k, v in ANC.get("founders", {}).items()}
    if any(not fd.get(n) for n in FOUNDER_SLOTS):
        return {}
    base = {n: fd[n] for n in FOUNDER_SLOTS}
    extra = ANC.get("hypotheses", {})
    out = {"gedcom": {"about": extra.get("gedcom", {}).get(
        "about", "GEDCOM birthplaces (work/genealogy/pedigree.json) mapped to proxy populations in family.json"),
        "founders": base}}
    for h, H in extra.items():
        if h != "gedcom":
            out[h] = dict(H, about=H.get("about", h),
                          founders={**base, **{int(k): v for k, v in H.get("founders", {}).items()}})
    return out


HYPOTHESES = _hypotheses()
SIM_NODES = {n: p for n, p in ((1, SELF), (2, DAD), (3, MOM), (7, GMA)) if p}
PROXY_NOTES = ANC.get("proxy_notes", [
    "Founders are drawn without replacement and removed from all source counts in each replicate"])


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


# ---------------------------------------------------------------- data
def load():
    c = F.load_cache()
    snps = c["snps"]                       # rsid chrom pos a1 a2 flipped in_family
    L = len(snps)
    G = np.empty((len(c["geno"]), L), dtype=np.uint8)
    for j, rec in enumerate(c["geno"]):
        G[j] = np.frombuffer(rec, dtype=np.uint8)
    c["geno"] = None
    configured = set(D.people())
    fam = {p: np.array([3 if v is None else v for v in col], dtype=np.uint8)
           for p, col in c["family"].items() if p in configured}
    # configured people with data who are not in the cache yet (e.g. a sibling's file added later)
    try:
        for d in sorted(configured):
            if d in fam or not os.path.isdir(os.path.join(F.RAW, d)):
                continue
            data, _ = F.load_person(d)
            if not data:
                continue
            g = np.full(L, 3, dtype=np.uint8)
            for i, (rs, ch, pos, a1, a2, fl, _inf) in enumerate(snps):
                v = data.get(rs)
                if not v or str(v[0]) != ch or str(v[1]) != pos:
                    continue
                gt = v[2].translate(F.COMP) if fl == "1" else v[2]
                if len(gt) == 2 and set(gt) <= {a1, a2}:
                    g[i] = gt.count(a1)
            if (g != 3).sum() > 10000:
                fam[d] = g
                log(f"added {d} from raw data ({int((g != 3).sum())} SNPs)")
    except Exception as e:  # pragma: no cover
        log("note: could not scan raw/ for new people:", e)
    chrom = np.array([int(s[1]) for s in snps])
    pos = np.array([int(s[2]) for s in snps])
    return snps, c["samples"], c["pops"], G, fam, chrom, pos, c.get("excluded", [])


def add_subpops(pops, samples):
    """Derived locality sub-groups (do not enter pbar, the equal-n ranking or source counts)."""
    for name, (base, keys) in SUBPOPS.items():
        pops[name] = [j for j in pops[base] if any(k in samples[j]["locality"] for k in keys)]
    return pops


def counts_of(Gc, idx):
    sub = Gc[np.asarray(idx)]
    called = sub != 3
    return np.where(called, sub, 0).sum(0, dtype=np.int64).astype(float), 2.0 * called.sum(0)


# ---------------------------------------------------------------- admixture
def admix(g, P, max_iter=300, tol=1e-9):
    """ML mixture proportions on the simplex; g: (L,) 0/1/2 floats, P: (L,K). Active-set Newton."""
    L, K = P.shape
    q = np.full(K, 1.0 / K)
    two_m_g = 2.0 - g

    def ll_of(q):
        f = np.clip(P @ q, 1e-12, 1 - 1e-12)
        return float(g @ np.log(f) + two_m_g @ np.log1p(-f))

    free = np.ones(K, bool)
    ll = ll_of(q)
    for _ in range(max_iter):
        f = np.clip(P @ q, 1e-12, 1 - 1e-12)
        a = g / f
        b = two_m_g / (1 - f)
        grad = P.T @ (a - b)
        w = a / f + b / (1 - f)
        Fi = np.where(free)[0]
        Pf = P[:, Fi]
        H = -(Pf * w[:, None]).T @ Pf
        try:
            Hi_g = np.linalg.solve(H, grad[Fi])
            Hi_1 = np.linalg.solve(H, np.ones(len(Fi)))
            d = (Hi_g.sum() / Hi_1.sum()) * Hi_1 - Hi_g
        except np.linalg.LinAlgError:
            d = grad[Fi] - grad[Fi].mean()
        step = 1.0
        qf = q[Fi]
        hit = (d < 0) & (qf + d < 0)
        if hit.any():
            step = float(np.min(-qf[hit] / d[hit]))
        while True:
            new = q.copy()
            new[Fi] = np.maximum(0.0, qf + step * d)
            new /= new.sum()
            ll_new = ll_of(new)
            if ll_new >= ll - 1e-9 or step < 1e-10:
                break
            step /= 2
        moved = float(np.abs(new - q).max())
        q, ll = new, ll_new
        free = q > 1e-12
        q[~free] = 0.0
        if moved < tol:
            f = np.clip(P @ q, 1e-12, 1 - 1e-12)
            grad = P.T @ (g / f - two_m_g / (1 - f))
            level = grad[free].mean()
            cand = np.where((~free) & (grad > level + 1e-6 * abs(level) + 1e-6))[0]
            if len(cand) == 0:
                break
            k = cand[np.argmax(grad[cand])]
            q *= (1 - 1e-4)
            q[k] = 1e-4
            free = q > 0
            ll = ll_of(q)
    return q / q.sum(), ll


class Panel:
    """Pooled counts per population on a SNP column subset; supports leaving people out."""

    def __init__(self, G, pops, cols):
        self.pops = pops
        self.Gc = G[:, cols]
        self.x, self.N = {}, {}
        for p, idx in pops.items():
            self.x[p], self.N[p] = counts_of(self.Gc, idx)
        fr = [(self.x[p] + .5) / (self.N[p] + 1) for p in pops
              if p not in NON_WEST_EURASIAN and p not in SUBPOPS]
        self.pbar = np.clip(np.mean(fr, axis=0), 0.01, 0.99)
        self.pop_of = {j: p for p, idx in pops.items() if p not in SUBPOPS for j in idx}

    def counts(self, members, drop=()):
        x = sum(self.x[p] for p in members)
        N = sum(self.N[p] for p in members)
        for j in drop:
            if self.pop_of[j] in members:
                gd = self.Gc[j]
                called = gd != 3
                x = x - np.where(called, gd, 0)
                N = N - 2.0 * called
        return x, N

    def freq(self, members, drop=()):
        x, N = self.counts(members, drop)
        return np.clip((x + LAMBDA * self.pbar) / (N + LAMBDA), EPS, 1 - EPS)

    def P(self, config, drop=()):
        return np.column_stack([self.freq(m, drop) for _n, m in CONFIGS[config]["sources"]])


def fit(panel, g, config, drop=(), jack_chrom=None, P=None):
    m = g != 3
    if P is None:
        P = panel.P(config, drop)
    q, ll = admix(g[m].astype(float), P[m])
    res = {"q": q, "ll": ll, "n_snps": int(m.sum())}
    if jack_chrom is not None:
        reps, ws = [], []
        for c in range(1, 23):
            nc = int((m & (jack_chrom == c)).sum())
            if nc == 0:
                continue
            mm = m & (jack_chrom != c)
            reps.append(admix(g[mm].astype(float), P[mm])[0])
            ws.append(nc)
        res["jk"] = np.array(reps)
        res["jk_w"] = np.array(ws, float)
        res["se"] = jk_se(res["jk"], res["jk_w"], q)
    return res


def jk_se(reps, sizes, full):
    """Weighted delete-one-group jackknife SE (Busing et al. 1999)."""
    n = sizes.sum()
    h = n / sizes
    g = len(sizes)
    est = g * full - ((1 - sizes / n)[:, None] * reps).sum(0)
    tau = h[:, None] * full - (h - 1)[:, None] * reps
    return np.sqrt((((tau - est) ** 2) / (h - 1)[:, None]).sum(0) / g)


# ---------------------------------------------------------------- nearest population (equal n)
def ll_table(x, N, pbar):
    al = x + PRIOR * pbar
    be = N - x + PRIOR * (1 - pbar)
    T = al + be
    den = np.log(T) + np.log(T + 1)
    return np.column_stack([np.log(be) + np.log(be + 1) - den,
                            math.log(2) + np.log(al) + np.log(be) - den,
                            np.log(al) + np.log(al + 1) - den]).astype(np.float32)


class EqualN:
    """Per-population predictive log-likelihood with all populations subsampled to NEQ people."""

    def __init__(self, Gc, pops, pbar, cache=True):
        rng = np.random.default_rng(SEED)
        self.Gc, self.pbar, self.pops = Gc, pbar, pops
        self.names = sorted(p for p in pops if len(pops[p]) >= NEQ and p not in SUBPOPS)
        self.subs = {p: [rng.choice(pops[p], NEQ, replace=False) for _ in range(R_NEQ)] for p in self.names}
        self.cache = cache
        self.tabs = {}

    def _tab(self, p, r, exclude=()):
        sub = self.subs[p][r]
        ex = [j for j in exclude if j in set(sub)]
        if ex:
            pool = [j for j in self.pops[p] if j not in set(exclude)]
            rng = np.random.default_rng(SEED + r + 7919 * int(ex[0]))
            sub = pool if len(pool) <= NEQ else rng.choice(pool, NEQ, replace=False)
            return ll_table(*counts_of(self.Gc, sub), self.pbar)
        if not self.cache:
            return ll_table(*counts_of(self.Gc, sub), self.pbar)
        if (p, r) not in self.tabs:
            self.tabs[(p, r)] = ll_table(*counts_of(self.Gc, sub), self.pbar)
        return self.tabs[(p, r)]

    def per_snp(self, g, exclude=()):
        called = g != 3
        gi = np.where(called, g, 0).astype(int)
        ar = np.arange(len(g))
        out = {}
        for p in self.names:
            acc = np.zeros(len(g))
            for r in range(R_NEQ):
                acc += self._tab(p, r, exclude)[ar, gi]
            out[p] = np.where(called, acc / R_NEQ, 0.0)
        return out


def ll_margins(llmap, blocks, called, top=None):
    tot = {p: float(v.sum()) for p, v in llmap.items()}
    best = max(tot, key=tot.get)
    ub, inv = np.unique(blocks[called], return_inverse=True)
    B = len(ub)
    bb = np.bincount(inv, weights=llmap[best][called], minlength=B)
    out = {}
    for p in sorted(tot, key=tot.get, reverse=True)[:top]:
        d = np.bincount(inv, weights=llmap[p][called], minlength=B) - bb
        loo = d.sum() - d
        se = math.sqrt((B - 1) / B * float(((loo - loo.mean()) ** 2).sum())) if p != best else 0.0
        out[p] = {"ll_minus_best": round(tot[p] - tot[best], 1), "se": round(se, 1),
                  "z": round((tot[p] - tot[best]) / se, 2) if se > 0 else None}
    return best, out, tot


# ---------------------------------------------------------------- pedigree simulation
class Meiosis:
    def __init__(self, chrom, pos, rng):
        self.rng = rng
        self.pos = pos
        self.spans = {}
        for c in range(1, 23):
            ix = np.where(chrom == c)[0]
            if len(ix):
                self.spans[c] = (ix[0], ix[-1] + 1, pos[ix[0]], pos[ix[-1]])

    def founder(self, g):
        het = g == 1
        r = self.rng.random(len(g)) < 0.5
        h0 = ((g == 2) | (het & r)).astype(np.int8)
        h1 = ((g == 2) | (het & ~r)).astype(np.int8)
        return h0, h1, g == 3

    def gamete(self, ind, morgan_bp):
        h0, h1, miss = ind
        out = np.empty_like(h0)
        for c, (a, b, p0, p1) in self.spans.items():
            ncross = self.rng.poisson((p1 - p0) / morgan_bp)    # 1 Morgan per morgan_bp
            cuts = np.sort(self.rng.uniform(p0, p1, ncross))
            seg = np.searchsorted(cuts, self.pos[a:b])
            use1 = ((seg + self.rng.integers(2)) % 2).astype(bool)
            out[a:b] = np.where(use1, h1[a:b], h0[a:b])
        return out, miss

    def child(self, fa, mo):
        g1, m1 = self.gamete(fa, 1.05e8)    # paternal ~0.95 cM/Mb, ~27 crossovers per sperm
        g2, m2 = self.gamete(mo, 0.67e8)    # maternal ~1.5 cM/Mb, ~42 crossovers per egg
        return g1, g2, m1 | m2


def simulate_family(panel, pops, founders, mei, rng):
    """Return ({node: genotype vector}, founder sample indices)."""
    chosen, used = {}, set()          # one global set: sub-groups overlap their parent population
    for n, p in sorted(founders.items()):
        pool = [j for j in pops[p] if j not in used]
        j = int(rng.choice(pool))
        used.add(j)
        chosen[n] = j
    ind = {n: mei.founder(panel.Gc[j]) for n, j in chosen.items()}
    for n in range(15, 0, -1):
        ind[n] = mei.child(ind[2 * n], ind[2 * n + 1])
    geno = {}
    for n in SIM_NODES:
        h0, h1, miss = ind[n]
        geno[n] = np.where(miss, 3, h0 + h1).astype(np.uint8)
    h0, _h1, miss = ind[3]                 # mom's paternal haplotype = grandpa's transmitted one
    geno["gp_hap"] = np.where(miss, 3, h0).astype(np.uint8)
    return geno, list(chosen.values())


# ---------------------------------------------------------------- QC and validation
def king_matrix(Gs, chunk=20000):
    """KING-robust kinship for all pairs of rows of Gs (uint8, 3 = missing)."""
    n = Gs.shape[0]
    hh = np.zeros((n, n)); opp = np.zeros((n, n)); hm = np.zeros((n, n))
    for a in range(0, Gs.shape[1], chunk):
        A = Gs[:, a:a + chunk]
        H = (A == 1).astype(np.float32)
        Z = (A == 0).astype(np.float32)
        T = (A == 2).astype(np.float32)
        M = (A != 3).astype(np.float32)
        hh += H @ H.T
        zt = Z @ T.T
        opp += zt + zt.T
        hm += H @ M.T                      # het in i where j is called
    den = hm + hm.T
    with np.errstate(divide="ignore", invalid="ignore"):
        k = np.where(den > 0, (hh - 2 * opp) / den, 0.0)
    return k, opp


def kinship_scan(G, pops, samples, excluded):
    """Within-population and within-region relatives (excluded samples included in the scan)."""
    base = {p: list(idx) for p, idx in pops.items() if p not in SUBPOPS}
    for j in excluded:
        base.setdefault(samples[j]["pop"], []).append(j)
    region = collections.defaultdict(list)
    for p, idx in base.items():
        region[F.POPS.get(p, "other")] += idx
    pairs = []
    for reg, idx in region.items():
        idx = sorted(idx)
        k, _ = king_matrix(G[idx])
        iu = np.triu_indices(len(idx), 1)
        for a, b in zip(*iu):
            if k[a, b] > KIN_DROP:
                i, j = idx[a], idx[b]
                pairs.append({"a": samples[i]["iid"], "b": samples[j]["iid"], "pop_a": samples[i]["pop"],
                              "pop_b": samples[j]["pop"], "kinship": rnd(k[a, b], 3),
                              "degree": "duplicate/MZ" if k[a, b] > 0.354 else
                                        "first" if k[a, b] > 0.177 else "second",
                              "ia": i, "ib": j})
    drop = []
    for pr in sorted(pairs, key=lambda d: -d["kinship"]):    # greedy: closest pairs first
        if pr["ia"] in excluded or pr["ib"] in excluded or pr["ia"] in drop or pr["ib"] in drop:
            continue
        na, nb = int(samples[pr["ia"]]["hosnps"] or 0), int(samples[pr["ib"]]["hosnps"] or 0)
        drop.append(pr["ib"] if nb < na else pr["ia"])
    for pr in pairs:
        pr["dropped"] = next((samples[j]["iid"] for j in (pr["ia"], pr["ib"]) if j in excluded or j in drop),
                             None)
        del pr["ia"], pr["ib"]
    return pairs, sorted(set(drop))


EURO_REGIONS = {"Italy", "Iberia", "WestEurope", "NorthEurope", "CentralMed", "Balkans", "EastEurope"}


def strand_calibration(G, pops, fam, snps):
    """Family dosage vs pooled European reference frequency (a strand flip would give slope ~ -1)."""
    idx = sorted(j for p, ix in pops.items() if p not in SUBPOPS and F.POPS.get(p) in EURO_REGIONS for j in ix)
    x = np.zeros(G.shape[1]); N = np.zeros(G.shape[1])
    for a in range(0, len(idx), 200):
        sub = G[idx[a:a + 200]]
        cl = sub != 3
        x += np.where(cl, sub, 0).sum(0); N += 2 * cl.sum(0)
    pf = np.where(N > 0, x / np.maximum(N, 1), np.nan)
    out = {"reference": f"pooled European/Mediterranean reference individuals (n={len(idx)})",
           "snps_flipped_to_match": sum(1 for s_ in snps if s_[5] == "1"),
           "palindromic_snps_kept": sum(1 for s_ in snps if {s_[3], s_[4]} in ({"A", "T"}, {"C", "G"})),
           "people": {}}
    for p, g in fam.items():
        m = (g != 3) & np.isfinite(pf)
        d, f = g[m] / 2.0, pf[m]
        slope = float(np.cov(d, f)[0, 1] / np.var(f, ddof=1))
        bins = np.minimum((f * 10).astype(int), 9)
        dev = [abs(float(d[bins == b].mean() - f[bins == b].mean())) for b in range(10) if (bins == b).any()]
        rare = (f < 0.01) | (f > 0.99)
        opp = int((((f < 0.01) & (g[m] == 2)) | ((f > 0.99) & (g[m] == 0))).sum())
        out["people"][p] = {"n_snps": int(m.sum()), "slope_dosage_on_ref_freq": rnd(slope, 3),
                            "r": rnd(np.corrcoef(d, f)[0, 1], 3),
                            "max_abs_bin_deviation_10bins": rnd(max(dev), 4),
                            "opposite_homozygotes_at_near_fixed": opp, "near_fixed_snps": int(rare.sum())}
    out["verdict"] = ("no strand/allele flips: slope ~ +1 for everyone"
                      if all(0.9 < v["slope_dosage_on_ref_freq"] < 1.1 for v in out["people"].values())
                      else "CHECK: slope far from +1 for someone")
    return out


def relationship_checks(fam):
    """Opposite homozygotes for every pair and Mendelian errors for parent-child trios."""
    ppl = list(fam)
    pairs = {}
    for i, a in enumerate(ppl):
        for b in ppl[i + 1:]:
            m = (fam[a] != 3) & (fam[b] != 3)
            o = int((((fam[a] == 0) & (fam[b] == 2)) | ((fam[a] == 2) & (fam[b] == 0)))[m].sum())
            pairs[f"{a}|{b}"] = {"opposite_homozygotes": o, "snps": int(m.sum()), "rate": float(f"{o / m.sum():.2g}")}
    trios = {}
    for c in D.children():
        if c in fam and MOM in fam and DAD in fam:
            g, mo, da = fam[c], fam[MOM], fam[DAD]
            m = (g != 3) & (mo != 3) & (da != 3)
            ok = np.zeros(len(g), bool)
            for x in (0, 1):
                for y in (0, 1):
                    pa = np.where(da == 2, 1, np.where(da == 0, 0, x))
                    pm = np.where(mo == 2, 1, np.where(mo == 0, 0, y))
                    ok |= (pa + pm) == g
            trios[f"{c}<-{MOM}+{DAD}"] = {"mendelian_errors": int((~ok & m).sum()), "snps": int(m.sum())}
    return {"note": "parent-child pairs expect ~0 opposite homozygotes (genotyping error only); unrelated "
                    "pairs ~6%", "pairs": pairs, "trios": trios}


GRANDPA_FILE = os.path.join(WORK, "grandpa", "grandpa_haplotype.tsv.gz")


def load_grandpa(snps):
    """Grandpa's transmitted haplotype on the HO SNPs: 1 = a1, 0 = a2, 3 = unknown."""
    if not os.path.exists(GRANDPA_FILE):
        return None
    import gzip
    where = {s_[0]: i for i, s_ in enumerate(snps)}
    h = np.full(len(snps), 3, dtype=np.uint8)
    with gzip.open(GRANDPA_FILE, "rt") as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            i = where.get(r["rsid"])
            al = r.get("grandpa_transmitted") or ""
            if i is None or len(al) != 1 or r.get("ambiguous_hethet", "0") not in ("0", ""):
                continue
            rs, ch, ps, a1, a2, fl = snps[i][:6]
            if r["chromosome"] != ch or r["position_grch37"] != ps:
                continue
            al = al.translate(F.COMP) if fl == "1" else al
            if al == a1:
                h[i] = 1
            elif al == a2:
                h[i] = 0
    return h


def pseudo_haploid(g, rng, mask=None):
    r = rng.random(len(g)) < 0.5
    h = np.where(g == 2, 1, np.where(g == 0, 0, np.where(g == 1, r.astype(np.uint8), 3))).astype(np.uint8)
    if mask is not None:
        h = np.where(mask, h, 3).astype(np.uint8)
    return h


def ll_table_hap(x, N, pbar):
    al = x + PRIOR * pbar
    be = N - x + PRIOR * (1 - pbar)
    T = al + be
    return np.column_stack([np.log(be / T), np.log(al / T)]).astype(np.float32)


class Pairwise:
    """Equal-n pairwise log-likelihood ratios with leave-out (targets/founders never in a subsample)."""

    def __init__(self, Gc, pops, pbar, pairs):
        rng = np.random.default_rng(SEED + 17)
        self.Gc, self.pops, self.pbar = Gc, pops, pbar
        self.pairs, self.subs = [], {}
        for a, b in pairs:
            m = min(len(pops[a]), len(pops[b])) - 1
            self.pairs.append((a, b, m))
            for p in (a, b):
                if (p, m) not in self.subs:
                    self.subs[(p, m)] = [np.sort(rng.choice(pops[p], m, replace=False)) for _ in range(R_PAIR)]
        self._cnt, self._tab = {}, {}
        self.ar = np.arange(Gc.shape[1])

    def _members(self, p, m, r, exclude):
        """Subsample r of population p with excluded people replaced (if the pool allows)."""
        sub = self.subs[(p, m)][r]
        ex = [j for j in exclude if j in sub]
        if not ex:
            return sub, [], []
        pool = np.setdiff1d(self.pops[p], np.concatenate([sub, np.asarray(list(exclude))]))
        rng = np.random.default_rng(SEED + 31 * r + 7919 * int(ex[0]))
        repl = [int(k) for k in rng.choice(pool, min(len(ex), len(pool)), replace=False)] if len(pool) else []
        return sub, ex, repl

    def size_after(self, p, m, r, exclude):
        _sub, ex, repl = self._members(p, m, r, exclude)
        return m - len(ex) + len(repl)

    def table(self, p, m, r, exclude=(), hap=False, memo=None, size=None):
        """Log-likelihood table; `size` trims the subsample so both sides of a contrast have equal n
        (unequal n biases the ratio toward the larger sample)."""
        mk = ll_table_hap if hap else ll_table
        sub, ex, repl = self._members(p, m, r, exclude)
        n_now = m - len(ex) + len(repl)
        trim = size is not None and size < n_now
        if (p, m, r) not in self._cnt:
            self._cnt[(p, m, r)] = counts_of(self.Gc, sub)
        if not ex and not trim:
            key = (p, m, r, hap)
            if key not in self._tab:
                self._tab[key] = mk(*self._cnt[(p, m, r)], self.pbar)
            return self._tab[key]
        mkey = (p, m, r, hap, size)
        if memo is not None and mkey in memo:
            return memo[mkey]
        x, N = self._cnt[(p, m, r)]
        drop = list(ex)
        if trim:
            keep = [int(j) for j in sub if j not in set(ex)] + repl
            rng = np.random.default_rng(SEED + 101 * r + n_now)
            drop += [int(j) for j in rng.choice(keep, n_now - size, replace=False)]
            repl = [k for k in repl if k not in set(drop)]
        for j, sign in [(j, -1) for j in drop if j in sub] + [(k, 1) for k in repl]:
            gj = self.Gc[j]
            cl = gj != 3
            x = x + sign * np.where(cl, gj, 0).astype(float)
            N = N + sign * 2.0 * cl
        t = mk(x, N, self.pbar)
        if memo is not None:
            memo[mkey] = t
        return t

    def per_snp(self, g, a, b, m, exclude=(), hap=False, reps=R_PAIR, memo=None):
        called = g != 3
        gi = np.where(called, g, 0).astype(int)
        acc = np.zeros(len(g))
        for r in range(reps):
            size = None
            if exclude:
                size = min(self.size_after(a, m, r, exclude), self.size_after(b, m, r, exclude))
            acc += (self.table(a, m, r, exclude, hap, memo, size)[self.ar, gi]
                    - self.table(b, m, r, exclude, hap, memo, size)[self.ar, gi])
        return np.where(called, acc / reps, 0.0)

    def value(self, g, a, b, m, **kw):
        v = self.per_snp(g, a, b, m, **kw)
        return float(v.sum()) / max(1, int((g != 3).sum())) * 1e4


def block_se(per, blocks, called):
    ub, inv = np.unique(blocks[called], return_inverse=True)
    B = len(ub)
    bs = np.bincount(inv, weights=per[called], minlength=B)
    loo = bs.sum() - bs
    return math.sqrt((B - 1) / B * float(((loo - loo.mean()) ** 2).sum())) / max(1, int(called.sum())) * 1e4


def pct(arr, v):
    return rnd((np.asarray(arr) <= v).mean(), 3)


# ---------------------------------------------------------------- helpers
def rnd(v, k=4):
    return round(float(v), k)


def comp_dict(names, q):
    return {n: rnd(x) for n, x in zip(names, q)}


def src_names(config):
    return [n for n, _m in CONFIGS[config]["sources"]]


def near_idx(config):
    return [i for i, n in enumerate(src_names(config)) if n not in FAR]


def mahal(v, cloud, loo_dist=True, min_n=MIN_CLOUD):
    """Distance of v to a cloud (rows) and the fraction of cloud members at least as far
    (each member measured against the cloud without itself). Clouds with fewer than MIN_CLOUD
    members are not used: their covariance is dominated by the ridge term."""
    if len(cloud) < min_n:
        return {"skipped": f"cloud has {len(cloud)} < {min_n} members; see pairwise contrasts"}
    mu = cloud.mean(0)
    k = cloud.shape[1]
    Ci = np.linalg.inv(np.cov(cloud.T) + np.eye(k) * RIDGE)
    d = float(np.sqrt((v - mu) @ Ci @ (v - mu)))
    dl = []
    if loo_dist and len(cloud) > 3:
        for i in range(len(cloud)):
            rest = np.delete(cloud, i, axis=0)
            mi = rest.mean(0)
            Cii = np.linalg.inv(np.cov(rest.T) + np.eye(k) * RIDGE)
            dl.append(float(np.sqrt((cloud[i] - mi) @ Cii @ (cloud[i] - mi))))
    dl = np.array(dl)
    f = float((dl >= d).mean()) if len(dl) else None
    return {"distance": round(d, 2),
            "frac_members_further": rnd(f, 3) if len(dl) else None,
            "frac_se": rnd(math.sqrt(max(f * (1 - f), 1 / len(dl)) / len(dl)), 3) if len(dl) else None,
            "members_median_distance": round(float(np.median(dl)), 2) if len(dl) else None}


def broad(cfg, q):
    c = dict(zip(src_names(cfg), q))
    ital = sum(v for k, v in c.items() if k in ITALIAN_SOURCES)
    ib = c.get("Iberia", 0) + c.get("Basque", 0)
    return {"Italian_total": rnd(ital), "Iberian_total": rnd(ib),
            "Balkan_Aegean": rnd(c.get("Balkan_Aegean", 0)), "NW_Europe": rnd(c.get("NW_Europe", 0)),
            "other": rnd(1 - ital - ib - c.get("Balkan_Aegean", 0) - c.get("NW_Europe", 0))}


# ---------------------------------------------------------------- display
SIM_NODE_FOR = {s: SELF for s in D.siblings()}     # a full sibling has the same expectation as self
NODE_OF_ROLE = {"self": 1, "sibling": 1, "father": 2, "mother": 3, "maternal_grandmother": 7}
# LOO populations used for the plain-language typicality wording (group name, adjective).
TYP_POPS = {"Italian_South": ("Southern Italians", "Southern Italian"), "TSI": ("Tuscans", "Tuscan"),
            "Italian_North": ("Northern Italians", "Northern Italian"), "Spanish": ("Spaniards", "Spanish")}
TYPICAL = 0.1           # share of reference people more unusual at/above which a profile counts as typical
MENDEL_MAX = 0.005      # trio Mendelian-error rate below which errors are genotyping error only
PC_OPP_MAX = 0.001      # parent-child opposite-homozygote rate at/below which it is genotyping error only
UNREL_OPP_MIN = 0.04    # opposite-homozygote rate at/above which a pair is not closely related (unrelated ~6%)
MID_Z_MAX = 2.0         # self vs midparent: largest |difference| / SE still counted as chance
MID_SE_FLOOR = 0.005    # floor on that SE (components pinned at 0 have a zero jackknife SE)
ITALIAN_POPS = set(S_ITALY + NC_ITALY + list(SUBPOPS) + ["Sardinian"])
IBERIAN_POPS = set(IBERIA + ["Basque"])


def _label(p):
    """Display label (dnalib label, with the maternal grandmother spelled out)."""
    return "Grandma (mom's mother)" if D.role(p) == "maternal_grandmother" else D.label(p)


def _short(p):
    """'me', 'mom', 'dad', 'grandma', 'brother', ... for running text."""
    return D.label(p).split(" (")[0].lower()


def _pron(p):
    """(subject, possessive) by sex, falling back on the role; neutral 'they/their' when unknown."""
    sx = D.sex(p) or {"mother": "F", "maternal_grandmother": "F", "father": "M"}.get(D.role(p))
    return {"M": ("he", "his"), "F": ("she", "her")}.get(sx, ("they", "their"))


def _join(words, last="and"):
    return ", ".join(words[:-1]) + f" {last} " + words[-1] if len(words) > 1 else "".join(words)


def _tree_mix(out, who):
    """(Iberian, Italian) share of a person's great-great-grandparent slots in the documented tree, or None."""
    fd = out.get("pedigree_simulation", {}).get("hypotheses", {}).get("gedcom", {}).get("founders")
    n = 6 if who == "grandpa" else NODE_OF_ROLE.get(D.role(who))
    if not fd or not n:
        return None
    fd = {int(k): v for k, v in fd.items()}
    k = 2 ** (4 - (n.bit_length() - 1))
    slots = [fd.get(i) for i in range(n * k, (n + 1) * k)]
    return sum(p in IBERIAN_POPS for p in slots) / k, sum(p in ITALIAN_POPS for p in slots) / k


def _pc(x):
    return f"{100 * x:.0f}%"


def _pair(out, key, who):
    pr = out.get("pairwise", {}).get("pairs", {}).get(key)
    if not pr:
        return None
    if who == "grandpa":
        g = pr.get("grandpa_haploid")
        if not g:
            return None
        return {"v": g["llr_per10k"], "se": g["se"], "pa": g["pct_among_a_pseudohaploid"],
                "pb": g["pct_among_b_pseudohaploid"], "na": pr["n_a"], "nb": pr["n_b"],
                "acc": g["loo_accuracy_pseudohaploid"]}
    f = pr["family"].get(who)
    if not f:
        return None
    return {"v": f["llr_per10k"], "se": f["se"], "pa": f["pct_among_a"], "pb": f["pct_among_b"],
            "na": pr["n_a"], "nb": pr["n_b"], "acc": pr["loo_accuracy"]}


def _num(c):
    return f"{c['v']:+.0f} \u00b1 {c['se']:.0f}"


def _sim_z(out, who, key, hyp="gedcom"):
    fl = out.get("family_logic", {})
    if who == "grandpa":
        d = fl.get("grandpa_transmitted_haplotype", {}).get("pairwise_vs_simulation", {}).get(hyp, {})
    else:
        node = SIM_NODE_FOR.get(who, who)
        d = fl.get("observed_vs_pedigree_simulations", {}).get(node, {}).get("pairwise_vs_simulation", {}).get(hyp, {})
        if who != node:          # observed value of the sibling, expectation of self
            c = _pair(out, key, who)
            e = d.get(key)
            if c and e:
                return {"observed": c["v"], "sim_mean": e["sim_mean"], "sim_sd": e["sim_sd"],
                        "z": round((c["v"] - e["sim_mean"]) / e["sim_sd"], 2)}
    return d.get(key)


def _fits(out, who, cfg="A_core"):
    node = SIM_NODE_FOR.get(who, who)
    return out.get("family_logic", {}).get("hypothesis_fits", {}).get("by_person", {}).get(node, {}).get(cfg, [])


def _placed(typ, mix):
    """True when the Italy/Spain wording applies: the overall profile is typical of at least one Italian/Iberian
    reference cloud, or the documented tree has Italian/Iberian founders; None when neither is known."""
    t = [typ[k] for k in TYP_POPS if k in (typ or {})]
    if not t and not mix:
        return None
    return bool(t and max(t) >= TYPICAL) or bool(mix and (mix[0] or mix[1]))


def _person_display(out, who, typ=None):
    S = []          # (statement, evidence) pairs
    nums = {}
    flags = {}
    it = _pair(out, "TSI|Spanish", who)
    itn = _pair(out, "Italian_North|Spanish", who)
    mix = _tree_mix(out, who)
    placed = _placed(typ, mix)
    if it and not placed:
        nums["Tuscan-vs-Spanish contrast"] = f"{_num(it)} (more Italian-like than {_pc(it['pb'])} of Spaniards)"
        if placed is False:
            flags["italy"] = "outside"
            S.append(("Not close to the Italian/Iberian reference groups this analysis tests.",
                      ("Overall-profile typicality (share of reference people more unusual): "
                       + ", ".join(f"{TYP_POPS[k][0]} {typ[k]:.2f}" for k in TYP_POPS if k in typ)
                       + f", below {TYPICAL} for all of them; " if any(k in (typ or {}) for k in TYP_POPS) else "")
                      + ("no Italian/Iberian founders in the documented tree" if mix else "no documented tree founders")
                      + ", so the Italy/Spain comparisons are not interpreted."))
    elif it:
        nums["Tuscan-vs-Spanish contrast"] = (f"{_num(it)} (more Italian-like than {_pc(it['pb'])} of Spaniards "
                                              f"and than {_pc(it['pa'])} of Tuscans)")
        if it["pb"] >= 0.99 and it["pa"] >= 0.05:
            flags["italy"] = "italian"
            S.append(("Clearly Italian, not Spanish.",
                      f"On the Tuscan-vs-Spanish comparison (people placed correctly {_pc(it['acc'])} of the time) "
                      f"this genome is more Italian-like than all {it['nb']} Spanish reference people."))
        elif it["pb"] >= 0.99:
            flags["italy"] = "between"
            mix = _tree_mix(out, who)
            why = ", as expected for a mixed Italian/Spanish family" if mix and mix[0] and mix[1] else ""
            S.append((f"In between Italian and Spanish{why}.",
                      f"More Italian-like than every one of the {it['nb']} Spaniards, but less Italian-like than "
                      f"{_pc(1 - it['pa'])} of the {it['na']} Tuscans."))
        elif it["pb"] >= 0.85:
            flags["italy"] = "spanish_shifted"
            S.append(("Mostly Spanish-like, with a clear pull toward Italy.",
                      f"More Italian-like than {_pc(it['pb'])} of the {it['nb']} Spanish reference people "
                      f"(Tuscan-vs-Spanish {_num(it)})."))
        else:
            flags["italy"] = "spanish"
            S.append(("Looks like a typical Spaniard on the Italy-vs-Spain comparison.",
                      f"{_pc(it['pb'])} of Spaniards are less Italian-like."))
    if itn:
        nums["Northern-Italian-vs-Spanish contrast"] = f"{_num(itn)} (above {_pc(itn['pb'])} of Spaniards)"
    # Sicily / Southern Italy
    sic = {k: _pair(out, f"Sicilian|{k}", who) for k in ("TSI", "Italian_North", "Italian_Central", "Spanish")}
    sic = {k: v for k, v in sic.items() if v}
    if sic:
        low = sum(1 for v in sic.values() if v["pa"] == 0)
        edge = sum(1 for v in sic.values() if v["pa"] <= 0.1)
        n = next(iter(sic.values()))["na"]
        nums["Sicilian contrasts (vs Tuscan / N Italian / Central Italian / Spanish)"] = "; ".join(
            f"{_num(v)}, {_pc(v['pa'])} of Sicilians lower" for v in sic.values())
        flags["sicily_low"] = low
        if flags.get("italy") in ("italian", "between") or who == "grandpa":
            if edge == len(sic) and low >= 1:
                a0, a1 = min(v["acc"] for v in sic.values()), max(v["acc"] for v in sic.values())
                S.append(("Not a typical Sicilian.",
                          f"On all {len(sic)} Sicily-vs-elsewhere comparisons (each places {_pc(a0)}-{_pc(a1)} of "
                          f"reference people in the right group) this genome is at the bottom of the {n} Sicilian reference people "
                          f"(less Sicilian-like than all of them on {low}, than all but one on the rest). "
                          "Partial Sicilian ancestry is still possible."))
            elif low == 0 and all(v["pa"] >= 0.1 for v in sic.values()):
                S.append(("Within the range of real Sicilians.",
                          f"Inside the Sicilian reference range on all {len(sic)} Sicily comparisons."))
    so = _pair(out, "Italian_South|Italian_North", who)
    sot = _pair(out, "Italian_South|TSI", who)
    if so:
        nums["Southern-vs-Northern Italian contrast"] = (f"{_num(so)} ({_pc(so['pa'])} of Southern Italians lower, "
                                                        f"{_pc(so['pb'])} of Northern Italians lower)")
        flags["south_pa"] = so["pa"]
    if sot:
        nums["Southern-Italian-vs-Tuscan contrast"] = (f"{_num(sot)} ({_pc(sot['pa'])} of Southern Italians lower, "
                                                      f"{_pc(sot['pb'])} of Tuscans lower)")
    if so and flags.get("italy") in ("italian",):
        if so["pa"] <= 0.1:
            S.append(("Sits with Central/Northern Italians rather than with Southern Italians as a whole.",
                      f"On Southern-vs-Northern Italian (correct {_pc(so['acc'])} of the time) the value {_num(so)} "
                      f"is lower than {_pc(1 - so['pa'])} of the {so['na']} Southern Italians and in the middle of "
                      f"Northern Italians ({_pc(so['pb'])} lower)."))
    cn = _pair(out, "Italian_Central|Italian_North", who)
    if cn:
        nums["Central-vs-Northern Italian contrast (poorly separable)"] = f"{_num(cn)} (accuracy {_pc(cn['acc'])})"
    # eastern Mediterranean / Balkan shift
    sg = _pair(out, "Spanish|Greek", who)
    tg = _pair(out, "TSI|Greek", who)
    if sg and flags.get("italy") in ("spanish_shifted", "spanish", "between"):
        nums["Spanish-vs-Greek contrast"] = f"{_num(sg)} (more Greek-ward than {_pc(1 - sg['pa'])} of Spaniards)"
        if sg["pa"] <= 0.02:
            flags["emed"] = True
            # Italian lines in the documented tree are a sufficient explanation; any further reading (e.g. which
            # line, extra Balkan ancestry) is optional family.json text "ancestry_ho.greek_pull"
            mix = _tree_mix(out, who)
            why = (["Italian ancestry alone can explain this"] if mix and mix[1] > 0 else []) + \
                [x for x in [D.text("ancestry_ho.greek_pull", "")] if x]
            why = "; ".join(why)
            S.append(("Has an Eastern-Mediterranean/Balkan-like pull that typical Spaniards lack.",
                      f"More Greek-ward than {_pc(1 - sg['pa'])} of the {sg['na']} Spaniards."
                      + (f" {why[0].upper()}{why[1:]}." if why else "")))
    if tg and (flags.get("italy") == "italian" or who == "grandpa"):
        nums["Tuscan-vs-Greek contrast"] = f"{_num(tg)} ({_pc(tg['pa'])} of Tuscans lower)"
    # pedigree checks
    tree = []
    fits = _fits(out, who)
    if fits and who != "grandpa":
        best = fits[0]
        good = [f for f in fits if f["not_distinguishable_from_best"]]
        bad = [f for f in fits if f["fit"] < 0.05]
        tree.append(("Family-tree simulations (admixture profile, main configuration): "
                     + "; ".join(f"{f['hypothesis']} {f['fit']:.2f}\u00b1{f['se']:.2f}" for f in fits)
                     + f". Hypotheses within {NO_RANK} of the best ({', '.join(f['hypothesis'] for f in good)}) "
                     "cannot be ranked against each other"
                     + (f"; rejected (fit < 0.05): {', '.join(f['hypothesis'] for f in bad)}." if bad else "."),))
    zs = {}
    for k in ("Sicilian|TSI", "Sicilian|Italian_North", "Italian_South|Italian_North",
              "Italian_South_CampCal|Italian_North", "TSI|Spanish", "Spanish|Greek",
              "Croatian|Italian_North", "Albanian|Italian_South"):
        z = _sim_z(out, who, k)
        if z:
            zs[k] = z
    if zs:
        tree.append(("Observed vs the documented tree (GEDCOM simulation, same contrasts): "
                     + "; ".join(f"{k} observed {z['observed']:+.0f}, expected {z['sim_mean']:+.0f}\u00b1"
                                 f"{z['sim_sd']:.0f} (z {z['z']:+.1f})" for k, z in zs.items()) + ".",))
        flags["gedcom_z"] = {k: z["z"] for k, z in zs.items()}
    return S, nums, tree, flags


def build_display(out):
    """Plain-language, calibrated per-person statements (reconciles the verification reviews)."""
    P = out.get("pairwise", {}).get("pairs", {})
    acc = {k: v["loo_accuracy"] for k, v in P.items()}
    disp = {"how_to_read": (
        "Each statement rests on comparisons that were first tested on real reference people with "
        "known origins (leave-one-out). Numbers are log-likelihood ratios per 10,000 SNPs ± a block-"
        "jackknife standard error; positive means closer to the first-named population. 'X% of Sicilians "
        "lower' places the person among real Sicilians. Ancestry percentages from the admixture model are "
        "not shown as headline numbers: they are biased (they under-count Southern Italian ancestry even in "
        "real Sicilians) and are only used through simulations of the family tree."),
        "what_can_and_cannot_be_told_apart": [], "people": {}}
    sep = disp["what_can_and_cannot_be_told_apart"]
    good = [k for k in ("Sicilian|TSI", "Sicilian|Italian_North", "Sicilian|Spanish", "Sicilian|Greek",
                        "Italian_South|Italian_North", "Italian_South|TSI", "TSI|Spanish", "Spanish|Greek",
                        "Italian_North|Spanish", "Croatian|Italian_North") if acc.get(k, 0) >= 0.85]
    weak = [k for k in ("Sicilian|Italian_South", "Italian_South_CampCal|Sicilian", "Italian_South_CampCal|Salento_Griko",
                        "Italian_Central|Italian_North", "Albanian|Greek", "IBS|Spanish") if k in acc and acc[k] < 0.85]
    sep.append("Separable for one person (share of reference people placed correctly): "
               + ", ".join(f"{k.replace('|', ' vs ')} {acc[k]:.2f}" for k in good) + ".")
    sep.append("NOT reliably separable: " + ", ".join(f"{k.replace('|', ' vs ')} {acc[k]:.2f}" for k in weak)
               + ". So Sicily vs mainland Southern Italy (Campania/Calabria/Apulia), Albanian vs Greek, and "
                 "Central vs Northern Italy cannot be resolved for an individual.")
    # remarks on earlier versions of the report (family.json text)
    sep.extend(D.text("ancestry_ho.separability_notes", []))
    fl = out.get("family_logic", {})
    fits_all = {w: {f["hypothesis"]: f for f in _fits(out, w)} for w in out.get("people", {})}

    def fitv(w, h):
        f = fits_all.get(w, {}).get(h)
        return f"{f['fit']:.2f}±{f['se']:.2f}" if f else "n/a"

    def fitnum(w, h):
        f = fits_all.get(w, {}).get(h)
        return f["fit"] if f else None

    def typ_of(ts):
        """{pop: share of (pseudo-haploid) LOO members more unusual} for the TYP_POPS clouds that were fitted."""
        return {k: ts[k]["frac_members_further"] for k in TYP_POPS if "frac_members_further" in ts.get(k, {})}

    def fill(tpl, who, node, typ):
        """Optional family.json note with computed fields: {fit:HYP}, {fit:HYP@role}, {z:PAIR}, {z:PAIR@HYP},
        {acc:PAIR}, {typical:POP:group name} (a 'Typical of ...' sentence from the overall-profile typicality)."""
        def one(m):
            kind, arg = m.groups()
            if kind == "fit":
                h, _, r = arg.partition("@")
                return fitv(SIM_NODE_FOR.get(D.person(r), D.person(r)) if r else node, h)
            if kind == "z":
                k, _, h = arg.partition("@")
                return f"{(_sim_z(out, who, k, h or 'gedcom') or {}).get('z', 0):+.1f}"
            if kind == "acc":
                return f"{acc.get(arg, 0):.2f}"
            pop, _, words = arg.partition(":")
            tn, her = typ.get(pop), _pron(who)[1]
            if tn is None:
                return ""
            if tn >= TYPICAL:
                return f"Typical of {words} in {her} overall profile (share of {words} more unusual {tn:.2f}). "
            return f"Not especially typical of {words} in {her} overall profile (share more unusual {tn:.2f}). "
        return re.sub(r"\{(fit|z|acc|typical):([^{}]+)\}", one, tpl)

    def headline(who, flags):
        it = flags.get("italy")
        parts = {"italian": "Italian (clearly not Spanish)",
                 "between": "A mix of Italian and Spanish, sitting between the two populations",
                 "spanish_shifted": "Mostly Spanish, with clear Italian ancestry"
                                    + (" and some Eastern-Mediterranean-like shift" if flags.get("emed") else ""),
                 "spanish": "Spanish-like",
                 "outside": "Not close to the Italian/Iberian reference groups this analysis tests"}.get(it, "")
        if it == "italian" and flags.get("south_pa", 1) <= 0.1:
            parts += "; closest to Central/Northern Italians"
        if flags.get("sicily_low") and it in ("italian", "between"):
            parts += "; not a typical Sicilian"
        g = fits_all.get(SIM_NODE_FOR.get(who, who), {}).get("gedcom")
        if g:
            parts += ("; consistent with the documented family tree" if g["fit"] >= 0.1
                      else "; does not match the documented family tree well")
        parts = parts.removeprefix("; ")
        return parts[:1].upper() + parts[1:] + "." if parts else ""

    ibs_open = []
    italy_flags = {}
    for who in out.get("people", {}):
        pr = out["people"][who]
        typ = pr.get("typicality_vs_loo_populations", {}).get("A_core", {})
        typ = {k: v["frac_members_further"] for k, v in typ.items() if "frac_members_further" in v}
        S, nums, tree, flags = _person_display(out, who, typ)
        italy_flags[who] = flags.get("italy")
        rl = D.role(who)
        b = pr.get("broad_by_config", {}).get("A_core", {})
        node = SIM_NODE_FOR.get(who, who)
        simb = (out.get("pedigree_simulation", {}).get("hypotheses", {}).get("gedcom", {}).get("people", {})
                .get(node, {}).get("A_core", {}).get("broad_mean", {}))
        if typ:
            nums["Overall-profile typicality (share of reference people more unusual; groups of 20+ only)"] = \
                ", ".join(f"{k} {v:.2f}" for k, v in typ.items())
        if b and "Italian_total_se" in b and simb:
            nums["Admixture 'Italian' total (biased low; compare with the tree's simulated value)"] = (
                f"{b['Italian_total']:.2f} ± {b['Italian_total_se']:.2f} vs {simb['Italian_total']:.2f} "
                "expected from the GEDCOM tree")
        notes = []
        # alternative trees described in family.json ("brief"): reported when this person's DNA rejects them
        for h, H in HYPOTHESES.items():
            f = fitnum(node, h)
            if H.get("brief") and f is not None and f < 0.05:
                g = fitnum(node, "gedcom")
                notes.append(f"{H['brief']} is rejected for this person (fit {fitv(node, h)}); the documented "
                             f"GEDCOM tree is {'not' if g is None or g >= 0.05 else 'rejected too'} "
                             f"(fit {fitv(node, 'gedcom')}).")
        z = flags.get("gedcom_z", {})
        if who in D.children():
            mid = fl.get("me_vs_midparent", {}).get("A_core", {})
            rel = out.get("validation", {}).get("family_relationships", {}).get("trios", {}).get(f"{who}<-{MOM}+{DAD}")
            if rel and rel["snps"]:
                if rel["mendelian_errors"] / rel["snps"] < MENDEL_MAX:
                    notes.append(f"Confirmed child of {_short(MOM)} and {_short(DAD)} ({rel['mendelian_errors']} "
                                 f"inconsistent genotypes out of {rel['snps']:,}, i.e. genotyping error only).")
                else:
                    notes.append(f"CHECK: {rel['mendelian_errors']:,} genotypes out of {rel['snps']:,} "
                                 f"({100 * rel['mendelian_errors'] / rel['snps']:.1f}%) are inconsistent with being "
                                 f"the child of {_short(MOM)} and {_short(DAD)}, more than genotyping error "
                                 f"(< {100 * MENDEL_MAX:.1f}%); check sample labels before reading the trio results.")
            # control trees for one parent's side ("control": {"role": ...}): self's profile vs that parent's
            for h, H in HYPOTHESES.items():
                cp = D.person((H.get("control") or {}).get("role", ""))
                if who != SELF or not cp or fitnum(SELF, h) is None or fitnum(cp, h) is None:
                    continue
                what, par = H["control"].get("as_parent", f"the '{h}' tree"), _short(cp)
                if fitnum(SELF, h) >= 0.05 and fitnum(cp, h) < 0.05:
                    notes.append(f"My own profile cannot tell the variants of {par}'s side apart (it even fits "
                                 f"{what}, {fitv(SELF, h)}); {par}'s own DNA rejects that ({fitv(cp, h)}), so "
                                 f"{_pron(cp)[1]} result is the decisive one.")
                else:
                    notes.append(f"Under {what}: my own profile fits {fitv(SELF, h)}, {par}'s own DNA "
                                 f"{fitv(cp, h)}.")
            if mid and who == SELF:
                # each component's difference in units of max(jackknife SE, simulated SD of the difference)
                dz = {k: abs(d) / max(mid.get("jackknife_se", {}).get(k, 0),
                                      mid.get("simulated_sd_of_difference", {}).get(k, 0), MID_SE_FLOOR)
                      for k, d in mid.get("me_minus_midparent", {}).items()}
                kz = max(dz, key=dz.get) if dz else None
                if kz is None or dz[kz] < MID_Z_MAX:
                    notes.append(f"The ancestry profile matches the average of {_short(MOM)} and {_short(DAD)} "
                                 f"(largest component difference {mid['max_abs_diff']:.2f}, within chance).")
                else:
                    notes.append(f"CHECK: the ancestry profile departs from the average of {_short(MOM)} and "
                                 f"{_short(DAD)} by {dz[kz]:.1f} SE ({kz} "
                                 f"{mid['me_minus_midparent'][kz]:+.2f}; largest component difference "
                                 f"{mid['max_abs_diff']:.2f}), more than chance; check sample labels.")
            if z.get("Italian_South|Italian_North", 0) < -2:
                pars = [q for q in (MOM, DAD) if q]
                par = {q: (_sim_z(out, q, "Italian_South|Italian_North") or {}).get("z") for q in pars}
                if len(pars) == 2 and all(par[q] is not None and par[q] <= -1 for q in pars):
                    fam = (f"Both parents lean the same way ({_short(MOM)} z {par[MOM]:+.1f}, {_short(DAD)} z "
                           f"{par[DAD]:+.1f}), so this is a family-wide Northern lean relative to the tree, not "
                           "something specific to this person")
                else:
                    fam = "Parents: " + ", ".join(f"{_short(q)} z {par[q]:+.1f}" if par[q] is not None
                                                  else f"{_short(q)} n/a" for q in pars)
                notes.append("Less 'mainland Southern Italian' than the tree predicts (Southern-vs-Northern z "
                             f"{z['Italian_South|Italian_North']:+.1f}; against Campania+Calabria only z "
                             f"{z.get('Italian_South_CampCal|Italian_North', 0):+.1f}). " + fam + "; Sicily contrasts "
                             f"vs the tree z {z.get('Sicilian|TSI', 0):+.1f}, "
                             f"{z.get('Sicilian|Italian_North', 0):+.1f}.")
        if rl == "mother":
            gp = fl.get("grandpa_transmitted_haplotype")
            zsn = z.get('Italian_South|Italian_North')
            g = fitnum(who, "gedcom")
            if g is not None:
                z1, z2 = z.get('Sicilian|TSI', 0), z.get('Sicilian|Italian_North', 0)
                notes.append(f"Her mix is {'consistent' if g >= 0.05 else 'not consistent'} with the documented "
                             f"tree (fit {fitv(who, 'gedcom')}), and Sicily contrasts "
                             f"{'match it' if max(abs(z1), abs(z2)) < 2 else 'deviate from it'} (z "
                             f"{z1:+.1f}, {z2:+.1f})"
                             + ("" if zsn is None else
                                f"; on the pooled Southern-vs-Northern contrast she is about {abs(zsn):.0f} SD "
                                + ("more Northern" if zsn < 0 else "more Southern") + " than the tree predicts")
                             + " (Campania+Calabria only: "
                             f"z {z.get('Italian_South_CampCal|Italian_North', 0):+.1f}).")
            if gp:
                fr = typ_of(gp["typicality_vs_pseudohaploid_loo"]["A_core"])
                hit = [k for k in fr if fr[k] >= TYPICAL]
                if len(hit) == 1 and len(fr) > 1:
                    t = hit[0]
                    txt = (f"Her {TYP_POPS[t][1]} ancestry shows up on her father's side: the half of her DNA she "
                           f"got from him has a {TYP_POPS[t][1].replace(' ', '-')}-type profile (typical of "
                           f"{TYP_POPS[t][0]}: {fr[t]:.2f}; atypical of "
                           + ", ".join(f"{TYP_POPS[k][0]} {v:.2f}" for k, v in fr.items() if k != t) + ")")
                else:
                    txt = ("The half of her DNA she got from her father, compared with reference people (share more "
                           "unusual: " + ", ".join(f"{TYP_POPS[k][0]} {v:.2f}" for k, v in fr.items()) + ")")
                dg = gp.get("grandpa_minus_grandma_A_core", {})
                if "NW_Europe" in dg.get("diff", {}):
                    d, e = dg["diff"]["NW_Europe"], dg["se"]["NW_Europe"]
                    num = f"({d:+.2f} ± {e:.2f})"
                    if d + 2 * e < 0:
                        txt += (f" and has much less NW-European-like ancestry than {_short(GMA)}'s {num}. "
                                "Her mother's side is Northern-leaning.")
                    elif d - 2 * e > 0:
                        txt += f" and has more NW-European-like ancestry than {_short(GMA)}'s {num}."
                    else:
                        txt += f"; NW-European-like ancestry similar to {_short(GMA)}'s {num}."
                else:
                    txt += "."
                notes.append(txt)
        # control tree for this person's own side ("control": {"role": <this role>})
        for h, H in HYPOTHESES.items():
            c = H.get("control") or {}
            f, g = fitnum(who, h), fitnum(who, "gedcom")
            if not rl or c.get("role") != rl or f is None:
                continue
            ka, kb = c.get("contrast", "Spanish|Greek").split("|")
            words = c.get("contrast_words", [f"{ka}-ward", f"{kb}-ward"])
            za = (_sim_z(out, who, f"{ka}|{kb}", h) or {}).get("z")
            summ = D.text(f"ancestry_ho.tree_summary.{rl}", "")
            notes.append(f"{c.get('tree', f'The {h} tree')} is {'rejected' if f < 0.05 else 'not rejected'} (fit "
                         f"{fitv(who, h)}"
                         + ("" if za is None else
                            f"; under it {_pron(who)[1]} {ka}-vs-{kb} value would be {abs(za):.1f} SD too "
                            + (words[1] if za < 0 else words[0]))  # positive LLR = closer to the first population
                         + f"). The documented {summ + ' ' if summ else ''}tree "
                         f"{'fits' if g is not None and g >= 0.05 else 'does not fit'} (fit {fitv(who, 'gedcom')}; "
                         f"same contrast z {z.get(f'{ka}|{kb}', 0):+.1f}).")
        # optional per-role interpretation from family.json text "ancestry_ho.notes.<role>" (computed fields filled in)
        for tpl in D.text(f"ancestry_ho.notes.{rl}", []):
            notes.append(fill(tpl, who, node, typ))
        ibs = _pair(out, "IBS|Spanish", who)
        # a clearly Italian person preferring the IBS sample is anomalous (for Spanish ancestry it is not)
        if ibs and ibs["pa"] >= 0.9 and flags.get("italy") in ("italian", "between"):
            sp = _pron(who)[0]
            notes.append(f"Open item: {sp} {'match' if sp == 'they' else 'matches'} the 1000 Genomes Spanish sample "
                         f"(IBS) better than the main Spanish panel more strongly than {_pc(ibs['pa'])} of IBS people "
                         f"do, although {sp} {'are' if sp == 'they' else 'is'} clearly Italian. Unexplained; possibly "
                         "a chip/platform effect. It does not change the Italian placement.")
            ibs_open.append(who)
        disp["people"][who] = {"label": _label(who), "headline": headline(who, flags),
                               "statements": [{"text": t, "evidence": e} for t, e in S],
                               "notes": notes, "family_tree_check": [t[0] for t in tree], "key_numbers": nums}
    gp = out.get("family_logic", {}).get("grandpa_transmitted_haplotype")
    if gp:
        ts = gp["typicality_vs_pseudohaploid_loo"]["A_core"]
        fr = typ_of(ts)
        S, nums, tree, flags = _person_display(out, "grandpa", fr)
        italy_flags["grandpa"] = flags.get("italy")
        hit = [k for k in fr if fr[k] >= TYPICAL]
        miss = [k for k in fr if k not in hit]
        gvs = gp.get("vs_simulation", {}).get("A_core", {}).get("gedcom", {})
        gzs = {k: d["z"] for k, d in gp.get("pairwise_vs_simulation", {}).get("gedcom", {}).items()}
        gfit = gvs.get("frac_members_further")
        gmatch = gfit is not None and gfit >= 0.05 and (not gzs or max(abs(x) for x in gzs.values()) < 1.5)
        stm = [{"text": t, "evidence": e} for t, e in S]
        sn, sc = _pair(out, "Italian_South|Italian_North", "grandpa"), _pair(out, "Sicilian|Italian_North", "grandpa")
        where = None
        if sn and sc and _placed(fr, _tree_mix(out, "grandpa")):
            where = ("with Southern Italians" if sn["pa"] >= 0.1 and sn["pb"] >= 0.9 else
                     "with Northern Italians" if sn["pa"] < 0.1 and sn["pb"] < 0.9 else
                     "between Southern and Northern Italians")
            stm.append(
                {"text": f"On direct comparisons he sits {where}"
                         + (", not as a typical Sicilian" if sc["pa"] <= 0.1 else "") + ".",
                 "evidence": (lambda a, c: f"Southern-vs-Northern Italian {_num(a)}: {_pc(a['pa'])} of Southern "
                              f"Italians and {_pc(a['pb'])} of Northern Italians lower; Sicilian-vs-Northern "
                              f"{_num(c)}: " + (f"below all {c['na']} Sicilians" if c['pa'] == 0 else
                                                f"{_pc(c['pa'])} of the {c['na']} Sicilians lower")
                              + f", {_pc(c['pb'])} of Northern Italians lower "
                              "(reference people reduced to one DNA copy, like him).")(sn, sc)})
        if fr:
            groups = [TYP_POPS[k][0] for k in hit]
            others = [TYP_POPS[k][0] for k in miss]
            stm.append(
                {"text": (f"Profile typical of {_join(groups)}" + (f", not of {_join(others, 'or')}." if others else ".")
                          if groups else f"Profile not typical of {_join(others, 'or')}."),
                 "evidence": "Share of real people (reduced to one DNA copy, like him) who are more unusual: "
                             + ", ".join(f"{TYP_POPS[k][0]} {v:.2f}" for k, v in fr.items()) + "."})
        if gfit is not None:
            stm.append(
                {"text": ("Matches simulations of his family-tree ancestry." if gmatch else
                          "Broadly consistent with simulations of his family-tree ancestry, with some contrasts "
                          "off."),
                 "evidence": f"Fit to the GEDCOM tree {gfit:.2f}; " + (
                     (lambda zs: "no pairwise contrasts simulated." if not zs else
                      f"all pairwise contrasts within {max(abs(x) for x in zs.values()):.1f} SD of the tree's "
                      "expectation." if max(abs(x) for x in zs.values()) < 1.5 else
                      "pairwise contrasts outside 1.5 SD of the tree's expectation: "
                      + ", ".join(f"{k} z {x:+.1f}" for k, x in zs.items() if abs(x) >= 1.5) + ".")(gzs))})
        head = (f"{TYP_POPS[hit[0]][1].replace(' ', '-')}-type profile" if len(hit) == 1 else
                f"Profile typical of {_join([TYP_POPS[k][0] for k in hit])}" if hit else
                "Profile not typical of any single reference group")
        if where == "between Southern and Northern Italians" and "Italian_South" in hit:
            head += " mixed with some Central/Northern Italian"
        if gfit is not None:
            gtree = D.text("ancestry_ho.grandpa_tree", "")
            head += (("; consistent with" if gfit >= 0.05 else "; not a good match to") + " his family-tree ancestry"
                     + (f" ({gtree})" if gtree else ""))
        disp["people"]["grandpa-maternal (inferred)"] = {
            "label": "Grandpa (mom's father): the half of mom's DNA that came from him, reconstructed from "
                     "mom, grandma and mom's children",
            "headline": head + ".",
            "statements": stm,
            "notes": [("Which Southern region (Sicily vs Campania/Calabria) cannot be told apart. "
                       if "Italian_South" in hit or where in ("with Southern Italians",
                                                              "between Southern and Northern Italians") else "")
                      + "The haplotype is a "
                      f"reconstruction ({gp['n_snps']:,} of {out['meta']['n_common_snps']:,} SNPs usable) and covers "
                      "only the half of grandpa's DNA that mom inherited."],
            "family_tree_check": [t[0] for t in tree], "key_numbers": nums}
    v = out.get("validation", {})
    fr = v.get("family_relationships", {}).get("pairs", {})

    def opp(a, b):
        return (fr.get(f"{a}|{b}") or fr.get(f"{b}|{a}") or {}).get("opposite_homozygotes")

    def opp_rate(a, b):
        r = fr.get(f"{a}|{b}") or fr.get(f"{b}|{a}") or {}
        return r.get("rate", r["opposite_homozygotes"] / r["snps"] if r.get("snps") else None)
    pc = [(a, b) for a, b in ((SELF, MOM), (SELF, DAD), (MOM, GMA)) if opp(a, b) is not None]
    ok = [(a, b) for a, b in pc if opp_rate(a, b) is not None and opp_rate(a, b) <= PC_OPP_MAX]
    bad = [(a, b) for a, b in pc if (a, b) not in ok]
    rel = ""
    if ok:
        rel = ("Relationships confirmed from the DNA: "
               + ", ".join(f"{_short(a)}-{_short(b)} {opp(a, b)}" for a, b in ok)
               + " conflicting genotypes (parent-child: error only)")
    if opp(MOM, DAD) is not None:
        r = opp_rate(MOM, DAD)
        rel += ((f" vs {opp(MOM, DAD)} for " if ok else f"{opp(MOM, DAD)} conflicting genotypes for ")
                + f"{_short(MOM)}-{_short(DAD)} " + ("(unrelated)" if r is not None and r >= UNREL_OPP_MIN else
                   f"(CHECK: rate {100 * (r or 0):.1f}%, lower than expected for unrelated people)"))
    rel += ". " if rel else ""
    if bad:
        rel += ("CHECK: " + ", ".join(f"{_short(a)}-{_short(b)} {opp(a, b)} conflicting genotypes "
                                      + ("(rate n/a)" if opp_rate(a, b) is None else f"({100 * opp_rate(a, b):.2f}%)")
                                      for a, b in bad)
                + f", more than genotyping error for parent and child (<= {100 * PC_OPP_MAX:.1f}%). ")
    disp["family_checks"] = rel + f"Strand check: {v.get('strand_calibration', {}).get('verdict', 'n/a')}."
    # Italian wording (South Asian explanation, Sicily open item) only when someone places as Italian or the
    # documented tree has Italian founders
    italian = (any(f in ("italian", "between") for f in italy_flags.values())
               or any((_tree_mix(out, w) or (0, 0))[1] > 0 for w in list(out.get("people", {})) + ["grandpa"]))
    sa = out.get("south_asian_test", {}).get("family", {})
    if sa:
        try:
            with open(os.path.join(RESULTS, "ancestry.json")) as f:
                old = json.load(f).get("people", {})
        except (OSError, ValueError):
            old = {}
        sh = sorted(100 * old[p]["broad"]["South Asian"] for p in sa if "South Asian" in old.get(p, {}).get("broad", {}))
        size = D.text("ancestry_ho.south_asian_range", "") or (
                "small" if not sh else f"{sh[0]:.0f}%" if f"{sh[0]:.0f}" == f"{sh[-1]:.0f}"
                else f"{sh[0]:.0f}-{sh[-1]:.0f}%")
        emed = ("the ordinary Eastern-Mediterranean share of Italians" if italian else
                "an ordinary Eastern-Mediterranean-like share")
        e1 = np.mean([x["Iranian_E1"] for x in sa.values()])
        e2 = np.mean([x["Iranian_E2"] for x in sa.values()])
        disp["south_asian"] = (
            f"The {size} 'South Asian' in the old 1000 Genomes analysis "
            + (f"was most likely {emed}, not South Asian ancestry. " if e2 < e1 else
               f"is only partly explained by {emed}. ")
            + "With an Iranian group standing in for "
            "the old Gujarati reference and only Western European alternatives, an 'Iranian' share appears ("
            + ", ".join(f"{k} {100 * x['Iranian_E1']:.0f}±{100 * x['Iranian_E1_se']:.0f}%" for k, x in sa.items())
            + "); once Southern Italian, Levant and Balkan sources are offered it "
            + ("drops to " if e2 < e1 else "is ")
            + ", ".join(f"{k} {100 * x['Iranian_E2']:.0f}±{100 * x['Iranian_E2_se']:.0f}%" for k, x in sa.items())
            + ". This is an analogy, not an exact re-test: the old per-person values are not reproduced.")
    # Earlier statements of this family's reports that were withdrawn or corrected (family.json text).
    disp["withdrawn_or_corrected_claims"] = list(D.text("ancestry_ho.withdrawn_claims", []))
    disp["open_items"] = ([("My" if p == SELF else D.label(p).split(' (')[0] + "'s")
                           + " unexplained preference for the IBS sample over the main Spanish panel."
                           for p in ibs_open]
                          + list(D.text("ancestry_ho.open_items", []))
                          + (["Sicily vs Campania/Calabria cannot be resolved with ~11 Sicilians and ~10 Campanian/"
                              "Calabrian reference people."] if italian else []))
    return disp


# ---------------------------------------------------------------- main
def main():
    t0 = time.time()
    snps, samples, pops, G, fam, chrom, pos, excluded = load()
    people = [p for p in D.people() if p in fam]
    log(f"loaded {G.shape[0]} reference samples x {G.shape[1]} SNPs; people: {people}")
    if SELF not in people:
        log("note: no genotypes for self in the cache; run ancestry_ho_fetch.py first. Nothing written.")
        return
    if not HYPOTHESES:
        log("note: no family.json ancestry founders (16..31 with ref_pop); pedigree-hypothesis tests skipped")
    # ---- 0. QC: relatives in the reference, strand calibration, family relationships
    rel_pairs, auto_drop = kinship_scan(G, pops, samples, set(excluded))
    if auto_drop:
        log("WARNING: new related reference pairs, dropping", [samples[j]["iid"] for j in auto_drop])
        pops = {p: [j for j in idx if j not in set(auto_drop)] for p, idx in pops.items()}
    pops = add_subpops(pops, samples)
    validation = {
        "reference_relatives": {
            "method": f"KING-robust kinship on all HO SNPs within each region; pairs > {KIN_DROP} listed",
            "pairs": rel_pairs,
            "excluded": {samples[j]["iid"]: F.EXCLUDE_RELATIVES.get(samples[j]["iid"], "auto: kinship scan")
                         for j in list(excluded) + auto_drop}},
        "strand_calibration": strand_calibration(G, pops, fam, snps),
        "family_relationships": relationship_checks({p: fam[p] for p in people}),
        "subpopulations": {n: {"of": SUBPOPS[n][0], "n": len(pops[n]),
                               "localities": sorted({samples[j]["locality"] for j in pops[n]})} for n in SUBPOPS},
        "italian_south_composition": dict(collections.Counter(
            "Campania" if any(k in samples[j]["locality"] for k in ("Campania", "Naples", "Salerno")) else
            "Calabria" if any(k in samples[j]["locality"] for k in ("Calabria", "Cosenza")) else
            "Salento_Griko (Lecce)" if "Lecce" in samples[j]["locality"] else samples[j]["locality"]
            for j in pops["Italian_South"]))}
    log("QC done:", len(rel_pairs), "related reference pairs;", validation["strand_calibration"]["verdict"])
    for cfg in CONFIGS.values():
        for _n, mem in cfg["sources"]:
            for p in mem:
                assert p in pops, p
    for H in HYPOTHESES.values():
        for p in H["founders"].values():
            assert p in pops, f"family.json ancestry founder ref_pop {p!r} is not a reference population"
    cols = np.where(np.all([fam[p] != 3 for p in people], axis=0))[0]
    log(f"common SNP set: {len(cols)}")
    panel = Panel(G, pops, cols)
    ch_c, pos_c = chrom[cols], pos[cols]
    blocks_c = ch_c * 1000 + pos_c // WIN_BP
    famc = {p: fam[p][cols] for p in people}
    eqn = EqualN(panel.Gc, pops, panel.pbar)
    pw = Pairwise(panel.Gc, pops, panel.pbar, PAIRS)
    gp = load_grandpa(snps)
    gpc = None
    if gp is not None and MOM in people and GMA in people:
        gpc = gp[cols]
        if int((gpc != 3).sum()) < 10000:
            gpc = None
    gp_mask = None if gpc is None else gpc != 3

    out = {"meta": {}, "configs": {}, "people": {}, "calibration": {}, "pedigree_simulation": {},
           "family_logic": {}, "south_asian_test": {}}

    # ---- 1. supervised admixture
    adm = collections.defaultdict(dict)
    Pfull = {cfg: panel.P(cfg) for cfg in CONFIGS}
    for cfg in CONFIGS:
        out["configs"][cfg] = {"about": CONFIGS[cfg]["about"],
                               "sources": {n: {"pops": m, "n": sum(len(pops[p]) for p in m)}
                                           for n, m in CONFIGS[cfg]["sources"]}}
        for p in people:
            adm[p][cfg] = fit(panel, famc[p], cfg, jack_chrom=ch_c, P=Pfull[cfg])
    log("family admixture done")
    me_full = None
    if SELF in fam:
        colsm = np.where(fam[SELF] != 3)[0]
        panel_m = Panel(G, pops, colsm)
        me_full = {cfg: fit(panel_m, fam[SELF][colsm], cfg, jack_chrom=chrom[colsm]) for cfg in MAIN_CONFIGS}
        eq_m = EqualN(panel_m.Gc, pops, panel_m.pbar, cache=False)
        g = fam[SELF][colsm]
        near_me_full = ll_margins(eq_m.per_snp(g), chrom[colsm] * 1000 + pos[colsm] // WIN_BP, g != 3, top=10)
        del panel_m, eq_m
        log("self full-SNP analysis done")

    # ---- 2. nearest population (equal n), plus the naive full-n ranking for comparison
    near, naive_rank = {}, {}
    for p in people:
        g = famc[p]
        called = g != 3
        near[p] = ll_margins(eqn.per_snp(g), blocks_c, called)
        gi = np.where(called, g, 0).astype(int)
        ar = np.arange(len(g))
        tot = {q: float(ll_table(panel.x[q], panel.N[q], panel.pbar)[ar, gi][called].sum()) for q in pops}
        naive_rank[p] = sorted(tot, key=tot.get, reverse=True)[:5]
    log("nearest population done")

    # ---- 3. LOO calibration
    loo = []
    for lp in LOO_POPS:
        for j in pops[lp]:
            gd = panel.Gc[j]
            row = {"pop": lp, "iid": samples[j]["iid"], "locality": samples[j]["locality"],
                   "n_snps": int((gd != 3).sum())}
            for cfg in CONFIGS:
                P = Pfull[cfg].copy()      # only the target's own source changes
                for k, (_n, mem) in enumerate(CONFIGS[cfg]["sources"]):
                    if lp in mem:
                        P[:, k] = panel.freq(mem, (j,))
                row[cfg] = fit(panel, gd, cfg, P=P)["q"]
            if gpc is not None:        # pseudo-haploid copy on grandpa's SNPs, for the grandpa comparison
                hh = pseudo_haploid(gd, np.random.default_rng(SEED + j), gp_mask)
                row["hap"] = hh
                for cfg in ("A_core", "D_sicily_split"):
                    P = Pfull[cfg].copy()
                    for k, (_n, mem) in enumerate(CONFIGS[cfg]["sources"]):
                        if lp in mem:
                            P[:, k] = panel.freq(mem, (j,))
                    g2 = np.where(hh == 3, 3, 2 * hh).astype(np.uint8)
                    row["hap_" + cfg] = fit(panel, g2, cfg, P=P)["q"]
            tot = {q: float(v.sum()) for q, v in eqn.per_snp(gd, exclude=(j,)).items()}
            row["best_pop"] = max(tot, key=tot.get)
            row["ll"] = tot
            loo.append(row)
        log(f"LOO {lp} ({len(pops[lp])}) done")

    calib = {}
    loo_cloud = {cfg: {} for cfg in CONFIGS}
    for lp in LOO_POPS:
        rows = [r for r in loo if r["pop"] == lp]
        n = len(rows)
        bp = collections.Counter(r["best_pop"] for r in rows)
        reg = POP_REGION.get(lp)
        c = {"n": n, "nearest_population_equal_n": {
            "assigned_to_own_pop": rnd(bp[lp] / n, 3),
            "assigned_to_own_region": rnd(sum(v for k, v in bp.items() if POP_REGION.get(k) == reg) / n, 3),
            "best_pop_counts": dict(bp.most_common()),
            "note": None if lp in eqn.names else f"n={len(pops[lp])} < {NEQ}: not a candidate population"},
            "admixture": {}}
        for cfg in CONFIGS:
            names = src_names(cfg)
            Q = np.array([r[cfg] for r in rows])
            loo_cloud[cfg][lp] = Q
            top = collections.Counter(names[int(np.argmax(q))] for q in Q)
            d = {"mean": comp_dict(names, Q.mean(0)),
                 "sd": comp_dict(names, Q.std(0, ddof=1) if n > 1 else np.zeros(len(names))),
                 "p10": comp_dict(names, np.percentile(Q, 10, axis=0)),
                 "p90": comp_dict(names, np.percentile(Q, 90, axis=0)),
                 "largest_component_counts": dict(top.most_common())}
            own = OWN_SOURCE.get(cfg, {}).get(lp)
            if own:
                d["largest_is_own_source"] = rnd(top[own] / n, 3)
                d["own_source_mean"] = rnd(Q[:, names.index(own)].mean())
            if cfg in MAIN_CONFIGS:
                d["broad_mean"] = {k: rnd(np.mean([broad(cfg, q)[k] for q in Q])) for k in broad(cfg, Q[0])}
            c["admixture"][cfg] = d
        calib[lp] = c
    # profile classifier: each LOO individual -> nearest LOO-population cloud (Mahalanobis),
    # leaving the individual out of its own cloud
    big = [lp for lp in LOO_POPS if len(pops[lp]) >= MIN_CLOUD]
    for cfg in ("A_core", "D_sicily_split", "B_extended"):
        ix = near_idx(cfg)
        for lp in big:
            Qo = loo_cloud[cfg][lp][:, ix]
            own_hit = reg_hit = 0
            for i in range(len(Qo)):
                best, bd = None, 1e9
                for lp2 in big:
                    cl = loo_cloud[cfg][lp2][:, ix]
                    if lp2 == lp:
                        cl = np.delete(cl, i, axis=0)
                    d = mahal(Qo[i], cl, loo_dist=False, min_n=MIN_CLOUD - 1)["distance"]
                    if d < bd:
                        best, bd = lp2, d
                own_hit += best == lp
                reg_hit += POP_REGION.get(best) == POP_REGION.get(lp)
            calib[lp]["admixture"][cfg]["profile_classifier_own_pop"] = rnd(own_hit / len(Qo), 3)
            calib[lp]["admixture"][cfg]["profile_classifier_own_region"] = rnd(reg_hit / len(Qo), 3)
    out["calibration"] = {
        "method": "leave-one-out: each individual removed from its own population/source counts (and from "
                  "the equal-n subsamples), analysed on the family's common SNP set",
        "caution": "nearest_population_equal_n.assigned_to_own_pop is a 65-way ranking with attractor "
                   "populations (IBS, Bulgarian); it understates separability. Use 'pairwise'.",
        "populations": calib}
    log("calibration summaries done")

    # ---- 2b. pairwise contrasts: LOO calibration, family, grandpa
    pair_out = {}
    loo_pair = {}
    for a, b, m in pw.pairs:
        key = f"{a}|{b}"
        vals, hvals = {a: [], b: []}, {a: [], b: []}
        for grp in (a, b):
            for j in pops[grp]:
                gd = panel.Gc[j]
                vals[grp].append(pw.value(gd, a, b, m, exclude=(j,)))
                if gpc is not None:
                    hh = pseudo_haploid(gd, np.random.default_rng(SEED + 3 * j), gp_mask)
                    hvals[grp].append(pw.value(hh, a, b, m, exclude=(j,), hap=True))
        xa, xb = np.array(vals[a]), np.array(vals[b])
        loo_pair[key] = (xa, xb)
        d = {"n_eq": m, "n_a": len(xa), "n_b": len(xb),
             "loo_accuracy": rnd(((xa > 0).sum() + (xb < 0).sum()) / (len(xa) + len(xb)), 3),
             "loo_accuracy_balanced": rnd(((xa > 0).mean() + (xb < 0).mean()) / 2, 3),
             "a_median": round(float(np.median(xa)), 1), "a_range": [round(float(xa.min()), 1), round(float(xa.max()), 1)],
             "b_median": round(float(np.median(xb)), 1), "b_range": [round(float(xb.min()), 1), round(float(xb.max()), 1)],
             "family": {}}
        for p in people:
            g = famc[p]
            per = pw.per_snp(g, a, b, m)
            called = g != 3
            v = float(per.sum()) / called.sum() * 1e4
            d["family"][p] = {"llr_per10k": round(v, 1), "se": round(block_se(per, blocks_c, called), 1),
                              "pct_among_a": pct(xa, v), "pct_among_b": pct(xb, v)}
        if gpc is not None:
            ha, hb = np.array(hvals[a]), np.array(hvals[b])
            per = pw.per_snp(gpc, a, b, m, hap=True)
            v = float(per.sum()) / gp_mask.sum() * 1e4
            d["grandpa_haploid"] = {
                "llr_per10k": round(v, 1), "se": round(block_se(per, blocks_c, gp_mask), 1),
                "pct_among_a_pseudohaploid": pct(ha, v), "pct_among_b_pseudohaploid": pct(hb, v),
                "a_median": round(float(np.median(ha)), 1), "b_median": round(float(np.median(hb)), 1),
                "loo_accuracy_pseudohaploid": rnd(((ha > 0).sum() + (hb < 0).sum()) / (len(ha) + len(hb)), 3)}
        pair_out[key] = d
        log(f"pair {key}: acc {d['loo_accuracy']}")
    out["pairwise"] = {
        "method": "per-SNP log P(g|A) - log P(g|B), beta-binomial predictive, A and B subsampled to equal n "
                  f"(n_eq), {R_PAIR} replicates, per 10k SNPs; positive = closer to A. LOO members never in "
                  "their own subsample; family SE = 5-Mb block jackknife; pct = fraction of LOO members of "
                  "A (or B) with a value <= the person's. Grandpa: haploid likelihood vs pseudo-haploid LOO.",
        "pairs": pair_out}

    # ---- 3b. pedigree simulations
    rng = np.random.default_rng(SEED)
    mei = Meiosis(ch_c, pos_c, rng)
    sim_cfgs = MAIN_CONFIGS + ["E2_sas_eastmed"]
    sim_who = list(SIM_NODES.values()) + (["grandpa_hap"] if gpc is not None else [])
    sims = {h: {who: {cfg: [] for cfg in sim_cfgs} for who in sim_who} for h in HYPOTHESES}
    sim_llr = {h: {who: {f"{a}|{b}": [] for a, b in SIM_PAIRS} for who in sim_who} for h in HYPOTHESES}
    pm = {(a, b): m for a, b, m in pw.pairs}
    for h, H in HYPOTHESES.items():
        for _ in range(N_SIM):
            geno, used = simulate_family(panel, pops, H["founders"], mei, rng)
            if gpc is not None:
                geno["gp2"] = np.where(gp_mask & (geno["gp_hap"] != 3), 2 * geno["gp_hap"], 3).astype(np.uint8)
                geno["gp1"] = np.where(gp_mask, geno["gp_hap"], 3).astype(np.uint8)
            for cfg in sim_cfgs:
                P = panel.P(cfg, tuple(used))
                for n_, who in SIM_NODES.items():
                    sims[h][who][cfg].append(fit(panel, geno[n_], cfg, P=P)["q"])
                if gpc is not None:
                    sims[h]["grandpa_hap"][cfg].append(fit(panel, geno["gp2"], cfg, P=P)["q"])
            memo = {}
            for a, b in SIM_PAIRS:
                m = pm[(a, b)]
                for n_, who in SIM_NODES.items():
                    sim_llr[h][who][f"{a}|{b}"].append(
                        pw.value(geno[n_], a, b, m, exclude=tuple(used), reps=R_PAIR_SIM, memo=memo))
                if gpc is not None:
                    sim_llr[h]["grandpa_hap"][f"{a}|{b}"].append(
                        pw.value(geno["gp1"], a, b, m, exclude=tuple(used), hap=True, reps=R_PAIR_SIM, memo=memo))
        log(f"simulations {h} done")
    simsum = {}
    for h in HYPOTHESES:
        simsum[h] = {"about": HYPOTHESES[h]["about"], "founders": HYPOTHESES[h]["founders"], "people": {}}
        for who in sim_who:
            simsum[h]["people"][who] = {"pairwise": {k: {"mean": round(float(np.mean(v)), 1),
                                                         "sd": round(float(np.std(v, ddof=1)), 1)}
                                                     for k, v in sim_llr[h][who].items()}}
            for cfg in sim_cfgs:
                Q = np.array(sims[h][who][cfg])
                names = src_names(cfg)
                d = {"mean": comp_dict(names, Q.mean(0)), "sd": comp_dict(names, Q.std(0, ddof=1))}
                if cfg in MAIN_CONFIGS:
                    d["broad_mean"] = {k: rnd(np.mean([broad(cfg, q)[k] for q in Q])) for k in broad(cfg, Q[0])}
                simsum[h]["people"][who][cfg] = d
    out["pedigree_simulation"] = {"method": f"see module docstring (3b); {N_SIM} simulated families per hypothesis; "
                                            f"pairwise contrasts with {R_PAIR_SIM} replicates, founders left out",
                                  "proxy_notes": PROXY_NOTES, "hypotheses": simsum}

    # ---- assemble per person
    for p in people:
        pr = {"n_snps_common": int((famc[p] != 3).sum()), "admixture": {}, "broad_by_config": {}}
        for cfg in CONFIGS:
            r = adm[p][cfg]
            names = src_names(cfg)
            pr["admixture"][cfg] = {"components": comp_dict(names, r["q"]), "se": comp_dict(names, r["se"])}
        for cfg in MAIN_CONFIGS:
            pr["broad_by_config"][cfg] = broad(cfg, adm[p][cfg]["q"])
            it_ix = [i for i, n in enumerate(src_names(cfg)) if n in ITALIAN_SOURCES]
            r = adm[p][cfg]
            pr["broad_by_config"][cfg]["Italian_total_se"] = rnd(
                jk_se(r["jk"][:, it_ix].sum(1, keepdims=True), r["jk_w"], np.array([r["q"][it_ix].sum()]))[0])
        it = [pr["broad_by_config"][c]["Italian_total"] for c in MAIN_CONFIGS]
        ib = [pr["broad_by_config"][c]["Iberian_total"] for c in MAIN_CONFIGS]
        pr["stability"] = {"Italian_total_range": [min(it), max(it)], "Iberian_total_range": [min(ib), max(ib)]}
        best, marg, tot = near[p]
        pr["nearest_population_equal_n"] = {
            "caution": "65-way ranking with attractor populations (e.g. 90/173 LOO Spaniards rank IBS first); "
                       "a diagnostic only, not the person's population. See pairwise.",
            "top": dict(list(marg.items())[:12]),
            "key_pops": {k: marg[k] for k in LOO_POPS if k in marg},
            "naive_full_n_top5": naive_rank[p]}
        typ = {}
        for cfg in ("A_core", "D_sicily_split", "C_coarse"):
            ix = near_idx(cfg)
            v = adm[p][cfg]["q"][ix]
            typ[cfg] = {lp: mahal(v, loo_cloud[cfg][lp][:, ix]) for lp in LOO_POPS}
        pr["typicality_vs_loo_populations"] = typ
        nsn = pr["n_snps_common"]
        cand = [q for q in LOO_POPS + ["IBS"] if q in eqn.names]
        lt = {}
        for lp in cand:
            if lp not in LOO_POPS:
                continue
            others = [q for q in cand if q != lp]
            mine = (tot[lp] - max(tot[q] for q in others)) / nsn * 1e4
            theirs = np.array([(r["ll"][lp] - max(r["ll"][q] for q in others)) / r["n_snps"] * 1e4
                               for r in loo if r["pop"] == lp])
            lt[lp] = {"margin_per10k": round(mine, 2), "loo_median": round(float(np.median(theirs)), 2),
                      "loo_p10": round(float(np.percentile(theirs, 10)), 2),
                      "frac_loo_members_lower": rnd((theirs <= mine).mean(), 3)}
        pr["typicality_nearest_pop_margin"] = lt
        pr["pairwise"] = {k: dict(v["family"][p], loo_accuracy=v["loo_accuracy"]) for k, v in pair_out.items()}
        out["people"][p] = pr
    if me_full is not None:
        out["people"][SELF]["full_snp_set"] = {
            "n_snps": int(len(colsm)),
            "admixture": {cfg: {"components": comp_dict(src_names(cfg), r["q"]),
                                "se": comp_dict(src_names(cfg), r["se"])} for cfg, r in me_full.items()},
            "nearest_population_equal_n": {"best": near_me_full[0], "top": near_me_full[1]}}

    # ---- 4. family logic
    fl = {}
    if all(x in people for x in (SELF, MOM, DAD)):
        mid = {}
        S = sims.get("gedcom")
        for cfg in MAIN_CONFIGS:
            names = src_names(cfg)
            qme, qm, qd = (adm[x][cfg]["q"] for x in (SELF, MOM, DAD))
            diff = qme - (qm + qd) / 2
            jd = adm[SELF][cfg]["jk"] - (adm[MOM][cfg]["jk"] + adm[DAD][cfg]["jk"]) / 2
            se = jk_se(jd, adm[SELF][cfg]["jk_w"], diff)
            mid[cfg] = {"me": comp_dict(names, qme), "midparent": comp_dict(names, (qm + qd) / 2),
                        "me_minus_midparent": comp_dict(names, diff), "jackknife_se": comp_dict(names, se)}
            if S:
                simd = np.array(S[SELF][cfg]) - (np.array(S[MOM][cfg]) + np.array(S[DAD][cfg])) / 2
                mid[cfg]["simulated_sd_of_difference"] = comp_dict(names, simd.std(0, ddof=1))
            mid[cfg]["max_abs_diff"] = rnd(np.abs(diff).max())
            if GMA in people:
                qg = adm[GMA][cfg]["q"]
                mid[cfg]["grandpa_maternal_inferred_2mom_minus_grandma"] = comp_dict(names, 2 * qm - qg)
                if S:
                    sg = 2 * np.array(S[MOM][cfg]) - np.array(S[GMA][cfg])
                    mid[cfg]["grandpa_inferred_simulated_gedcom_mean"] = comp_dict(names, sg.mean(0))
                    mid[cfg]["grandpa_inferred_simulated_gedcom_sd"] = comp_dict(names, sg.std(0, ddof=1))
        fl["me_vs_midparent"] = mid
    if HYPOTHESES:
        ped = {}
        for who in SIM_NODES.values():
            if who not in people:
                continue
            ped[who] = {}
            for cfg in ("A_core", "D_sicily_split", "C_coarse", "B_extended"):
                ix = near_idx(cfg)
                names = src_names(cfg)
                v = adm[who][cfg]["q"]
                ped[who][cfg] = {}
                for h in HYPOTHESES:
                    Q = np.array(sims[h][who][cfg])
                    z = (v - Q.mean(0)) / np.maximum(Q.std(0, ddof=1), 0.01)
                    if cfg == "A_core":
                        pz = {}
                        for k, arr in sim_llr[h][who].items():
                            o = pair_out[k]["family"][who]
                            sd = math.sqrt(float(np.var(arr, ddof=1)))
                            pz[k] = {"observed": o["llr_per10k"], "sim_mean": round(float(np.mean(arr)), 1),
                                     "sim_sd": round(sd, 1), "z": round((o["llr_per10k"] - float(np.mean(arr))) / sd, 2)}
                        ped[who].setdefault("pairwise_vs_simulation", {})[h] = pz
                    ped[who][cfg][h] = {
                        "simulated_mean": comp_dict(names, Q.mean(0)),
                        "observed_minus_sim_mean": comp_dict(names, v - Q.mean(0)),
                        "z": {n: round(float(x), 2) for n, x in zip(names, z) if n not in FAR},
                        "mahalanobis": mahal(v[ix], Q[:, ix])}
        fl["observed_vs_pedigree_simulations"] = ped
        fits = {}
        for who in ped:
            fits[who] = {}
            for cfg in ("A_core", "C_coarse", "D_sicily_split"):
                r_ = sorted(((h, ped[who][cfg][h]["mahalanobis"]) for h in HYPOTHESES
                             if "frac_members_further" in ped[who][cfg][h]["mahalanobis"]),
                            key=lambda t: -t[1]["frac_members_further"])
                if not r_:
                    continue
                top = r_[0][1]["frac_members_further"]
                fits[who][cfg] = [{"hypothesis": h, "fit": d["frac_members_further"], "se": d["frac_se"],
                                   "not_distinguishable_from_best": top - d["frac_members_further"] < NO_RANK}
                                  for h, d in r_]
        fl["hypothesis_fits"] = {"note": f"fit = fraction of simulated relatives at least as far from the simulated "
                                         f"mean (Mahalanobis); fits within {NO_RANK} of the best are not ranked",
                                 "by_person": fits}
    if gpc is not None:
        gpres = {"n_snps": int(gp_mask.sum()), "source": "work/grandpa/grandpa_haplotype.tsv.gz "
                 "(transmitted allele known, no ambiguous het/het), analysed as a haploid genome",
                 "admixture": {}, "typicality_vs_pseudohaploid_loo": {}, "vs_simulation": {}}
        g2 = np.where(gp_mask, 2 * gpc, 3).astype(np.uint8)
        for cfg in MAIN_CONFIGS:
            r = fit(panel, g2, cfg, jack_chrom=ch_c, P=Pfull[cfg])
            names = src_names(cfg)
            gpres["admixture"][cfg] = {"components": comp_dict(names, r["q"]), "se": comp_dict(names, r["se"])}
            if GMA in people and cfg == "A_core":
                dlt = r["q"] - adm[GMA][cfg]["q"]
                jd = r["jk"] - adm[GMA][cfg]["jk"]
                gpres["grandpa_minus_grandma_A_core"] = {"diff": comp_dict(names, dlt),
                                                         "se": comp_dict(names, jk_se(jd, r["jk_w"], dlt))}
            if cfg in ("A_core", "D_sicily_split"):
                ix = near_idx(cfg)
                gpres["typicality_vs_pseudohaploid_loo"][cfg] = {
                    lp: mahal(r["q"][ix], np.array([x["hap_" + cfg] for x in loo if x["pop"] == lp])[:, ix])
                    for lp in LOO_POPS}
            if cfg in sim_cfgs and HYPOTHESES:
                ix = near_idx(cfg)
                gpres["vs_simulation"][cfg] = {h: mahal(r["q"][ix], np.array(sims[h]["grandpa_hap"][cfg])[:, ix])
                                               for h in HYPOTHESES}
        gpres["pairwise_vs_simulation"] = {
            h: {k: {"observed": pair_out[k]["grandpa_haploid"]["llr_per10k"],
                    "sim_mean": round(float(np.mean(v)), 1), "sim_sd": round(float(np.std(v, ddof=1)), 1),
                    "z": round((pair_out[k]["grandpa_haploid"]["llr_per10k"] - float(np.mean(v)))
                               / float(np.std(v, ddof=1)), 2)}
                for k, v in sim_llr[h]["grandpa_hap"].items()} for h in HYPOTHESES}
        fl["grandpa_transmitted_haplotype"] = gpres
    out["family_logic"] = fl

    # ---- 5. South Asian test
    sas = {"method": CONFIGS["E1_sas_like"]["about"] + " | " + CONFIGS["E2_sas_eastmed"]["about"],
           "family": {}, "loo_reference": {}, "simulated_gedcom": {}}
    n1, n2 = src_names("E1_sas_like"), src_names("E2_sas_eastmed")
    for p in people:
        e1, e2 = adm[p]["E1_sas_like"], adm[p]["E2_sas_eastmed"]
        c2 = dict(zip(n2, e2["q"]))
        sas["family"][p] = {
            "Iranian_E1": rnd(e1["q"][n1.index("Iranian")]), "Iranian_E1_se": rnd(e1["se"][n1.index("Iranian")]),
            "Iranian_E2": rnd(c2["Iranian"]), "Iranian_E2_se": rnd(e2["se"][n2.index("Iranian")]),
            "E2_eastern_total": rnd(c2["Balkan_Aegean"] + c2["Levant"] + c2["Anatolia_Caucasus"] + c2["Iranian"]),
            "E2_components": comp_dict(n2, e2["q"])}
    for lp in LOO_POPS:
        rows = [r for r in loo if r["pop"] == lp]
        a1 = np.array([r["E1_sas_like"][n1.index("Iranian")] for r in rows])
        a2 = np.array([r["E2_sas_eastmed"][n2.index("Iranian")] for r in rows])
        sas["loo_reference"][lp] = {"Iranian_E1_mean": rnd(a1.mean()),
                                    "Iranian_E1_p10_p90": [rnd(np.percentile(a1, 10)), rnd(np.percentile(a1, 90))],
                                    "Iranian_E2_mean": rnd(a2.mean())}
    for who in SIM_NODES.values() if HYPOTHESES else ():
        Q = np.array(sims["gedcom"][who]["E2_sas_eastmed"])
        sas["simulated_gedcom"][who] = {"Iranian_E2_mean": rnd(Q[:, n2.index("Iranian")].mean()),
                                        "Levant_E2_mean": rnd(Q[:, n2.index("Levant")].mean())}
    out["south_asian_test"] = sas

    out["meta"] = {
        "script": "scripts/ancestry_ho.py",
        "reference": "AADR v66.p1 HO (Mallick et al. 2024, doi:10.7910/DVN/FFIDCW); cache work/ancestry_ho/",
        "people": people, "n_common_snps": int(len(cols)), "lambda": LAMBDA, "predictive_prior": PRIOR,
        "equal_n": NEQ, "equal_n_replicates": R_NEQ, "equal_n_candidates": len(eqn.names),
        "jackknife_window_bp": WIN_BP, "n_sim": N_SIM,
        "loo_populations": {lp: len(pops[lp]) for lp in LOO_POPS},
        "subpopulations": {n: len(pops[n]) for n in SUBPOPS}, "pairwise_replicates": R_PAIR,
        "min_cloud": MIN_CLOUD, "excluded_reference": [samples[j]["iid"] for j in list(excluded) + auto_drop],
        "runtime_s": round(time.time() - t0, 1)}

    out["validation"] = validation
    out["display"] = build_display(out)
    os.makedirs(RESULTS, exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=1)
    mpath = os.path.join(F.OUT, "manifest.json")
    try:
        with open(mpath) as f:
            man = json.load(f)
        man["validation_downstream"] = {
            "by": "scripts/ancestry_ho.py", "reference_relatives": validation["reference_relatives"],
            "strand_calibration": validation["strand_calibration"],
            "family_relationships": validation["family_relationships"]}
        with open(mpath, "w") as f:
            json.dump(man, f, indent=1)
    except OSError as e:
        log("could not update manifest:", e)
    with open(OUT_TSV, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["section", "person_or_pop", "config", "component", "value", "se_or_sd"])
        for p, pr in out["people"].items():
            for cfg, d in pr["admixture"].items():
                for k, v in d["components"].items():
                    w.writerow(["admixture", p, cfg, k, v, d["se"][k]])
            for k, d in pr["nearest_population_equal_n"]["top"].items():
                w.writerow(["nearest_pop_ll_minus_best", p, "", k, d["ll_minus_best"], d["se"]])
        for lp, c in calib.items():
            for cfg, d in c["admixture"].items():
                for k, v in d["mean"].items():
                    w.writerow(["loo_mean", lp, cfg, k, v, d["sd"][k]])
            w.writerow(["loo_nearest_own_pop_rate", lp, "", "",
                        c["nearest_population_equal_n"]["assigned_to_own_pop"], ""])
        for h, H in simsum.items():
            for who, d in H["people"].items():
                for cfg, s in d.items():
                    if cfg == "pairwise":
                        continue
                    for k, v in s["mean"].items():
                        w.writerow([f"sim_{h}", who, cfg, k, v, s["sd"][k]])
    with open(os.path.join(F.OUT, "loo_rows.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        hdr = ["pop", "iid", "locality", "best_pop_equal_n"]
        for cfg in CONFIGS:
            hdr += [f"{cfg}:{n}" for n in src_names(cfg)]
        w.writerow(hdr)
        for r in loo:
            row = [r["pop"], r["iid"], r["locality"], r["best_pop"]]
            for cfg in CONFIGS:
                row += [f"{x:.4f}" for x in r[cfg]]
            w.writerow(row)
    log(f"wrote {OUT_JSON} and {OUT_TSV} in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    if "--display-only" in sys.argv:
        with open(OUT_JSON) as f:
            res = json.load(f)
        res["display"] = build_display(res)
        with open(OUT_JSON, "w") as f:
            json.dump(res, f, indent=1)
        print(json.dumps(res["display"], indent=1))
    else:
        main()
