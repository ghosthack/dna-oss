"""Loading and normalizing consumer DNA raw files (23andMe v5, AncestryDNA v2).

Genotypes are stored as sorted allele strings ("AG", "CC", "A" for hemizygous,
"DI" for indels). No-calls are dropped. Everything is GRCh37, plus strand.
"""
import glob
import json
import os
import shutil
import subprocess
import zipfile

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# The data folder (raw/ in, work/ and results/ out, plus family.json) lives outside the code: set DNA_DATA.
DATA = os.path.abspath(os.path.expanduser(os.environ.get("DNA_DATA", CODE)))
RAW = os.path.join(DATA, "raw")
WORK = os.path.join(DATA, "work")
RESULTS = os.path.join(DATA, "results")
# Third-party tools and reference files fetched by fetch_refs.sh (never committed).
TOOLS = os.path.abspath(os.path.expanduser(os.environ.get("DNA_TOOLS", os.path.join(CODE, "tools"))))
BCFTOOLS = os.environ.get("BCFTOOLS") or shutil.which("bcftools") or "/opt/homebrew/bin/bcftools"

CHROMS = [str(c) for c in range(1, 23)] + ["X", "Y", "MT"]
ANCESTRY_CHROM = {"23": "X", "24": "Y", "26": "MT"}  # 25 (PAR) is dropped

# Approximate GRCh37 chromosome lengths, used for drawing.
CHROM_LEN = {
    "1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276, "5": 180915260,
    "6": 171115067, "7": 159138663, "8": 146364022, "9": 141213431, "10": 135534747,
    "11": 135006516, "12": 133851895, "13": 115169878, "14": 107349540, "15": 102531392,
    "16": 90354753, "17": 81195210, "18": 78077248, "19": 59128983, "20": 63025520,
    "21": 48129895, "22": 51304566, "X": 155270560,
}


def _lines(path):
    if path.endswith(".zip"):
        z = zipfile.ZipFile(path)
        name = next(n for n in z.namelist() if n.endswith(".txt"))
        fh = z.open(name)
    else:
        fh = open(path, "rb")
    with fh:
        for line in fh:
            if line[:1] != b"#":
                yield line.decode().strip().split("\t")


def load_23andme(path):
    """Return {rsid: (chrom, pos, gt)}."""
    out = {}
    for f in _lines(path):
        if len(f) != 4 or f[0] == "rsid":
            continue
        rs, ch, pos, gt = f
        if "-" in gt or not gt:
            continue
        out[rs] = (ch, int(pos), "".join(sorted(gt)))
    return out


def load_ancestry(path):
    out = {}
    for f in _lines(path):
        if len(f) != 5 or f[0] == "rsid":
            continue
        rs, ch, pos, a1, a2 = f
        ch = ANCESTRY_CHROM.get(ch, ch)
        if ch not in CHROMS or "0" in (a1, a2):
            continue
        out[rs] = (ch, int(pos), "".join(sorted(a1 + a2)))
    return out


def merge(primary, secondary):
    """Union of two chips for one person; discordant shared SNPs become no-calls.

    Returns (merged, stats).
    """
    merged = dict(secondary)
    merged.update(primary)
    shared = primary.keys() & secondary.keys()
    discordant = [rs for rs in shared if _norm(primary[rs]) != _norm(secondary[rs])]
    for rs in discordant:
        del merged[rs]
    return merged, {"shared": len(shared), "discordant": len(discordant),
                    "only_primary": len(primary) - len(shared),
                    "only_secondary": len(secondary) - len(shared),
                    "total": len(merged)}


def _norm(rec):
    ch, _, gt = rec
    # Hemizygous X/Y/MT may be "A" on one chip and "AA" on the other.
    if ch in ("X", "Y", "MT") and len(gt) == 2 and gt[0] == gt[1]:
        return gt[0]
    return gt


def find(person, pattern):
    hits = sorted(glob.glob(os.path.join(RAW, person, pattern)))
    return hits[0] if hits else None


def has_genotypes(person):
    """True once a person's raw genotype file is in raw/<person>/ (a folder with only report archives doesn't count)."""
    return bool(find(person, "genome_*.zip") or find(person, "genome_*.txt") or find(person, "dna-data-*.zip")
                or find(person, "AncestryDNA*.txt"))


def load_person(person):
    """Load a person's best available measured genotypes (merging chips if present)."""
    t23 = find(person, "genome_*.zip") or find(person, "genome_*.txt")
    anc = find(person, "dna-data-*.zip") or find(person, "AncestryDNA*.txt")
    if not t23 and not anc:
        return None, None
    if t23 and anc:
        return merge(load_23andme(t23), load_ancestry(anc))
    data = load_23andme(t23) if t23 else load_ancestry(anc)
    return data, {"total": len(data)}


# GRCh38 positions (Ensembl) and plus-strand REF/ALT of single-marker reports that the imputed BCFs can fill when a
# chip lacks the probe (older chips such as 23andMe v4 lack several). DPYD HapB3, GAA IVS1 and GJB2 K15T are not in the panel.
IMPUTED38 = {
    "rs1799963": ("chr11", 46739505, "G", "A"),      # F2 prothrombin (stored under an internal ID on v5)
    "rs116855232": ("chr13", 48045719, "C", "T"),    # NUDT15 *3
    "rs4349859": ("chr6", 31398010, "G", "A"),       # HLA-B27 tag
    "rs116488202": ("chr6", 31377139, "C", "T"),     # HLA-B27 tag
    "rs73885319": ("chr22", 36265860, "A", "G"),     # APOL1 G1
    "rs4713586": ("chr6", 32691805, "A", "G"),       # HLA-DQ2.2 / DQ4 caveat
}


def imputed_genotypes(person, rsids):
    """{rsid: genotype} from the imputed GRCh38 BCFs, only where both haplotype dosages are confident (<0.1 or >0.9)."""
    out = {}
    for rs in rsids:
        if rs not in IMPUTED38:
            continue
        chrom, pos, ref, alt = IMPUTED38[rs]
        bcf = glob.glob(os.path.join(WORK, "imputed", person, "*", f"{chrom}.bcf"))
        if not bcf:
            continue
        proc = subprocess.run([BCFTOOLS, "query", "-r", f"{chrom}:{pos}",
                               "-f", "%REF\t%ALT\t[%HDS]\n", bcf[0]], capture_output=True, text=True)
        if proc.returncode:
            raise RuntimeError(f"bcftools query failed on {bcf[0]}: {proc.stderr.strip()}")
        for line in proc.stdout.splitlines():
            r, a, hds = line.split("\t")
            try:
                h = [float(x) for x in hds.split(",")]
            except ValueError:  # "." = no dosage at this site
                continue
            if (r, a) == (ref, alt) and all(x < 0.1 or x > 0.9 for x in h):
                out[rs] = "".join(sorted(alt if x > 0.5 else ref for x in h))
    return out


# ---------------------------------------------------------------- 1000 Genomes reference panel (public downloads)
KG_BASE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/"
KG_PANEL = (KG_BASE + "working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.{chrom}.filtered."
            "SNV_INDEL_SV_phased_panel.vcf.gz")
KG_PED = KG_BASE + "20130606_g1k_3202_samples_ped_population.txt"


def download(url, dest, timeout=60, check=None):
    """Fetch url to dest via a temporary file, so an interrupted or failed download never leaves a partial dest.
    `check(path)` may raise to reject the downloaded file."""
    import urllib.request
    tmp = f"{dest}.part"
    try:
        with urllib.request.urlopen(url, timeout=timeout) as r, open(tmp, "wb") as fh:
            shutil.copyfileobj(r, fh)
        if check:
            check(tmp)
        os.replace(tmp, dest)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def append_gzip(path, text):
    """Append text to a gzip cache as a new member, atomically: a run killed mid-write leaves the old file intact."""
    import gzip
    tmp = f"{path}.tmp"
    if os.path.exists(path):
        shutil.copyfile(path, tmp)
    with gzip.open(tmp, "at") as fh:
        fh.write(text)
    os.replace(tmp, path)


def write_atomic(path, text):
    tmp = f"{path}.tmp"
    with open(tmp, "w") as fh:
        fh.write(text)
    os.replace(tmp, path)


def kg_sample_files(ref_dir):
    """Make sure ref_dir has the 1000 Genomes sample order (samples.txt) and population table (ped.txt)."""
    os.makedirs(ref_dir, exist_ok=True)
    samples_path = os.path.join(ref_dir, "samples.txt")
    if not os.path.exists(samples_path):
        names = subprocess.run([BCFTOOLS, "query", "-l", KG_PANEL.format(chrom="chr22")], capture_output=True,
                               text=True, cwd=ref_dir, check=True).stdout
        if len(names.split()) < 3000:
            raise RuntimeError(f"1000 Genomes sample list looks truncated ({len(names.split())} samples)")
        write_atomic(samples_path, names)
    ped_path = os.path.join(ref_dir, "ped.txt")
    if not os.path.exists(ped_path):
        def check(tmp):
            with open(tmp) as fh:
                n = sum(1 for _ in fh)
            if n < 3000:
                raise RuntimeError(f"1000 Genomes population table looks truncated ({n} lines)")
        download(KG_PED, ped_path, check=check)
    return samples_path, ped_path


# ---------------------------------------------------------------- family configuration (DATA/family.json)
# Each person is a folder in raw/ with a role. Supported shape: one "self" (the main viewer), any number of
# "sibling"s (full siblings of self), "mother", "father" and the mother's mother ("maternal_grandmother").
# Everyone except self is optional; analyses that need a missing person are skipped.
ROLES = ["self", "sibling", "mother", "father", "maternal_grandmother"]
ROLE_LABEL = {"self": "Me", "mother": "Mom", "father": "Dad", "maternal_grandmother": "Grandma (maternal)"}
DEFAULT_FAMILY = {"people": {
    "me": {"role": "self"}, "brother": {"role": "sibling", "sex": "M"}, "sister": {"role": "sibling", "sex": "F"},
    "mom": {"role": "mother", "sex": "F"}, "dad": {"role": "father", "sex": "M"},
    "grandma-maternal": {"role": "maternal_grandmother", "sex": "F"}}}


def load_family():
    """family.json from the data folder; without one, the conventional folder names above are used."""
    path = os.path.join(DATA, "family.json")
    fam = json.load(open(path, encoding="utf-8")) if os.path.exists(path) else DEFAULT_FAMILY
    fam = dict(fam)
    ppl = fam.get("people", {})
    if not os.path.exists(path):  # defaults: keep only folders that exist
        ppl = {p: v for p, v in ppl.items() if os.path.isdir(os.path.join(RAW, p))}
    # Role order, then file order within a role, so "self" comes first and siblings follow.
    fam["people"] = dict(sorted(ppl.items(), key=lambda kv: ROLES.index(kv[1]["role"])))
    return fam


FAMILY = load_family()


def people():
    """All configured person ids (raw/ folder names): self, siblings, mother, father, maternal grandmother."""
    return list(FAMILY["people"])


def person(role):
    """The id of the single person with this role (self, mother, father, maternal_grandmother), or None."""
    return next((p for p, v in FAMILY["people"].items() if v["role"] == role), None)


def siblings():
    return [p for p, v in FAMILY["people"].items() if v["role"] == "sibling"]


def children():
    """The mother's and father's children in the data: self first, then siblings."""
    return [p for p in (person("self"), *siblings()) if p]


def role(p):
    return FAMILY["people"].get(p, {}).get("role")


def sex(p):
    """"M", "F" or None (unknown) from family.json."""
    return FAMILY["people"].get(p, {}).get("sex")


def name(p):
    """A person's own first name if family.json gives one, else their role label (used where "Me" is ambiguous)."""
    v = FAMILY["people"].get(p, {})
    return v.get("name") or label(p)


def label(p, viewer=None):
    """Column label from the viewer's point of view: "Me", "Mom", "Dad", "Brother", "Sister", "Grandma (maternal)"."""
    viewer = viewer or person("self")
    if p == viewer:
        return "Me"
    r = role(p)
    if r in ("self", "sibling"):
        return {"M": "Brother", "F": "Sister"}.get(sex(p), "Sibling")
    return ROLE_LABEL.get(r, p)


def text(key, default=""):
    """Optional family-specific wording from family.json "text" (genealogy context, medical context, ...)."""
    return FAMILY.get("text", {}).get(key, default)


# Shown at the top of every generated page. Inline styles so it works with each page's own CSS, light or dark.
DISCLAIMER_HTML = (
    '<div role="note" style="border:1.5px solid #c9a227;background:rgba(201,162,39,.12);border-radius:6px;'
    'padding:8px 12px;margin:0 0 16px;font-size:14px;line-height:1.45">'
    "<b>For learning and research only. Not medical advice.</b> These results come from consumer genotyping chips "
    "analysed with open-source research code. They are not a diagnosis and have not been clinically validated. "
    "Anything health-related, including medication notes, needs review by a clinician (a doctor or genetic "
    "counselor) and confirmation with a clinical test before any decision about care.</div>")


def is_autosome(ch):
    return ch.isdigit()


def is_het(gt):
    return len(gt) == 2 and gt[0] != gt[1]
