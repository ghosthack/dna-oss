"""Scan the files git would commit for private family data before it leaves the machine.

Run before every commit (scripts/install_hooks.sh installs a pre-commit hook). Exits 1 if anything is found.

  python3 scripts/leak_check.py [--repo DIR] [--data DIR] [--staged] [--no-genotypes]

What is checked (every file from git ls-files --cached --others --exclude-standard, or the index with --staged):
  denylist     every string value in DATA/family.json, split into tokens: people names and ids, the virtual
               grandfather's name, words of 4+ letters from labels, places and text (prose is filtered against a
               dictionary of common words), 4-digit years (1800-2030) in text and ancestry values, the GEDCOM path
               and its basename, page filenames, and each line of an optional DATA/leak_denylist.txt. Accent- and
               case-insensitive, whole words. File paths are checked too.
  generic      email addresses, home-directory paths (/Users/<name>/, /home/<name>/), consumer-DNA download
               filenames that embed a person's name, GEDCOM cross-reference ids outside test files.
  fingerprint  values measured for this family in DATA/results: star-allele diplotypes (on a line with their gene,
               or anywhere when distinctive), mt and Y haplogroup labels (very generic ones such as "R1b" are allowed),
               and "rsNNN" written next to a two-letter genotype that equals a family member's measured genotype.

DATA is $DNA_DATA, else the repo root if it holds a family.json. Without a family.json only the generic checks run.

The report never prints a private value: only the category, file:line and a masked hint (first letter and length);
genotype hits print the rsID only. Paths that contain a denylisted word are masked as well.

Allowing a hit: put "leak-ok" anywhere on the line, or add a regex to .leakallow in the repo root (one per line,
matched against the flagged text; "path:<regex>" skips matching files). .leakallow is committed, so never put a
private value in it: private exceptions go in DATA/leak_allow.txt (same format, never committed).
"""
import argparse
import csv
import fnmatch
import json
import os
import re
import subprocess
import sys
import unicodedata

HERE = os.path.dirname(os.path.abspath(__file__))
MARK = "leak-ok"

# Words that may appear in family.json labels but are too generic to flag: population, panel and method names that
# the analysis code itself uses as identifiers or table labels (reference groups, AADR subgroups, file and field
# names); plurals ("Italians") are covered by generic_term(). Keep it to words the scripts use: a word that only
# appears in family.json or in prose stays on the denylist.
GENERIC_TERMS = set("""
italian italy european europe spanish spaniards spain iberia iberian sicily sicilian greek balkan asian
american albanian tuscan campania campanian calabria calabrian griko salento lecce campcal naples salerno
bulgarian croatian slavic levant mediterranean british finnish utah ceph
aadr gedcom ancestrydna mahalanobis slco neanderthal sprime haplogroup
""".split())
# 4-digit years (birth, emigration, record dates) in family.json free text and ancestry values. Round centuries
# ("1900s", "~2000 BCE") and BCE/BC dates are too common to identify anyone and are skipped.
YEAR = re.compile(r"(?<!\d)(18\d\d|19\d\d|20[0-2]\d|2030)(?!\d)(?!\s*BCE?\b)")

STOPWORDS = set("""
self sibling siblings mother father brother sister maternal paternal grandmother grandfather grandma grandpa mom dad
family parent parents child children unknown assumed region regions label place places founder founders county
province state country city town village parish island river valley north south east west central northern southern
eastern western upper lower near from with born died line lineage none null true false yes also this that their
them they have were been about into only some more most other than then when where which while ancestry genealogy
text name names page pages index html json people person role default note notes reference population proxy
hypothesis hypotheses
""".split())

# Person ids (raw/ folder names) that are conventional and not personal.
GENERIC_IDS = {"me", "self", "brother", "sister", "sibling", "mom", "dad", "mother", "father", "grandma",
               "grandma-maternal", "grandmother", "maternal-grandmother", "index"}
# Placeholder user names that are fine in documentation.
GENERIC_USERS = {"user", "username", "name", "you", "yourname", "your-name", "me", "example", "shared", "runner",
                 "USER", "USERNAME", "Shared"}
GENERIC_EMAIL = re.compile(r"(@example\.(com|org|net)$|^git@|^noreply@|@users\.noreply\.github\.com$)", re.I)
# Haplogroup labels too broad to identify anyone ("H", "U5", "R1b", "J1c", "HV").
GENERIC_HAPLO = re.compile(r"^[A-Z]{1,2}\d?[a-z]?$")

EMAIL = re.compile(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9-]+(?:\.[A-Za-z0-9-]+)*\.[A-Za-z]{2,}")
HOME = re.compile(r"/(?:Users|home)/([A-Za-z0-9._-]+)/")
DOWNLOAD = re.compile(r"genome_(?!\*)([A-Za-z][A-Za-z-]*(?:_[A-Za-z][A-Za-z-]*)*?)_v\d+_"
                      r"|imputed_genotype_data_[A-Za-z0-9_-]*?_([A-Z][A-Za-z-]+)\b")
DOWNLOAD_PLACEHOLDER = re.compile(r"^(name|yourname|first(name)?(_last(name)?)?|john_doe|jane_doe|x+)$", re.I)
XREF = re.compile(r"@[A-Z]{1,2}\d+@")
# rsID followed (after a little punctuation) by two alleles: "rs123 AG", "rs123: A/G", "rs123 (A;G)", "rs123='TT'".
RS_GT = re.compile(r"(?<![A-Za-z0-9])(rs\d{2,})(?![0-9])[\s:=,(\[\"'>-]{0,4}([ACGTDI])\s?[/;|]?\s?([ACGTDI])(?![A-Za-z0-9])")


def norm(s):
    """Accent- and case-insensitive form used for denylist matching."""
    s = unicodedata.normalize("NFKD", s)
    return "".join(c for c in s if not unicodedata.combining(c)).casefold()


def mask(s):
    s = s.strip()
    return f"{s[:1]}…[{len(s)}]" if s else "…"


def word_pattern(tok):
    """Whole-word regex for a normalized token (letters/digits on either side break the match, "_" does not)."""
    e = re.escape(tok)
    if re.match(r"\w", tok):
        e = r"(?<![^\W_])" + e
    if re.search(r"\w$", tok):
        e = e + r"(?![^\W_])"
    return e


# ---------------------------------------------------------------- private inputs (DATA)
def data_dir(arg):
    if arg:
        return os.path.abspath(os.path.expanduser(arg))
    if os.environ.get("DNA_DATA"):
        return os.path.abspath(os.path.expanduser(os.environ["DNA_DATA"]))
    # Same default as dnalib.DATA: the repo itself, used only if someone put their family.json there.
    d = os.path.dirname(HERE)
    return d if os.path.exists(os.path.join(d, "family.json")) else None


def common_words():
    """Lower-case dictionary words (proper nouns are capitalized there, so they are kept as private)."""
    words = set(STOPWORDS)
    for path in ("/usr/share/dict/words", "/usr/dict/words"):
        if os.path.exists(path):
            with open(path, encoding="utf-8", errors="ignore") as fh:
                words.update(w.strip() for w in fh if w[:1].islower())
            break
    return words


def generic_term(n):
    """True for a GENERIC_TERMS word or its plural."""
    return n in GENERIC_TERMS or n.endswith("s") and n[:-1] in GENERIC_TERMS


def denylist(data):
    """{normalized token: kind} from family.json and leak_denylist.txt."""
    toks = {}
    common = None

    def add(tok, kind, digits=False):
        t = norm(tok.strip())
        if len(t) >= 2 and (digits or not t.isdigit()):
            toks.setdefault(t, kind)

    def add_years(s):
        # birth and record years: the only all-digit tokens kept (round centuries are too common to identify anyone)
        for y in YEAR.findall(s):
            if not y.endswith("00"):
                add(y, "year", digits=True)

    def add_name(s, kind):
        add(s, kind)
        for w in re.findall(r"[^\W\d_]{2,}", s):
            if norm(w) not in STOPWORDS:
                add(w, kind)

    def is_common(n):
        nonlocal common
        if common is None:
            common = common_words()
        # also accept simple inflections ("variants", "carried", "testing") of dictionary words
        return n in STOPWORDS or any(n.endswith(suf) and n[:len(n) - len(suf)] + add_back in common
                                     for suf, add_back in (("", ""), ("s", ""), ("es", ""), ("ed", ""), ("ed", "e"),
                                                           ("ing", ""), ("ing", "e"), ("ies", "y"), ("'s", "")))

    def add_words(s, kind):
        # Free text is mostly ordinary English: only capitalized non-dictionary words (names, places) are private.
        # Short labels (founder labels, places, surnames) keep every non-dictionary word.
        prose = kind == "text" or len(s.split()) >= 6
        for w in re.findall(r"[^\W\d_]{4,}", s):
            n = norm(w)
            if is_common(n) or (prose and not w[:1].isupper()) or generic_term(n):
                continue
            add(w, kind)

    def walk(v, path):
        if isinstance(v, dict):
            for k, x in v.items():
                walk(x, path + (str(k),))
        elif isinstance(v, list):
            for x in v:
                walk(x, path)
        elif isinstance(v, str):
            if len(path) == 3 and path[0] == "people" and path[2] in ("role", "sex"):
                return
            if (path[:1] == ("people",) and path[-1] == "name") or path == ("maternal_grandfather", "name"):
                add_name(v, "name")
            elif path == ("genealogy", "gedcom"):
                add(v, "gedcom")
                base = os.path.basename(v)
                add(base, "gedcom")
                stem = os.path.splitext(base)[0]
                if len(stem) >= 4:
                    add(stem, "gedcom")
            elif path == ("genealogy", "self_id"):
                add(v, "gedcom")
            elif path[:1] == ("pages",):
                stem = os.path.splitext(os.path.basename(v))[0]
                if norm(stem) not in GENERIC_IDS:
                    add_name(stem, "filename")
            elif path[:2] == ("ancestry", "founders") and path[-1] == "label" or \
                    path[:2] == ("ancestry", "places") and path[-1] in ("place", "surnames"):
                add_words(v, "ancestry")  # surnames and birthplaces
                add_years(v)
            elif path[:1] in (("text",), ("ancestry",)):
                # free text and model settings are mostly generic science vocabulary (populations, panels);
                # only words that look like names (capitalized, non-dictionary, not a population term) count
                add_words(v, "text")
                add_years(v)
            else:
                add_words(v, path[0] if path else "value")  # only the top-level key: deeper keys can be ids

    fam_path = os.path.join(data, "family.json")
    fam = None
    if os.path.exists(fam_path):
        try:
            with open(fam_path, encoding="utf-8") as fh:
                fam = json.load(fh)
        except (OSError, ValueError) as e:
            print(f"leak_check: cannot parse family.json ({type(e).__name__}); denylist skipped", file=sys.stderr)
    if fam:
        walk(fam, ())
        for pid in fam.get("people", {}):
            if norm(pid) not in GENERIC_IDS:
                add_name(pid, "person id")
        for viewer in fam.get("pages", {}):
            if norm(viewer) not in GENERIC_IDS:
                add_name(viewer, "person id")
    raw = os.path.join(data, "raw")
    if os.path.isdir(raw):  # folder names are person ids even without family.json
        for pid in os.listdir(raw):
            if os.path.isdir(os.path.join(raw, pid)) and norm(pid) not in GENERIC_IDS:
                add_name(pid, "person id")
    extra = os.path.join(data, "leak_denylist.txt")
    if os.path.exists(extra):
        with open(extra, encoding="utf-8") as fh:
            for line in fh:
                if line.strip() and not line.lstrip().startswith("#"):
                    add(line.strip(), "denylist file")
    return toks, fam


def read_tsv(data, name):
    path = os.path.join(data, "results", name)
    if not os.path.exists(path):
        return []
    with open(path, encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def star_alleles(dip):
    return re.findall(r"\*\d+[A-Za-z]*", dip)


def fingerprints(data):
    """Diplotypes [(gene, variant, anywhere)] and haplogroup labels [(kind, label)] measured in this family."""
    dips = {}
    for r in read_tsv(data, "pharmacogenomics.tsv"):
        gene, dip = (r.get("gene") or "").strip(), (r.get("diplotype") or "").strip()
        # "ambiguous (variants at *3 ×2, *6)": the alleles are the fingerprint; otherwise drop the "(…)" note
        dip = re.sub(r"\s*×\d+", "", dip[11:-1]).replace(", ", "+") if dip.startswith("ambiguous (variants at ") \
            else dip.split(" (")[0].strip()
        dip = dip.replace("variants at ", "")
        if not gene or not dip or dip in ("*1/*1", "?", "NA", "-"):
            continue
        stars = star_alleles(dip)
        # Distinctive: no reference allele, a sub-allele letter suffix, or an ambiguous "+" call.
        anywhere = bool(stars) and (all(s != "*1" for s in stars) or any(s[-1].isalpha() for s in stars) or "+" in dip)
        variants = {dip}
        if "/" in dip:
            a, b = dip.split("/", 1)
            variants.add(f"{b}/{a}")
        for v in variants:
            dips[(gene, v)] = dips.get((gene, v), False) or anywhere
    haplos = set()
    for r in read_tsv(data, "mt_haplogroups.tsv"):
        haplos.add(("mt haplogroup", (r.get("Haplogroup") or "").strip()))
    for r in read_tsv(data, "y_haplogroups.tsv"):
        haplos.add(("Y haplogroup", (r.get("haplogroup") or "").strip()))
        haplos.add(("Y haplogroup", (r.get("isogg_2016") or "").strip()))
    haplos = sorted((k, h) for k, h in haplos if h and not GENERIC_HAPLO.match(h))
    return [(g, v, a) for (g, v), a in sorted(dips.items())], haplos


def measured(data, rsids):
    """{rsid: {genotype, ...}} over everyone in DATA/raw, for the given rsIDs only. Never printed."""
    os.environ["DNA_DATA"] = data
    sys.path.insert(0, HERE)
    try:
        import dnalib
    except Exception as e:  # noqa: BLE001 - report and carry on with the other checks
        print(f"leak_check: cannot import dnalib ({type(e).__name__}); genotype check skipped", file=sys.stderr)
        return None
    ids = dnalib.people()
    if os.path.isdir(dnalib.RAW):
        ids += [p for p in sorted(os.listdir(dnalib.RAW)) if p not in ids and os.path.isdir(os.path.join(dnalib.RAW, p))]
    out = {}
    for p in ids:
        try:
            d, _ = dnalib.load_person(p)
        except Exception as e:  # noqa: BLE001
            print(f"leak_check: cannot read raw data for one person ({type(e).__name__}); skipped", file=sys.stderr)
            continue
        for rs in rsids & (d or {}).keys():
            gt = d[rs][2]
            out.setdefault(rs, set()).add(gt * 2 if len(gt) == 1 else gt)  # hemizygous "A" reads as "AA"
    return out


# ---------------------------------------------------------------- repo files
def git(repo, *args, data=None):
    return subprocess.run(["git", "-C", repo, *args], capture_output=True, input=data, check=True).stdout


def ignore_rules(repo):
    rules = []
    path = os.path.join(repo, ".gitignore")
    if os.path.exists(path):
        with open(path, encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line.strip() or line.startswith("#"):
                    continue
                neg = line.startswith("!")
                pat = line[1:] if neg else line
                rules.append((neg, pat.endswith("/"), pat.startswith("/") or "/" in pat.rstrip("/"), pat.strip("/")))
    return rules


def ignored(rel, is_dir, rules):
    hit = False
    for neg, dir_only, anchored, pat in rules:
        if dir_only and not is_dir:
            continue
        ok = fnmatch.fnmatch(rel, pat) if anchored else fnmatch.fnmatch(os.path.basename(rel), pat)
        if ok:
            hit = not neg
    return hit


def walk_files(repo):
    """Fallback without git: every file under repo, honoring the top-level .gitignore."""
    rules = ignore_rules(repo)
    for root, dirs, files in os.walk(repo):
        relroot = os.path.relpath(root, repo)
        relroot = "" if relroot == "." else relroot
        keep = []
        for d in sorted(dirs):
            rel = os.path.join(relroot, d)
            if d == ".git":
                continue
            if os.path.islink(os.path.join(root, d)):
                if not ignored(rel, False, rules):
                    files.append(d)  # git stores a symlink as a file holding its target
            elif not ignored(rel, True, rules):
                keep.append(d)
        dirs[:] = keep
        for f in sorted(files):
            rel = os.path.join(relroot, f)
            if not ignored(rel, False, rules):
                yield rel


def read_worktree(repo, rel):
    path = os.path.join(repo, rel)
    if os.path.islink(path):
        return os.readlink(path).encode()
    try:
        with open(path, "rb") as fh:
            return fh.read()
    except OSError:
        return None


def repo_files(repo, staged):
    """[(relative path, bytes or None)] of everything git would commit."""
    try:
        if staged:
            out = []
            entries = [e for e in git(repo, "ls-files", "-s", "-z").split(b"\0") if e]
            blobs = []
            for e in entries:
                meta, rel = e.split(b"\t", 1)
                mode, sha, _ = meta.split(b" ")
                if mode != b"160000":  # submodules carry no content here
                    blobs.append((rel.decode("utf-8", "replace"), sha))
            if blobs:
                raw = git(repo, "cat-file", "--batch", data=b"".join(s + b"\n" for _, s in blobs))
                pos = 0
                for rel, _ in blobs:
                    nl = raw.index(b"\n", pos)
                    size = int(raw[pos:nl].split(b" ")[2])
                    out.append((rel, raw[nl + 1:nl + 1 + size]))
                    pos = nl + 1 + size + 1
            return out
        rels = git(repo, "ls-files", "-z", "--cached", "--others", "--exclude-standard").split(b"\0")
        rels = sorted({r.decode("utf-8", "replace") for r in rels if r})
    except (OSError, subprocess.CalledProcessError):
        if staged:
            sys.exit("leak_check: --staged needs a git repository")
        rels = list(walk_files(repo))
    return [(r, read_worktree(repo, r)) for r in rels]


def is_test(rel):
    parts = rel.replace("\\", "/").split("/")
    base = parts[-1]
    return any(p in ("test", "tests", "testdata", "fixtures") for p in parts[:-1]) or base.startswith("test_") \
        or re.search(r"_test\.\w+$", base) is not None


def load_allow(*paths):
    text_rx, path_rx = [], []
    for path in paths:
        if not path or not os.path.exists(path):
            continue
        with open(path, encoding="utf-8") as fh:
            for n, line in enumerate(fh, 1):
                line = line.rstrip("\n")
                if not line.strip() or line.startswith("#"):
                    continue
                try:
                    if line.startswith("path:"):
                        path_rx.append(re.compile(line[5:]))
                    else:
                        text_rx.append(re.compile(line, re.I))
                except re.error:
                    print(f"leak_check: bad regex in {os.path.basename(path)}:{n}; ignored", file=sys.stderr)
    return text_rx, path_rx


# ---------------------------------------------------------------- scan
def main():
    ap = argparse.ArgumentParser(description="Fail if files to be committed contain private family data.")
    ap.add_argument("--repo", help="repository to scan (default: the git repo of the current directory)")
    ap.add_argument("--data", help="private data folder (default: $DNA_DATA, else the repo root if it has family.json)")
    ap.add_argument("--staged", action="store_true", help="scan the staged (index) content instead of the work tree")
    ap.add_argument("--no-genotypes", action="store_true", help="skip loading raw genotype files")
    args = ap.parse_args()

    repo = args.repo
    if not repo:
        try:
            repo = git(os.getcwd(), "rev-parse", "--show-toplevel").decode().strip()
        except (OSError, subprocess.CalledProcessError):
            repo = os.path.dirname(HERE)
    repo = os.path.abspath(repo)
    data = data_dir(args.data)
    if data and not os.path.isdir(data):
        print("leak_check: data folder not found; only generic checks run", file=sys.stderr)
        data = None

    toks, dips, haplos = {}, [], []
    if data:
        toks, fam = denylist(data)
        if fam is None:
            print("leak_check: no usable family.json in the data folder; family denylist skipped", file=sys.stderr)
        dips, haplos = fingerprints(data)
    elif not args.data and not os.environ.get("DNA_DATA"):
        print("leak_check: no data folder (set DNA_DATA); only generic checks run", file=sys.stderr)
    deny_rx = re.compile("|".join(word_pattern(t) for t in sorted(toks, key=len, reverse=True))) if toks else None
    dip_rx = [(re.compile(r"(?<![A-Za-z0-9])" + re.escape(g) + r"(?![A-Za-z0-9])"),
               re.compile(r"(?<![\w*])" + re.escape(v) + r"(?![\w])"), g, v, a) for g, v, a in dips]
    haplo_rx = [(k, h, re.compile(r"(?<![\w'-])" + re.escape(h) + r"(?![\w'+-])")) for k, h in haplos]
    allow_text, allow_path = load_allow(os.path.join(repo, ".leakallow"),
                                        data and os.path.join(data, "leak_allow.txt"))

    def masked_path(rel):
        if not deny_rx:
            return rel
        return "/".join(f"<{mask(c)}>" if deny_rx.search(norm(c)) else c for c in rel.split("/"))

    hits = []          # (category, rel, line, hint)
    rs_candidates = []  # (rel, line, rsid, genotype): resolved after the scan

    def allowed(subject, rel, line=""):
        # A rule matches the flagged value or its whole line. .leakallow is committed, so its own rules never
        # excuse a hit inside it.
        return rel != ".leakallow" and any(r.search(subject) or r.search(line) for r in allow_text)

    def hit(cat, rel, n, subject, hint, line=""):
        if not allowed(subject, rel, line):
            hits.append((cat, rel, n, hint))

    files = repo_files(repo, args.staged)
    for rel, content in files:
        if any(r.search(rel) for r in allow_path) and rel != ".leakallow":
            continue
        if deny_rx:
            for m in deny_rx.finditer(norm(rel)):
                hit("denylist", rel, 0, m.group(0), f"{toks[m.group(0)]} in file name {mask(m.group(0))}")
        if content is None or b"\0" in content[:8192]:
            continue  # missing or binary
        test = is_test(rel)
        for n, line in enumerate(content.decode("utf-8", "replace").splitlines(), 1):
            if MARK in line:
                continue
            if deny_rx:
                for m in deny_rx.finditer(norm(line)):
                    hit("denylist", rel, n, m.group(0), f"{toks[m.group(0)]} {mask(m.group(0))}", line)
            for m in EMAIL.finditer(line):
                if not GENERIC_EMAIL.search(m.group(0)):
                    hit("generic", rel, n, m.group(0), f"email {mask(m.group(0))}", line)
            for m in HOME.finditer(line):
                if m.group(1) not in GENERIC_USERS and not m.group(1).startswith("$"):
                    hit("generic", rel, n, m.group(0), f"home path /…/{mask(m.group(1))}/", line)
            for m in DOWNLOAD.finditer(line):
                who = m.group(1) or m.group(2)
                if not DOWNLOAD_PLACEHOLDER.match(who):
                    hit("generic", rel, n, m.group(0), f"download filename with a name {mask(who)}", line)
            if not test:
                for m in XREF.finditer(line):
                    hit("generic", rel, n, m.group(0), f"GEDCOM xref {mask(m.group(0))}", line)
            for gene_rx, var_rx, gene, v, anywhere in dip_rx:
                if var_rx.search(line) and (anywhere or gene_rx.search(line)):
                    hit("fingerprint", rel, n, v, f"{gene} diplotype {mask(v)}", line)
            for kind, h, rx in haplo_rx:
                if rx.search(line):
                    hit("fingerprint", rel, n, h, f"{kind} {mask(h)}", line)
            for m in RS_GT.finditer(line):
                rs_candidates.append((rel, n, m.group(1), "".join(sorted(m.group(2) + m.group(3))), line))

    if rs_candidates and data and not args.no_genotypes:
        gts = measured(data, {c[2] for c in rs_candidates})
        for rel, n, rs, gt, line in rs_candidates if gts is not None else []:
            if gt in gts.get(rs, ()):
                hit("fingerprint", rel, n, rs, f"{rs} with a family member's measured genotype", line)

    hits = sorted(set(hits), key=lambda h: (h[1], h[2], h[0], h[3]))
    if not hits:
        print(f"leak_check: {len(files)} files clean")
        return 0
    print(f"leak_check: {len(hits)} possible leak(s) in {len(files)} files (values masked):")
    for cat, rel, n, hint in hits:
        print(f"  {cat:11s} {masked_path(rel)}:{n}  {hint}")
    print("Fix them, or mark a false positive with 'leak-ok' on the line or a regex in .leakallow "
          "(private exceptions: DATA/leak_allow.txt).")
    return 1


if __name__ == "__main__":
    sys.exit(main())
