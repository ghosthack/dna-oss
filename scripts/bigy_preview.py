"""What the chip already says below R-P311, as input for the landing page's Big Y-700 preview.

Counts derived/ancestral calls per major R1b-P311 branch using yhaplo's ISOGG 2016 SNP list and work/y.genos.txt.
Only men placed at or below R-P311 in results/y_haplogroups.tsv are counted; with none, nothing is written.
Usage: python3 -I scripts/bigy_preview.py  → results/bigy_preview.json
"""
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dnalib import RESULTS, WORK, children, person

GENOS = os.path.join(WORK, "y.genos.txt")
ISOGG = os.path.join(WORK, "yhaplo", "isogg.snps.cleaned.2016.01.04.txt")
BRANCHES = {  # ISOGG 2016 code → common name, most specific first when matching
    "R1b1a2a1a2a": "DF27", "R1b1a2a1a2b": "U152", "R1b1a2a1a2c": "L21", "R1b1a2a1a2d": "L238",
    "R1b1a2a1a2e": "DF19", "R1b1a2a1a2f": "DF99", "R1b1a2a1a1": "U106",
}
KEY_SNPS = {"DF27": 21380200, "P312": 22157311}  # GRCh37
P311_ISOGG = "R1b1a2a1a"  # ISOGG 2016 code of R-P311
OUT = os.path.join(RESULTS, "bigy_preview.json")


def p311_men():
    """Person ids that yhaplo placed at or below R-P311 (results/y_haplogroups.tsv), plus father-line men whose
    sparser chip stops just above it (e.g. R-L51) when another man of the same line reaches R-P311."""
    path = os.path.join(RESULTS, "y_haplogroups.tsv")
    if not os.path.exists(path):
        return set()
    rows = [l.rstrip("\n").split("\t") for l in open(path)][1:]
    codes = {r[0]: r[2] for r in rows if len(r) > 2}
    men = {p for p, c in codes.items() if c.startswith(P311_ISOGG)}
    line = [p for p in (person("father"), *children()) if p in codes]  # father and sons share one Y
    if men & set(line):
        men |= {p for p in line if P311_ISOGG.startswith(codes[p]) and len(codes[p]) >= len(P311_ISOGG) - 1}
    return men


def main():
    men = p311_men()
    if not men or not os.path.exists(GENOS):
        print("no man at or below R-P311: no Big Y preview")
        if os.path.exists(OUT):
            os.remove(OUT)  # don't leave a preview from an earlier run
        return
    lines = open(GENOS).read().split("\n")
    idx = {int(p): i for i, p in enumerate(lines[0].split("\t")[1:]) if p}
    people = {l.split("\t")[0]: l.split("\t")[1:] for l in lines[1:] if l and l.split("\t")[0] in men}
    out = {"branches": {}, "key_snps": {}}
    for person, calls in people.items():
        counts = {name: {"derived": 0, "ancestral": 0} for name in BRANCHES.values()}
        for line in open(ISOGG):
            m = re.match(r"(\S+)\s+(R1b1a2a1a\S+)\s+(\d+)\s+([ACGT])->([ACGT])", line)
            if not m or int(m.group(3)) not in idx:
                continue
            branch = next((n for code, n in sorted(BRANCHES.items(), key=lambda kv: -len(kv[0]))
                           if m.group(2).startswith(code)), None)
            if branch:
                g = calls[idx[int(m.group(3))]]
                if g == m.group(5):
                    counts[branch]["derived"] += 1
                elif g == m.group(4):
                    counts[branch]["ancestral"] += 1
        out["branches"][person] = counts
        out["key_snps"][person] = {n: (calls[idx[p]] if p in idx else "not on chip") for n, p in KEY_SNPS.items()}
    json.dump(out, open(OUT, "w"), indent=1)
    print(OUT)


if __name__ == "__main__":
    main()
