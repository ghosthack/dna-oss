"""Chip-testable hereditary cancer markers, for the doctor pages.

Colorectal: MUTYH Y179C and G396D (biallelic = MUTYH-associated polyposis; one copy = carrier) and APC I1307K
(low-penetrance CRC risk allele, ~1.5-2x, mainly Ashkenazi; not a polyposis variant). Breast/ovarian: the three
Ashkenazi BRCA founder variants 23andMe tests (i-probes, I/D-coded). Plus-strand alt alleles on GRCh37.
A negative result excludes almost nothing: Lynch syndrome and the many
other BRCA1/2, PALB2, APC and MUTYH variants are not reliably tested by arrays.
Usage: python3 -I scripts/hereditary_markers.py  → results/hereditary_markers.tsv
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dnalib import RESULTS, load_person, people

MARKERS = {  # probe → (variant allele, name, panel)
    "rs34612342": ("C", "MUTYH Y179C", "colorectal"),
    "rs36053993": ("T", "MUTYH G396D", "colorectal"),
    "rs1801155": ("A", "APC I1307K", "colorectal"),
    "i4000377": ("D", "BRCA1 185delAG", "breast-ovarian"),
    "i4000378": ("I", "BRCA1 5382insC", "breast-ovarian"),
    "i4000379": ("D", "BRCA2 6174delT", "breast-ovarian"),
}


def copies(gt, alt, probe, person):
    """Variant allele count as text; "" (= not tested) for a no-call or an encoding we don't recognise.

    i-probes for indels are expected as D/I letters. If a chip version coded one with bases instead, counting "D"
    would silently report 0 copies, so such a genotype is treated as untested and reported."""
    if not gt:
        return ""
    expected = set("DI") if alt in "DI" else set("ACGT")
    if not set(gt) <= expected:
        print(f"WARNING: {person} {probe} genotype {gt!r} is not in the expected {'/'.join(sorted(expected))} "
              "coding; reported as not tested", file=sys.stderr)
        return ""
    return str(gt.count(alt))


def main():
    rows = []
    for p in sorted(people()):
        data, _ = load_person(p)
        if not data:
            continue
        for probe, (alt, name, panel) in MARKERS.items():
            g = data.get(probe)
            rows.append({"person": p, "panel": panel, "marker": name, "probe": probe,
                         "alt_copies": copies(g and g[2], alt, probe, p)})
    path = os.path.join(RESULTS, "hereditary_markers.tsv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, ["person", "panel", "marker", "probe", "alt_copies"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(path, len(rows), "rows")


if __name__ == "__main__":
    main()
