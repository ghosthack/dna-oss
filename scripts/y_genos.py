"""Write chrY genotypes of the males as a yhaplo sample-major text file.

Usage: python3 -I scripts/y_genos.py
Writes work/y.genos.txt (then run: .venv/bin/yhaplo -i work/y.genos.txt -o work/yhaplo -hp -hpd).
Males are the people with sex "M" in family.json (fathers are male by role); where the sex is not given,
the chip-based sex from results/summary.json (analyze.py) is used, so this agrees with the pages.
If there are no males, no y.genos.txt is written (and a stale one is removed).
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze import is_male
from dnalib import RESULTS, WORK, load_person, people, role, sex

ROLE_SEX = {"father": "M", "mother": "F", "maternal_grandmother": "F"}


def qc_sex():
    """{person: "M"/"F"} as analyze.py recorded it in summary.json ({} if not there)."""
    try:
        with open(os.path.join(RESULTS, "summary.json")) as fh:
            return {p: q.get("sex") for p, q in json.load(fh).get("qc", {}).items()}
    except (OSError, ValueError, AttributeError):
        return {}


def y_calls(data):
    calls = {}
    for ch, pos, gt in data.values():
        if ch == "Y" and len(set(gt)) == 1 and gt[0] in "ACGT":
            calls[pos] = gt[0]
    return calls


def main():
    males, chip = {}, qc_sex()
    path = os.path.join(WORK, "y.genos.txt")
    for p in people():
        s = sex(p) or ROLE_SEX.get(role(p)) or chip.get(p)
        if s == "F":
            continue
        data, _ = load_person(p)
        if not data:
            continue
        if s is None and not is_male(data):  # no summary.json yet: same rule analyze.py uses
            continue
        males[p] = y_calls(data)
    if not males:
        if os.path.exists(path):
            os.remove(path)  # don't leave genotypes from an earlier run for yhaplo
        print("no males with genotype data: no y.genos.txt, no Y haplogroups")
        return
    positions = sorted(set().union(*males.values()))
    with open(path, "w") as fh:
        fh.write("ID\t" + "\t".join(map(str, positions)) + "\n")
        for p, c in males.items():
            fh.write(p + "\t" + "\t".join(c.get(pos, ".") for pos in positions) + "\n")
    for p, c in males.items():
        print(f"{p}: {len(c)} Y sites")


if __name__ == "__main__":
    main()
