"""Parse the family GEDCOM and list the direct ancestors of one person with birth places.

Usage: python3 -I scripts/pedigree.py [path.ged] [root_id]
The GEDCOM path and the root person's id default to family.json "genealogy": {"gedcom", "self_id"}
(a relative path is taken from the data folder). Without a GEDCOM nothing is written.
Writes work/genealogy/pedigree.json (Ahnentafel-numbered ancestors).
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dnalib import DATA, FAMILY, WORK  # noqa: E402

GEN = FAMILY.get("genealogy", {})
GED = sys.argv[1] if len(sys.argv) > 1 else GEN.get("gedcom") and os.path.join(DATA, os.path.expanduser(GEN["gedcom"]))
ROOT_ID = sys.argv[2] if len(sys.argv) > 2 else GEN.get("self_id", "@I0@")
OUT = os.path.join(WORK, "genealogy", "pedigree.json")


def parse(path):
    indi, fam, cur, sub = {}, {}, None, None
    for raw in open(path, encoding="utf-8-sig"):
        parts = raw.rstrip("\r\n").split(" ", 2)
        if len(parts) < 2:
            continue
        level, tag = int(parts[0]), parts[1]
        val = parts[2] if len(parts) > 2 else ""
        if level == 0:
            sub = None
            if val == "INDI":
                cur = indi.setdefault(tag, {"id": tag, "famc": []})
            elif val == "FAM":
                cur = fam.setdefault(tag, {"id": tag})
            else:
                cur = None
            continue
        if cur is None:
            continue
        if level == 1:
            sub = tag
            if tag == "NAME" and "name" not in cur:
                cur["name"] = val.replace("/", "").strip()
            elif tag == "SEX":
                cur["sex"] = val
            elif tag == "FAMC":
                cur["famc"].append(val)
            elif tag in ("HUSB", "WIFE"):
                cur[tag.lower()] = val
        elif level == 2 and sub in ("BIRT", "DEAT") and tag in ("DATE", "PLAC"):
            cur[f"{sub.lower()}_{tag.lower()}"] = val
    return indi, fam


def main():
    if not GED or not os.path.exists(GED):
        print("no GEDCOM (family.json genealogy.gedcom or first argument); pedigree skipped")
        return
    indi, fam = parse(GED)
    out, queue = {}, [(1, ROOT_ID)]
    while queue:
        n, pid = queue.pop(0)
        p = indi.get(pid)
        if not p:
            continue
        out[n] = {k: p.get(k) for k in ("id", "name", "sex", "birt_date", "birt_plac", "deat_date", "deat_plac")}
        for fid in p["famc"][:1]:
            f = fam.get(fid, {})
            if f.get("husb"):
                queue.append((2 * n, f["husb"]))
            if f.get("wife"):
                queue.append((2 * n + 1, f["wife"]))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(out, open(OUT, "w"), indent=1, ensure_ascii=False)
    for n in sorted(out):
        gen = n.bit_length() - 1
        p = out[n]
        print(f"{'  ' * gen}{n:>3} g{gen} {p['name'] or '?'} | b. {p['birt_date'] or '?'}, {p['birt_plac'] or '?'}")


if __name__ == "__main__":
    main()
