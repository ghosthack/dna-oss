"""Write a multi-sample chrM VCF (vs rCRS) from the chip data, for Haplogrep 3 --chip.

Usage: python3 -I scripts/mt_vcf.py
Writes work/mt.vcf and prints per-person variant counts and pairwise differences.
"""
import itertools
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dnalib import TOOLS, WORK, load_person, people as family_people


def rcrs():
    with open(os.path.join(TOOLS, "ref", "rCRS.fasta")) as fh:
        return "".join(l.strip() for l in fh if not l.startswith(">"))


def mt_calls(data):
    """{position: base}; positions with conflicting probes or indels are dropped."""
    calls, bad = {}, set()
    for ch, pos, gt in data.values():
        if ch != "MT" or len(set(gt)) != 1 or gt[0] not in "ACGT":
            continue
        if calls.get(pos, gt[0]) != gt[0]:
            bad.add(pos)
        calls[pos] = gt[0]
    return {p: b for p, b in calls.items() if p not in bad}


def main():
    ref = rcrs()
    people = {}
    for p in family_people():
        data, _ = load_person(p)
        if data:
            people[p] = mt_calls(data)
    out = os.path.join(WORK, "mt.vcf")
    if not people:
        if os.path.exists(out):
            os.remove(out)  # don't leave samples from an earlier run for Haplogrep
        print("no genotype data: no mt.vcf written")
        return
    positions = sorted(set().union(*people.values()))
    with open(out, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n##contig=<ID=chrM,length=16569>\n"
                 '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(people) + "\n")
        for pos in positions:
            r = ref[pos - 1]
            alts = sorted({c[pos] for c in people.values() if pos in c and c[pos] != r})
            gts = []
            for c in people.values():
                b = c.get(pos)
                gts.append("." if b is None else "0" if b == r else str(alts.index(b) + 1))
            fh.write(f"chrM\t{pos}\t.\t{r}\t{','.join(alts) or '.'}\t.\tPASS\t.\tGT\t" + "\t".join(gts) + "\n")
    for p, c in people.items():
        var = sorted(pos for pos, b in c.items() if b != ref[pos - 1])
        print(f"{p:17s} sites={len(c)} non-rCRS={len(var)}")
    for a, b in itertools.combinations(people, 2):
        shared = people[a].keys() & people[b].keys()
        diff = sum(people[a][p] != people[b][p] for p in shared)
        print(f"  {a} vs {b}: {diff} differences / {len(shared)} shared sites")


if __name__ == "__main__":
    main()
