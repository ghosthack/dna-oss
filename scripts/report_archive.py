"""23andMe's own displayed results, from archived report pages in raw/<person>/23andme_reports/.

Reads the "Reports Summary" PDFs printed from each account (all sections, incl. Ancestry), the
"All Health & Traits Reports" page archive (section → report → displayed result) and the Sleep Paralysis report (a proprietary 23andMe polygenic score we cannot recompute: its variants and weights
are unpublished). These are 23andMe's labels, kept separate from our own calls so the pages can compare them.

Usage: python3 -I scripts/report_archive.py  → results/23andme_reports.json
"""
import glob
import json
import os
import re
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dnalib import RAW, RESULTS

SECTIONS = ("Health Predisposition", "Carrier Status", "Wellness", "Traits")
PDFTOTEXT = os.environ.get("PDFTOTEXT") or shutil.which("pdftotext") or "/opt/local/bin/pdftotext"
NOT_A_RESULT = ("Upgrade to view report", "Tasks Required", "See report")


def parse_all_reports(path):
    out, section = {}, None
    for line in open(path, encoding="utf-8"):
        m = re.match(r"## (.+)", line)
        if m:
            section = m.group(1).strip() if m.group(1).strip() in SECTIONS else None
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if section and line.startswith("|") and len(cells) >= 2 and cells[0] not in ("Report", "---"):
            out.setdefault(section, {})[cells[0]] = cells[1]
    return out


def parse_summary_pdf(path):
    """23andMe "Reports Summary" PDF → {section: {report: result}}; Ancestry rows keep their indent-free labels."""
    text = subprocess.run([PDFTOTEXT, "-layout", path, "-"], capture_output=True, text=True, check=True).stdout
    out, section, indents = {}, None, []
    for line in text.replace("\u2011", "-").splitlines():  # PDF uses non-breaking hyphens in some names
        m = re.match(r"\s*(Health Predisposition|Carrier Status|Wellness|Ancestry|Traits) Reports\s+\d+\+ reports", line)
        if m:
            section = m.group(1)
            continue
        m = re.match(r"(\s*)(\S.{0,90}?\S)\s{3,}(\S(?:.{0,80}\S)?)\s*$", line)
        if section and m and not m.group(3).startswith(NOT_A_RESULT):
            out.setdefault(section, {})[m.group(2)] = m.group(3)
            if section == "Ancestry" and m.group(3).endswith("%"):
                indents.append((m.group(2), m.group(3), len(m.group(1))))
    levels = sorted({i for _, _, i in indents})
    out["_ancestry_tree"] = [[label, pct, levels.index(i)] for label, pct, i in indents]
    return out


def parse_sleep_paralysis(path):
    text = open(path, encoding="utf-8").read()
    cat = re.search(r"\*\*Result category: ([^.*]+)", text)
    pct = re.search(r"(\d+)% of people with genetics like yours have experienced sleep paralysis by their (\d+)s", text)
    if not cat:
        return None
    return {"category": cat.group(1).strip(), "percent": int(pct.group(1)) if pct else None,
            "by_age": f"{pct.group(2)}s" if pct else None}


def main():
    out = {}
    for d in sorted(glob.glob(os.path.join(RAW, "*", "23andme_reports"))):
        person = os.path.basename(os.path.dirname(d))
        rec, sources = {"reports": {}}, []
        # Oldest first, so newer captures overwrite: page archive (.md), then printed summaries (.pdf), by date.
        files = [(re.search(r"(\d{4}-\d{2}-\d{2})", os.path.basename(f)).group(1), f)
                 for f in glob.glob(os.path.join(d, "23andMe_all_health_traits_reports_*.md"))
                 + glob.glob(os.path.join(d, "23andMe_reports_summary_*.pdf"))]
        for _, f in sorted(files):
            parsed = parse_summary_pdf(f) if f.endswith(".pdf") else parse_all_reports(f)
            tree = parsed.pop("_ancestry_tree", None)
            if tree:
                rec["ancestry_tree"] = tree
            for section, reps in parsed.items():
                rec["reports"].setdefault(section, {}).update(reps)
            sources.append(os.path.basename(f))
        rec["source"] = ", ".join(sources)
        sp = [p for p in glob.glob(os.path.join(d, "*Sleep_Paralysis*.md")) if parse_sleep_paralysis(p)]
        if sp:
            rec["sleep_paralysis"] = parse_sleep_paralysis(max(sp, key=os.path.getmtime))
        elif "Sleep Paralysis" in rec["reports"].get("Wellness", {}):  # summary label only, no percentage
            rec["sleep_paralysis"] = {"category": rec["reports"]["Wellness"]["Sleep Paralysis"], "percent": None,
                                      "by_age": None}
        if not rec["reports"]:
            del rec["reports"], rec["source"]
        if rec:
            out[person] = rec
    json.dump(out, open(os.path.join(RESULTS, "23andme_reports.json"), "w"), ensure_ascii=False, indent=1)
    print(os.path.join(RESULTS, "23andme_reports.json"))


if __name__ == "__main__":
    main()
