"""Developer helper (not part of the pipeline): look up which plus-strand allele causes which amino-acid change
(Ensembl VEP), with retry/backoff. Used when writing the allele tables in scripts/.

Usage: python3 -I scripts/ensembl.py rs1800562 rs1799945 ...
Prints, per variant: gene, plus-strand alt allele, amino-acid change and position, gene strand;
for non-coding variants: alleles, minor allele and MAF.
"""
import json
import sys
import time
import urllib.error
import urllib.request

API = "https://rest.ensembl.org"


def get(path):
    for attempt in range(8):
        try:
            req = urllib.request.Request(API + path, headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code != 429 and e.code < 500:  # 400/404: a bad rsID; retrying won't help
                raise LookupError(f"Ensembl {e.code} for {path}") from e
            wait = e.headers.get("Retry-After")  # Ensembl allows ~15 requests/s and says how long to back off
            time.sleep(float(wait) if wait and wait.replace(".", "", 1).isdigit() else 1.5 * (attempt + 1))
        except (urllib.error.URLError, json.JSONDecodeError, TimeoutError, OSError):
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"Ensembl request failed: {path}")


def main():
    for rs in sys.argv[1:]:
        try:
            res = get(f"/vep/human/id/{rs}?canonical=1")
        except LookupError as e:
            print(f"{rs} | not found ({e})", flush=True)
            continue
        if not res:
            print(f"{rs} | no VEP result", flush=True)
            continue
        vep = res[0]
        changes = sorted({(t["gene_symbol"], t["variant_allele"], t["amino_acids"], t.get("protein_start"), t["strand"])
                          for t in vep.get("transcript_consequences", [])
                          if t.get("canonical") == 1 and "amino_acids" in t})
        if changes:
            desc = "; ".join(f"{g} alt(+)={a} {aa.split('/')[0]}{pos}{aa.split('/')[-1]} gene_strand={s}"
                             for g, a, aa, pos, s in changes)
        else:
            var = get(f"/variation/human/{rs}")
            desc = f"non-coding; minor={var.get('minor_allele')} MAF={var.get('MAF')}"
        print(f"{rs} {vep['allele_string']} | {desc}", flush=True)
        time.sleep(0.3)


if __name__ == "__main__":
    main()
