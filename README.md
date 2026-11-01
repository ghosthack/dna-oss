# dna-oss: family DNA analysis from consumer chip data

Local, offline-first analysis of 23andMe and AncestryDNA raw data for a family:
you, your full siblings, your parents and your mother's mother. It turns the raw
downloads into a set of static HTML pages: one "Family DNA" page per child, a
"what's different" page per sibling, and a page per person to show a doctor.

**Your genetic data never leaves your machine.** The scripts only download public
reference data: 1000 Genomes slices, the AADR Human Origins panel, PGS Catalog
scoring files, Ensembl lookups by rsID, and Haplogrep/yhaplo trees.

> **For learning and research only. Not medical advice.** See [DISCLAIMER.md](DISCLAIMER.md).
> Consumer chips miss most disease variants and have false positives. A clinician must review
> anything health-related and confirm it with a clinical test. Every generated page carries this notice.

## What it does

| Area | Details |
|---|---|
| QC & relatedness | Chip merge (23andMe + AncestryDNA), Mendelian error checks, KING-style kinship for every pair |
| Family inheritance | Phasing through the trio, maternal crossover maps, grandparent painting (which of mom's half came from grandma vs grandpa), sibling sharing |
| Virtual grandpa | Reconstructs the half of the mother's father's genome that he passed to the mother, from mother + grandmother + children |
| Lineages | mtDNA (Haplogrep 3) and Y (yhaplo) haplogroups, Big Y-700 preview for R1b-P311 men |
| Pharmacogenomics | CPIC genes (CYP2C19, CYP2C9, VKORC1, SLCO1B1, TPMT, NUDT15, DPYD, HLA tags, ...) |
| Health & traits | APOE, carrier screens, hereditary markers (BRCA founders, MUTYH, APC I1307K), wellness, taste, nutrition, pigmentation (HIrisPlex-S) |
| Polygenic scores | PGS Catalog scores on 23andMe's imputed data, ranked against 1000 Genomes, with a mid-parent sanity check |
| Ancestry | 1000 Genomes admixture, fine-scale AADR Human Origins analysis with optional family-tree hypotheses, Neanderthal segments |
| 23andMe comparison | Parses archived 23andMe report pages and shows where our calls agree |

### Scope

Some modules were designed for **Southern European ancestry** and say less for other families:

- Fine-scale ancestry (`ancestry_ho.py`) compares Italian and Iberian reference groups; its verdict rows
  (e.g. Italian vs Spanish) assume that background. The pedigree hypotheses work with any AADR populations.
- The Big Y-700 preview only covers the R1b-P311 branch. Other Y lines get their haplogroup but no preview.
- Neanderthal segments are compared against Sprime calls from a Spanish (IBS) sample.

Contributions that generalize these are welcome.

## Quick start (demo family)

```bash
./fetch_refs.sh                       # third-party tools + reference files into tools/ and .venv/
python3 -I scripts/make_demo.py       # writes a synthetic 5-person family to demo-data/
DNA_DATA=$PWD/demo-data ./run_all.sh --offline
open demo-data/results/index.html
```

The demo family is simulated: random founders, then realistic recombination
and inheritance. It exercises every family-structure step. Reference-panel
steps (ancestry, PGS, Neanderthal) find little overlap with its synthetic
background SNPs.

## Your own family

1. Make a **data folder outside this repo** (e.g. `~/dna-data`) with one folder per person under `raw/`.
   Folder names are your choice; filenames inside don't matter.

   ```
   ~/dna-data/
     family.json
     raw/
       me/    23andMe Raw Data zip (+ Imputed and Phased zips), AncestryDNA zip (optional)
       sis/   ...
       mom/   ...
       dad/   ...
       gma/   ...
   ```

   23andMe: Settings → 23andMe Data. **Raw Data** is required. **Imputed Genotype Data**
   is needed for polygenic scores and fills markers missing from older chips.
   **Phased Genotype Data** is an optional cross-check. AncestryDNA: Settings → DNA → Download DNA data.

2. Write `family.json` (full reference: [docs/family-json.md](docs/family-json.md)). Only `people` is required:

   ```json
   {"people": {
     "me":  {"role": "self", "name": "Alex", "sex": "M"},
     "sis": {"role": "sibling", "name": "Maya", "sex": "F"},
     "mom": {"role": "mother", "sex": "F"},
     "dad": {"role": "father", "sex": "M"},
     "gma": {"role": "maternal_grandmother", "sex": "F"}
   }}
   ```

   Supported roles: `self` plus any of `sibling` (0..n), `mother`, `father` and
   `maternal_grandmother`. Steps that need someone who is missing are skipped.

3. Run:

   ```bash
   DNA_DATA=~/dna-data ./run_all.sh            # everything (fine-scale ancestry takes ~30 min)
   DNA_DATA=~/dna-data ./run_all.sh --quick    # skip polygenic scores and fine-scale ancestry
   ```

   Pages go to `~/dna-data/results/` (start at `index.html`). Intermediate files go to `~/dna-data/work/`.

Optional extras in the data folder:
- `raw/family_notes.json`: observed eye/hair/skin colours, relatives' traits and family medical history ([schema](docs/family-json.md#rawfamily_notesjson)).
- `raw/<person>/23andme_reports/`: saved 23andMe report pages, for the "23andMe vs ours" comparison.
- A GEDCOM file referenced from `family.json`, for family-tree ancestry hypotheses.

## Requirements

macOS or Linux with `python3` (3.12+), `bcftools`, `java` (11+), `curl`, `git` and `unzip`.
`./fetch_refs.sh --check` lists what is missing, with install hints.

## Privacy

- Keep the data folder out of version control. This repo's `.gitignore` excludes `raw/`, `work/`, `results/`, `family.json` and `demo-data/` in case you put data here anyway.
- `scripts/leak_check.py` scans everything git would commit. It looks for names and places from your `family.json`, personal paths, download filenames, and your actual results (diplotypes, haplogroups, rsID+genotype pairs). It never prints the private values. Install it as a pre-commit hook with `scripts/install_hooks.sh`.
- The generated HTML pages contain health information. Treat them like medical records.

## Layout

- `scripts/`: analysis code (`dnalib.py` holds shared loading, paths and the family config)
- `run_all.sh`: runs every step in order
- `fetch_refs.sh`: downloads third-party tools and references
- `docs/`: configuration reference

## License

MIT, see [LICENSE](LICENSE). Third-party tools and data keep their own licenses: [THIRD_PARTY.md](THIRD_PARTY.md).
