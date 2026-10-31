#!/bin/bash
# Rerun every analysis and rebuild results/index.html. Safe to rerun; downloads are cached.
# Usage: ./run_all.sh            (all steps)
#        ./run_all.sh --no-pgs   (skip polygenic scores)
#        ./run_all.sh --quick    (skip polygenic scores and the ~30-min fine-scale ancestry)
#        ./run_all.sh --offline  (also skip every step that downloads reference panels: ancestry, Neanderthal; for the demo)
set -euo pipefail
CODE="$(cd "$(dirname "$0")" && pwd)"
export DNA_DATA="${DNA_DATA:-$CODE}"            # folder with raw/ and family.json; work/ and results/ are written there
export DNA_TOOLS="${DNA_TOOLS:-$CODE/tools}"     # filled by ./fetch_refs.sh
VENV="${DNA_VENV:-$CODE/.venv}"
cd "$DNA_DATA"

PY="python3 -I"
S="$CODE/scripts"
BCFTOOLS="${BCFTOOLS:-$(command -v bcftools || echo /opt/homebrew/bin/bcftools)}"
mkdir -p work results
skip_pgs=false; skip_fine=false; offline=false
[[ "${1:-}" == "--no-pgs" ]] && skip_pgs=true
[[ "${1:-}" == "--quick" ]] && skip_pgs=true && skip_fine=true
[[ "${1:-}" == "--offline" ]] && skip_pgs=true && skip_fine=true && offline=true

step() { printf '\n== %s\n' "$*"; }

# Unpack and index any imputed downloads not yet extracted (e.g. a newly added relative's).
for zip in raw/*/imputed_genotype_data_*.zip; do
  [[ -e "$zip" ]] || continue
  person=$(basename "$(dirname "$zip")")
  if ! ls work/imputed/"$person"/*/chr1.bcf >/dev/null 2>&1; then
    step "Unpacking imputed data for $person"
    mkdir -p work/imputed/"$person"
    unzip -q -o "$zip" -d work/imputed/"$person"
    for f in work/imputed/"$person"/*/*.bcf; do "$BCFTOOLS" index -f "$f"; done
  fi
done

step "QC, relatedness, grandparent painting";   $PY "$S/analyze.py"
step "Pharmacogenomics";                         $PY "$S/pgx.py"
step "APOE";                                     $PY "$S/apoe.py"
step "Wellness, taste, nutrition, conditions";   $PY "$S/wellness.py"

step "mtDNA haplogroups"
$PY "$S/mt_vcf.py"
if [[ -s work/mt.vcf ]]; then
  (cd "$DNA_TOOLS/haplogrep3" && java -jar haplogrep3.jar classify --in "$DNA_DATA/work/mt.vcf" \
     --out "$DNA_DATA/work/mt_haplogroups.txt" --tree phylotree-rcrs@17.3 --chip --extend-report --hits 3 >/dev/null)
  awk -F'\t' 'NR==1 || $3=="\"1\""' work/mt_haplogroups.txt | cut -f1,2,4,6,7 | tr -d '"' > results/mt_haplogroups.tsv
fi

step "Y haplogroups"
$PY "$S/y_genos.py"
printf "person\thaplogroup\tisogg_2016\n" > results/y_haplogroups.tsv
if [[ -s work/y.genos.txt && $(wc -l < work/y.genos.txt) -gt 1 ]]; then  # header + one row per man
  "$VENV/bin/yhaplo" -i work/y.genos.txt -o work/yhaplo -hp -hpd >/dev/null
  awk '{printf "%s\t%s\t%s\n",$1,$2,$4}' work/yhaplo/haplogroups.y.txt >> results/y_haplogroups.tsv
fi
$PY "$S/bigy_preview.py" >/dev/null

if ! $skip_pgs; then
  step "Polygenic scores (may fetch new 1000 Genomes positions)"; $PY "$S/pgs.py"
fi

# Optional analyses: run if present.
for s in grandpa grandpa_cells pigmentation ancestry neanderthal hereditary_markers report_archive sibling_diffs doctor; do
  if $offline && [[ $s == ancestry || $s == neanderthal ]]; then continue; fi
  if [[ -f "$S/$s.py" ]]; then step "$s"; $PY "$S/$s.py"; fi
done

if ! $skip_fine && [[ -f "$S/ancestry_ho.py" ]]; then
  step "Fine-scale ancestry (AADR panel, ~30 min)"; "$VENV/bin/python" -I "$S/ancestry_ho.py"
fi

if [[ ! -f work/map/basemap.json && -f work/map/countries-50m.json ]]; then
  step "Ancestry basemap"; $PY "$S/basemap.py"
fi
step "Landing page"; $PY "$S/landing.py"
echo "Done → results/index.html"
