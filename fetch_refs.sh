#!/bin/bash
# Download the third-party tools and public reference files the pipeline expects, and create the
# Python venv. Idempotent: anything already present is skipped, so it is safe to rerun.
# Usage: ./fetch_refs.sh              (everything)
#        ./fetch_refs.sh --check      (only check system prerequisites)
#        ./fetch_refs.sh --no-venv    (skip the Python venv)
#
# Where things go (same variables as run_all.sh):
#   $DNA_TOOLS (default <repo>/tools)  haplogrep3 + mtDNA tree, ref/rCRS.fasta, pgs/*.txt.gz
#   $DNA_VENV  (default <repo>/.venv)  yhaplo, numpy (requirements.txt)
#   $DNA_DATA  (default <repo>)        work/map/countries-50m.json (basemap.py) and
#                                      work/neanderthal/{sprime,hg19ToHg38.over.chain.gz} (neanderthal.py)
# Not fetched here (the scripts download and cache these themselves on first run):
#   - 1000 Genomes high-coverage GRCh38 panel: remote bcftools range queries by pgs.py, ancestry.py
#     and neanderthal.py, cached in work/pgs/ref, work/ancestry/ref, work/neanderthal/ref (needs a
#     bcftools built with libcurl, as Homebrew's and most distro packages are).
#   - AADR v66.p1 Human Origins (Harvard Dataverse): ancestry_ho_fetch.py fetches only the needed
#     records by HTTP range request (~330 MB) into work/ancestry_ho/.
#   - Ensembl REST lookups (scripts/ensembl.py) are live queries, not cached files.
# Nothing about the family is sent anywhere: every request here is a plain download of a public file.
set -euo pipefail
CODE="$(cd "$(dirname "$0")" && pwd)"
TOOLS="${DNA_TOOLS:-$CODE/tools}"
VENV="${DNA_VENV:-$CODE/.venv}"
DATA="${DNA_DATA:-$CODE}"
PYTHON="${PYTHON:-python3}"

check_only=false; no_venv=false
for arg in "$@"; do
  case "$arg" in
    --check) check_only=true ;;
    --no-venv) no_venv=true ;;
    *) echo "unknown option: $arg" >&2; exit 2 ;;
  esac
done

# ---------- pinned versions ----------
HAPLOGREP_VERSION=3.3.2
HAPLOGREP_URL="https://github.com/genepi/haplogrep3/releases/download/v$HAPLOGREP_VERSION/haplogrep3-$HAPLOGREP_VERSION-linux.tar.gz"
HAPLOGREP_SHA256=957165143ddd7d6ad9f34cc3af56b784fc236081e031c99bea4a574d6ad70959  # the linux build is plain Java; works on macOS
MT_TREE=phylotree-rcrs@17.3            # must match --tree in run_all.sh
RCRS_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012920.1&rettype=fasta&retmode=text"
WORLD_ATLAS_URL="https://cdn.jsdelivr.net/npm/world-atlas@2.0.2/countries-50m.json"
CHAIN_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
PGS_URL="https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{id}/ScoringFiles/Harmonized/{id}_hmPOS_GRCh38.txt.gz"
# Polygenic scores: pgs.py scores every tools/pgs/*_hmPOS_GRCh38.txt.gz it finds, so this list is the
# score set. Override with PGS_IDS="PGS000001 PGS000002 ..." (add or drop scores freely). Optional extras the
# pages know how to describe: PGS000657 (heel bone density), PGS001811 (colon polyps), PGS001024 (hemorrhoids).
PGS_IDS="${PGS_IDS:-PGS000004 PGS000035 PGS000137 PGS000746 PGS000885
  PGS002247 PGS002257 PGS002728 PGS002729 PGS003784 PGS003850 PGS003998 PGS004115 PGS004119
  PGS004125 PGS004326 PGS004609}"
# Sprime introgressed-segment calls (Browning et al. 2018, Mendeley Data y7hyt83vxr v1), one tarball per
# 1000 Genomes population. neanderthal.py reads the population set by its POP constant (IBS).
SPRIME_POPS="${SPRIME_POPS:-IBS}"
MENDELEY="https://data.mendeley.com/public-files/datasets/y7hyt83vxr/files"
sprime_file() {  # population -> "<mendeley file id> <sha256>"
  case "$1" in
    BEB) echo "72c6314c-10ad-4626-8ad5-d4d846aa65bf 2c04a24eac05346bc4ac271ece9451035237d9f8df3abe5f641f0f99b8c5d415" ;;
    CDX) echo "f01566ae-0a9b-4847-b312-1af059fae3d1 ac7336a0cb0461a979db6d0a381957c39621112f4298a55f97ba47a04666e923" ;;
    CEU) echo "46dcb65f-47dc-4cd4-84e6-f0ffb5a189ef f695e484db15d1974d0fb6bd15a3e8d769a872f2a3803908d9197a057b9472f0" ;;
    CHB) echo "6625276a-bcc8-4ac3-be58-e7b15c786bae 345885d49379477ec48bb470bb8489bb339f24898db09676fcc10840e7337c05" ;;
    CHS) echo "0b77dfa9-dc7b-44b7-b7db-f2dde8a87104 c43b76dee120c1b60a9982a058e9c9ca26136c1bc7a8d349390421820133ccb4" ;;
    CLM) echo "03547a83-a35a-4827-a94d-c106361922c3 958f9223eace3506b581371dd4a5259f93a4fc4cafb73751a4859530b251244a" ;;
    FIN) echo "3ec4a2b4-9e65-40f6-89b1-6da39bdf7754 447e0234213167cd50718b02354bd711d5bbd8ee74e4228d1a700bee11caaa7c" ;;
    GBR) echo "a225099b-8404-498e-a8a8-53dc19cd409f e3492663fbf8390c7fb6aaaea43ecc686314c65b20265506050f403cba26b8a3" ;;
    GIH) echo "f129f344-b21c-4e39-bd54-1cf66ca26292 e1c26d457bf04cab9e548b28446b7733d7ad7a61e690761ccc3ee39b1e701e48" ;;
    IBS) echo "d93b3e25-142f-485e-9761-de5c7a43625f cc0907d2864ab6a67f632efa0d9317e2aa9fed1e612d98bda87fcceab9348604" ;;
    ITU) echo "117b0f85-6cca-44bd-8bf5-f10a14e8a242 d87ed8e0690dae34f5323face6b480fe3e2824ebdba9d0b0d3c3fb7542fbab3c" ;;
    JPT) echo "f6c6ce7a-77fa-4c2a-bd82-3d10996932b7 ef92f06f7e576d6b5c5777be9192542c00a736675e374a54da8d03941f89f376" ;;
    KHV) echo "964a196f-2472-4a30-a6be-9c31cc120710 aa0f04475e2cba5bf54d5571b2b615d64b199019319c70f2d5220f98efb20038" ;;
    MXL) echo "cc28e670-c93f-44b1-ad46-f149957b67da 81ab87c396c1e5db8e145e2104e985b22fc9b7be6fd531f20170efb2e7788b59" ;;
    Papuans) echo "7af62f7d-b3d8-4ffa-b982-e4680c1370eb 850402fffdf6a2805e9345664406bcf222e27ed9da3575a6babeec7248137820" ;;
    PEL) echo "d92b175b-2ea8-4553-87b1-d52fa7742a72 458476a50aada2f85975834908d2967a0ea1ddaf944605fafa3b87c1c08ec030" ;;
    PJL) echo "55746b93-b90e-48d6-94fc-9cdcd281d287 578a7abbf2fed26c54c10fa9b8f840d8e08df714a5b6e1018a8ae8d4967e8f52" ;;
    PUR) echo "4b677041-353a-459c-8537-a6a468f7b15e 9c6613d11f60d91b0dda1a17c0018aa941f5c1ec5a7186191d160ac1a0da894a" ;;
    STU) echo "1dff818f-ffbd-446d-8b80-15983ba430e5 642602f5708b9d41e42da4c288ec26fbd082826f8373010e3cb8b02635de55bd" ;;
    TSI) echo "67163b17-220c-4e36-8184-2f085f30e73b 6596f7fde651092c4e2d43e5917e3c17b50ab519b49caea59136214a4f1a5723" ;;
    *) return 1 ;;
  esac
}

step() { printf '\n== %s\n' "$*"; }
skip() { printf '   have %s\n' "${1#"$CODE/"}"; }

# ---------- prerequisites (never installed here; hints only) ----------
hint() {
  case "$1" in
    curl)     echo "macOS: brew install curl          Debian/Ubuntu: sudo apt install curl" ;;
    python3)  echo "macOS: brew install python@3.12   Debian/Ubuntu: sudo apt install python3 python3-venv (3.12 or newer)" ;;
    git)      echo "macOS: xcode-select --install     Debian/Ubuntu: sudo apt install git" ;;
    java)     echo "macOS: brew install openjdk       Debian/Ubuntu: sudo apt install default-jre (Java 11 or newer)" ;;
    bcftools) echo "macOS: brew install bcftools      Debian/Ubuntu: sudo apt install bcftools" ;;
    unzip)    echo "macOS: preinstalled               Debian/Ubuntu: sudo apt install unzip" ;;
    *)        echo "install $1 with your package manager" ;;
  esac
}
missing_required=(); missing_later=()
for t in curl tar gzip "$PYTHON"; do command -v "$t" >/dev/null || missing_required+=("$t"); done
$no_venv || command -v git >/dev/null || missing_required+=(git)
# Needed by run_all.sh, not by this script (java is also needed here to install the mtDNA tree).
java_ok=false; command -v java >/dev/null && java -version >/dev/null 2>&1 && java_ok=true  # macOS has a /usr/bin/java stub
$java_ok || missing_later+=(java)
BCFTOOLS="${BCFTOOLS:-$(command -v bcftools || true)}"
[[ -n "$BCFTOOLS" ]] || missing_later+=(bcftools)
command -v unzip >/dev/null || missing_later+=(unzip)
if command -v "$PYTHON" >/dev/null && ! "$PYTHON" -c 'import sys; sys.exit(sys.version_info < (3, 12))'; then
  echo "warning: $PYTHON is $("$PYTHON" -V 2>&1); the venv (numpy 2.5) needs 3.12+. Set PYTHON=/path/to/python3.12" >&2
  $no_venv || missing_required+=("python3.12+")
fi
for t in ${missing_required[@]+"${missing_required[@]}"} ${missing_later[@]+"${missing_later[@]}"}; do
  printf 'missing: %-10s %s\n' "$t" "$(hint "${t%%.*}")" >&2
done
if $check_only; then
  [[ ${#missing_required[@]} -eq 0 && ${#missing_later[@]} -eq 0 ]] && echo "All prerequisites found."
  exit $(( ${#missing_required[@]} + ${#missing_later[@]} > 0 ))
fi
[[ ${#missing_required[@]} -eq 0 ]] || { echo "Install the missing tools above and rerun." >&2; exit 1; }

sha256() { if command -v sha256sum >/dev/null; then sha256sum "$1"; else shasum -a 256 "$1"; fi | cut -d' ' -f1; }

# Pinned checksums for files whose URL has no version in it (refs.sha256: "<sha256>  <key>" per line).
PINS="$CODE/refs.sha256"
pinned() { [[ -f "$PINS" ]] && awk -v k="$1" '$2 == k {print $1; exit}' "$PINS"; return 0; }

# fetch URL DEST [SHA256 | pin:KEY]: download to DEST.part, verify, then move into place. Skips an existing DEST.
# "pin:KEY" looks the checksum up in refs.sha256; a file with no pin yet is accepted and its checksum printed.
fetch() {
  local url="$1" dest="$2" want="${3:-}" key=""
  if [[ "$want" == pin:* ]]; then key="${want#pin:}"; want="$(pinned "$key")"; fi
  if [[ -s "$dest" ]]; then skip "$dest"; return 0; fi
  mkdir -p "$(dirname "$dest")"
  printf '   get  %s\n' "${dest#"$CODE/"}"
  curl -fSL --retry 3 --retry-delay 5 --connect-timeout 30 -s -o "$dest.part" "$url"
  local got; got="$(sha256 "$dest.part")"
  if [[ -n "$want" && "$got" != "$want" ]]; then
    rm -f "$dest.part"
    echo "checksum mismatch for $url (got $got, expected $want)." >&2
    [[ -n "$key" ]] && echo "   If the upstream file was deliberately updated, update '$key' in refs.sha256." >&2
    return 1
  fi
  [[ -n "$key" && -z "$want" ]] && printf '   unpinned: %s  %s\n' "$got" "$key"
  mv "$dest.part" "$dest"
}

failed=()

# ---------- haplogrep3 (mtDNA haplogroups) ----------
step "haplogrep3 $HAPLOGREP_VERSION + $MT_TREE"
HG="$TOOLS/haplogrep3"
if [[ -f "$HG/haplogrep3.jar" ]]; then
  skip "$HG/haplogrep3.jar"
else
  mkdir -p "$TOOLS"; tmp="$(mktemp -d "$TOOLS/.haplogrep3.XXXXXX")"
  fetch "$HAPLOGREP_URL" "$tmp/haplogrep3.tar.gz" "$HAPLOGREP_SHA256"
  mkdir -p "$tmp/x" && tar -xzf "$tmp/haplogrep3.tar.gz" -C "$tmp/x"
  jar="$(find "$tmp/x" -name haplogrep3.jar | head -1)"
  [[ -n "$jar" ]] || { echo "haplogrep3.jar not found in $HAPLOGREP_URL" >&2; exit 1; }
  mkdir -p "$HG" && cp -R "$(dirname "$jar")"/. "$HG"/
  rm -rf "$tmp"
fi
tree_name="${MT_TREE%@*}"; tree_ver="${MT_TREE#*@}"
if [[ -f "$HG/trees/$tree_name/$tree_ver/tree.xml" ]]; then
  skip "$HG/trees/$tree_name/$tree_ver"
elif $java_ok; then
  # haplogrep3 resolves the tree from its repository index (haplogrep3.yaml) and unpacks it into trees/.
  (cd "$HG" && java -jar haplogrep3.jar install-tree "$MT_TREE") || failed+=("$MT_TREE")
else
  echo "   java not found: rerun after installing Java to install $MT_TREE" >&2; failed+=("$MT_TREE")
fi

# ---------- revised Cambridge Reference Sequence (mt_vcf.py) ----------
step "rCRS (NC_012920.1)"
if fetch "$RCRS_URL" "$TOOLS/ref/rCRS.fasta" pin:rCRS.fasta; then
  len="$(grep -v '^>' "$TOOLS/ref/rCRS.fasta" | tr -d '\n\r ' | wc -c | tr -d ' ')"
  [[ "$len" == 16569 ]] || { echo "   rCRS.fasta has $len bases, expected 16569; removed, rerun to retry" >&2
                             rm -f "$TOOLS/ref/rCRS.fasta"; failed+=(rCRS); }
else failed+=(rCRS); fi

# ---------- PGS Catalog harmonized scoring files (pgs.py) ----------
step "PGS Catalog scoring files (GRCh38 harmonized)"
for id in $PGS_IDS; do
  dest="$TOOLS/pgs/${id}_hmPOS_GRCh38.txt.gz"
  fetch "${PGS_URL//\{id\}/$id}" "$dest" "pin:${id}_hmPOS_GRCh38.txt.gz" || { failed+=("$id"); continue; }
  gzip -t "$dest" 2>/dev/null || { echo "   $id is corrupt; removed, rerun to retry" >&2; rm -f "$dest"; failed+=("$id"); }
done

# ---------- basemap source (basemap.py) ----------
step "world-atlas 2.0.2 countries-50m (Natural Earth)"
fetch "$WORLD_ATLAS_URL" "$DATA/work/map/countries-50m.json" pin:countries-50m.json || failed+=(world-atlas)

# ---------- Neanderthal references (neanderthal.py) ----------
step "Sprime segments ($SPRIME_POPS) + UCSC hg19->hg38 chain"
SP="$DATA/work/neanderthal/sprime"
for pop in $SPRIME_POPS; do
  read -r fid want < <(sprime_file "$pop" || echo "- -")
  if [[ "$fid" == "-" ]]; then echo "   unknown Sprime population: $pop" >&2; failed+=("sprime-$pop"); continue; fi
  if compgen -G "$SP/mendeley_data/$pop.chr*.ND_match" >/dev/null; then skip "$SP/mendeley_data/$pop.chr*"; continue; fi
  fetch "$MENDELEY/$fid/file_downloaded" "$SP/${pop}_sprime_results.tar.gz" "$want" || { failed+=("sprime-$pop"); continue; }
  tar -xzf "$SP/${pop}_sprime_results.tar.gz" -C "$SP"
  compgen -G "$SP/mendeley_data/$pop.chr*.ND_match" >/dev/null \
    || { echo "   $pop tarball did not unpack to mendeley_data/$pop.chr*.ND_match" >&2; failed+=("sprime-$pop"); }
done
fetch "$CHAIN_URL" "$DATA/work/neanderthal/hg19ToHg38.over.chain.gz" pin:hg19ToHg38.over.chain.gz || failed+=(liftover-chain)

# ---------- Python venv (yhaplo for the Y step, numpy for ancestry_ho.py) ----------
if ! $no_venv; then
  step "Python venv ($VENV)"
  req="$CODE/requirements.txt"; stamp="$VENV/.requirements.sha256"
  if [[ -x "$VENV/bin/python" && -f "$stamp" && "$(cat "$stamp")" == "$(sha256 "$req")" ]]; then
    skip "$VENV (requirements.txt unchanged)"
  else
    [[ -x "$VENV/bin/python" ]] || "$PYTHON" -m venv "$VENV"
    "$VENV/bin/python" -m pip install -q --upgrade pip
    if "$VENV/bin/python" -m pip install -q -r "$req"; then sha256 "$req" > "$stamp"; else failed+=(venv); fi
  fi
fi

if [[ ${#failed[@]} -gt 0 ]]; then
  printf '\nIncomplete: %s. Rerun ./fetch_refs.sh to retry (finished downloads are kept).\n' "${failed[*]}" >&2
  exit 1
fi
[[ ${#missing_later[@]} -eq 0 ]] || printf '\nReferences ready, but run_all.sh still needs: %s (see hints above).\n' "${missing_later[*]}"
echo "Done. Next: ./run_all.sh"
