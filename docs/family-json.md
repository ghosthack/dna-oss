# Configuration reference

Two optional-but-useful files live in your **data folder** (`$DNA_DATA`), never in this repo.

## `family.json`

Only `people` is required. Every other key is optional; when one is missing, the section that uses it is hidden or
gets neutral wording.

### `people` (required)

```json
"people": {
  "<folder id>": {"role": "self | sibling | mother | father | maternal_grandmother",
                  "name": "First name (optional, used on pages and in titles)",
                  "sex": "M | F (optional; inferred from chrY calls when absent)"}
}
```

- The key is the person's folder name under `raw/`.
- Exactly one `self`. Any number of `sibling`s (full siblings of self). At most one each of the other roles.
- Without a `family.json`, the folder names `me`, `brother`, `sister`, `mom`, `dad` and `grandma-maternal` are recognized.

### Family structure extras

| Key | Meaning |
|---|---|
| `maternal_grandfather.name` | Name shown for the virtual grandfather reconstruction |
| `maternal_grandfather.red_haired_daughter` | `true` if the mother has a red-haired sister; explains MC1R alleles traced to him |
| `pages` | `{"<child id>": "file.html"}` page filenames. Default: self → `index.html`, each sibling → `<id>.html` |
| `genealogy.gedcom`, `genealogy.self_id` | GEDCOM file (absolute, or relative to the data folder) and self's xref (e.g. `@I1@`; default: `@I0@`, leak-ok), for `scripts/pedigree.py` |

### `ancestry` (family-tree ancestry, optional)

Founders are self's 16 great-great-grandparents, numbered 16–31 as in an Ahnentafel
(16–23 father's side, 24–31 mother's side). Without these keys the generic ancestry
analyses still run and the family-tree parts are skipped.

```json
"ancestry": {
  "regions": ["Region A", "Region B"],
  "founders": {"16": {"region": "Region A", "label": "Surname (place)", "assumed": false, "ref_pop": "<AADR population>"}},
  "places": [{"lat": 0.0, "lon": 0.0, "place": "Town", "surnames": "Surname", "region": "Region A",
              "dx": 9, "dy": -4, "anchor": "start"}],
  "map": {"bbox": [lon_min, lat_min, lon_max, lat_max], "countries": ["Country"], "px_per_degree": 22},
  "hypotheses": {"<name>": {"about": "What this alternative tree assumes", "founders": {"23": "<AADR population>"}}},
  "proxy_notes": ["Caveats about the proxy populations"],
  "pgs_local_pop": "GBR"
}
```

- `regions`, `founders.*.region/label/assumed`, `places` and `map` drive the family-tree bars and map on the landing page.
- `founders.*.ref_pop`, `hypotheses` and `proxy_notes` drive the pedigree simulations in `scripts/ancestry_ho.py`. Population names must exist in the AADR Human Origins panel or in the subgroups defined in that script.
- `pgs_local_pop` adds a 1000 Genomes population as an extra comparison group for polygenic scores.

### `text` (optional wording)

Free-text snippets that personalize the pages: family history context, Big Y guesses,
a surname for the Y line, notes on ancestry results. Keys are `<script>.<what>`.
Some keys also take a per-person form, `<script>.<what>.<person id>`, or a
`{"<viewer id>": "..."}` dict where the wording depends on whose page it is.
Pages read fine without any of these. The keys the scripts read are listed with
`grep -o 'text(f\?"[^"]*' scripts/*.py`.

## `raw/family_notes.json`

Facts the chips can't measure. Everything is optional.

```json
{
  "observed":  {"<id>": {"eye": "brown", "hair": "dark brown", "hair_texture": "wavy", "skin": "olive"}},
  "relatives": {"aunt_maternal": "red hair"},
  "medical_history": {"<id>": {"high_cholesterol": "...", "type_2_diabetes": "...", "glaucoma": "...", "prostate_cancer": "..."}},
  "genealogy": {"source": "...", "maternal": "...", "paternal": "..."}
}
```

`observed` is compared with the pigmentation predictions. `medical_history` adds
family-history context to the doctor pages, e.g. hereditary cancer guidance when a
relative has a `*_cancer` entry. Keys are free-form snake_case condition names;
common ones get a display label, others are shown with underscores as spaces. Write free-text values from
self's point of view ("lighter than mine"); sibling pages rewrite "mine" automatically.
