# Third-party tools and data

None of these are included in this repository. `fetch_refs.sh` and the scripts download them from their
original sources, and each keeps its own license. Please cite the sources if you publish results.

| Component | Used for | License / terms |
|---|---|---|
| [Haplogrep 3](https://github.com/genepi/haplogrep3) 3.3.2 + PhyloTree 17 (rCRS) | mtDNA haplogroups | MIT; cite Schönherr et al. 2023 and van Oven & Kayser 2009 |
| [yhaplo](https://github.com/23andMe/yhaplo) | Y haplogroups | 23andMe BSD-style license (keep notice; credit 23andMe; no endorsement) |
| [PGS Catalog](https://www.pgscatalog.org) harmonized scoring files | Polygenic scores | EMBL-EBI terms; each score keeps its publication's license (e.g. PGS000657 is CC BY-NC-ND 4.0). Cite the Catalog and each score |
| [1000 Genomes high-coverage (NYGC)](https://www.internationalgenome.org) | Reference frequencies and percentiles | IGSR: free reuse, no restrictions |
| [AADR v66 Human Origins](https://doi.org/10.7910/DVN/FFIDCW) | Fine-scale ancestry | CC0; cite Mallick et al. 2024 |
| [Sprime archaic segments](https://data.mendeley.com/datasets/y7hyt83vxr/1) (Browning et al. 2018) | Neanderthal segments | CC BY 4.0 |
| UCSC hg19ToHg38 liftover chain | Coordinate conversion | Freely available |
| NCBI rCRS (NC_012920.1) | mtDNA reference | Public domain (NCBI) |
| [world-atlas](https://github.com/topojson/world-atlas) 2.0.2 / Natural Earth | Ancestry map | ISC / public domain |
| [Ensembl REST](https://rest.ensembl.org) | rsID positions and alleles | Free data, Apache-2.0 software |
| numpy, pandas and yhaplo's dependencies | Python runtime | BSD/MIT/Apache-style |

The pharmacogenomic, trait and pigmentation rules in `scripts/` are re-implemented from published guidance (CPIC,
HIrisPlex-S and the papers cited in each script). They are not copied from proprietary tools.
