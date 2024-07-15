# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

### TODO 
Resolve naming discrepancy between `somatic_variants_qc.sh` in QC repo and in nextflow pipeline which calls this repo `qc_somatic_variants.sh`

## [0.5.0] Dockerisation/Pipelining
### Added 
- Created elements for Dockerising analysis (tidyverse base image)
- Caputured R dependencies
- Added PERL 3.8.0 to project dockerfile
- Moved `somatic_variants_qc.sh` script into this analysis method to sit with the scripts it calls.
### Changed 
- Modified the font used in `maketileplot_from_maf.R` to use sans rather than Arial - avoid liscencing issues and having to install Arial fonts into the container for the time being. 


## [0.4.0] Publishable Unit 4
### Added
- Add `source_me.sh`

### Changed
- Refactor to remove hard coded paths and to use a source_me.sh

## [0.3.0] Publishable Unit 3
### Changed
- Nothing

## [0.2.0] Publishable Unit 2
### Changed
- Nothing

## [0.1.0] Publishable Unit 1
### Added
- Add `CHANGELOG.md` and `README.md`
- Add initial R script `maf2dndscv.R`, `maf2xlsx.R`, `maketileplot_from_maf.R`,
`plot_vaf_vs_depth_from_maf.R` and `tsv2xlsx.R`.
- Add initial shell script `combine_mafs_and_qc.sh` and `make_variant_release.sh`
- Add initial Perl script `summarise_ascat_estimates.pl`, `reformat_vcf2maf.pl` and `select_vcfsites_from_maf.pl`
