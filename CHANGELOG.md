# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [0.6.3] - 2025-02-20
### Added 
- Externalised with license, readme.

## [0.6.2] - 2025-02-13
### Fixed
- Fixed `reformat_vcf2maf.pl` but that only affects gnomAD filtering of keep/keepPA variants for genes in --transcript 
- See tag 0.5.2

## [0.6.1] - 2024-10-31
### Fixed
- Fixed `maketileplot_from_maf.R` error with edge cases where all genes and all samples are mutated

## [0.6.0] - 2024-10-31
### Added
- Add functionality to `plot_vaf_vs_depth_from_maf.R` to allow plotting of AF vs Depth from germline `VAF_norm` and `n_depth`
- Add functionality to `maketileplot_from_maf.R` to use a genelist with format `Hugo_Symbol/Gene` for non-human data sets

## [0.5.5] - 2024-10-23
### Added
- Add `--canonical` and `--exclude_noncoding` options to use with `--pass` option in `reformat_vcf2maf.pl`
### Fixed
- Fixed bug in `plot_vaf_vs_depth_from_maf.pl` that occurred when MAF has samples that are not in --samplelist

## [0.5.4] - 2024-10-01
### Fixed
- Fixed `reformat_vcf2maf.pl` filtering of genes/transcripts in the --transcripts list, if used.
- Previous incorrect behaviour: keep/voi status incorrect due to not checking for PASS, protein-coding
- Corrected behaviour: all transcript in --transcripts list are checked as PASS and protein-coding

## [0.5.3] - 2024-09-24
### Fixed
- Fixed `reformat_vcf2maf.pl` GATK germline VAF- and DP-based filtering. Does not affect the caveman/pindel filtering.

## [0.5.2] - 2024-08-16
### Added
- Added an option to `reformat_vcf2maf.pl` to use a list of alternative transcripts for variant consequence annotation.

## [0.5.1] Housekeeping
### Changed
- Ran styling on all R files. Kept deliberately seperate from 0.5.0 changes

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
