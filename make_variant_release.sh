#!/bin/bash

# This script converts MAFs to xlsx, creates a README file
# and copies over the relevant file to create a variant release

PROJECTDIR=$1
STUDY=$2
VERSION=$3
RELEASE=$4

# Check scripts directory for maf2xlsx.R script

if [[ -z "$PROJECTDIR" || -z "$STUDY" || -z "$VERSION" || -z "$RELEASE" ]]; then
    echo -e "\nUsage: $0 PROJECTDIR STUDY_PREFIX VERSION RELEASE\n"
    exit 1
elif [[ ! -e $PROJECTDIR/scripts/MAF/maf2xlsx.R ]]; then
    echo "Cannot find required script $PROJECTDIR/scripts/MAF/maf2xlsx.R"
    exit 1
fi


##### Process one tumour per patient #####

for cohort in all independent onePerPatient; do 
    REL_DIR=""
    if [[ -d "$PROJECTDIR/analysis/variants_combined/$VERSION/$cohort" ]]; then 
        if [[ $cohort == "all" ]]; then
            REL_DIR="all_tumours"
            NAME="allTum"
        elif [[ $cohort == "independent" ]]; then
            REL_DIR="independent_tumours"
            NAME="indepTum"
        elif [[ $cohort == "onePerPatient" ]]; then
            REL_DIR="one_tumour_per_patient"
            NAME="oneTum"
        fi
        echo "Found $PROJECTDIR/analysis/variants_combined/$VERSION/$cohort"
    
        mkdir -p  $RELEASE/$REL_DIR/matched_samples/QC_keepPA
        mkdir -p  $RELEASE/$REL_DIR/all_samples/QC_keepPA
        mkdir -p  $RELEASE/$REL_DIR/matched_samples/QC_keep
        mkdir -p  $RELEASE/$REL_DIR/all_samples/QC_keep
        
        # keep-PA, filtered, matched tumours only
        cp -v $PROJECTDIR/analysis/variants_combined/$VERSION/$cohort/keepPA_vaf_size_filt_matched_caveman_pindel_${cohort}.maf $RELEASE/$REL_DIR/matched_samples/$STUDY-filtered_mutations_matched_${NAME}_keepPA.maf
        
            rsync -a --exclude '*maf' --exclude 'Rplots.pdf' $PROJECTDIR/analysis/variants_combined/$VERSION/$cohort/plots_keepPA_vaf_size_filt_matched/* $RELEASE/$REL_DIR/matched_samples/QC_keepPA/
        
        Rscript $PROJECTDIR/scripts/MAF/maf2xlsx.R $RELEASE/$REL_DIR/matched_samples/$STUDY-filtered_mutations_matched_${NAME}_keepPA.maf
        
        cp -v $PROJECTDIR/analysis/variants_combined/$VERSION/$cohort/sample_list_matched.tsv $RELEASE/$REL_DIR/matched_samples/
        
        # keep, filtered, matched tumours only
        cp -v $PROJECTDIR/analysis/variants_combined/$VERSION/$cohort/keep_vaf_size_filt_matched_caveman_pindel_${cohort}.maf $RELEASE/$REL_DIR/matched_samples/$STUDY-filtered_mutations_matched_${NAME}_keep.maf
        
        cp -v $PROJECTDIR/analysis/variants_combined/$VERSION/$cohort/plots_keep_vaf_size_filt_matched/mutations_per_Mb.tsv $RELEASE/$REL_DIR/matched_samples/
        
        rsync -a --exclude '*maf' --exclude 'Rplots.pdf' --exclude mutations_per_Mb.tsv $PROJECTDIR/analysis/variants_combined/$VERSION/$cohort/plots_keep_vaf_size_filt_matched/* $RELEASE/$REL_DIR/matched_samples/QC_keep/
        
        Rscript $PROJECTDIR/scripts/MAF/maf2xlsx.R $RELEASE/$REL_DIR/matched_samples/$STUDY-filtered_mutations_matched_${NAME}_keep.maf
        
        cp -v $PROJECTDIR/analysis/variants_combined/$VERSION/$cohort/sample_list_matched.tsv $RELEASE/$REL_DIR/matched_samples/
        
         
        # keep-PA, filtered, matched and unmatched tumours
        cp -v $PROJECTDIR/analysis/variants_combined/$VERSION/$cohort/keepPA_vaf_size_filt_caveman_pindel_${cohort}.maf $RELEASE/$REL_DIR/all_samples/$STUDY-filtered_mutations_all_${NAME}_keepPA.maf
         
        rsync -a --exclude '*maf' --exclude 'Rplots.pdf' $PROJECTDIR/analysis/variants_combined/$VERSION/$cohort/plots_keepPA_vaf_size_filt/* $RELEASE/$REL_DIR/all_samples/QC_keepPA
        
        Rscript $PROJECTDIR/scripts/MAF/maf2xlsx.R $RELEASE/$REL_DIR/all_samples/$STUDY-filtered_mutations_all_${NAME}_keepPA.maf
        
        cp -v $PROJECTDIR/analysis/variants_combined/$VERSION/$cohort/sample_list.tsv $RELEASE/$REL_DIR/all_samples/
         
        # keep, filtered, matched and unmatched tumours
        cp -v $PROJECTDIR/analysis/variants_combined/$VERSION/$cohort/keep_vaf_size_filt_caveman_pindel_${cohort}.maf $RELEASE/$REL_DIR/all_samples/$STUDY-filtered_mutations_all_${NAME}_keep.maf
         
        cp -v $PROJECTDIR/analysis/variants_combined/$VERSION/$cohort/plots_keep_vaf_size_filt/mutations_per_Mb.tsv $RELEASE/$REL_DIR/all_samples/
        
        rsync -a --exclude '*maf' --exclude 'Rplots.pdf' --exclude mutations_per_Mb.tsv $PROJECTDIR/analysis/variants_combined/$VERSION/$cohort/plots_keep_vaf_size_filt/* $RELEASE/$REL_DIR/all_samples/QC_keep/
        
        Rscript $PROJECTDIR/scripts/MAF/maf2xlsx.R $RELEASE/$REL_DIR/all_samples/$STUDY-filtered_mutations_all_${NAME}_keep.maf
        
        cp -v $PROJECTDIR/analysis/variants_combined/$VERSION/$cohort/sample_list.tsv $RELEASE/$REL_DIR/all_samples/
    fi
done


##### Create a README.txt file for the release #####

date > $RELEASE/README.txt

cat >> $RELEASE/README.txt << END

Directories and files:

files.log
- a log file of files copied and coverted to xlsx format

NOTE: The DOMAINS column in the .xlsx files may be trunctated to 1000 characters.
(Excel has a ~32,000 character limit, but here we limit it to 1000. Use the original MAF
as to view all DOMAINS informtion)

# Main directories

one_tumour_per_patient
- variants from tumours when considering only one tumour per patient 

independent_tumours (if applicable)
- variants from all tumours, including multiple independent tumours from the same patient.
- does NOT include samples taken from the SAME tumour unless specified by Ingrid/Louise.

all_tumours
- variants from all tumours, including duplicates (samples from the same tumours) 

# Sub-directories

all_samples
- variants from all matched t/n pairs and unmatched tumours
- *keep.maf *keep.xlsx - all variants in coding and splice sites regions (including synonymous variants)
- *keepPA.maf *keepPA.xslx - protein-altering and splice site altering variants only

matched_samples
- variants from matched t/n pairs only
- *keep.maf *keep.xlsx - all variants in coding and splice sites regions (including synonymous variants)
- *keepPA.maf *keepPA.xslx - protein-altering and splice site altering variants only

QC_keep
QC_keepPA
- Quality control plots showing VAF vs depth for variant calls, variant call types,
top 12 recurrent sites and genes. FOR QC PURPOSES ONLY!

# Other files

sample_list.tsv 
- tumour/normal pairs used for analysis, including unmatched and insilico normal PDv38is_wes_v2

sample_list_matched.tsv
- matched tumour/normal pairs used for analysis

mutations_per_megabase.tsv
- mutation rate by sample, calculated for 'keep' MAFs only


For information on MAF files used for DERMATLAS:

https://drive.google.com/drive/u/1/folders/1CjKEPyBFIm3hwbKi1ja4xybMyF6uOMOC

------------------------------------------------------------
END

tree $RELEASE >> $RELEASE/README.txt



