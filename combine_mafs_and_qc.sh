#!/bin/bash

#OUTDIR=$1
#PROJECTDIR=$2
#RELEASE=$3
#LIST=$4

#RSCRIPT=/software/team113/dermatlas/R/R-4.2.2/bin/Rscript
#RSCRIPT2="module load R/4.1.3; Rscript"
#RSCRIPT2="/software/R-4.1.3/bin/Rscript"

module load R/4.1.3
RSCRIPT=Rscript

usage () {
	echo -e "\nCombine Caveman and Pindel call from a specific release.\n\nUsage $0: [options]
	Required options:
	-p   Full path to project directory 
	-o   Relative or absolute path to output directory
	-r   Name of release directory, e.g. 'release_v1'
	-l   Sample list type [all|independent|onePerPatient]\n"
	exit 1;
}


while getopts ":p:o:r:l:h" flag; do
	case "${flag}" in
		p) PROJECTDIR=$OPTARG
		   ;;
		o) OUTDIR=$OPTARG
		   ;;
		r) RELEASE=$OPTARG
		   ;;
		l) LIST=$OPTARG
		   ;;
		: )
		   echo "Invalid option: $OPTARG requires an argument" 1>&2
		   ;;	
		h | *) 
		   usage
		   exit 0
		   ;;
    esac
done
shift $((OPTIND -1))


if [[ -z "$PROJECTDIR" || -z "$OUTDIR" || -z "$RELEASE" || -z "$LIST" ]]; then
	usage
	#echo -e "\nUsage: somatic_variants_qc.sh /path/to/output/dir  /path/to/projectdir/\n"

	exit 0
fi


if [[ ! -d $OUTDIR ]]; then
	echo "Directory $OUTDIR does not exist"
	exit
else
	cd $OUTDIR
	echo "Current directory is $OUTDIR"
fi


if [[ ! -d $PROJECTDIR/analysis/caveman/$RELEASE/ ]]; then
	echo "Directory $PROJECTDIR/analysis/caveman/$RELEASE/ does not exist"
	exit
elif [[ ! -d $PROJECTDIR/analysis/pindel/$RELEASE/ ]]; then
	echo "Directory $PROJECTDIR/analysis/pindel/$RELEASE/ does not exist"
	exit
fi

# Check that directories and files exist

echo "Looking for script directory"

SCRIPTDIR=$PROJECTDIR/scripts/MAF

if [[ ! -d $SCRIPTDIR ]]; then
	echo "Directory $SCRIPTDIR does not exist"
	exit
elif [[ ! -e $SCRIPTDIR/plot_vaf_vs_depth_from_maf.R ]]; then
	echo "File $SCRIPTDIR/plot_vaf_vs_depth_from_maf.R does not exist"
	exit
elif [[ ! -e $SCRIPTDIR/maketileplot_from_maf.R ]]; then
	echo "$SCRIPTDIR/maketileplot_from_maf.R does not exist"
	exit
fi

echo "Found script directory"


# Make QC plots

#  Usage: Rscript plot.R [options]
#
#  Options:
#    --file, -f character  Variants file [required]
#    --genefile, -g character List of genes to include. 
#    --colname, -c character Column in variants file with gene name matching the gene list. Required with --genelist
#    --suffix, -s character Suffix to be added to output files when using a custom gene list. Required with --genelist
#    --noindels            Ingore indels (only SNVs and MNVs will be considered)
#    --indels_only, -i     Plot indels only
#    --protein_alt, -p     Exclude variants with main consequence synonymous, stop/start retained, splice_region
#
#  Options for AF vs Depth plots:
#    --width, -w numeric   Width (inches) of AF vs Depth by sample plots
#    --height numeric      Height (inches) of AF vs Depth by sample plots
#    --ncol integer        Number of columns for AF vs Depth plot
# 
#

# Check that files exist


for type in keep_matched keepPA_matched keep keepPA; do
	caveman_file=$PROJECTDIR/analysis/caveman/$RELEASE/QC_plots_${LIST}/${type}_caveman_$LIST.maf
	pindel_file=$PROJECTDIR/analysis/pindel/$RELEASE/QC_plots_${LIST}/${type}_pindel_$LIST.maf
	for file in $caveman_file $pindel_file; do 
		if [[ ! -e $file ]]; then
			echo "No such file $file"
			exit
		fi
	done
done

for type in keep_vaf_filt_matched keepPA_vaf_filt_matched keep_vaf_filt keepPA_vaf_filt keep_vaf_size_filt_matched keepPA_vaf_size_filt_matched keep_vaf_size_ filt keepPA_vaf_size_filt ; do
	pindel_file=$PROJECTDIR/analysis/pindel/$RELEASE/QC_plots_${LIST}/${type}_pindel_$LIST.maf
	if [[ ! -e $file ]]; then
		echo "No such file $file"
		exit
	fi
done

list_matched=$PROJECTDIR/analysis/caveman/$RELEASE/QC_plots_${LIST}/sample_list_matched.tsv
list_all=$PROJECTDIR/analysis/caveman/$RELEASE/QC_plots_${LIST}/sample_list.tsv

for list in $list_matched $list_all; do
	if [[ ! -e $list ]]; then
		echo "No such file $list"
		exit
	fi
done

# Combine Caveman and Pindel MAFs


for type in keep_matched keepPA_matched keep keepPA \
	keep_vaf_filt_matched keepPA_vaf_filt_matched \
	keep_vaf_filt keepPA_vaf_filt \
	keep_vaf_size_filt_matched keepPA_vaf_size_filt_matched \
	keep_vaf_size_filt keepPA_vaf_size_filt ; do
	mkdir -p plots_combined_${type}
	cd plots_combined_${type}
	echo $PWD
	echo "Processing $type"
	out_maf=caveman_pindel_${type}.maf

	caveman_file=$PROJECTDIR/analysis/caveman/$RELEASE/QC_plots_${LIST}/${type}_caveman_$LIST.maf
	pindel_file=$PROJECTDIR/analysis/pindel/$RELEASE/QC_plots_${LIST}/${type}_pindel_$LIST.maf
	if [[ "${type}"  == *vaf_filt* || "${type}" == *vaf_size_filt* ]]; then
		orig_type=`echo $type | sed 's/_vaf_filt//; s/_vaf_size_filt//'`
		caveman_file=$PROJECTDIR/analysis/caveman/$RELEASE/QC_plots_${LIST}/${orig_type}_caveman_$LIST.maf
		pindel_file=$PROJECTDIR/analysis/pindel/$RELEASE/QC_plots_${LIST}/${type}_pindel_$LIST.maf
	fi
	echo -e "Files:\n$caveman_file\n$pindel_file"
	cp $caveman_file $out_maf
	grep -v ^Hugo $pindel_file >> $out_maf
    echo "Checking wc"
    wc -l $caveman_file
    wc -l $pindel_file
    wc -l $out_maf

	if [[ "${type}"  == *matched* ]]; then
		echo "Copying $list_matched"
		cp $list_matched .
		cut -f 1 $list_matched > sample_list_matched_tum.tsv
	else 
		echo "Copying $list_all"
		cp $list_all .
		cut -f 1 sample_list.tsv > sample_list_tum.tsv
	fi
	# Make QC plots
	echo "Making QC plots"
	$RSCRIPT $SCRIPTDIR/plot_vaf_vs_depth_from_maf.R --file $out_maf --width 8 --height 10 -ncol 7
	echo "Running $RSCRIPT $SCRIPTDIR/plot_vaf_vs_depth_from_maf.R --file $out_maf --width 8 --height 10 -ncol 7"
	# Make tile plot of top genes
	if [[ ! -e top_recurrently_mutated_genes.tsv ]]; then
		echo "No such file $PWD/top_recurrently_mutated_genes.tsv"
		exit
	fi
	cut -f 3 top_recurrently_mutated_genes.tsv | sort -u | grep -v Hugo_ > top_genes.list
	echo "Plotting tile plot $type"

	if [[ "${type}"  == *matched* ]]; then
		$RSCRIPT $SCRIPTDIR/maketileplot_from_maf.R -a $out_maf  -s sample_list_matched_tum.tsv -g top_genes.list --sortbyfrequency -w 8 -t 5
	else
		$RSCRIPT $SCRIPTDIR/maketileplot_from_maf.R -a $out_maf  -s sample_list_tum.tsv -g top_genes.list --sortbyfrequency -w 8 -t 5
	fi

	cd ..
done


