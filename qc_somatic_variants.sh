#!/bin/bash

usage () {
	echo -e "\nUsage $0: [options]
	-l   File with a list of VCFs
	-m   Output MAF file suffix
	-s   Path to script directory
	-b   Genome build for the Build column in the MAF file, e.g. GRCh38
	-a   Name of the column with SNP allele frequencies, e.g. gnomAD_AF
	-f   Filter to use (filter1 or filter2)
	-p   Skip MAF creation and run plotting only, if MAFs already exist [optional]

	NOTE: Use 'filter1' to indicate filtering by VAF >= 0.1 or 'filter2' to filter by 'ONP', VAF and REF and ALT lengths (indels)\n"
	exit 1;
}


while getopts "l:m:s:b:a:f:ph" flag; do
	case "${flag}" in
		l) vcflist=$OPTARG
		   ;;
		m) maf_file=$OPTARG
		   ;;
		s) SCRIPTDIR=$OPTARG
		   ;;
		b) BUILD=$OPTARG
		   ;;
		a) AF_COL=$OPTARG
		   ;;
		f) which=$OPTARG
		   ;;
		p) skip_maf=1
		   ;;
		h | *) 
		   usage
		   exit 0
		   ;;
#		: )
#		   echo "Invalid option: $OPTARG requires an argument" 1>&2
#		   ;;	
    esac
done
shift $((OPTIND -1))


# Check input parameters

if [[ -z "$SCRIPTDIR" || -z "$vcflist" ]]; then
	usage
	exit 1
fi

if [[ ! -z "$which" && ! ($which == "filter1" || $which == "filter2") ]]; then
	echo "The filter option -f should be 'filter1' to indicate filtering by VAF >= 0.1 or 'filter2' to filter by 'ONP', VAF and REF and ALT lengths (indels)"
	usage
	exit 1
fi

# Check that directories and files exist

echo "Checking script directory and input files"

SCRIPTDIR=$SCRIPTDIR/MAF

if [[ ! -d $SCRIPTDIR ]]; then
	echo "Directory $SCRIPTDIR does not exist"
	exit 1
elif [[ ! -e $SCRIPTDIR/reformat_vcf2maf.pl ]]; then
	echo "File $SCRIPTDIR/reformat_vcf2maf.pl does not exist"
	exit 1
elif [[ ! -e $SCRIPTDIR/plot_vaf_vs_depth_from_maf.R ]]; then
	echo "File $SCRIPTDIR/plot_vaf_vs_depth_from_maf.R does not exist"
	exit 1
elif [[ ! -e $SCRIPTDIR/maketileplot_from_maf.R ]]; then
	echo "$SCRIPTDIR/maketileplot_from_maf.R does not exist"
	exit 1
elif [[ ! -e $vcflist ]]; then 
	echo "File $vcflist does not exist"
	exit 1
fi

echo "Output directory is current directory: $PWD"
if [[ ! -z "$skip_maf" ]]; then
	echo "Skipping vcf2maf; making plots only"
fi

# Check that each VCF exists

for vcf in `cat $vcflist`; do
	if [[ ! -e $vcf ]]; then
		echo "NOT FOUND: $vcf"	
		exit 1
	fi
done


# Reformat VCFs to MAF
#
#	Usage: reformatVCF2TSV.pl --vcflist [file] 
#	
#	where file is a file with a list of VCF files (full or relative paths). File formats
#	parsed are MuTect (v1), Strelka2, cgpCaVEMan and cgpPindel and may be a combination
#	of any of these. VCF type will be determined by header information.
#
#	Optional:
#
#		--af_col    Column name with population AFs to be used when filtering variants to 'keep'
#		            and  variants of interest (along with variant consequences. Population AF 
#		            cutoffs 0.01 (common SNPs) and 0.001 used for 'keep' and 'voi', respectively).
#		--build     NCBI Genome build [default: $BUILD]
#
#	Output options:
#		--pass      Print out PASS calls only (based only on PASS filter)
#		--voi_only  Print out variants of interest only (not compatible with --keep)
#		--keep      Print out 'keep' and 'keep-PA' variants only (not compatible with --voi_only or --keepPA)
#		--keepPA    Print out 'keep-PA' variants only (not compatible with --voi_only or --keep)
#
#		--tum_vaf   Print out variants with tumour VAF above this value [default: not used].

#		            Works with any of the 3 output options above.


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

if [[ -z "$skip_maf" ]]; then
	$SCRIPTDIR/reformat_vcf2maf.pl --vcflist $vcflist --build $BUILD --pass --sample_list sample_list.tsv --dbsnp_filter > pass_${maf_file}
	$SCRIPTDIR/reformat_vcf2maf.pl --vcflist $vcflist --build $BUILD --voi_only --dbsnp_filter > voi_${maf_file}

	head -n1 pass_${maf_file} > keep_${maf_file}
	awk 'BEGIN{FS="\t"}{if($28~/keep/){print}}' pass_${maf_file} >> keep_${maf_file}

	head -n1 keep_${maf_file} > keepPA_${maf_file}
	awk 'BEGIN{FS="\t"}{if($28~/keep-PA/){print}}' keep_${maf_file} >> keepPA_${maf_file}

	# Remove tumours run with in silico normal

	grep -v PDv38is_wes_v2 pass_${maf_file} > pass_matched_${maf_file}
	grep -v PDv38is_wes_v2 keep_${maf_file} > keep_matched_${maf_file}
	grep -v PDv38is_wes_v2 keepPA_${maf_file} > keepPA_matched_${maf_file}
	grep -v PDv38is_wes_v2 voi_${maf_file} > voi_matched_${maf_file}

	grep -v wes_v2 sample_list.tsv > sample_list_matched.tsv

	cut -f 1 sample_list.tsv > sample_list_tum.tsv
	cut -f 1 sample_list_matched.tsv > sample_list_matched_tum.tsv

	# Filter VCFs
	##--vcflist test.list  --build $BUILD  --keepPA --tum_vaf 0.1 --vaf_filter_type indel --indel_filter

	if [[ $which == "filter1" || $which == "filter2" ]]; then
		$SCRIPTDIR/reformat_vcf2maf.pl --vcflist $vcflist --build $BUILD --pass --tum_vaf 0.1 --vaf_filter_type indel --af_col $AF_COL --dbsnp_filter > pass_vaf_filt_${maf_file}
		$SCRIPTDIR/reformat_vcf2maf.pl --vcflist $vcflist --build $BUILD --voi_only --tum_vaf 0.1 --vaf_filter_type indel --af_col $AF_COL --dbsnp_filter > voi_vaf_filt_${maf_file}
		head -n1 pass_vaf_filt_${maf_file} > keep_vaf_filt_${maf_file}
		awk 'BEGIN{FS="\t"}{if($28~/keep/){print}}' pass_vaf_filt_${maf_file} >> keep_vaf_filt_${maf_file}

		head -n1 keep_vaf_filt_${maf_file} > keepPA_vaf_filt_${maf_file}
		awk 'BEGIN{FS="\t"}{if($28~/keep-PA/){print}}' keep_vaf_filt_${maf_file} >> keepPA_vaf_filt_${maf_file}

		grep -v PDv38is_wes_v2 pass_vaf_filt_${maf_file} > pass_vaf_filt_matched_${maf_file}
		grep -v PDv38is_wes_v2 keep_vaf_filt_${maf_file} > keep_vaf_filt_matched_${maf_file}
		grep -v PDv38is_wes_v2 keepPA_vaf_filt_${maf_file} > keepPA_vaf_filt_matched_${maf_file}
		grep -v PDv38is_wes_v2 voi_vaf_filt_${maf_file} > voi_vaf_filt_matched_${maf_file}

	fi

	if [[ $which == "filter2" ]]; then
		$SCRIPTDIR/reformat_vcf2maf.pl --vcflist $vcflist --build $BUILD --pass --tum_vaf 0.1 --vaf_filter_type indel --af_col $AF_COL --indel_filter --dbsnp_filter > pass_vaf_size_filt_${maf_file}
		$SCRIPTDIR/reformat_vcf2maf.pl --vcflist $vcflist --build $BUILD --voi_only --tum_vaf 0.1 --vaf_filter_type indel --af_col $AF_COL --indel_filter --dbsnp_filter > voi_vaf_size_filt_${maf_file}

		head -n1 pass_vaf_size_filt_${maf_file} > keep_vaf_size_filt_${maf_file}
		awk 'BEGIN{FS="\t"}{if($28~/keep/){print}}' pass_vaf_size_filt_${maf_file} >> keep_vaf_size_filt_${maf_file}

		head -n1 keep_vaf_size_filt_${maf_file} > keepPA_vaf_size_filt_${maf_file}
		awk 'BEGIN{FS="\t"}{if($28~/keep-PA/){print}}' keep_vaf_size_filt_${maf_file} >> keepPA_vaf_size_filt_${maf_file}

		grep -v PDv38is_wes_v2 pass_vaf_size_filt_${maf_file} > pass_vaf_size_filt_matched_${maf_file}
		grep -v PDv38is_wes_v2 keep_vaf_size_filt_${maf_file} > keep_vaf_size_filt_matched_${maf_file}
		grep -v PDv38is_wes_v2 keepPA_vaf_size_filt_${maf_file} > keepPA_vaf_size_filt_matched_${maf_file}
		grep -v PDv38is_wes_v2 voi_vaf_size_filt_${maf_file} > voi_vaf_size_filt_matched_${maf_file}

	fi
fi

# Make plots

#echo "Loading perl-5.38.0"
#module load perl-5.38.0
#which perl
#
#export PERL5LIB=/software/team113/dermatlas/perl5/5.38.0/perl5/


for type in keep keepPA voi; do
	mkdir -p plots_${type}
	cd plots_${type}
	sample_col=`head -n1 ../${type}_${maf_file} | sed 's/\t/\n/g' | grep -n Barcode | cut -f 1 -d ":"`
	plot_height=`cut -f $sample_col ../${type}_${maf_file} | grep -v Barcode | sort -u |wc -l`
	check=$(echo $plot_height / 5 | bc -l | perl -ne 's/\S+\.(\S+)/$1/;print')
	if [[ "$plot_height" -le 5 ]]; then
		plot_height=1
	elif [[ "$check" -gt 0 ]]; then
		let "plot_height=1+$(echo $plot_height / 5 | bc)"
	else
		plot_height=$plot_height
	fi
	plot_height=$(echo $plot_height*1.5 | bc -l)
	Rscript $SCRIPTDIR/plot_vaf_vs_depth_from_maf.R --file ../${type}_${maf_file} --samplefile ../sample_list.tsv --width 10 --height $plot_height -ncol 5
	cut -f 3 top_recurrently_mutated_genes.tsv |sort -u | grep -v Hugo_ > top_genes.list
	echo "Plotting tile plot $type"
	Rscript $SCRIPTDIR/maketileplot_from_maf.R -a ../${type}_${maf_file} -s ../sample_list_tum.tsv  -g top_genes.list --sortbyfrequency -w 8 -t 5 
	echo "Rscript $SCRIPTDIR/maketileplot_from_maf.R -a ../${type}_${maf_file} -s ../sample_list_tum.tsv  -g top_genes.list --sortbyfrequency -w 8 -t 5"
	cd ..
done

for type in keep_matched keepPA_matched voi_matched; do
	mkdir -p plots_${type}
	cd plots_${type}
	sample_col=`head -n1 ../${type}_${maf_file} | sed 's/\t/\n/g' | grep -n Barcode | cut -f 1 -d ":"`
	plot_height=`cut -f $sample_col ../${type}_${maf_file} | grep -v Barcode | sort -u |wc -l`
	check=$(echo $plot_height / 5 | bc -l | perl -ne 's/\S+\.(\S+)/$1/;print')
	if [[ "$plot_height" -le 5 ]]; then
		plot_height=1
	elif [[ "$check" -gt 0 ]]; then
		let plot_height=1+$(echo $plot_height / 5 | bc)
	else
		plot_height=$plot_height
	fi
	plot_height=$(echo $plot_height*1.5 | bc -l)
	Rscript $SCRIPTDIR/plot_vaf_vs_depth_from_maf.R --file ../${type}_${maf_file} --samplefile ../sample_list_matched.tsv --width 10 --height $plot_height -ncol 5
	cut -f 3 top_recurrently_mutated_genes.tsv | sort -u | grep -v Hugo_ > top_genes.list
	echo "Plotting tile plot $type"
	Rscript $SCRIPTDIR/maketileplot_from_maf.R -a ../${type}_${maf_file} -s ../sample_list_matched_tum.tsv -g top_genes.list --sortbyfrequency -w 8 -t 5
	echo "Rscript $SCRIPTDIR/maketileplot_from_maf.R -a ../${type}_${maf_file} -s ../sample_list_matched_tum.tsv -g top_genes.list --sortbyfrequency -w 8 -t 5"
	cd ..
done

if [[ $which == "filter1" || $which == "filter2" ]]; then
	type_array=("keep_vaf_filt" "keepPA_vaf_filt" "voi_vaf_filt")
	if [[ $which == "filter2" ]]; then
		type_array+=("keep_vaf_size_filt" "keepPA_vaf_size_filt" "voi_vaf_size_filt")
	fi
	for type in ${type_array[*]}; do
		mkdir -p plots_${type}
		cd plots_${type}
		sample_col=`head -n1 ../${type}_${maf_file} | sed 's/\t/\n/g' | grep -n Barcode | cut -f 1 -d ":"`
		plot_height=`cut -f $sample_col ../${type}_${maf_file} | grep -v Barcode | sort -u |wc -l`
		check=$(echo $plot_height / 5 | bc -l | perl -ne 's/\S+\.(\S+)/$1/;print')
		if [[ "$plot_height" -le 5 ]]; then
			plot_height=1
		elif [[ "$check" -gt 0 ]]; then
			let plot_height=1+$(echo $plot_height / 5 | bc)
		else
			plot_height=$plot_height
		fi
		plot_height=$(echo $plot_height*1.5 | bc -l)
		Rscript $SCRIPTDIR/plot_vaf_vs_depth_from_maf.R --file ../${type}_${maf_file} --samplefile ../sample_list.tsv --width 10 --height $plot_height -ncol 5
		cut -f 3 top_recurrently_mutated_genes.tsv |sort -u | grep -v Hugo_ > top_genes.list
		echo "Plotting tile plot $type"
		Rscript $SCRIPTDIR/maketileplot_from_maf.R -a ../${type}_${maf_file} -s ../sample_list_tum.tsv  -g top_genes.list --sortbyfrequency -w 8 -t 5 
		echo "Rscript $SCRIPTDIR/maketileplot_from_maf.R -a ../${type}_${maf_file} -s ../sample_list_tum.tsv  -g top_genes.list --sortbyfrequency -w 8 -t 5"
		cd ..
	done

	type_matched_array=("keep_vaf_filt_matched" "keepPA_vaf_filt_matched" "voi_vaf_filt_matched")
	if [[ $which == "filter2" ]]; then
		type_matched_array+=("keep_vaf_size_filt_matched" "keepPA_vaf_size_filt_matched" "voi_vaf_size_filt_matched")
	fi
	for type in ${type_matched_array[*]}; do
		mkdir -p plots_${type}
		cd plots_${type}
		sample_col=`head -n1 ../${type}_${maf_file} | sed 's/\t/\n/g' | grep -n Barcode | cut -f 1 -d ":"`
		plot_height=`cut -f $sample_col ../${type}_${maf_file} | grep -v Barcode | sort -u |wc -l`
		check=$(echo $plot_height / 5 | bc -l | perl -ne 's/\S+\.(\S+)/$1/;print')
		if [[ "$plot_height" -le 5 ]]; then
			plot_height=1
		elif [[ "$check" -gt 0 ]]; then
			let plot_height=1+$(echo $plot_height / 5 | bc)
		else
			plot_height=$plot_height
		fi
		plot_height=$(echo $plot_height*1.5 | bc -l)
		Rscript $SCRIPTDIR/plot_vaf_vs_depth_from_maf.R --file ../${type}_${maf_file} --samplefile ../sample_list_matched.tsv --width 10 --height $plot_height -ncol 5
		cut -f 3 top_recurrently_mutated_genes.tsv | sort -u | grep -v Hugo_ > top_genes.list
		echo "Plotting tile plot $type"
		Rscript $SCRIPTDIR/maketileplot_from_maf.R -a ../${type}_${maf_file} -s ../sample_list_matched_tum.tsv -g top_genes.list --sortbyfrequency -w 8 -t 5
		echo "Rscript $SCRIPTDIR/maketileplot_from_maf.R -a ../${type}_${maf_file} -s ../sample_list_matched_tum.tsv -g top_genes.list --sortbyfrequency -w 8 -t 5"
		cd ..
	done
fi

