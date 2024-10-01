#!/usr/bin/env perl

use strict;
use warnings;
no warnings 'experimental::smartmatch';
use List::MoreUtils qw(firstidx uniq);
use List::Util qw(max);
use File::Basename;
use Getopt::Long;
#use match::simple;

my $vcflist;
my $pass;
my $pop_af = 'N/A';
my $voi_only;
my $keep_only;
my $keepPA_only;
my $build = "GRCh38";
my $tum_vaf;
my $vaf_filter_type;
my $indel_filter;
my $common_snp_filter;
my $keep_multi;
my $cutoff;
my $print_header = 0;
my @main_csq_header;
my $sample_list;
my @sample_list;
my $transcript_list;
my %alt_transcripts;
my $help;
my $more_help;


GetOptions ("vcflist=s"     => \$vcflist,
			"pass"          => \$pass,
			"keep"          => \$keep_only,
			"keepPA"        => \$keepPA_only,
			"af_col=s"      => \$pop_af,
			"voi_only"      => \$voi_only,
			"build=s"       => \$build,
			"tum_vaf=f"     => \$cutoff,
			"vaf_filter_type=s" => \$vaf_filter_type,
			"indel_filter"  => \$indel_filter,
			"dbsnp_filter"  => \$common_snp_filter,
			"keep_multi"    => \$keep_multi,
			"sample_list=s" => \$sample_list,
			"transcripts=s" => \$transcript_list,
			"help"          => \$help,
			"more_help"    => \$more_help);

if ($more_help) {
	print(usage("long"));
	exit;
} elsif (!$vcflist || $help) {
	print(usage("short"));
	exit;
}

if (($keep_only && $voi_only) || ($keep_only && $keepPA_only) || ($voi_only && $keepPA_only)) {
	print STDERR "$keep_only, $voi_only, $keepPA_only\n";
	die "Choose EITHER --keep or --keepPA or --voi_only\n";
}

my $check_af = $pop_af && $pop_af ne 'N/A' ? 1 : 0;

if ($vaf_filter_type) {
	if ($vaf_filter_type ne 'all' && $vaf_filter_type ne 'snv' && $vaf_filter_type ne 'indel') {
		die "--vaf_filter_type must be one of: 'all', 'snv' or 'indel'\n";
		exit;
	} elsif (!$cutoff) {
		die "Must use --tum_vaf with --variant_filter_type\n";
	}
}
if ($cutoff) {
	if ($cutoff > 1 || $cutoff < 0) {
		die "--tum_vaf must a value ranging from 0 to 1\n";
	}
	elsif (!$vaf_filter_type) {
		die "Must use --vaf_filter_type with --tum_vaf\n";
	}
}

print STDERR "OPTIONS SELECTED:\n"; 
print STDERR "--vcflist: $vcflist\n";
print STDERR "--build: $build\n";
print STDERR "--af_col: $pop_af\n" if $pop_af;
print STDERR "--pass: $pass\n" if $pass;
print STDERR "--voi_only: $voi_only\n" if $voi_only;
print STDERR "--keep: $keep_only\n" if $keep_only;
print STDERR "--keepPA: $keepPA_only\n" if $keepPA_only;
print STDERR "--sample_list: $sample_list\n" if $sample_list;
print STDERR "--tum_vaf: $cutoff\n" if $cutoff;
print STDERR "--indel_filter: $indel_filter\n" if $indel_filter;
print STDERR "--dbsnp_filter: $common_snp_filter\n" if $common_snp_filter;
print STDERR "--vaf_filter_type: $vaf_filter_type\n" if $vaf_filter_type;
print STDERR "--keep_multi: $keep_multi\n" if $keep_multi;
print STDERR "--transcripts: $transcript_list\n" if $transcript_list;

if ($check_af == 0) {
	print STDERR "Option --af_col not used; 'voi' and 'keep' will not consider population frequencies\n";
}

if ($keep_multi) {
	print STDERR "Keep and split multi-allelic sites\n"
} else {
	print STDERR "Skipping multi-allelic sites with ',' in the ALT column\n";
}

my @vcflist = `cat $vcflist`;
chomp @vcflist;

if ($transcript_list) {
	open TRANS, "<$transcript_list" or die "Can't open $transcript_list\n";
	while (<TRANS>) {
		chomp;
		my ($gene, $trans_id) = split(/\t/, $_);
		die "Error in transcript file; gene and transcript ID not found\n" if !$gene || !$trans_id;
		$alt_transcripts{$gene} = $trans_id;
	}
	close TRANS;
	print STDERR "Using transcripts in $transcript_list\n";
}

# Get the TSV column headers, excluding CSQ headers
my @header = get_column_headers();

# Get all Ensembl conseqeuneces, ordered by severity
my @all_conseq = get_conseq_order();

# VEP consequences to exclude for 'keep' and 'keep-PA'
my %exclude_cons = (
	'splice_region_variant' => 1,
	'start_retained_variant' => 1,
	'stop_retained_variant' => 1,
	'synonymous_variant' => 1,
	'mature_miRNA_variant' => 1,
	'5_prime_UTR_variant' => 1,
	'3_prime_UTR_variant' => 1,
	'non_coding_transcript_exon_variant	' => 1,
	'intron_variant' => 1,
	'NMD_transcript_variant' => 1,
	'non_coding_transcript_variant' => 1,
	'upstream_gene_variant' => 1,
	'downstream_gene_variant' => 1,
	'intergenic_variant' => 1
);

my %keep_cons = (
	'start_retained_variant' => 1,
	'stop_retained_variant' => 1,
	'synonymous_variant' => 1,
);

# Convert VEP SO term to maftools-style terms
my %so2maf = %{convert_so2maf()};

my %csq_check_headers;

my @samples;

# Process each VCF
foreach my $vcf (@vcflist) {
	my @csq_header;
	my $mutect = 0;
	my $mutect2 = 0;
	my $strelka = 0;
	my $caveman = 0;
	my $pindel = 0;
	my $gatk_hc = 0;

	print STDERR "Processing $vcf\n";
	if ($vcf =~ /\.gz/) {
		print STDERR "Inferring bgzip file\n";
		open F, "zcat $vcf |" or die "Can't open file $vcf\n";
	} else {
		open F, "cat $vcf |" or die "Can't open file $vcf\n";
	}

	my ($tumour, $normal, $tum_index, $norm_index);
	while (<F>) {
		my $line = $_;
		chomp $line;
		# Deal with header information, infer caller, and check sample names, if present
		if ($line =~ /^#/) {
			if ($line =~ /##INFO=<ID=CSQ.+Format: (\S+)">/) {
				my $format = $1;
				@csq_header = split (/\|/, $format);
				if (!@csq_header) {
					die "ERROR: Found CSQ in header but no CSQ format\n";
				}
				# change some headers to make compatible for maftools
				@csq_header = fix_csq_header(\@csq_header);
				# Print out the header only once
				if (!@main_csq_header) {
					my @headers_to_check = ("BIOTYPE", "CANONICAL", "Consequence", "HGVSp", $pop_af, "Hugo_Symbol", "SOURCE", "Existing_variation", "Transcript_ID");
					%csq_check_headers = parse_csq_header(\@csq_header, @headers_to_check);
					foreach my $colname ("BIOTYPE", "CANONICAL", "Consequence", $pop_af, "Transcript_ID") {
						if (! defined($csq_check_headers{$colname})) {
							die "ERROR: CSQ in header is missing \"$colname\"\n" if $colname ne $pop_af || (($colname eq $pop_af) && $check_af);
						}
					}
					@main_csq_header = @csq_header;
				} elsif (! @main_csq_header ~~ @csq_header ) {
					die "VCF CSQ headers in $vcf are not consistent\n";
				}
				#my @headers_to_check = ("BIOTYPE", "CANONICAL", "Consequence", "HGVSp", $pop_af, "Hugo_Symbol", "SOURCE", "Existing_variation", "Transcript_ID");
#				#%csq_check_headers = parse_csq_header(\@main_csq_header, @headers_to_check);
				#%csq_check_headers = parse_csq_header(\@csq_header, @headers_to_check);
				if ($print_header == 0) {
					foreach my $remove ("Hugo_Symbol", "Existing_variation", "SOURCE") {
						$main_csq_header[$csq_check_headers{$remove}] = "REMOVE" if $csq_check_headers{$remove};
					}
					@main_csq_header = grep {$_ ne 'REMOVE'} @main_csq_header; 
					# Print out main header
					print join("\t", @header, @main_csq_header) . "\n";
					$print_header = 1;
				}
			} elsif ($line =~ /##GATKCommandLine/ && $line =~ /tumor/ && $line !~ /Mutect2/) {
				$normal = $1 if $line =~ /normal_sample_name=(\S+)/;
				$tumour = $1 if $line =~ /tumor_sample_name=(\S+)/;
				if (!$normal || !$tumour) {
					die "No tumour or normal sample found in Mutect file header.\n";
				}
				$mutect = 1;
				print STDERR "Inferring Mutect_v1 format\n";
			} elsif ($line =~ /##GATKCommandLine=<ID=Mutect2/) {
				$mutect2 = 1;
				print STDERR "Inferring Mutect2 format\n";
				print STDERR "Multi-allelic sites will be skipped for Mutect2 calls\n";
			} elsif ($line =~ /##GATKCommandLine/ && $line !~ /Mutect2/) {
				$gatk_hc = 1;
				print STDERR "Inferring GATK germline format\n";
			} elsif ($line =~ /##cmdline=\S+Strelka/) {
				if ($line =~ /--normalBam=\S+\/([^\/]+).sample.dupmarked.bam/) {
					$normal = $1;
				}
				if ($line =~ /--tumorBam=\S+\/([^\/]+).sample.dupmarked.bam/) {
					$tumour = $1;
				}
				if (!$normal || !$tumour) {
					die "No tumour or normal sample found in Strelka file header.\n";
				}
				$strelka = 1;
				print STDERR "Inferring Strelka format\n";
			} elsif ($line =~ /##cavemanVersion/) {
				$caveman = 1;
				print STDERR "Inferring CaVEMan format\n";
			} elsif ($line =~ /##source.+pindel/) {
				$pindel = 1;
				print STDERR "Inferring Pindel format\n";
			} elsif ($line =~ /^##SAMPLE=<ID=NORMAL.+SampleName=([^,>]+)/) {
				$normal = $1;
			} elsif ($line =~ /^##SAMPLE=<ID=TUMOUR.+SampleName=([^,>]+)/) {
				$tumour = $1;
			} elsif ($line =~ /^##normal_sample=(\S+)/) {
				$normal = $1;
			} elsif ($line =~ /^##tumor_sample=(\S+)/) {
				$tumour = $1;
			} elsif ($line =~ /#CHR/) {
				my @colnames = split(/\t/, $line);
				if ($mutect == 1 || $mutect2 == 1) {
					if (!$normal || !$tumour) {
						die "Can't find tumour and normal sample names in VCF header.";
					} elsif ($colnames[9] eq $normal && $colnames[10] eq $tumour) {
						$norm_index = 9;
						$tum_index = 10;
					} elsif ($colnames[9] eq $tumour && $colnames[10] eq $normal) {
						$norm_index = 10;
						$tum_index = 9;
					} else {
						die "Sample names $tumour and $normal do not match FORMAT headers $colnames[9] and $colnames[10]\n";
					}
				} elsif ($strelka == 1 || $pindel == 1 || $caveman == 1) {
					if (!$normal || !$tumour) {
						die "Can't find tumour and normal sample names in VCF header.";
					}
					$norm_index = 9;
					$tum_index = 10;
				} elsif ($gatk_hc == 1) {
					@samples = @colnames[9..$#colnames];
					chomp @samples;
					$norm_index = 9;
					#$tumour = "NA";
				} else {
					die "Cannot find GATK, Mutect, Strelka, Pindel or CaVEMan VCF headers\n";
				}
				# Write a list of t and n pairs
				if ($sample_list) {
					if ($gatk_hc != 1) {
						push @sample_list, "$tumour\t$normal\n";
						#print SAM "$tumour\t$normal\n";
					} elsif ($gatk_hc == 1) {
						push @sample_list, join("\n", @samples);
						#print SAM join("\n", @samples);
					}
				}
			}
			next;
		}
		# Process the VCF line
		# If file is GATK-HC germline, split the lines up by sample
		# make parsing easier
		my @lines;
		my @gatk_hc_samples;
		if ($gatk_hc == 1) {
			my @cols = split(/\t/, $line);
			foreach my $sample_index (0..$#samples) {
				my $index = 9 + $sample_index;
				my $sample_format = $cols[$index];
				next if $sample_format =~ /^0\S0/ || $sample_format =~ /^\S\S\.:/ || $sample_format =~ /^\.\S\S:/;
				push @lines, join("\t", @cols[0..8], $sample_format);
				push @gatk_hc_samples, $samples[$sample_index];
			}
		} else {
			push @lines, $line;
		}
		foreach my $vcfline_num (0..$#lines) {
			my @c = split(/\t/, $lines[$vcfline_num]);
			my ($chrom, $pos, $ids, $ref, $alt, $qual, $passfilt) = @c[0..6];
			my ($info, $csq) = ($1, $2) if $c[7] =~ /(\S+)CSQ=([^;\s]+)/;
			# This check handles the rare case where Pindel outputs REF==ALT,
			# in which case VEP does not provide an CSQ
			if (!$csq) {
				print STDERR "WARNING: SKIPPING LINE - No CSQ for $line\n";
				next;
			}

			# Skip if FILTER ne PASS and $pass defined
			if ($c[6] ne 'PASS' && $pass) {
				print STDERR "Skipping $line\n";
				next;
			}

			# Skip if --dbsnp_filter and found COMMON
			if ($gatk_hc != 1 && $common_snp_filter && $c[7] =~ /dbSNP=COMMON/) {
				print STDERR "Skipping common snp $line\n";
				next;
			}

			# Skip if --keep_multi and ALT has multiple alleles
			# Or format is Mutect2 (all multi-allelic sites are flagged anyway
			if ($alt =~ /,/ && (! $keep_multi || $mutect2 == 1)) {
				print STDERR "Skipping multi-allelic site $line\n";
				next;
			}

			# Iterate through the different alternative alleles
			my @alts = split(/,/, $alt);
			foreach my $alt_num (0..$#alts) {
				next if $alts[$alt_num] eq "*";
				$alt = $alts[$alt_num];

				# Calculate the start and end postitions, ref and alt allleles, variant type for maftools format
				#my ($pos_start, $pos_end, $ref_maf, $alt_maf) = calc_positions($pos, $type, $len, $ref, $alt);
				my ($pos_start, $pos_end, $ref_maf, $alt_maf, $type, $len) = calc_positions($pos, $ref, $alt);

				# Parse tumour and normal format columns
				my $vaf;
				my $vaf_norm;
				my $tum_counts;
				my $norm_counts;
				my $dp_tum;
				my $dp_norm;

				my @format = split(/:/, $c[8]);

				my @tum_format = split (/:/, $c[$tum_index]) if $gatk_hc != 1;
				my @norm_format = split (/:/, $c[$norm_index]);

				# Parse the FORMAT column
				if ($mutect == 1) {
					($vaf, $vaf_norm, $tum_counts, $norm_counts, $dp_tum, $dp_norm)  = mutect_stats(\@format, \@tum_format, \@norm_format);
				} elsif ($mutect2 == 1) {
					#($vaf, $vaf_norm, $tum_counts, $norm_counts, $dp_tum, $dp_norm)  = mutect2_stats(\@format, \@tum_format, \@norm_format);
					#($vaf_norm, $norm_counts, $dp_norm)  = gatkhc_stats(\@format, \@norm_format);
					($vaf, $vaf_norm, $tum_counts, $norm_counts, $dp_tum, $dp_norm)  = mutect2_stats(\@format, \@tum_format, \@norm_format);
					if ($vaf_norm eq 'NA') {
						print STDERR "Skipping non-ref multiallelic : $line\n";
						next;
					}
				} elsif ($strelka == 1) {
					($vaf, $vaf_norm, $tum_counts, $norm_counts, $dp_tum, $dp_norm)  = strelka_stats(\@format, \@tum_format, \@norm_format);
				} elsif ($caveman == 1) {
					($vaf, $vaf_norm, $tum_counts, $norm_counts, $dp_tum, $dp_norm)  = caveman_stats(\@format, \@tum_format, \@norm_format, $ref, $alt);
				} elsif ($pindel == 1) {
					($vaf, $vaf_norm, $tum_counts, $norm_counts, $dp_tum, $dp_norm)  = pindel_stats(\@format, \@tum_format, \@norm_format, $ref, $alt);
				} elsif ($gatk_hc == 1) {
					($vaf_norm, $norm_counts, $dp_norm)  = gatkhc_stats(\@format, \@norm_format, $alt_num);
					if ($vaf_norm eq 'NA') {
						print STDERR "Skipping non-ref multiallelic : $line\n";
						next;
					}
					$tum_counts = "NA\tNA";
					$vaf  = "NA";
					$dp_tum = "NA";
					$normal = $gatk_hc_samples[$vcfline_num];
					$tumour = $gatk_hc_samples[$vcfline_num];
				}
			
				# Apply a hard tumour VCF cutoff if provided by user
				if ($gatk_hc != 1 && $cutoff && $vaf < $cutoff) {
					if ($vaf_filter_type eq 'all' || ($vaf_filter_type eq 'indel' && ($strelka == 1 || $pindel == 1)) || ($vaf_filter_type eq 'snv' && ($mutect == 1 || $caveman == 1))) {
						print STDERR "Skipping : low tumour VAF in $line\n";
						next;
					}
				} elsif ($gatk_hc == 1 && $cutoff && $vaf_norm < $cutoff) {
					my $type;
					if (length($ref) == length($alt)) {
						$type = 'snv'; # filter snv and mnv the same way
					} else {
						$type = 'indel';
					}
					if ($vaf_filter_type eq 'all' || ($vaf_filter_type eq 'indel' && $type eq 'indel') || ($vaf_filter_type eq 'snv' && $type eq 'snv')) {
						print STDERR "Skipping : low tumour VAF in $line\n";
						next;
					}
				}
				# Additional filtering for indels
				# Keep if length of indel is less than 25 bp but not a complex indel with both alleles > 10, or ONP
				# if not, keep if tum VAF > 0.25 and depth in T and N are both >= 20. Subtract 1 to length to account for anchor base
				if ($indel_filter && ($strelka == 1 || $pindel == 1 || $gatk_hc == 1)) {
					if (! (length($ref) - 1 <= 25 && length($alt) - 1 <= 25 && !(length($ref) - 1 > 10 && length($alt) - 1 > 10) && $type ne 'ONP')) {
						if ($gatk_hc != 1 && $vaf > 0.25 && $dp_tum >= 20 && $dp_norm >= 20 && $type ne 'ONP') {
							# do nothing
						} elsif ($gatk_hc == 1 && $vaf_norm > 0.25 && $dp_norm >= 20 && $type ne 'ONP') {
							# do nothing
						} else {
							print STDERR "Skipping indel: $line\n";
							next;
						}
					}
				}

				my @newlines;
				my %genes_found;

				# Parse each consequence in the CSQ section.
				# Print one line per consequence (per gene, per allele)
				my @csq = split /\,/, $csq;
				foreach my $c (@csq) {
					# split with -1 to include empty fields
					my @cs = split(/\|/, $c, -1);

					# Check choose the consequences that match the alt allele;
					# This is used to get the consequences for multi-allelic sites
					# (NB: mutect v1, pindel and strelka do not have multi-allelic sites
					# so no consequences should be skipped. 
					if ($cs[0] ne $alt_maf) {
						print STDERR "Alt is $alt and found $cs[0]. Skipping $c\n";
						next;
					}
					foreach my $i (0..$#cs) {
						$cs[$i] = '-' if !$cs[$i];
					}
					# Parse out most severe conseq if multiple types, 
					# flag if variant of interest, pass, check gnomAD or
					# other vaf, check biotype, main consequence, canonical status
					# in order to flag consequence lines to keep
					my ($voi, $keep, $main_csq) = flag_voi(\@cs, $passfilt, \%csq_check_headers, $check_af);

					# Get the maftools consequence term
					if (! $so2maf{$main_csq}) {
						die "Cannot find Maftools term for $main_csq\n";
					}
					my $class_maf = $so2maf{$main_csq};

					# Fix frameshift mutation and protein-altering annotations
					if ($class_maf eq 'Frame_Shift') {
						if ($type eq 'INS') {
							$class_maf .= "_Ins";
						} elsif ($type eq 'DEL') {
							$class_maf .= "_Del";
						} else {
							die "Found Frame_Shift but $type does not match $class_maf\n";
						}
					} elsif ($class_maf eq "protein_altering_variant") {
						my $len_diff = length($ref) - length($alt);
						$class_maf = $len_diff % 3 == 0 ? "In_Frame" : "Frame_Shift";
						if ($type  eq 'DEL') {
							$class_maf .= "_Del";
						} else {
							$class_maf .= "_Ins";
						}
					}

					# Get the protein change from HGVSp
					my $prot_change = get_protein_change($cs[$csq_check_headers{HGVSp}], $main_csq);
					# Convert HEX "%3D" to "=" in case --no_escape was not used when running VEP
					$cs[$csq_check_headers{HGVSp}] =~ s/%3D/=/;
					# Get SYMBOL from CSQ section so it can be printed with maftools required columns
					my $symbol = $cs[$csq_check_headers{Hugo_Symbol}];
					# 
					# Remove SYMBOL from @cs array before appending
					# and remove unused columns
					foreach my $remove ("Hugo_Symbol", "Existing_variation", "SOURCE") {
						$cs[$csq_check_headers{$remove}] = "REMOVE" if $csq_check_headers{$remove};
					}
					@cs = grep {$_ ne 'REMOVE'} @cs; 
					my @columns = (
						$symbol,
						'Sanger',
						$build,
						$chrom,
						$pos_start,
						$pos_end,
						$class_maf,
						$type,
						$ref_maf,
						$alt_maf,
						$tumour,
						$prot_change,
						$pos,
						$ids,
						$qual,
						$passfilt,
						$tumour,
						$normal,
						$vaf,
						$tum_counts,
						$dp_tum,
						$vaf_norm,
						$norm_counts,
						$dp_norm,
						$voi,
						$keep,
						$main_csq,
						$ref,
						$alt,
						@cs
					);
					# Print out line, depending on options
					if ($voi_only) {
						#print join("\t", @columns) . "\n" if $voi eq 'yes';
						$genes_found{$symbol}++;
						push @newlines, join("\t", @columns) . "\n" if $voi eq 'yes';
					} elsif ($keepPA_only) {
						#print join("\t", @columns) . "\n" if $keep =~ /keep/;
						$genes_found{$symbol}++;
						push @newlines, join("\t", @columns) . "\n" if $keep =~ /keep-PA/;
					} elsif ($keep_only) {
						#print join("\t", @columns) . "\n" if $keep =~ /keep/;
						$genes_found{$symbol}++;
						push @newlines, join("\t", @columns) . "\n" if $keep =~ /keep/;
					} else {
						$genes_found{$symbol}++;
						push @newlines, join("\t", @columns) . "\n";
					}
				}
				# Re-order variants affecting more than one gene, to prioritize
				# when using maftools oncoplot
		#		if (keys %genes_found > 1) {
					@newlines = reorder_by_conseq(\@newlines, \@all_conseq, 28, 0);
		#		}
				print @newlines if @newlines;
			}
		}
	}
	close F;
}

if ($sample_list) {
	@sample_list = uniq @sample_list;
	open SAM, ">$sample_list" || die "Can't open $sample_list\n";
	print SAM @sample_list;
	close SAM;
}

###############
# Subroutines #
###############

# Defined column headers

sub get_column_headers {
	my @header = (
		'Hugo_Symbol',
		'Center',
		'NCBI_Build',
		'Chromosome',
		'Start_Position',
		'End_Position',
		'Variant_Classification',
		'Variant_Type',
		'Reference_Allele',
		'Tumor_Seq_Allele2',
		'Tumor_Sample_Barcode',
		'HGVSp_Short',
		'POS_VCF',
		'ID',
		'QUAL',
		'FILTER',
		'TUMOUR',
		'NORMAL',
		'VAF_tum',
		't_ref_count',
		't_alt_count',
		't_depth',
		'VAF_norm',
		'n_ref_count',
		'n_alt_count',
		'n_depth',
		'Variant_of_interest',
		'Keep',
		'Main_consequence_VEP',
		'REF_VEP',
		'ALT_VEP'
	);

	return @header;

}

# Infer variant type for maftools

sub get_variant_type {
	my ($ref, $alt) = @_;
	my $type;
	my $len;
	$ref =~ s/-//;
	$alt =~ s/-//;
	if ($ref && $alt && length($alt) == length($ref)) {
		if (length($ref) == 1) {
			$type = 'SNP';
			$len = 1;
		} elsif (length($ref) == 2) {
			$type = 'DNP';
			$len = 2;
		} elsif (length($ref) == 3) {
			$type = 'TNP';
			$len = 3;
		} else {
			$type = 'ONP';
			$len = length($ref);
		}
	} elsif (!$alt || length($ref) > length($alt)) {
		$type = "DEL";
		$len = length($ref) - length($alt);
	} elsif (!$ref || length($ref) < length($alt)) {
		$type = "INS";
		$len = length($alt) - length($ref);
	}

	return ($type, $len);
}

# Parse the CSQ header names and get indexes

sub parse_csq_header {
	my ($csq, @names) = @_;
	my @csq = @$csq;
	my %c_idx;
	foreach my $i (0..$#csq) {
		foreach my $colname (@names) {
			if ($csq[$i] eq $colname) {
				$c_idx{$colname} = $i;
			}
		}
	}

	return %c_idx;
}

# Change the CSQ headers to maftools headers

sub fix_csq_header {
	my $header = shift @_;
	my %replace = (
		'SYMBOL' => 'Hugo_Symbol',
		#'Gene' => 'Ensembl_Gene_ID',
		'STRAND' => 'TRANSCRIPT_STRAND',  # Transcript strand; MAF 'stand' is variant strand, always '+'
		'EXON' => 'Exon_Number',
		'INTRON' => 'Intron_Number',
		'Feature' => 'Transcript_ID',
	);
	foreach my $name (@$header) {
		$name = $replace{$name} if exists $replace{$name};
	}

	return @$header;
}

# Reorder CSQ by consequences. In cases where a variant affects 2 genes,
# maftools keeps the first non-duplicate in the order they are found, based on: 
# Chromosome, Start_Position, Tumor_Sample_Barcode, Reference_Allele, Tumor_Seq_Allele2

sub reorder_by_conseq {
	my ($csq_lines, $all_conseq, $main_csq_idx, $gene_idx) = @_;
	my @sorted;
	my %vars;
	foreach my $csq_line (@$csq_lines) {
		my @cols = split(/\t/, $csq_line);
		push @{$vars{$cols[$main_csq_idx]}{$cols[$gene_idx]}}, $csq_line;
	}
	# sort by main consequence then by gene
	foreach my $cons (@$all_conseq) {
		my @genes = sort keys %{$vars{$cons}};
		if (@genes > 1) {
			foreach my $gene (sort keys %{$vars{$cons}}) {
				push @sorted, @{$vars{$cons}{$gene}};
			}
		} elsif (@genes == 1) {
			push @sorted, @{$vars{$cons}{$genes[0]}};
		}
	}

	return @sorted;
}

# Flag variants of interest and variants to keep

sub flag_voi {
	my ($csq, $flag, $csq_idx, $check_af) = @_;
	my @csq = @$csq;
	my %csq_idx = %$csq_idx;
	my $main_csq;
	my $voi_af_fail = 0;
	my $all_csq = $csq[$csq_idx{Consequence}];
	my $highest_af = 0;
	
	# From custom annotations, if multiple AFs, VEP combines them with "&" 
	if ($check_af) {
		if ($csq[$csq_idx{$pop_af}] eq '-' || $csq[$csq_idx{$pop_af}] !~ /\&/) {
			$highest_af = $csq[$csq_idx{$pop_af}];
			$highest_af =~ s/\S+://;
		} else {
			$highest_af = choose_af($csq[$csq_idx{$pop_af}]);
		}
	}
	# 'splice_region_variant' doesn't have much meaning;
	# could be synonymous, intron, etc. Listed first if 
	# with 'synonymous'; use the next consequence
	if ($all_csq =~ /^splice_region_variant\&(\S+)/) {
		$all_csq = $1;
	}
	if ($all_csq =~ /^(\w+)/) {
		$main_csq = $1;
	} else {
		$main_csq = $all_csq;
	}

	# Determine variants of interest, variants to keep and main consequence
	# If a list of non-canonical transcripts is provided, first check for those
	if ($transcript_list && exists($alt_transcripts{$csq[$csq_idx{Hugo_Symbol}]})) {
		if ($alt_transcripts{$csq[$csq_idx{Hugo_Symbol}]} ne $csq[$csq_idx{Transcript_ID}]) {
			return ('no', 'no', $main_csq)
		} elsif (($csq[$csq_idx{BIOTYPE}] ne 'protein_coding') || $flag ne 'PASS') {
			return ('no', 'no', $main_csq);
		}
	} elsif (($csq[$csq_idx{BIOTYPE}] ne 'protein_coding') || ($csq[$csq_idx{CANONICAL}] ne 'YES') || $flag ne 'PASS') {
		return ('no', 'no', $main_csq);
	} elsif ($check_af && $highest_af ne '.' && $highest_af ne '-' && $highest_af >= 0.01) {  # common SNP cutoff
		return('no', 'no', $main_csq);
	}
	if ($check_af && $highest_af ne '.' && $highest_af ne '-' && $highest_af >= 0.001) {  # Variant of interest cutoff
		$voi_af_fail = 1;
	}
	if ($exclude_cons{$main_csq}) {
		if ($keep_cons{$main_csq}) {
			return('no', 'keep', $main_csq);
		} else {
			return ('no', 'no', $main_csq);
		}
	} else {
		if ($check_af && $voi_af_fail == 1) {
			return ('no', 'keep-PA', $main_csq);
		} else {
			return ('yes', 'keep-PA', $main_csq);
		}
	}
}

# Split the values in pop_af column and get the max

sub choose_af {
	my $afs = shift @_;
	my @afs = split('&', $afs);
	@afs = map { s/\S+://; $_ } @afs;
	return eval(max @afs);
}

# VAF and allele counts from Mutect

sub mutect_stats {
	my ($format, $tum_format, $norm_format) = @_;
	my $vaf_index = firstidx { $_ eq 'FA' } @$format; 
	my $ad_index = firstidx { $_ eq 'AD' } @$format; 
	my $dp_index = firstidx { $_ eq 'DP' } @$format; 
	if (! defined($dp_index)) {
		die "Can't find DP in " . join(':', @$format) . "\n";
	}	
	if (! (defined($vaf_index) && defined($ad_index))) {
		die "Can't find FA and/or AD in " . join(':', @$format) . "\n";
	}
	my $vaf = sprintf("%.3f", $tum_format->[$vaf_index]);
	my $vaf_norm = sprintf("%.3f", $norm_format->[$vaf_index]);
	my $tum_counts = join("\t", split(',', $tum_format->[$ad_index]));
	my $norm_counts = join("\t", split(',', $norm_format->[$ad_index]));
	if (! (defined($vaf) && defined($vaf_norm) && defined($tum_counts) && defined($norm_counts))) {
		die join(' ', "Can't find tumour or normal allele counts and/or AFs in" , join(':', @$tum_format) , join(':', @$norm_format)) . "\n";
	}

	return ($vaf, $vaf_norm, $tum_counts, $norm_counts, $tum_format->[$dp_index], $norm_format->[$dp_index]);
}

sub mutect2_stats_old {
	my ($format, $tum_format, $norm_format) = @_;
#	my $vaf_index = firstidx { $_ eq 'AF' } @$format; 
	my $ad_index = firstidx { $_ eq 'AD' } @$format; 
	my $dp_index = firstidx { $_ eq 'DP' } @$format; 
	if (! defined($dp_index)) {
		die "Can't find DP in " . join(':', @$format) . "\n";
	}	
#	if (! (defined($vaf_index) && defined($ad_index))) {
	if (! defined($ad_index)) {
		die "Can't find AD in " . join(':', @$format) . "\n";
	}
#	my $vaf = sprintf("%.3f", $tum_format->[$vaf_index]);
#	my $vaf_norm = sprintf("%.3f", $norm_format->[$vaf_index]);
	my $dp = $tum_format->[$dp_index];
	my $dp_norm = $norm_format->[$dp_index];
	my $tum_counts = join("\t", split(',', $tum_format->[$ad_index]));
	my $norm_counts = join("\t", split(',', $norm_format->[$ad_index]));
	my $vaf = sprintf("%.3f", $tum_counts / $dp);
	my $vaf_norm = sprintf("%.3f", $norm_counts / $dp_norm);
	if (! (defined($vaf) && defined($vaf_norm) && defined($tum_counts) && defined($norm_counts))) {
		die join(' ', "Can't find tumour or normal allele counts and/or AFs in" , join(':', @$tum_format) , join(':', @$norm_format)) . "\n";
	}

	#return ($vaf, $vaf_norm, $tum_counts, $norm_counts, $tum_format->[$dp_index], $norm_format->[$dp_index]);
	return ($vaf, $vaf_norm, $tum_counts, $norm_counts, $dp, $dp_norm);
}

sub mutect2_stats {
	# GT:AD:DP:GQ:PL
	#my ($format, $norm_format) = @_;
	my ($format, $tum_format, $norm_format) = @_;
#	my $vaf_index = firstidx { $_ eq 'FA' } @$format; 
	my $ad_index = firstidx { $_ eq 'AD' } @$format; 
	my $dp_index = firstidx { $_ eq 'DP' } @$format; 
	my $gt_index = firstidx { $_ eq 'GT' } @$format; 
	if (! defined($dp_index)) {
		die "Can't find DP in " . join(':', @$format) . "\n";
	}	
	if (! (defined($ad_index))) {
		die "Can't find AD in " . join(':', @$format) . "\n";
	}
	# GATK HC may have multiallelic sites
	my @gt_tum = split(//, $tum_format->[$gt_index]);
	my @gt_norm = split(//, $norm_format->[$gt_index]);
	# Ignore GTs that are 1/2, 1/3, etc
	if (($gt_tum[0] != $gt_tum[2]) && ($gt_tum[0] > 0 && $gt_tum[2] > 0)) {
		return ("NA", "NA", "NA", "NA", "NA", "NA");
	}
	my @ad_tum = split(/,/, $tum_format->[$ad_index]);
	my $ad_tum = $ad_tum[$gt_tum[2]]; # first allele is the ref allele
	my @ad_norm = split(/,/, $norm_format->[$ad_index]);
	my $ad_norm = $ad_norm[$gt_tum[2]]; # use the index from the tum to get the counts for the normal
#	print STDERR join("\n", @ad_tum) . "\n";
#	print STDERR "$ad_tum[$gt_tum[2]], $gt_tum[2], $ad_tum\n";
#	print STDERR "$dp_index, $norm_format->[$dp_index]\n";
#	print STDERR join(" ", @$norm_format) . "\n";
	if ($norm_format->[$dp_index] == 0 || $tum_format->[$dp_index] == 0) {
		# sometimes a GT is given but allele counts are 0
		return (0, 0, 0, 0, 0, 0);
	}
	my $vaf_norm = sprintf("%.3f", $ad_norm / ($norm_format->[$dp_index]));
	my $vaf_tum = sprintf("%.3f", $ad_tum / ($tum_format->[$dp_index]));
	my @norm_counts = split(',', $norm_format->[$ad_index]);
	my $norm_counts = join("\t", $norm_counts[$gt_norm[0]], $norm_counts[$gt_tum[2]]); # use the index from the tum to get the counts for the normal
	my @tum_counts = split(',', $tum_format->[$ad_index]);
	my $tum_counts = join("\t", $tum_counts[$gt_tum[0]], $tum_counts[$gt_tum[2]]);
#	print STDERR "INDEX $ad_index $norm_format->[$ad_index]\n";
#	print STDERR "NORM " . join(" ", @norm_counts) . "\n";
#	print STDERR "TUM " . join(" ", @tum_counts) . "\n";
	if ((! (defined($vaf_norm) && defined($norm_counts))) || (! (defined($vaf_tum) && defined($tum_counts)))) {
		die join(' ', "Can't find tumour or normal allele counts and/or AFs in" , join(':', @$norm_format), join(':', @$tum_format)) . "\n";
	}

	#return ($vaf_norm, $norm_counts, $norm_format->[$dp_index]);
	return ($vaf_tum, $vaf_norm, $tum_counts, $norm_counts, $tum_format->[$dp_index], $norm_format->[$dp_index]);
}
# VAF and allele counts from Strelka

sub strelka_stats {
	my ($format, $tum_format, $norm_format) = @_;
	my ($altcount, $refcount, $altcount_norm, $refcount_norm);
	# Strelka format: use Tier 1 counts in TAR (ref plus other non-supporting reads)
	# and TIR (indel suppporting reads)
	# DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50:BCN50
	my $refcount_idx = firstidx { $_ eq 'TAR' } @$format;	
	my $altcount_idx = firstidx { $_ eq 'TIR' } @$format;	
	my $dp_index = firstidx { $_ eq 'DP' } @$format; 
	if (! defined($dp_index)) {
		die "Can't find DP in " . join(':', @$format) . "\n";
	}	
	if (! (defined($altcount_idx) && defined($refcount_idx))) {
		die "Can't find TAR and/or TIR in " . join(':', @$format) . "\n";
	}
	if ($tum_format->[$altcount_idx] =~ /(\S+),\S+/) {
		$altcount = $1;
	}
	if ($tum_format->[$refcount_idx] =~ /(\S+),\S+/) {
		$refcount = $1;
	}
	if ($norm_format->[$altcount_idx] =~ /(\S+),\S+/) {
		$altcount_norm = $1;
	}
	if ($norm_format->[$refcount_idx] =~ /(\S+),\S+/) {
		$refcount_norm = $1;
	}
	if (! (defined($refcount) && defined($altcount) && defined($refcount_norm) && defined($altcount_norm))) {
		die join(' ', "Can't find tumour or normal TAR/TIR allele counts in ", 
		join(':', @$format) , join(':', @$tum_format) , join(':', @$norm_format)) . "\n";
	}
	my $vaf = $refcount == 0 && $altcount == 0 ? '0' : sprintf("%.3f", $altcount/($refcount + $altcount));
	my $vaf_norm = $refcount_norm == 0 && $altcount_norm == 0 ? '0' : sprintf("%.3f", $altcount_norm/($refcount_norm + $altcount_norm));
	my $tum_counts = join("\t", $refcount, $altcount);
	my $norm_counts = join("\t", $refcount_norm, $altcount_norm);

	return ($vaf, $vaf_norm, $tum_counts, $norm_counts, $tum_format->[$dp_index], $norm_format->[$dp_index]);
}

# Get average values in the VCF FORMAT column for MNVs, based on specific indexes

sub get_stat {
	my ($indexes, $array, $var_size) = @_;
	my $tot = 0;
	foreach my $i (@$indexes) {
		$tot += $array->[$i];
	}

	return eval($tot/$var_size);
}

# VAF and allele counts from CaVEMan

sub caveman_stats {
	my ($format, $tum_format, $norm_format, $ref, $alt) = @_;
	my ($altcount, $refcount, $altcount_norm, $refcount_norm);
	my @ref_idx;
	my @alt_idx;
	my @format = @$format;
	my @ref = split('', $ref);
	my @alt = split('', $alt);

	# SNVs: GT:FAZ:FCZ:FGZ:FTZ:RAZ:RCZ:RGZ:RTZ:PM
	# MNVs: GT_1:FAZ_1:FCZ_1:FGZ_1:FTZ_1:RAZ_1:RCZ_1:RGZ_1:RTZ_1:PM_1:GT_2:FAZ_2:FCZ_2:FGZ_2:FTZ_2:RAZ_2:RCZ_2:RGZ_2:RTZ_2:PM_2
	my @all_idx = grep { $format[$_] =~ /Z/ } 0..$#format;
	my @pm_idx = grep { $format[$_] =~ /PM/ } 0..$#format;
	my $var_size = @pm_idx;
	if ($var_size == 1) {
		@ref_idx = grep { $format[$_] =~ /${ref}Z/ } 0..$#format;
		@alt_idx = grep { $format[$_] =~ /${alt}Z/ } 0..$#format;
	} else {
		my $base_num = 0;
		foreach my $i (0..$#ref) {
		#foreach my $refbase (@ref) {
			$base_num++;
			my $refbase = $ref[$i];
			my $altbase = $alt[$i];
			push @ref_idx, grep { $format[$_] =~ /${refbase}Z_$base_num/ } 0..$#format;
		#}
		#foreach my $altbase (@alt) {
			push @alt_idx, grep { $format[$_] =~ /${altbase}Z_$base_num/ } 0..$#format;
		}
	}
	if (!@alt_idx || !@ref_idx) {
		die "Can't find allele counts in " . join(':', @format) . "\n";
	}
	my $vaf = sprintf("%.3f", get_stat(\@pm_idx, $tum_format, $var_size));
	my $vaf_norm = sprintf("%.3f", get_stat(\@pm_idx, $norm_format, $var_size));
	my $tum_ref = int(0.5 + get_stat(\@ref_idx, $tum_format, $var_size));
	my $tum_alt = int(0.5 + get_stat(\@alt_idx, $tum_format, $var_size));
	my $norm_ref = int(0.5 + get_stat(\@ref_idx, $norm_format, $var_size));
	my $norm_alt = int(0.5 + get_stat(\@alt_idx, $norm_format, $var_size));
	my $tum_depth = int(0.5 + get_stat(\@all_idx, $tum_format, $var_size));;
	my $norm_depth = int(0.5 + get_stat(\@all_idx, $norm_format, $var_size));
	my $tum_counts = join("\t", $tum_ref, $tum_alt);
	my $norm_counts = join("\t", $norm_ref, $norm_alt);

	if (! (defined($vaf) && defined($vaf_norm) && defined($tum_counts) && defined($norm_counts))) {
		die join(' ', "Can't find tumour or normal allele counts and/or AFs in" , join(':', @$tum_format) , join(':', @$norm_format)) . "\n";
	}

	return ($vaf, $vaf_norm, $tum_counts, $norm_counts, $tum_depth, $norm_depth);
}

# VAF and allele counts from Pindel

sub pindel_stats {
	my ($format, $tum_format, $norm_format, $ref, $alt) = @_;
	my ($altcount, $refcount, $altcount_norm, $refcount_norm);

	# Format GT:PP:NP:PB:NB:PD:ND:PR:NR:PU:NU:FD:FC
	my $FC_index = firstidx { $_ eq 'FC' } @$format; 
	my $FD_index = firstidx { $_ eq 'FD' } @$format;
##	my $vaf = sprintf("%.10g", $tum_format->[$FC_index]/$tum_format->[$FD_index]);
##	my $vaf_norm = sprintf("%.10g", $norm_format->[$FC_index]/$norm_format->[$FD_index]);
	my $vaf = sprintf("%.3f", $tum_format->[$FC_index]/$tum_format->[$FD_index]);
	my $vaf_norm = sprintf("%.3f", $norm_format->[$FC_index]/$norm_format->[$FD_index]);
	my $tum_ref = $tum_format->[$FD_index] - $tum_format->[$FC_index];
	my $norm_ref = $norm_format->[$FD_index] - $norm_format->[$FC_index];
	my $tum_counts = join("\t", $tum_ref, $tum_format->[$FC_index]);
	my $norm_counts = join("\t", $norm_ref, $norm_format->[$FC_index]);

	return ($vaf, $vaf_norm, $tum_counts, $norm_counts, $tum_format->[$FD_index], $norm_format->[$FD_index]);
}


sub gatkhc_stats {
	# GT:AD:DP:GQ:PL
	my ($format, $norm_format, $alt_num) = @_;
	$alt_num++;
#	my $vaf_index = firstidx { $_ eq 'FA' } @$format; 
	my $ad_index = firstidx { $_ eq 'AD' } @$format; 
	my $dp_index = firstidx { $_ eq 'DP' } @$format; 
	my $gt_index = firstidx { $_ eq 'GT' } @$format; 
	if (! defined($dp_index)) {
		die "Can't find DP in " . join(':', @$format) . "\n";
	}	
	if (! (defined($ad_index))) {
		die "Can't find AD in " . join(':', @$format) . "\n";
	}
	# GATK HC may have multiallelic sites
	my @gt = split(//, $norm_format->[$gt_index]);
	# Check to see if the alt allele indexes are the same
	if ($alt_num != $gt[2]) {
		return ("NA", "NA", "NA");
	}
	# Ignore GTs that are 1/2, 1/3, etc
	if (($gt[0] != $gt[2]) && ($gt[0] > 0 && $gt[2] > 0)) {
		return ("NA", "NA", "NA");
	}
	my @ad = split(/,/, $norm_format->[$ad_index]);
	my $ad = $ad[$gt[2]]; # first allele is the ref allele
#	print STDERR join("\n", @ad) . "\n";
#	print STDERR "$ad[$gt[2]], $gt[2], $ad\n";
#	print STDERR "$dp_index, $norm_format->[$dp_index]\n";
#	print STDERR join(" ", @$norm_format) . "\n";
	if ($ad == 0 || $norm_format->[$dp_index] == 0) {
		# sometimes a GT is given but allele counts are 0
		return (0, 0, 0);
	}
	my $vaf_norm = sprintf("%.3f", $ad / ($norm_format->[$dp_index]));
	my @norm_counts = split(',', $norm_format->[$ad_index]);
	my $norm_counts = join("\t", $norm_counts[$gt[0]], $norm_counts[$gt[2]]);
	if (! (defined($vaf_norm) && defined($norm_counts))) {
		die join(' ', "Can't find tumour or normal allele counts and/or AFs in" , join(':', @$norm_format)) . "\n";
	}

	return ($vaf_norm, $norm_counts, $norm_format->[$dp_index]);
}

# Order of Ensembl conseqeuences as of VEP v105

sub get_conseq_order {
	my @order = (
		"transcript_ablation",
		"splice_acceptor_variant",
		"splice_donor_variant",
		"stop_gained",
		"frameshift_variant",
		"stop_lost",
		"start_lost",
		"transcript_amplification",
		"inframe_insertion",
		"inframe_deletion",
		"missense_variant",
		"protein_altering_variant",
		"splice_region_variant",
		"incomplete_terminal_codon_variant",
		"start_retained_variant",
		"stop_retained_variant",
		"synonymous_variant",
		"coding_sequence_variant",
		"mature_miRNA_variant",
		"5_prime_UTR_variant",
		"3_prime_UTR_variant",
		"non_coding_transcript_exon_variant",
		"intron_variant",
		"NMD_transcript_variant",
		"non_coding_transcript_variant",
		"upstream_gene_variant",
		"downstream_gene_variant",
		"TFBS_ablation",
		"TFBS_amplification",
		"TF_binding_site_variant",
		"regulatory_region_ablation",
		"regulatory_region_amplification",
		"feature_elongation",
		"regulatory_region_variant",
		"feature_truncation",
		"intergenic_variant"
	);

	return @order;
}

# Convert VEP consequences to maftools format.
# Based on maftools vcf2maf.pl with exceptions below

sub convert_so2maf {
	my %so_terms = (
		"transcript_ablation" => "Splice_Site", # vcf2maf.pl uses "Splice site"; next highest consequence
		"splice_acceptor_variant" => "Splice_Site",
		"splice_donor_variant" => "Splice_Site",
		"stop_gained" => "Nonsense_Mutation",
		"frameshift_variant" => "Frame_Shift",
		"stop_lost" => "Nonstop_Mutation",
		"start_lost" => "Translation_Start_Site",
		"transcript_amplification" => "Unknown", # vcf2maf uses "Intron"
		"inframe_insertion" => "In_Frame_Ins",
		"inframe_deletion" => "In_Frame_Del",
		"missense_variant" => "Missense_Mutation",
		"protein_altering_variant" => "protein_altering_variant", # vcf2maf indcludes these with frameshifts or inframe ins/del
		"splice_region_variant" => "Splice_Region",
		"incomplete_terminal_codon_variant" => "Silent", # vcf2maf uses "Silent"
		"start_retained_variant" => "Silent", # not inf vcf2maf but stop_retained is "Silent"
		"stop_retained_variant" => "Silent",
		"synonymous_variant" => "Silent",
		"coding_sequence_variant" => "Targeted_Region", # vcf2maf uses "Missense"
		"mature_miRNA_variant" => "RNA",
		"5_prime_UTR_variant" => "5'UTR",
		"3_prime_UTR_variant" => "3'UTR",
		"non_coding_transcript_exon_variant" => "RNA",
		"intron_variant" => "Intron",
		"NMD_transcript_variant" => "Silent", # vcf2maf uses "Silent"
		"non_coding_transcript_variant" => "RNA",
		"upstream_gene_variant" => "5'Flank",
		"downstream_gene_variant" => "3'Flank",
		"TFBS_ablation" => "Targeted_Region",
		"TFBS_amplification" => "Targeted_Region",
		"TF_binding_site_variant" => "Targeted_Region",
		"regulatory_region_ablation" => "Targeted_Region",
		"regulatory_region_amplification" => "Targeted_Region",
		"feature_elongation" => "Targeted_Region",
		"regulatory_region_variant" => "Targeted_Region",
		"feature_truncation" => "Targeted_Region",
		"intergenic_variant" => "IGR"
    );

    return(\%so_terms);
}

# Variant start and end genomic postions for maftools

sub calc_positions {
	#my ($pos, $type, $len, $ref, $alt) = @_;
	my ($pos, $ref, $alt) = @_;
	# In VCF, indel first base is the same.
	# DNP, TNP, ONP may be indels (mostly ONPs)
	# Remove anchor base and adjust coordinate
	if (substr($ref, 0, 1) eq substr($alt, 0, 1)) { 
#		print STDERR "Found anchor base $ref $alt\n";
		$ref = substr($ref, 1);
		$alt = substr($alt, 1);
		$ref = '-' if !$ref;	
		$alt = '-' if !$alt;	
		$pos++
	}
	my ($type, $len) = get_variant_type($ref, $alt);
	if ($type eq 'SNP') {
		return ($pos, $pos, $ref, $alt, $type, $ref);
	} elsif ($type eq'DNP') {
		return ($pos, eval($pos + 1), $ref, $alt, $type, $ref);
	} elsif ($type eq 'TNP') {
		return ($pos, eval($pos + 2), $ref, $alt, $type, $ref);
	} elsif ($type eq 'ONP') {
		return ($pos, eval($pos + $len - 1), $ref, $alt, $type, $ref);
	} elsif ($type eq 'DEL') {
		#return (eval($pos + 1), eval($pos + $len), $ref, $alt, $type, $ref);
		return ($pos, eval($pos + $len - 1), $ref, $alt, $type, $ref);
	} elsif ($type eq 'INS') {
		return (eval($pos - 1), $pos, $ref, $alt, $type, $ref);
	} else {
		die "Error: cannot calculate start and stop positions for $pos, $type, $len $ref $alt $type $ref\n";
	}
}


# Get the 1-letter protein change from HGVSp

sub get_protein_change {
	my ($hgvsp, $class) = @_;

#	# Borrowed from vcf2maf.pl https://github.com/mskcc/vcf2maf/blob/main/vcf2maf.pl
	my %aa3to1 = qw( Ala A Arg R Asn N Asp D Asx B Cys C Glu E Gln Q Glx Z Gly G His H Ile I Leu L
    Lys K Met M Phe F Pro P Ser S Thr T Trp W Tyr Y Val V Xxx X Ter * );

	$hgvsp =~ s/%3D/=/;
	$hgvsp =~ s/\S+:(p\.\S+)/$1/;
	foreach my $long (sort keys %aa3to1) {
		$hgvsp =~ s/$long/$aa3to1{$long}/g;
	}

	return $hgvsp;
}


# Usage message

sub usage {

	my $which = shift(@_);

	my $usage = <<END;

    Convert a VCF file to tab-delimited MAF file, compatible with Maftools.

    Usage: reformat_vcf2maf.pl --vcflist [file] [options]
    
    where file is a file with a list of VCF files (full or relative paths). File formats
    parsed are MuTect (v1), Strelka2, Mutect2, GATK gemline, cgpCaVEMan and cgpPindel 
    and may be a combination of any of these. VCF type will be determined by header information.

    General options:

        --build [string]      NCBI Genome build for the 'NCBI_Build' MAF column. [default: GRCh38] 
        --sample_list [file]  Print to file a list of tumour and normal samples from the VCF list.
        --help                Prints help to screen.
        --more_help           Prints extended help to screen.

    Optional variant filters:

        --tum_vaf [vaf]       Print out variants with tumour VAF above this value. [default: not used]
                                If the VCF is detected as germline, --tum_vaf applies to the normal VAF.
                                Value is a fraction between 0 and 1 (e.g. 0.1, not a percentage).
                                Compatible with --pass/--voi_only/--keep/--keepPA.

        --vaf_filter_type [all|snv|indel]    Apply --tum_vaf cutoff to this type of variant only.
                                               The snv includes mnvs from Caveman/Mutect.

        --indel_filter        Apply indel filtering based on REF and ALT length, VAF and T/N depth; remove ONPs.
        --dbsnp_filter        Exclude variants flagged with "dbSNP=COMMON" in the VCF INFO column.
        --keep_multi          Keep and split multi-allelic sites to one per line. (ALT allele has at least one comma ",")

    Optional consequence filters:

        Applying none of the 4 options below generates a MAF with all variants and in all transcripts.
        Otherwise, select only 1 of 4 options below:

        --pass                Print out PASS calls only.
        --keep                Print out 'keep' and 'keep-PA' variants only.
        --keepPA              Print out 'keep-PA' variants only.
        --voi_only            Print out 'variants of interest' only.

            **** See extended help for definitions ****

        --af_col [string]     Column name with population AFs to be used when filtering by consequence.
                                Applies to --keep/--keepPA/--af_col. AF cutoffs are 0.01 for --keepPA/--keep and 0.001 for --voi.

        --transcripts [file]  A file with a list of "Hugo_Symbol[TAB]Transcript_ID" (no header), indicating the
                                non-canonical transcript to use when selecting consequence annotations. See --more_help. 

END

	my $extra = <<END;
    Example command line:

        # Print out protein-altering variants in canonical transcripts, unless it has gnomAD_AF >= 0.01,
        # and use alternative transcripts in alt_transcripts.tsv. Remove indels with VAF < 0.1.
        # Remove variants annotated with "dbSNP=COMMON" and apply the internal indel filter.

        reformat_vcf2maf.pl --vcflist vcf.list --keepPA --tum_vaf 0.1 --vaf_filter_type indel \
            --indel_filter --af_col gnomAD_AF --transcripts alt_transcripts.tsv --dbsnp_filter > out.nmaf


    IMPORTANT NOTES:

    * A variant is designated in the 'keep' column as 'keep' if the MAF if it is PASS is the VCF, has population
           AF < 0.01, the variant consequence is in a coding exon or splice site, and the variant is in a canonical
           transcript. Includes synonymous variants.

    * A variant is designated in the 'keep' column as 'keep-PA' based on the same criteria above, except the variant
           consequence is protein-altering (PA). i.e. synonymous variants are excluded.

    * A variant is designated in the 'voi' column as 'yes' based on the same criteria as 'keep-PA' variants, except
           the population AF < 0.001.

    * If the --pass option is used without any other options, all PASS variants in all transcripts and all consequence
           types are printed to the MAF, except multi-allelic sites (unless --keep_multi is used).

    * The --pass flag can be used with any of the variant or consequence filter options.

    * If --af_col is not used, the criteria for the 'voi' and 'keep' columns will exclude population AFs.

    * The --transcripts option is used to specify a specific Ensembl transcript to use. For example, if an old Ensembl
           version incorrectly assigns the canonical transcript, or if there are different Ensembl genes/transcripts 
           IDs associated with the same HUGO symbol and you only want to pick one gene/transcript. Default is to use
           the canonical transcript as indicated in the VCF VEP annotation.

    * If --transcripts is used, the 'CANONICAL' MAF column will not be altered, but the 'voi' and 'keep' columns
           will be correct, i.e. for genes in the transcripts list, the Ensembl canonical transcript 
           will be set to 'voi' = 'no' and 'keep' = 'no', while the transcripts in the --transcripts file
           will be evaluated as described above and the status set accordingly.

END

	if ($which eq 'long') {
		print STDERR "LONG\n";
		$usage = $usage . $extra;
	}
	return $usage;
}

