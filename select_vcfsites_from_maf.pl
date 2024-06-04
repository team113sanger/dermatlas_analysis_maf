#!/usr/bin/env perl

use strict;
use warnings;

# Look for MAF columns in DERMALTAS MAF
# Chromosome	POS_VCF	REF_VEP	ALT_VEP

my $help = <<END;

	Usage: $0 file.maf file.vcf sampleID

END

my $maf = shift @ARGV or die $help;
my $vcf = shift @ARGV or die $help;
my $sample = shift @ARGV or die $help;

my @find = ("Chromosome", "POS_VCF", "REF_VEP", "ALT_VEP", "TUMOUR");
my %sites;
my %index;

open MAF, "cat $maf | head -n1 |" or die "Can't open $maf\n";
while (<MAF>) {
	my $line = $_;
	chomp $line;
	my @colnames = split(/\t/, $line);
	my $i = 0;
	foreach my $name (@colnames) {
		if (grep{ /^$name$/ } @find) {
			$index{$name} = $i;
		}
		$i++;
	}
}
close MAF;

# Get site, ref and alt allleles from MAF

open MAF, "<$maf" or die "Can't open $maf\n";
while (<MAF>) {
	next if /^#/;
	my $line = $_;
	chomp $line;
	my @col = split(/\t/, $line);
	next unless $col[$index{TUMOUR}] eq $sample;
	my ($chr, $pos, $ref, $alt) = ($col[$index{Chromosome}], $col[$index{POS_VCF}], $col[$index{REF_VEP}], $col[$index{ALT_VEP}]);
	$sites{$chr}{$pos}{$ref}{$alt} = 1;
	print STDERR "$chr $pos $ref $alt\n";
	
}
close MAF;

if ($vcf =~ /\.gz$/) {
	open VCF, "zcat $vcf |" or die "Can't open $vcf\n";
} else {
	open VCF, "<$vcf" or die "Can't open $vcf\n";
}

print STDERR "Reading VCF\n";

while (<VCF>) {
	my $line = $_;
	if ($line =~ /^#/) {
		print $line;
		next;
	} else {
		my @col = split(/\t/, $line);
		my ($chr, $pos, $ref, $alt) = @col[0, 1, 3, 4];
		if ($sites{$chr}{$pos}{$ref}{$alt} && $sites{$chr}{$pos}{$ref}{$alt} == 1) {
			print $line;
		}
	}
}
close VCF;

