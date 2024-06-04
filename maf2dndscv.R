.libPaths(c("/software/team113/dermatlas/R/R-4.2.2/lib/R/library/"))

library(dplyr)
library(stringr)
library(GetoptLong)

cov <- NULL

# Get command line paramaters

outfile <- "dndscv_variants.tsv"

spec = "

Convert a MAF to dndscv input format. Merge by patient if required

  Usage: Rscript maf2dndscv [options]

  Options:
    <infile=s> MAF file [required]
    <outfile=s> Output file [default: dndscv_variants.tsv]
    <by_patient=s> 'y' or 'n' [required]

"

GetoptLong(spec, template_control = list(opt_width = 21))

if (! file.exists(infile)) {
	stop(paste("File", infile, "does not exist"))
} else if (!(by_patient %in% c("y", "n"))) {
	stop("Option --by_patient must be 'y' or 'n'")
} else {
	check_dir <- dirname(outfile)
	if (! dir.exists(check_dir)) {
		stop(paste("Directory for output file", outfile, "does not exist"))
	}
}

print(paste("Output file name is", outfile))	

# Read in MAF and format into dndscv input format

maf <- read.table(infile, header = T, sep = "\t", stringsAsFactors = F)

if (by_patient == "y") {
	maf$Tumor_Sample_Barcode <- str_sub(maf$Tumor_Sample_Barcode, end = -2)
}

dndscv_table <- maf %>% filter(!(Variant_Type %in% c("DNP", "TNP", "ONP"))) %>%
        select(Tumor_Sample_Barcode, Chromosome, POS_VCF, REF_VEP, ALT_VEP) %>%
        rename(sampleID = Tumor_Sample_Barcode, chr = Chromosome, pos = POS_VCF, ref = REF_VEP, mut = ALT_VEP) %>%
        mutate(chr = str_replace(chr, "chr", "")) %>% distinct()


write.table(dndscv_table, file = outfile, sep = "\t", row.names = F, quote = F)
