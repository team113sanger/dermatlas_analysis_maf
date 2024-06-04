library(writexl)
library(stringr)
library(dplyr)

args <- commandArgs(trailingOnly = T)

if (is.na(args[1])) {
	stop("Input a TSV file")
}

infile <- args[1]

if (! file.exists(infile)) {
	stop(paste("File", infile, "does not exist"))
} 

outfile <- paste0(sub(".tsv|.txt", "", infile) ,".xlsx")
print(paste("Output file is", outfile))

# In cases where different rows have different columns, need to define the number of columns and use fill=T

max <- max(count.fields(infile, sep = "\t"))
#print(max)

tsv_file <- read.table(infile, header = F, sep = "\t", quote = "", stringsAsFactors = F, comment.char = "", check.names = F, col.names = paste0("V",seq_len(max)), fill = T)

#head(tsv_file)
write_xlsx(tsv_file, path = outfile, col_names = F, format_headers = F)



