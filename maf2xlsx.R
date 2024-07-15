library(writexl)
library(stringr)
library(dplyr)

args <- commandArgs(trailingOnly = T)

if (is.na(args[1])) {
  stop("Input a MAF file")
}

infile <- args[1]

if (!file.exists(infile)) {
  stop(paste("File", infile, "does not exist"))
}

maf <- read.table(infile, header = T, sep = "\t", quote = "", stringsAsFactors = F, comment.char = "", check.names = F)

# outfile <- paste0(sub(".maf", "", basename(infile)) ,".xlsx")
outfile <- paste0(sub(".maf", "", infile), ".xlsx")

## print(outfile)

# Check for maximum string length

maf2 <- maf %>% mutate(DOMAINS = str_sub(DOMAINS, 1, 1000))
## print(head(maf2))
print("Limiting DOMAINS column to 1000 characters")

write_xlsx(maf2, path = outfile, col_names = T, format_headers = F)

print(paste("Output file is", outfile))
