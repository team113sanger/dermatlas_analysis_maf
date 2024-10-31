suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(plyr))
suppressMessages(library(gridExtra))
suppressMessages(library(grid))
suppressMessages(library(extrafont))
suppressMessages(library(tidyverse))
suppressMessages(library(cowplot))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))

suppressMessages(extrafont::loadfonts())

# showtext_auto(TRUE)
# font_paths()
# font_add("Arial", "Arial.ttf")
## font_add("Arial_Bold", "Arial_Bold.ttf")
# font_add("Arial_Italic", "Arial_Italic.ttf")
sessionInfo()

args <- commandArgs(trailingOnly = TRUE)

option_list <- list(
  make_option(c("-a", "--fileA"),
    type = "character", default = NULL,
    help = "dataset file name", metavar = "character"
  ),
  make_option(c("-b", "--fileB"),
    type = "character", default = NULL,
    help = "dataset file name", metavar = "character"
  ),
  make_option(c("-c", "--cutoff"),
    type = "integer", default = NULL,
    help = "Min. number of samples with mutation [default = %default]", metavar = "number"
  ),
  make_option(c("-s", "--samples"),
    type = "character", default = NULL,
    help = "File of sample names for fileA [default = %default]", metavar = "character"
  ),
  make_option(c("-g", "--genelist"),
    type = "character", default = NULL,
    help = "File of genes to plot [default = %default]", metavar = "character"
  ),
  make_option("--geneid_list",
    action = "store_true", default = FALSE,
    help = "List in --genelist is in format 'Hugo_Symbol/Gene'"
  ),
  make_option(c("-f", "--filter_col"),
    type = "character", default = NULL,
    help = "Column name and value to filter on (colname, value) [default = %default]", metavar = "character"
  ),
  make_option(c("-m", "--sortbymutations"),
    action = "store_true", default = FALSE,
    help = "Sort genes by number of mutations"
  ),
  make_option(c("-n", "--sortbysamplemut"),
    action = "store_true", default = FALSE,
    help = "Sort samples by number of mutations"
  ),
  make_option(c("-l", "--sortbygenelist"),
    action = "store_true", default = FALSE,
    help = "Sort genes by order given in genelist file"
  ),
  make_option(c("-e", "--sortbysamplelist"),
    action = "store_true", default = FALSE,
    help = "Sort sample by order given in sample file"
  ),
  make_option(c("--sortbyfrequency"),
    action = "store_true", default = FALSE,
    help = "Order by groups of genes and samples"
  ),
  make_option(c("--group"),
    type = "character", default = NULL,
    help = "Group samples by specific genes", metavar = "character"
  ),
  make_option(c("-p", "--prefix"),
    type = "character", default = "gene_tileplot",
    help = "Output file prefix (prefix.pdf) [default = %default]"
  ),
  make_option(c("-w", "--width"),
    type = "integer", default = 8,
    help = "PDF output width (inches) [default = %default]", metavar = "number"
  ),
  make_option(c("-t", "--height"),
    type = "integer", default = 8,
    help = "PDF output height (inches) [default = %default]", metavar = "number"
  ),
  make_option(c("-x", "--ncol"),
    type = "integer", default = 5,
    help = "Number of columns in the legend [default = %default]", metavar = "number"
  ),
  make_option(c("-r", "--nrow"),
    type = "integer", default = 2,
    help = "Number of rows in the legend [default = %default]", metavar = "number"
  ),
  make_option(c("--legendpos"),
    type = "character", default = "top",
    help = "ggplot2 legend.position; comma-separated list accepted if plotA and B have different positions [default = %default]", metavar = "number"
  ),
  make_option(c("-v", "--vertical"),
    action = "store_true", default = FALSE,
    help = "If fileA and fileB, plot A and B vertically [default = %default]"
  ),
  make_option(c("-A", "--relsizeA"),
    type = "integer", default = 1,
    help = "If fileA and fileB, relative width of plotA for plot_grid function (or height if --vertical) [default = %default]", metavar = "number"
  ),
  make_option(c("-B", "--relsizeB"),
    type = "integer", default = 1,
    help = "If fileA and fileB, relative width of plotB for plot_grid function (or height if --vertical) [default = %default]", metavar = "number"
  ),
  make_option(c("--vline"),
    type = "character", default = NULL,
    help = "Comma-separated list of x-axis coordinates to place a vertical line, e.g. 2.5 will place a line between 2nd and 3rd sample [default = %default]"
  ),
  make_option(c("--vlinetype"),
    type = "character", default = NULL,
    help = "Comma-separated list of geom_line types, to go with --vline, e.g. 1,2 where 1=solid, 2=dashed [default = %default]"
  ),
  make_option(c("--vlineplot"),
    type = "character", default = NULL,
    help = "Comma-separated list indicating which plot to draw the geom_vline(s), 1 or 2, to go with --vline, e.g. 1,2 if the first line is for plot A, 2nd for B [default = %default]"
  ),
  make_option(c("--hline"),
    type = "character", default = NULL,
    help = "Comma-separated list of x-axis coordinates to place a horizontal line, e.g. 2.5 will place a line between 2nd and 3rd sample [default = %default]"
  ),
  make_option(c("--hlinetype"),
    type = "character", default = NULL,
    help = "Comma-separated list of geom_line types, to go with --hline, e.g. 1,2 where 1=solid, 2=dashed [default = %default]"
  ),
  make_option(c("--hlineplot"),
    type = "character", default = NULL,
    help = "Comma-separated list indicating which plot to draw the geom_hline(s), 1 or 2, to go with --hline, e.g. 1,2 if the first line is for plot A, 2nd for B [default = %default]"
  ),
  make_option(c("--sizex"),
    type = "character", default = 8,
    help = "Font size for x-axis text [default = %default]"
  ),
  make_option(c("--sizey"),
    type = "character", default = 8,
    help = "Font size for x-axis text [default = %default]"
  ),
  make_option(c("--frequency"),
    type = "character", default = NULL,
    help = "Tab-separated file containing gene and frequency to plot in a heatmap [default = %default]"
  ),
  make_option(c("--sortbysidebar"),
    action = "store_true", default = FALSE,
    help = "Sort genes by the order provided in the frequency file [default = %default]"
  ),
  make_option(c("--label"),
    type = "character", default = NULL,
    help = "Add A and B labels to plots [default = NULL]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

#################
# CHECK OPTIONS #
#################

if (length(args) == 0) {
  stop("At least one argument must be supplied (--fileA).", call. = FALSE)
}

if (isTRUE(opt$sortbymutations) && isTRUE(opt$sortbygenelist)) {
  stop("Cannot sort by number of mutations AND gene list")
}

if (isTRUE(opt$geneid_list) && is.null(opt$genelist)) {
  stop("Option --geneid_list must be used with --genelist")
}

#############
# FUNCTIONS #
#############

###########################################
# Parse geom_hline and geom_vline options #
###########################################

parse_geomlines <- function(pos, type, plot) {
  rows <- length(strsplit(pos, ",")[[1]])
  line.df <- data.frame(matrix(ncol = 3, nrow = rows))
  colnames(line.df) <- c("plot", "pos", "type")

  line.df$pos <- as.numeric(strsplit(pos, ",")[[1]])
  line.df$type <- strsplit(type, ",")[[1]]
  line.df$plot <- strsplit(plot, ",")[[1]]
  return(line.df)
}


############################################
# Rename/merge different Effect categories #
############################################

shorten_consequence <- function(df) {
  print("Shortening Effects")

  data <- data.frame(lapply(data, function(x) {
    gsub("missense_variant", "Missense", x)
  }))
  data <- data.frame(lapply(data, function(x) {
    gsub("stop_gained", "Nonsense", x)
  }))
  data <- data.frame(lapply(data, function(x) {
    gsub("splice_acceptor_variant", "Splice site", x)
  }))
  data <- data.frame(lapply(data, function(x) {
    gsub("splice_donor_variant", "Splice site", x)
  }))
  data <- data.frame(lapply(data, function(x) {
    gsub("frameshift_variant", "Frameshift", x)
  }))
  data <- data.frame(lapply(data, function(x) {
    gsub("inframe_deletion", "Inframe indel", x)
  }))
  data <- data.frame(lapply(data, function(x) {
    gsub("inframe_insertion", "Inframe indel", x)
  }))
  data <- data.frame(lapply(data, function(x) {
    gsub("synonymous_variant", "Synonymous", x)
  }))
  data <- data.frame(lapply(data, function(x) {
    gsub("stop_retained_variant", "Synonymous", x)
  }))
  data <- data.frame(lapply(data, function(x) {
    gsub("start_retained_variant", "Synonymous", x)
  }))
  data <- data.frame(lapply(data, function(x) {
    gsub("coding_sequence_variant", "Other", x)
  }))
  data <- data.frame(lapply(data, function(x) {
    gsub("protein_altering_variant", "Other", x)
  }))
  data <- data.frame(lapply(data, function(x) {
    gsub("stop_lost", "Start/stop codon lost", x)
  }))
  data <- data.frame(lapply(data, function(x) {
    gsub("start_lost", "Start/stop codon lost", x)
  }))
  data <- data.frame(lapply(data, function(x) {
    gsub("incomplete_terminal_codon_variant", "Other", x)
  }))
  data <- data.frame(lapply(data, function(x) {
    gsub("tanscript_ablation", "Other", x)
  }))
  # data <- data.frame(lapply(data, function(x)  { gsub("_", " ", x)}))
  return(data)
}


##########################################
# Plot specific samples and sample order #
##########################################

filter_by_sample <- function(data, samples, sample_col) {
  print(paste("Plotting samples in sample list", samples))
  samplelist <- read.table(samples, header = F)
  # data <- as.data.frame(data %>% filter(Sample %in% samplelist$V1)) %>% droplevels(except = c(1,2))
  # data <- as.data.frame(data %>% filter(!!as.name(sample_col) %in% samplelist$V1)) %>% droplevels(except = c(1,2))
  data <- as.data.frame(data %>% filter(!!as.name(sample_col) %in% samplelist$V1))
  data[, sample_col] <- factor(data[, sample_col], levels = samplelist$V1)
  return(data)
}


######################################
# Arrange sample levels for plotting #
######################################

sort_by_sample <- function(data, samples, sample_col) {
  samplelist <- read.table(samples, header = F)
  print("Sorting samples by sample list")
  data[, sample_col] <- factor(data[, sample_col], levels = samplelist$V1)
  return(data)
}

##########################################
# Order by most frequently mutated genes #
##########################################

order_by_freq2 <- function(data, gene_col, sample_col, cons_col) {
  replace_vals <- function(x) {
    replace(x, x > 1, 1)
  }
  print("Sorting by gene and sample")
  # Make a temporary copy of the data frame and rename columns
  data_tmp <- data %>% rename(Gene_tmp = !!as.name(gene_col), Sample_tmp = !!as.name(sample_col), Effect = !!as.name(cons_col))

  # Get the count of mutated samples per gene and sort to get the gene order
  gene_counts <- data_tmp %>%
    select(Gene_tmp, Sample_tmp) %>%
    distinct() %>%
    count(Gene_tmp, sort = T)
  gene_order <- gene_counts[, "Gene_tmp"]

  # Now get the sample order
  # Create a matrix, transpose and sort by sample

  data_mat <- as.matrix(data_tmp)
  data_mat <- data_tmp[order(match(data_tmp[, "Gene_tmp"], gene_order)), ]
  data_tmp_cast <- as.data.frame(dcast(data_tmp, Sample_tmp ~ Gene_tmp, fun.aggregate = length, value.var = "Effect", drop = FALSE)) %>%
    mutate_if(is.numeric, replace_vals) %>%
    rowwise() %>%
    mutate(total = sum(c_across(where(is.numeric)))) %>%
    arrange(desc(total))
  sample_order <- data_tmp_cast$Sample_tmp

  # Now order the factors in the original data frame
  data[, sample_col] <- factor(data[, sample_col], levels = sample_order)
  data[, gene_col] <- factor(data[, gene_col], levels = gene_order)
  data <- data %>% arrange(factor(!!as.name(gene_col), levels = gene_order))
  return(data)
}

order_by_freq <- function(data, gene_col, sample_col, cons_col) {
  print("Sorting by gene and sample")
  # Make a temporary copy of the data frame and rename columns
  data_tmp <- data %>% rename(Gene_tmp = !!as.name(gene_col), Sample_tmp = !!as.name(sample_col), Effect = !!as.name(cons_col))

  # Get the count of mutated samples per gene and sort to get the gene order
  gene_counts <- data_tmp %>%
    select(Gene_tmp, Sample_tmp) %>%
    distinct() %>%
    count(Gene_tmp, sort = T)
  gene_order <- gene_counts[, "Gene_tmp"]
  data_tmp$Gene_tmp <- factor(data_tmp$Gene_tmp, levels = gene_order)

  # Now get the sample order
  # transpose and sort by sample

  # data_mat <- as.matrix(data_tmp)
  data_tmp_dcast <- dcast(data_tmp, Sample_tmp ~ Gene_tmp, fun.aggregate = length, value.var = "Effect", drop = F)
  val <- data_tmp_dcast[, -1, drop = F]
  sam <- data_tmp_dcast[, 1]
  val <- data.frame(apply(val, 2, function(x) replace(x, x > 1, 1)))
  merged <- cbind(sam, val)
  # sorted_sam <- merged[do.call(order, as.list(-merged[,2:ncol(merged), with=F])),]$sam
  sorted_sam <- as.character(merged[do.call(order, as.list(-merged[, 2:ncol(merged)])), ]$sam)
  merged$sam <- factor(merged$sam, levels = sorted_sam)
  merged <- merged %>% arrange(match(sam, levels(merged$sam)))

  # Now order the factors in the original data frame
  data[, sample_col] <- factor(data[, sample_col], levels = sorted_sam)
  data[, gene_col] <- factor(data[, gene_col], levels = gene_order)
  # data <- data %>% arrange(factor(!!as.name(gene_col), levels = gene_order))
  return(data)
}


#########################
# Apply cutoff if given #
#########################

# to account for multiple mutations per gene per sample, first group by
# gene and sample

filter_by_minmutations <- function(data, cutoff, gene_col, sample_col) {
  print(paste("Plotting genes mutated in at least", cutoff, "samples"))
  # data <- as.data.frame(data %>% add_count(!!as.name(gene_col)) %>% filter(n >= cutoff) %>% select(-n))
  gene_counts <- data %>%
    select(!!as.name(gene_col), !!as.name(sample_col)) %>%
    distinct() %>%
    add_count(!!as.name(gene_col)) %>%
    filter(n >= cutoff)
  data <- data %>% filter(!!as.name(gene_col) %in% gene_counts[, gene_col])
  data[, gene_col] <- factor(data[, gene_col], levels = unique(data[, gene_col]))
  return(data)
}


#######################
# Plot specific genes #
#######################

filter_by_genelist <- function(data, genelist, gene_col) {
  print(paste("Plotting genes in gene list ", genelist))
  genes <- read.table(genelist, header = F)
  data <- as.data.frame(data %>% filter(!!as.name(gene_col) %in% genes$V1))
  data[, gene_col] <- factor(data[, gene_col], levels = genes$V1)
  return(data)
}


##########################
# Filter by column value #
##########################

filter_by_column <- function(data, col_filt) {
  col <- strsplit(col_filt, split = ",")
  colname <- col[[1]][1]
  colvalue <- col[[1]][2]
  print(paste("Filtering by value", colvalue, "in column", colname))
  data <- data %>% filter(!!as.name(colname) == colvalue)
  return(data)
}


############################
# Plot specific gene order #
############################

sort_by_genelist <- function(data, genelist, gene_col) {
  genes <- read.table(genelist, header = F)
  print(paste("Sorting genes by gene list", genelist))
  ## data$Gene <- factor(data$(!!as.name(gene_col)), levels = genes$V1)
  data[, gene_col] <- factor(data[, gene_col], levels = genes$V1)
  return(data)
}


#####################################
# Plot genes by number of mutations #
#####################################

sort_genes_by_mutations <- function(data, gene_col, sample_col) {
  print("Sorting genes by number of samples with one or more mutations")
  gene_counts <- data %>%
    select(!!as.name(gene_col), !!as.name(sample_col)) %>%
    distinct() %>%
    count(!!as.name(gene_col), sort = T)
  data[, gene_col] <- factor(data[, gene_col], levels = gene_counts[, gene_col])
  return(data)
}


#######################################
# Plot samples by number of mutations #
#######################################

sort_samples_by_mutations <- function(data, gene_col, sample_col) {
  print("Sorting samples by number of mutated genes")
  sample_counts <- data %>%
    select(!!as.name(gene_col), !!as.name(sample_col)) %>%
    distinct() %>%
    count(!!as.name(sample_col), sort = T, .drop = F)
  data[, sample_col] <- factor(data[, sample_col], levels = sample_counts[, sample_col])
  return(data)
}


#######################################################
# Plot samples grouped by mutations in specific genes #
#######################################################

group_samples_by_gene <- function(data, gene_col, sample_col, group) {
  print(paste("Grouping samples by gene(s) in group", group))
  group_list <- unlist(strsplit(group, ","))
  sample_order <- vector()
  # gene_match <- data
  for (gene in group_list) {
    sample_match <- data %>%
      filter(!!as.name(gene_col) == gene) %>%
      select(!!as.name(sample_col)) %>%
      droplevels()
    # 			filter(!(!!as.name(sample_col) %in% sample_order))
    sample_list <- setdiff(levels(sample_match[, sample_col]), sample_order)
    sample_order <- append(sample_order, sample_list)
    print("sample_order")
    print(str(sample_order))
  }
  old_order <- levels(data[, sample_col])
  new_order <- append(sample_order, old_order[!old_order %in% sample_order])
  data[, sample_col] <- factor(data[, sample_col], levels = new_order)
  return(data)
}


########################################################
# Output lines from all genes/samples used in the plot #
########################################################

write_tsv <- function(data, prefix) {
  write.table(data, file = paste(prefix, ".tsv", sep = ""), row.names = F, quote = F, sep = "\t")
}

#####################################################
# Output list of the genes/samples used in the plot #
#####################################################

write_txt <- function(data, prefix) {
  write.table(data, file = paste(prefix, ".txt", sep = ""), row.names = F, col.names = F, quote = F, sep = "\t")
}

######################
# Output plot to PDF #
######################

write_pdf <- function(plotlist, prefix = "gene_table", width = 5, height = 5, relsizeA = 1, relsizeB = 1, ncol = 1, direction = "horizontal", label = opt$label) {
  if (length(plotlist) == 1) {
    print("Plotting 1 plot")
    plot <- plot_grid(plotlist = plotlist)
    save_plot(filename = paste(prefix, ".pdf", sep = ""), plot, base_width = width, base_height = height)
  } else {
    print("Plotting 2 plots")
    if (direction == "vertical") {
      if (is.null(label)) {
        plot <- plot_grid(plotlist = plotlist, ncol = 1, rel_heights = c(relsizeA, relsizeB), align = "v", axis = "rl")
      } else {
        plot <- plot_grid(plotlist = plotlist, ncol = 1, rel_heights = c(relsizeA, relsizeB), align = "v", axis = "rl", labels = c("A", "B"))
      }
    } else {
      if (is.null(label)) {
        plot <- plot_grid(plotlist = plotlist, ncol = 2, rel_widths = c(relsizeA, relsizeB), align = "h", axis = "bt")
      } else {
        plot <- plot_grid(plotlist = plotlist, ncol = 2, rel_widths = c(relsizeA, relsizeB), align = "h", axis = "bt", labels = c("A", "B"))
      }
    }
    save_plot(filename = paste(prefix, ".pdf", sep = ""), plot, base_width = width, base_height = height)
  }
}


####################################
# Prepare legend colours and order #
####################################

make_mutation_plot <- function(data, legendpos = "top", nrow = 2, ncol = 5, vlines = NULL, hlines = NULL, gene_col = gene_col, sample_col = sample_col, cons_col = cons_col) {
  # rename columns for plotting
  data <- data %>%
    select(!!as.name(gene_col), !!as.name(cons_col), !!as.name(sample_col)) %>%
    rename(Effect = !!as.name(cons_col), Gene = !!as.name(gene_col), Sample = !!as.name(sample_col))

  # plot colours

  # E41A1C=red
  # 377EB8=blue
  # 4DAF4A=green
  ## 984EA3=purple
  ## F781BF=magenta
  ## FFFF33=yellow
  # FF7F00=orange

  # colour list

  cols <- c("darkblue", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#F781BF", "#FF7F00", "lightblue1", "grey20", "grey40", "grey60", "grey80", "lightblue3")

  # Order for the legend key

  order <- c(
    "Missense",
    "Nonsense",
    "Splice site",
    "Start/stop codon lost",
    # 		"Start lost",
    # 		"Stop lost",
    "Frameshift",
    "Inframe indel",
    "Synonymous",
    "Other"
  )

  # Check for Effect that is not  in the order list

  diff <- setdiff(data$Effect, order)
  if (length(diff) > 0) {
    order <- append(order, data$Effect)
  }

  # Get colours for the legend

  cols <- cols[1:length(order)]
  legenddf <- data.frame(order, cols)

  # Use only the Effect values present in the data
  legendtext2 <- legenddf[legenddf$order %in% data$Effect, ]

  # 	# now make the Effect levels the same, for plotting
  data$Effect <- factor(data$Effect, levels = legendtext2$order)

  #####################################################################################
  # Data wrangling - allows genes with multiple mutations to be plotted as half tiles #
  #####################################################################################

  # Get distinct gene/sample combos

  singles <- data %>%
    add_count(Gene, Sample) %>%
    filter(n == 1 | n > 2) %>%
    mutate(Effect = replace(as.character(Effect), n > 2, "More than 2 mutations")) %>%
    distinct(Gene, Sample, .keep_all = TRUE) %>%
    select(-n)
  doubles <- data %>%
    add_count(Gene, Sample) %>%
    filter(n == 2) %>%
    select(-n)
  # Get genes mutated more than 2 times per sample
  if (nrow(singles %>% filter(Effect == "More than 2 mutations")) > 0) {
    legendtext2 <- rbind(legendtext2, c(order = "More than 2 mutations", cols = "black"))
  }

  # now make the Effect levels the same, for plotting
  # 	data$Effect <- factor(data$Effect, levels = legendtext2$order)
  singles$Effect <- factor(singles$Effect, levels = legendtext2$order)

  # For genes with 2 mutations, draw half squares
  doubles <- doubles %>%
    group_by(Gene, Sample) %>%
    mutate(id = row_number())
  data.first <- doubles %>% filter(id == "1")
  data.second <- doubles %>% filter(id == "2")
  # Get the x and y axis lables and max. break value

  labels.y <- rev(levels(data$Gene))
  labels.x <- levels(data$Sample)

  breaks.y <- length(labels.y)
  breaks.x <- length(labels.x)

  # Get genes/samples with no mutations to fill in tiles - empty tiles have no outline
  # Need to use coordinates on the axes instead of characters in order to plot
  # half tiles properly:

  # Get all genes and all samples

  gene.index <- rev(levels(data$Gene))
  sample.index <- levels(data$Sample)

  # Assign index to gene

  data.all <- data
  data.all$Gene.index <- match(data$Gene, gene.index)
  data.all$Sample.index <- match(data$Sample, sample.index)

  ## 	data$Tumor_Sample_Barcode_idx <- match(data$Tumor_Sample_Barcode, levels(data$Tumor_Sample_Barcode))
  ## 	data.first$Tumor_Sample_Barcode_idx <- match(data.first$Tumor_Sample_Barcode, levels(data$Tumor_Sample_Barcode)) + 0.25
  ## 	data.second$Tumor_Sample_Barcode_idx <- match(data.first$Tumor_Sample_Barcode, levels(data$Tumor_Sample_Barcode)) - 0.25

  # Assign NA to missing data
  data.all.tmp <- data.all
  data.all.tmp$Gene.index <- factor(data.all.tmp$Gene.index, levels = c(1:length(levels(data$Gene))))
  data.all.tmp$Sample.index <- factor(data.all.tmp$Sample.index, levels = c(1:length(levels(data$Sample))))
  # 	missing.check <- dcast(data.all.tmp, Gene.index ~ Sample.index, value.var = "Effect", drop = FALSE)
  # 	missing <- melt(dcast(data.all.tmp, Gene.index ~ Sample.index, value.var = "Effect", drop = FALSE), id.var = "Gene.index") %>%
  # missing <- melt(dcast(data.all.tmp, Gene ~ Sample.index, value.var = "Effect", drop = FALSE), id.var = "Gene")
  # 	filter(is.na(value)))

  # 	str(data.all.tmp)
  missing <- melt(dcast(data.all.tmp, Gene ~ Sample.index, value.var = "Effect", fun.aggregate = length, drop = FALSE), id.var = "Gene") %>%
    filter(value == 0) %>%
    select(-value) %>%
    rename(Sample.index = variable) %>%
    mutate(Sample.index = as.numeric(Sample.index))
  missing["Effect"] <- NA

  # Now, adjust the sample coordinates for genes with multiple mutations; need to shift half-tiles by 0.25

  data.first$Gene.index <- match(data.first$Gene, gene.index)
  data.second$Gene.index <- match(data.second$Gene, gene.index)
  data.first$Sample.index <- match(data.first$Sample, sample.index) + 0.25
  data.second$Sample.index <- match(data.second$Sample, sample.index) - 0.25

  # Get the indexes for the genes/samples with only 1 mutation per gene
  # and more than 2

  singles$Gene.index <- match(singles$Gene, gene.index)
  singles$Sample.index <- match(singles$Sample, sample.index)
  # 	multiples$Gene.index <- match(multiples$Gene, gene.index)
  # 	multiples$Sample.index <- match(multiples$Sample, sample.index)
  # 	multiples$Gene.index <- factor(multiples$Gene.index, levels = c(1:length(levels(data$Gene))))

  # plot each layer
  plot <- ggplot(singles, aes(x = Sample.index, y = Gene)) +
    geom_tile(aes(fill = Effect, width = 1, height = 1), colour = "white", linewidth = 0.2, alpha = 0.75) +
    geom_tile(data = missing, aes(fill = Effect, width = 1, height = 1), colour = "white", linewidth = 0.2, alpha = 0.75) +
    theme_bw() +
    theme(
      axis.line = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(),
      axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "grey90"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = sizex, family = "mono", color = "black", angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = sizey, family = "sans", face = "italic", color = "black"),
      legend.position = legendpos,
      legend.margin = margin(b = 0, unit = "cm"),
      # below: box around legend
      legend.background = element_blank(),
      legend.key.size = unit(0.25, "cm"), legend.title = element_text(family = "sans", size = 9, face = "bold"), legend.text = element_text(size = 9)
    ) +
    scale_y_discrete(expand = c(0, 0), limits = gene.index) +
    ## 			#scale_y_discrete(breaks = c(1:breaks.y), expand = c(0, 0), labels = labels.y, limits = c(1:breaks.y)) +
    # scale_y_discrete(expand = c(0, 0), breaks = labels.y, labels = labels.y) +
    scale_x_continuous(expand = c(0, 0), breaks = c(1:breaks.x), labels = sample.index) +
    # scale_x_continuous(breaks = c(1:breaks.x), labels = labels.x ,limits = c(1:breaks.x), expand = c(0, 0)) +
    # scale_fill_manual(values = cols, na.value = "grey90", na.translate = FALSE, labels = order, name = "") +
    scale_fill_manual(values = legendtext2$cols, na.value = "grey90", na.translate = FALSE, labels = legendtext2$order, name = "", drop = F) +
    guides(fill = guide_legend(override.aes = list(colour = NULL), nrow = nrow, ncol = ncol, byrow = T))

  print("Checking data.first")
  if (nrow(as.data.frame(data.first)) > 0) {
    plot <- plot + geom_tile(data = data.first, aes(fill = Effect, width = 0.5, height = 1), colour = "white", linewidth = 0.2, alpha = 0.75)
  }
  print("Checking data.second")
  if (nrow(as.data.frame(data.second)) > 0) {
    plot <- plot + geom_tile(data = data.second, aes(fill = Effect, width = 0.5, height = 1), colour = "white", linewidth = 0.2, alpha = 0.75)
  }
  # 	print("Checking multiples")
  # 	if(nrow(as.data.frame(multiples)) > 1) {
  # 		plot <- plot + geom_tile(data = multiples, aes(x = Sample.index, y = Gene, colour = "black"), inherit.aes = FALSE, colour = "white",size = 0.005, alpha = 0.75)
  # 	}
  if (!is.null(vlines) && nrow(vlines) > 0) {
    for (i in 1:nrow(vlines)) {
      plot <- plot + geom_vline(xintercept = vlines[i, 1], linetype = vlines[i, 2], color = "grey10", size = 0.2)
    }
  }
  if (!is.null(hlines) && nrow(hlines) > 0) {
    for (i in 1:nrow(hlines)) {
      plot <- plot + geom_hline(yintercept = hlines[i, 1], linetype = hlines[i, 2], color = "grey10", size = 0.2)
    }
  }
  return(plot)
}


make_frequency_plot <- function(data, gene_col, perc_col, group) {
  max_scale <- round_any(max(data[, perc_col], na.rm = T), 10, f = ceiling)
  scale_breaks <- seq(from = 0, to = max_scale, by = 10)
  max_scale <- max_scale + 1

  plot <- ggplot(data, aes(x = !!as.name(group), y = !!as.name(gene_col), fill = !!as.name(perc_col))) +
    geom_tile(aes(width = 0.5, height = 1), colour = "white", linewidth = 0.2, alpha = 0.75) +
    theme_bw() +
    theme(
      axis.line = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(),
      axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.border = element_blank(), panel.background = element_rect(fill = "white"),
      # uncomment below to add the x-axis labels from discrete scale
      # axis.text.x = element_text(face = "bold", size = 7, angle = 90, hjust=0.5, vjust=0.5),
      axis.text.x = element_blank(),
      legend.position = "right",
      plot.margin = margin(l = -0.1, unit = "cm"),
      legend.margin = margin(l = 0, unit = "cm")
    ) +
    # scale_y_discrete(breaks=c(1:length(names$V1)),labels=names$V1,limits=rev(levels(names$V1))) +
    # scale_y_discrete(breaks=levels(datf1$Gene),expand = c(0, 0),labels=names$V1,limits=rev(levels(datf1$Gene))) +
    # scale_fill_continuous(low="yellow", high="", na.value="white",
    # scale_fill_viridis_c(direction = -1, na.value="white", limits = c(0,21), breaks = c(0,5,10,15,20),
    scale_fill_viridis_c(
      direction = -1, na.value = "white", limits = c(0, max_scale), breaks = scale_breaks,
      guide = guide_colorbar(
        title = "% Cases",
        title.theme = element_text(size = 7),
        barwidth = 0.6, barheight = 5,
        label.theme = element_text(size = 7)
      )
    ) +
    # add sample names
    # scale_x_discrete(breaks=c(1),labels=c("Human"),expand = c(0, 0))
    scale_x_discrete(expand = c(0, 0))

  return(plot)
}


########
# MAIN #
########

# Parse sample and gene lists

samplefiles <- NULL
genefiles <- NULL
label <- NULL

sizex <- 8
sizey <- 8

if (!is.null(opt$sizex)) {
  sizex <- opt$sizex
}

if (!is.null(opt$sizey)) {
  sizey <- opt$sizey
}

if (!is.null(opt$samples)) {
  samplefiles <- as.list(strsplit(opt$samples, ",")[[1]])
}

if (!is.null(opt$genelist)) {
  genefiles <- as.list(strsplit(opt$genelist, ",")[[1]])
}

# Draw plots horizontally or vertically
direction <- "horizontal"
if (isTRUE(opt$vertical)) {
  direction <- "vertical"
}

# Where to draw geom_vline and geom_hline
vline.df <- data.frame(matrix(ncol = 3, nrow = 0))
hline.df <- data.frame(matrix(ncol = 3, nrow = 0))

if (!is.null(opt$vline)) {
  vline.df <- parse_geomlines(opt$vline, opt$vlinetype, opt$vlineplot)
}

if (!is.null(opt$hline)) {
  hline.df <- parse_geomlines(opt$hline, opt$hlinetype, opt$hlineplot)
}

# Legend position
legendpos.list <- as.list(strsplit(opt$legendpos, ",")[[1]])

# list of input files
filelist <- c(opt$fileA)
if (!is.null(opt$fileB)) {
  filelist <- append(filelist, opt$fileB)
}
filelist

# list of plots
plotlist <- vector(mode = "list", length = 1)

file_count <- 0
# Filter everything first, then sort; dplyr doesn't retain factor levels (??)

for (file in filelist) {
  file_count <- file_count + 1

  # Read input file

  data <- read.table(file, header = T, sep = "\t", check.names = F, quote = "", fill = FALSE)
  # colnames(data)<-c("Gene", "Hugo_Symbol", "Main_consequence_VEP", "Tumor_Sample_Barcode")
  # colnames(data)<-c("Gene", "Effect", "Sample")

  data <- shorten_consequence(data)

  # If gene list is format HugoSymbol/ENSGID
  if (isTRUE(opt$geneid_list)) {
    data <- data %>%
      mutate(Hugo_Symbol_Orig = Hugo_Symbol) %>%
      unite("Hugo_Symbol", c(Hugo_Symbol_Orig, Gene), sep = "/", remove = F)
  } else { 
  # if no gene symbol, use Ensembl ID in Gene column
    data <- data %>% mutate(Hugo_Symbol = ifelse(Hugo_Symbol == "-", Gene, Hugo_Symbol))
  }
  data$Gene <- factor(data$Gene, levels = unique(data$Gene))
  data$Hugo_Symbol <- factor(data$Hugo_Symbol, levels = unique(data$Hugo_Symbol))
  data$Tumor_Sample_Barcode <- factor(data$Tumor_Sample_Barcode, levels = unique(data$Tumor_Sample_Barcode))
  data$Main_consequence_VEP <- factor(data$Main_consequence_VEP, levels = unique(data$Main_consequence_VEP))

  if (!is.null(opt$samples)) {
    data <- filter_by_sample(data, samplefiles[[file_count]], "Tumor_Sample_Barcode")
  }

  if (!is.null(opt$genelist)) {
    data <- filter_by_genelist(data, genefiles[[file_count]], "Hugo_Symbol")
  }
  if (!is.null(opt$filter_col)) {
    data <- filter_by_column(data, opt$filter_col)
  }

  if (!is.null(opt$cutoff)) {
    data <- filter_by_minmutations(data, opt$cutoff, "Hugo_Symbol", "Tumor_Sample_Barcode")
  }

  # sort genes by number of samples with one or more mutations
  if (isTRUE(opt$sortbymutations)) {
    data <- sort_genes_by_mutations(data, "Hugo_Symbol", "Tumor_Sample_Barcode")
  }

  if (isTRUE(opt$sortbysample)) {
    data <- sort_by_sample(data, samplefiles[[file_count]], "Tumor_Sample_Barcode")
  }

  if (isTRUE(opt$sortbygenelist)) {
    data <- sort_by_genelist(data, genefiles[[file_count]], "Hugo_Symbol")
  }

  # sort samples by number of mutated genes
  if (isTRUE(opt$sortbysamplemut)) {
    data <- sort_samples_by_mutations(data, "Hugo_Symbol", "Tumor_Sample_Barcode")
  }
  if (!is.null(opt$group)) {
    data <- group_samples_by_gene(data, "Hugo_Symbol", "Tumor_Sample_Barcode", opt$group)
  }

  if (isTRUE(opt$sortbyfrequency)) {
    data <- order_by_freq(data, "Hugo_Symbol", "Tumor_Sample_Barcode", "Main_consequence_VEP")
  }

  if (!is.null(opt$frequency)) {
    if (isTRUE(opt$sortbysidebar)) {
      data.freq <- read.table(opt$frequency, sep = "\t", na.strings = "NA", header = T, check.names = F)
      colnames(data.freq) <- c("Hugo_Symbol", "Percent")
      data$Hugo_Symbol <- factor(data$Hugo_Symbol, levels = data.freq$Hugo_Symbol)
    }
  }
  # 	data <- data %>% group_by(as.character(Hugo_Symbol))
  print("Finished filtering/sorting")

  prefix <- NULL
  if (!is.null(opt$fileB) && file_count == 1) {
    prefix <- paste(opt$prefix, "-A", sep = "")
  } else if (!is.null(opt$fileB) && file_count == 2) {
    prefix <- paste(opt$prefix, "-B", sep = "")
  } else {
    prefix <- opt$prefix
  }
  print(paste("Writing gene list to ", prefix, ".tsv", sep = ""))
  write_tsv(data, prefix)

  write_txt(levels(data$Tumor_Sample_Barcode), paste0(prefix, "_sample_order"))
  write_txt(levels(data$Hugo_Symbol), paste0(prefix, "_gene_order"))

  # make plot and add to list of plots
  vlines <- NULL
  hlines <- NULL
  if (nrow(vline.df) > 0) {
    vlines <- as.data.frame(vline.df %>% filter(plot == file_count) %>% select(-plot))
  }
  if (nrow(hline.df) > 0) {
    hlines <- as.data.frame(hline.df %>% filter(plot == file_count) %>% select(-plot))
  }
  if (file_count <= length(legendpos.list)) {
    legendpos <- legendpos.list[[file_count]]
  } else {
    legendpos <- opt$legendpos
  }
  print("Making plot")
  plotlist[[file_count]] <- make_mutation_plot(data, nrow = opt$nrow, ncol = opt$ncol, legendpos = legendpos, vlines = vlines, hlines = hlines, gene_col = "Hugo_Symbol", sample_col = "Tumor_Sample_Barcode", cons_col = "Main_consequence_VEP")


  if (!is.null(opt$frequency)) {
    data.freq <- read.table(opt$frequency, sep = "\t", na.strings = "NA", header = T, check.names = F)
    colnames(data.freq) <- c("Hugo_Symbol", "Percent")
    data.freq["Group"] <- "Human"
    if (!isTRUE(opt$sortbysidebar)) {
      data.freq <- data %>%
        select(Tumor_Sample_Barcode, Hugo_Symbol) %>%
        left_join(x = ., y = data.freq, by = "Hugo_Symbol")
      # Re-order the levels for plotting
      data.freq$Hugo_Symbol <- factor(data.freq$Hugo_Symbol, levels = rev(levels(data$Hugo_Symbol)))
    } else {
      data.freq$Hugo_Symbol <- factor(data.freq$Hugo_Symbol, levels = rev(data.freq$Hugo_Symbol))
      data$Hugo_Symbol <- factor(data$Hugo_Symbol, levels = rev(data.freq$Hugo_Symbol))
    }
    # 		# Re-order the levels for plotting
    # 		data.freq$Hugo_Symbol <- factor(data.freq$Hugo_Symbol, levels = rev(levels(data$Hugo_Symbol)))

    plotlist[[2]] <- make_frequency_plot(data.freq, "Hugo_Symbol", "Percent", "Group")
  }
}

print(paste("Writing plot to ", opt$prefix, ".pdf", sep = ""))
# if (is.null(opt$fileB)) {
# 	print("making 1 plot")
# 	write_pdf(plotlist, prefix = opt$prefix, width = opt$width, height = opt$height, label = label)
# } else {
# 	ncol = 2
# 	if (direction == "vertical") {
# 		ncol = 1
# 	}
# 	write_pdf(plotlist, prefix = opt$prefix, width = opt$width, height = opt$height, relsizeA = opt$relsizeA, relsizeB = opt$relsizeB, ncol = ncol, direction = direction, label = label)
# }
warnings()

if (length(plotlist) == 2) {
  if (direction == "horizontal") {
    ncol_plot <- 2
  } else if (direction == "vertical") {
    ncol_plot <- 1
  }
  write_pdf(plotlist, prefix = opt$prefix, width = opt$width, height = opt$height, relsizeA = opt$relsizeA, relsizeB = opt$relsizeB, ncol = ncol_plot, direction = direction, label = opt$label)
} else {
  # print(plotlist[[1]])
  write_pdf(plotlist, prefix = opt$prefix, width = opt$width, height = opt$height, label = opt$label)
}



warnings()
