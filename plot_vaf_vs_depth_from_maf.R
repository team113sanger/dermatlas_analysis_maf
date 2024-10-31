suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(tidyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(stringr))
suppressMessages(library(GetoptLong))

genefile <- NULL
colname <- NULL
suffix <- NULL
samplefile <- NULL

type <- "pdf"
ncol_r <- 4
width_r <- 10
height_r <- 7
ncol <- 4
width <- 10
height <- 7

mut_order <- c("C>A/G>T", "C>G/G>C", "C>T/G>A", "T>A/A>T", "T>C/A>G", "T>G/A>C", "MNV", "Indel")
colours <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#F781BF", "#999999", "#525252")
## colours <- colorRampPalette(brewer.pal(8, "Dark2"))(16)

spec <- "

Create VAF vs. depth plots for the top 12 recurrent genes and sites found in
the variants file, or a subset of genes specified from a list of genes.

  Usage: Rscript plot.R [options]

  Options:
    <file=s> Variants file [required]
    <genefile=s> List of genes to include.
    <samplefile=s> List of samples to include.
    <colname=s> Column in variants file with gene name matching the gene list. Required with --genelist
    <suffix=s> Suffix to be added to output files when using a custom gene list. Required with --genelist
    <noindels> Ingore indels (only SNVs and MNVs will be considered)
    <indels_only> Plot indels only
    <protein_alt> Exclude variants with main consequence synonymous, stop/start retained, splice_region
    <germline> Use n_depth and VAF_norm for AF vs Depth plots

  Options for AF vs Depth plots:
    <width=f> Width (inches) of AF vs Depth by sample plots
    <height=f> Height (inches) of AF vs Depth by sample plots
    <ncol=i> Number of columns for AF vs Depth plot

"

GetoptLong(spec, template_control = list(opt_width = 21))

# Functions -------------------------------------------------------------------

# Get the top recurrently mutated sites

top_sites <- function(df) {
  top <- df %>%
    select(Gene, Hugo_Symbol, Chromosome, POS_VCF, Reference_Allele, Tumor_Seq_Allele2) %>%
    add_count(Gene, Hugo_Symbol, Chromosome, POS_VCF, Reference_Allele, Tumor_Seq_Allele2) %>%
    # arrange(-n, Chromosome, POS_VCF) %>%
    arrange(-n, Gene, POS_VCF) %>%
    distinct() %>%
    unite("Gene_Name_Site", c(Hugo_Symbol, Gene, POS_VCF, Tumor_Seq_Allele2), sep = "/", remove = F) %>%
    unite("Gene_Name", c(Hugo_Symbol, Gene), sep = "/", remove = F) %>%
    slice(1:12) %>%
    data.frame()

  return(top)
}


# Get the top recurrently mutated genes

top_genes <- function(df) {
  top <- df %>%
    select(Gene, Hugo_Symbol, Tumor_Sample_Barcode) %>%
    distinct() %>%
    select(Gene, Hugo_Symbol) %>%
    add_count(Gene, Hugo_Symbol) %>%
    arrange(-n, Gene, Hugo_Symbol) %>%
    distinct() %>%
    select(Gene, Hugo_Symbol) %>%
    unite("Gene_Name", c(Hugo_Symbol, Gene), sep = "/", remove = F) %>%
    slice(1:12) %>%
    data.frame()

  return(top)
}

# Summarize mutation types
# Cannot plot every type of indel of MNV; classyify as MNV or indel instead

process_variants <- function(df) {
  df <- df %>%
    filter(str_detect(Keep, "keep")) %>%
    select(Gene, Hugo_Symbol, Chromosome, POS_VCF, Reference_Allele, Tumor_Seq_Allele2, Variant_Type, VAF_tum, t_depth, VAF_norm, n_depth, Tumor_Sample_Barcode, Main_consequence_VEP, Consequence, Keep) %>%
    unite("Gene_Name", c(Hugo_Symbol, Gene), sep = "/", remove = F) %>%
    mutate(reflen = str_length(Reference_Allele), altlen = str_length(Tumor_Seq_Allele2)) %>%
    mutate(Reference_Allele_c = ifelse(Variant_Type == "SNP",
      case_when(
        Reference_Allele == "A" ~ "T",
        Reference_Allele == "T" ~ "A",
        Reference_Allele == "C" ~ "G",
        Reference_Allele == "G" ~ "C"
      ), Reference_Allele
    )) %>%
    mutate(Tumor_Seq_Allele2_c = ifelse(Variant_Type == "SNP",
      case_when(
        Tumor_Seq_Allele2 == "A" ~ "T",
        Tumor_Seq_Allele2 == "T" ~ "A",
        Tumor_Seq_Allele2 == "C" ~ "G",
        Tumor_Seq_Allele2 == "G" ~ "C"
      ), Tumor_Seq_Allele2
    )) %>%
    mutate(Mutation = case_when(
      (Reference_Allele == "-" | Tumor_Seq_Allele2 == "-") ~ "Indel",
      reflen == 1 & altlen == 1 & (Reference_Allele == "C" | Reference_Allele == "T") ~ paste0(Reference_Allele, ">", Tumor_Seq_Allele2, "/", Reference_Allele_c, ">", Tumor_Seq_Allele2_c),
      reflen == 1 & altlen == 1 & (Reference_Allele == "G" | Reference_Allele == "A") ~ paste0(Reference_Allele_c, ">", Tumor_Seq_Allele2_c, "/", Reference_Allele, ">", Tumor_Seq_Allele2),
      reflen == altlen ~ "MNV",
      reflen != altlen ~ "Indel"
    )) %>%
    mutate(Direction = case_when(
      (Reference_Allele == "-" | Tumor_Seq_Allele2 == "-") ~ "NA",
      reflen == 1 & altlen == 1 & (Reference_Allele == "C" | Reference_Allele == "T") ~ "Forward",
      reflen == 1 & altlen == 1 & (Reference_Allele == "G" | Reference_Allele == "A") ~ "Reverse",
      TRUE ~ "NA"
    )) %>%
    # 		mutate(Type = paste0(Reference_Allele, ">", Tumor_Seq_Allele2)) %>%
    select(-reflen, -altlen, -Reference_Allele_c, -Tumor_Seq_Allele2_c, -Keep)

  return(df)
}


# Create a ggplot object for VAF vs depth, by sample

make_af_vs_d_plots <- function(df, ncol, type, depth_col, vaf_col) {
  if (type == "snv") {
    colours <- brewer.pal(n = 8, name = "Dark2")
  } else {
    colours <- "black"
  }
  plot <- ggplot(df, aes(x = .data[[depth_col]], y = .data[[vaf_col]], colour = Mutation)) +
    geom_point(shape = 1, size = 2, alpha = 0.60, stroke = 1) +
    theme_bw() +
    facet_wrap(~Tumor_Sample_Barcode, ncol = ncol, drop = F) +
    scale_y_log10(breaks = c(0.01, 0.05, 0.1, 0.2, 0.4, 0.8)) +
    scale_x_log10(breaks = c(5, 10, 20, 40, 100, 300)) +
    ylab("Variant Allele Fraction") +
    xlab("Read Depth") +
    # scale_fill_manual(values = colours, drop = F) +
    scale_color_manual(values = colours, drop = F) +
    theme(
      strip.text = element_text(size = 8),
      legend.text = element_text(size = 9),
      axis.title = element_text(size = 10),
      axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 9)
    ) +
    geom_hline(aes(yintercept = .05), linetype = "dashed", colour = "grey50") +
    geom_hline(aes(yintercept = .10), linetype = "dashed", colour = "grey50") +
    geom_hline(aes(yintercept = .20), linetype = "dashed", colour = "grey50") +
    geom_vline(aes(xintercept = 10), linetype = "dashed", colour = "grey50") +
    geom_vline(aes(xintercept = 20), linetype = "dashed", colour = "grey50") +
    geom_vline(aes(xintercept = 50), linetype = "dashed", colour = "grey50") +
    geom_vline(aes(xintercept = 100), linetype = "dashed", colour = "grey50") +
    geom_vline(aes(xintercept = 300), linetype = "dashed", colour = "grey50")
  # theme(strip.text = element_text(size=7), legend.text = element_text(size=7), axis.title = element_text(size=8), axis.text = element_text(size=6))
  # facet_wrap(~Sample, ncol=4, scales = "free_y") + scale_y_log10()
  return(plot)
}


# Create a ggplot object for VAF vs depth, by gene

make_plots <- function(df, ncol, depth_col, vaf_col, labels = NULL) {
  plot <- ggplot(df, aes(x = .data[[depth_col]], y = .data[[vaf_col]], colour = Mutation)) +
    geom_point(shape = 1, size = 2, alpha = 0.70, stroke = 2) +
    theme_bw()

  if (is.null(labels)) {
    plot <- plot + facet_wrap(~Gene_Name, ncol = ncol, drop = F)
  } else {
    plot <- plot + facet_wrap(~Gene_Name_Site, labeller = as_labeller(labels), ncol = ncol, drop = F)
  }

  plot <- plot +
    scale_y_log10(breaks = c(0.01, 0.05, 0.1, 0.2, 0.4, 0.8)) +
    scale_x_log10(breaks = c(5, 10, 20, 40, 100, 300)) +
    ylab("Variant Allele Fraction") +
    xlab("Read Depth") +
    # scale_fill_manual(values = colours, drop = F) +
    scale_color_manual(values = colours, drop = F) +
    theme(
      strip.text = element_text(size = 8),
      legend.text = element_text(size = 9),
      axis.title = element_text(size = 10),
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9)
    ) +
    geom_hline(aes(yintercept = .05), linetype = "dashed", colour = "grey50") +
    geom_hline(aes(yintercept = .10), linetype = "dashed", colour = "grey50") +
    geom_hline(aes(yintercept = .20), linetype = "dashed", colour = "grey50") +
    geom_vline(aes(xintercept = 10), linetype = "dashed", colour = "grey50") +
    geom_vline(aes(xintercept = 20), linetype = "dashed", colour = "grey50") +
    geom_vline(aes(xintercept = 50), linetype = "dashed", colour = "grey50") +
    geom_vline(aes(xintercept = 100), linetype = "dashed", colour = "grey50") +
    geom_vline(aes(xintercept = 300), linetype = "dashed", colour = "grey50")
  # theme(strip.text = element_text(size=7), legend.text = element_text(size=7), axis.title = element_text(size=8), axis.text = element_text(size=6))
  # facet_wrap(~Sample, ncol=4, scales = "free_y") + scale_y_log10()

  return(plot)
}


# Create a ggplot object for VAF vs depth, by sample

make_af_vs_d_barplots <- function(df, ncol) {
  data <- df %>%
    select(Tumor_Sample_Barcode, Mutation, Direction) %>%
    group_by(Tumor_Sample_Barcode, Mutation, Direction) %>%
    mutate(Count = n()) %>%
    distinct() %>%
    data.frame()
  data$Direction <- factor(data$Direction, levels = c("Reverse", "Forward"))
  colnames(data) <- c("Sample", "Mutation", "Direction", "Count")
  plot <- ggplot(data, aes(x = Mutation, y = Count, fill = Direction)) +
    geom_bar(stat = "identity") +
    facet_wrap(~Sample, ncol = ncol, scales = "free_y", drop = F) +
    labs(x = "", y = "Number of mutations") +
    scale_fill_manual(values = c("#56B4E9", "#66CC99")) +
    theme(
      legend.position = "none",
      axis.text = element_text(color = "black", size = 10),
      axis.text.x = element_text(angle = 90, hjust = 1, size = 9, face = "bold", family = "mono"),
      axis.ticks.y = element_blank(),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 7)
    )
}


# Create a ggplot object for VAF vs depth, by sample, proportion

make_af_vs_d_barplots_perc <- function(df, ncol) {
  data <- df %>%
    select(Tumor_Sample_Barcode, Mutation, Direction) %>%
    group_by(Tumor_Sample_Barcode, Mutation, Direction) %>%
    mutate(Count = n()) %>%
    distinct() %>%
    data.frame()
  data$Direction <- factor(data$Direction, levels = c("Reverse", "Forward"))
  colnames(data) <- c("Sample", "Mutation", "Direction", "Count")
  plot <- ggplot(data, aes(x = Mutation, y = Count, fill = Direction)) +
    geom_bar(stat = "identity", position = "fill") +
    facet_wrap(~Sample, ncol = ncol, scales = "free_y", drop = F) +
    labs(x = "", y = "Proportion of mutations") +
    scale_fill_manual(values = c("#56B4E9", "#66CC99")) +
    theme(
      legend.position = "none",
      axis.text = element_text(color = "black", size = 10),
      axis.text.x = element_text(angle = 90, hjust = 1, size = 9, face = "bold", family = "mono"),
      axis.ticks.y = element_blank(),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 7)
    )
}

# Make plot name based on user options

make_file_name <- function(prefix, type, suffix, noindels, indels_only, protein_alt) {
  outfile <- paste(prefix, type, sep = "_")
  # tsvfile <- paste0("top_recurrently_mutated_", type)
  if (!is.null(suffix)) {
    outfile <- paste(outfile, suffix, sep = "_")
  }
  if (type != "samples") {
    if (isTRUE(noindels)) {
      outfile <- paste(outfile, "noIndels", sep = "_")
    } else if (isTRUE(indels_only)) {
      outfile <- paste(outfile, "indelsOnly", sep = "_")
    }
    if (isTRUE(protein_alt)) {
      outfile <- paste(outfile, "protAlt", sep = "_")
    }
  }

  return(outfile)
}


# Print out plot

print_plot <- function(plot, plot_name, width = width, height = height) {
  pdf(file = paste0(plot_name, ".pdf"), width = width, height = height)
  print(plot)
  dev.off()
}

# Write out TSV file

write_var_table <- function(df, filename) {
  write.table(df, file = paste0(filename, ".tsv"), row.names = F, sep = "\t", quote = F)
}

# Main ------------------------------------------------------------------------

# Read in gene list file if specified

if (!is.null(genefile)) {
  if (is.null(colname) || is.null(suffix)) {
    stop("Specify --colname and --suffix with --genefile")
  }
  print(paste("Reading gene list", genefile))
  genelist <- read.table(genefile, header = F, sep = "\t", stringsAsFactors = F, quote = "")
} else if (isTRUE(noindels) && isTRUE(indels_only)) {
  stop("Cannot select both --noIndels and --indels_only")
}

# Read in sample list file if specified

if (!is.null(samplefile)) {
  print(paste("Reading sample list", samplefile))
  samplelist <- read.table(samplefile, header = F, sep = "\t", stringsAsFactors = F, quote = "")
}

# Read variants file

df <- read.table(file, header = T, sep = "\t", stringsAsFactors = F, quote = "", fill = FALSE)

# Keep original calls to plot VAF vs depth

df_all <- df

all_samples <- ""
if (!is.null(samplefile)) {
  all_samples <- sort(samplelist$V1)
  df <- df %>% filter(Tumor_Sample_Barcode %in% samplelist$V1)
  df_all <- df_all %>% filter(Tumor_Sample_Barcode %in% samplelist$V1)
} else {
  all_samples <- sort(unique(df_all$Tumor_Sample_Barcode))
}

# Process the variants

df_all <- process_variants(df_all)
df_snv <- df_all %>% filter(!(Variant_Type %in% c("DEL", "INS")))
df_indel <- df_all %>% filter(Variant_Type %in% c("DEL", "INS"))
df_snv_only <- df_all %>% filter(Variant_Type == "SNP")

df_snv$Tumor_Sample_Barcode <- factor(df_snv$Tumor_Sample_Barcode, levels = all_samples)
df_indel$Tumor_Sample_Barcode <- factor(df_indel$Tumor_Sample_Barcode, levels = all_samples)
df_snv_only$Tumor_Sample_Barcode <- factor(df_snv_only$Tumor_Sample_Barcode, levels = all_samples)

# Make a VAF vs Read Depth plot for SNV/MNV

if (nrow(df_snv) > 0) {
  if (isTRUE(germline)) {
    snv_plot <- make_af_vs_d_plots(df_snv, ncol, "snv", "n_depth", "VAF_norm")
  } else {
    snv_plot <- make_af_vs_d_plots(df_snv, ncol, "snv", "t_depth", "VAF_tum")
  }
  snv_filename <- make_file_name("AF_vs_depth_snv", "samples", suffix, noindels, indels_only, protein_alt)
  print_plot(snv_plot, snv_filename, width, height)
} else {
  print("No SNV/MNV found. Skipping AF_vs_depth_svn plot")
}

# Make a VAF vs Read Depth plot for indels

if (nrow(df_indel) > 0) {
  if (isTRUE(germline)) {
    indel_plot <- make_af_vs_d_plots(df_indel, ncol, "indel", "n_depth", "VAF_norm")
  } else {
    indel_plot <- make_af_vs_d_plots(df_indel, ncol, "indel", "t_depth", "VAF_tum")
  }
  indel_filename <- make_file_name("AF_vs_depth_indel", "samples", suffix, noindels, indels_only, protein_alt)
  print_plot(indel_plot, indel_filename, width, height)
} else {
  print("No indels found. Skipping AF_vs_depth_indel plot")
}

# Make a bar plot for SNV types
if (nrow(df_snv_only) > 0) {
  snv_barplot <- make_af_vs_d_barplots(df_snv_only, ncol)
  snv_bar_filename <- make_file_name("mutation_types_barplot", "samples", suffix, noindels, indels_only, protein_alt)
  print_plot(snv_barplot, snv_bar_filename, width, height)
  snv_barplot_perc <- make_af_vs_d_barplots_perc(df_snv_only, ncol)
  snv_barperc_filename <- make_file_name("mutation_types_barplot_proportion", "samples", suffix, noindels, indels_only, protein_alt)
  print_plot(snv_barplot_perc, snv_barperc_filename, width, height)
} else {
  print("No SNVs found. Skipping mutation types barplot")
}


# Plot recurrent genes
# Set levels for Chrom

# chr_order <- df %>%
# 	filter(! df$Chromosome %in% c('X', 'Y')) %>%
# 	select(Chromosome) %>%
# 	unique() %>%
# 	arrange()
#
#
# chr_order <- c(as.character(chr_order$Chromosome), "X", "Y")
# df$Chromosome <- factor(df$Chromosome, levels = chr_order)

# Get subset of genes, if gene list provided

if (exists("genelist")) {
  print(paste("Checking genelist"))
  df <- subset(df, df[, colname] %in% genelist$V1)
}

# Get subset of protein-altering variants, if specified

if (isTRUE(protein_alt)) {
  print("Using protein-altering variants only")
  df <- df %>% filter(Keep == "keep-PA")
}

# Process the variants

df <- process_variants(df)

# Exclude indels, if specified

if (isTRUE(noindels)) {
  print("Excluding indels")
  df <- df %>% filter(Mutation != "Indel")
} else if (isTRUE(indels_only)) {
  print("Including indels only")
  df <- df %>% filter(Mutation == "Indel")
}

# Get the top recurrent sites

top_12_sites <- top_sites(df)
# print("top_12_sites")
# write.table(file="test.tsv", top_12_sites, quote=F, sep="\t")
# Get the original rows for the top sites

sites_df <- df %>%
  semi_join(top_12_sites, by = c("Gene", "Chromosome", "POS_VCF", "Reference_Allele", "Tumor_Seq_Allele2")) %>%
  unite("Gene_Name_Site", c(Hugo_Symbol, Gene, POS_VCF, Tumor_Seq_Allele2), sep = "/", remove = F)
# sites_df
# Get the gene order for plotting

# gene_order <- unique(top_12_sites$Gene_Name)
gene_order <- top_12_sites$Gene_Name_Site

# Set the gene order and mutation levels for plotting
sites_df$Gene_Name_Site <- factor(sites_df$Gene_Name_Site, levels = gene_order)
sites_df$Mutation <- factor(sites_df$Mutation, levels = mut_order)

# Make a VAF vs Read Depth plot
facet_labels <- top_12_sites$Gene_Name
names(facet_labels) <- gene_order
# print(facet_labels)

if(isTRUE(germline)) {
  sites_plot <- make_plots(sites_df, ncol_r, "n_depth", "VAF_norm", facet_labels)
} else {
  sites_plot <- make_plots(sites_df, ncol_r, "t_depth", "VAF_tum", facet_labels)
}

sites_filename <- make_file_name("AF_vs_depth_recurrent", "sites", suffix, noindels, indels_only, protein_alt)
print_plot(sites_plot, sites_filename, width_r, height_r)

# Print out the correspoding TSV file
sites_tsv_filename <- make_file_name("top_recurrently_mutated", "sites", suffix, noindels, indels_only, protein_alt)
write_var_table(sites_df, sites_tsv_filename)

# Get the top recurrently mutated genes
top_12_genes <- top_genes(df)

print("Top recurrently mutated genes")
# head(top_12_genes, 20)

# Get the original rows for the top genes

genes_df <- df %>% semi_join(top_12_genes, by = c("Gene", "Hugo_Symbol"))

# Get the gene order for plotting

# gene_order2 <- unique(top_12_genes$Gene_Name)
gene_order2 <- top_12_genes$Gene_Name

# Set the gene order and mutation levels for plotting

genes_df$Gene_Name <- factor(genes_df$Gene_Name, levels = gene_order2)
genes_df$Mutation <- factor(genes_df$Mutation, levels = mut_order)

# Make a VAF vs Read Depth plot
if(isTRUE(germline)) {
  genes_plot <- make_plots(genes_df, ncol_r, "n_depth", "VAF_norm")
} else {
  genes_plot <- make_plots(genes_df, ncol_r, "t_depth", "VAF_tum")
}

genes_filename <- make_file_name("AF_vs_depth_recurrent", "genes", suffix, noindels, indels_only, protein_alt)
print_plot(genes_plot, genes_filename, width_r, height_r)

# Print out the correspoding TSV file
genes_tsv_filename <- make_file_name("top_recurrently_mutated", "genes", suffix, noindels, indels_only, protein_alt)
write_var_table(genes_df, genes_tsv_filename)

warnings()
# sessionInfo()
