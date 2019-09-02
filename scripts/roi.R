#!/usr/bin/env Rscript

# Script will extract regions of interest from a gtf file a reformat the output in tab seperated format:
# 
# Given a gtf file, this script will extract a region of interest
#
#  USAGE:
#
#       Rscript roi.R \
#         --roi <GENOMIC_FEATURE_OF_INTEREST> \
#         --ann <PATH_TO_ANNOTATION_FILE> \
#         --output <PATH_TO_OUTPUT_FILE>
#
# ARGUMENTS:
#
#       roi: A genomic feature that will capture the region of interest
#            in from a gtf file. Must be one of: 
#                        CDS, exon, five_prime_utr, gene, start_codon, 
#                        stop_codon, three_primer_utr, transcript, or intron
#       ann: The path to the annotation file
#       output: The path to the output ROI file
#
####### parse command line arguments ------------------------------------
# loading the library seems a bit cleaner here
library("optparse")

option.list = list(
  make_option(c("-r", "--roi"), type = "character", default = NULL,
              help = "The genomic region of interest. (CDS, exon, five_prime_utr, gene, start_codon, stop_codon, three_primer_utr, transcript, or intron)"),
  make_option(c("-a", "--ann"), type = "character", default = NULL,
              help = "The input annotation file (e.g. genome.gtf, annotation.bed)"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "The path to the output ROI file")
)

opt.parser <- OptionParser(option_list = option.list)
opt <- parse_args(opt.parser)

# Collect the opts
roi <- opt$roi
gtffile <- opt$ann
outfile <- opt$output

#### main ---------------------------------------------------------------

# Useful tutorial:
# https://bioconductor.statistik.tu-dortmund.de/packages/3.3/bioc/vignettes/ensembldb/inst/doc/ensembldb.html

# Load dependencies
pacs...man <- c("dplyr","tibble","stringr", "data.table","GenomicRanges","GenomicFeatures","Biostrings","BSgenome","AnnotationHub","stringdist",'tidyr','biomaRt','org.Gg.eg.db','ggplot2')

# Load the packages
lapply(pacs...man, FUN = function(X) {
        do.call("library", list(X)) 
})

# Extract the apprpriate GTF file from ensembl
# ftp://ftp.ensembl.org/pub/release-87/gtf/gallus_gallus/Gallus_gallus.Gallus_gallus-5.0.87.gtf.gz

# To extract introns, follow instructions from here: https://www.biostars.org/p/13290/

#setwd("/Users/Alec/Documents/Bioinformatics/MDV_Project/p0100_music/data")

# Get today's date and user name
date <- format(Sys.Date(), "%Y%m%d")
user <- "steep"

# Load the gtf file
?read.table()
df <- read.table(gtffile, sep = '\t', header = FALSE, comment.char="#")
# Label the columns
set.seed(123)
df <- as_tibble(df)
dim(df)
names(df) <- c('CHROM','SOURCE','FEATURE','START','STOP','SCORE','STRAND','FRAME','ATTRIBUTE')

# Select necessary columns (learn to spell)
df <- df %>% dplyr::select(CHROM,FEATURE,START,STOP,ATTRIBUTE)

# Grab the region of interest (e.g. exons, introns)
df <- dplyr::filter(df, FEATURE == roi)

# Split the ATTRIBUTE column
df <- df %>% separate(ATTRIBUTE, 
                      c("GENE_ID", "GENE_VERSION","TX_ID","TX_VERSION", "EXON_NUMBER"),
                     extra = "drop",
                      sep = ";")

# Remove the '*_* ' from the columns
df$GENE_ID <- str_replace(df$GENE_ID, 'gene_id ','')
df$GENE_VERSION <- str_replace(df$GENE_VERSION, 'gene_version ','')
df$TX_ID <- str_replace(df$TX_ID, 'transcript_id ','')
df$TX_VERSION <- str_replace(df$TX_VERSION, 'transcript_version ','')
df$EXON_NUMBER <- str_replace(df$EXON_NUMBER, 'exon_number ','')

# Add the gene symbol
df$SYMBOL <- mapIds(org.Gg.eg.db, df$GENE_ID, "SYMBOL", "ENSEMBL")

# Select the appropriate columns
df <- df %>% dplyr::select(CHROM,START,STOP,GENE_ID)

# Dedup the df
df <- df[!duplicated(df),]

# Save the df to a file
write.table(df, file = paste0(outfile), sep ='\t', 
            quote = FALSE, row.names = FALSE, col.names = FALSE)



