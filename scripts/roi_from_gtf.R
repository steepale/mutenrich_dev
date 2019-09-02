#!/usr/bin/env Rscript

# Script will extract exons of genes from a gtf file a reformat the output in tab seperated format:
# Ensure the gtf file is in the current directory, output file will be saved in current directory

# Usage:
# Rscript exons roi_from_gtf.R Gallus_gallus.Gallus_gallus-5.0.87.gtf

# output example:
#1	30	83	ENSGALG00000045540
#1	125	317	ENSGALG00000045540
#1	1677	1769	ENSGALG00000045540
#1	1977	2024	ENSGALG00000045540
#1	2051	2157	ENSGALG00000045540
#1	2352	2458	ENSGALG00000045540
#1	2548	2665	ENSGALG00000045540

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

#setwd("/Users/Alec/Documents/Bioinformatics/MDV_Project/p0100_music/data")
args = commandArgs(trailingOnly = TRUE)
region <- args[1]
gtffile <- args[2]
outfile <- args[3]


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
df <- dplyr::filter(df, FEATURE == region)

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



