###############################################################
# Script to process amplicon sequencing data with DADA2       #
# Based on DADA2 tutorial by Benjamin Callahan                #
# Data: Mice - Antibiotic experience                          #
# By: ArrietaLab - University of Calgary                      #
# Author: Emily Mercer                                        # 
# Dates: 3-5 April 2023                                       #
# Location: Oswaldo-Cruz Institute - Rio de Janeiro (Brazil)  #
###############################################################

## 1.1 Install Packages ####

# Install DADA2 package

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager") # Installs BiocManager, enabling you to download packages from the BioConductor repository (similar to CRAN)
#library(BiocManager) # Loads BiocManager package

#BiocManager::install("dada2") # Install DADA2 package from the BioConductor repository (only required once)
library(dada2); packageVersion("dada2") # Loads DADA2 package (Callahan et al., 2016)

# For any installation issues, visit: https://benjjneb.github.io/dada2/dada-installation.html

# Install tidyverse 

#install.packages("tidyverse") # Installs tidyverse package (only required once)
library(tidyverse); packageVersion("tidyverse") # Loads tidyverse package for tidier coding (Wickham et al., 2019)

## 1.2 Set Up R Environment ####

# Prepare R environment
rm(list = ls(all = TRUE)) # Clears your global environment

# Set working directory
setwd("~/Downloads/data/mouse_demultiplexed/")
getwd() # Check that your working directory is correctly set

# Set path to folder where sequencing files are located
path <- "~/Downloads/data/mouse_demultiplexed/" # Creates path object
list.files(path, pattern = "fastq") # Lists fastq sequencing files at path location

# Save forward and reverse sequence file names into separate objects for downstream processing
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE)) # Forward sequences
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE)) # Reverse sequences

## 1.3 Inspect Sequence Quality ####

# Inspect sequence quality
plotQualityProfile(fnFs[1:4]) # Forward sequences quality for the first 4 sequences
plotQualityProfile(fnRs[1:4]) # Reverse sequences quality for the first 4 sequences

plotQualityProfile(fnFs) # Forward sequences quality for all forward sequences (may take a while to run)
plotQualityProfile(fnRs) # Reverse sequences quality for all forward sequences (may take a while to run)

# Overall sequence quality of all sequences
plotQualityProfile(fnFs, aggregate = TRUE) # Forward sequences quality aggregated in one plot
plotQualityProfile(fnRs, aggregate = TRUE) # Forward sequences quality aggregated in one plot

# Save sequence plots (optional)
jpeg(file = "Forward Sequences 1-4.jpeg", width = 2000, height = 2000, units = "px", res = 300)
plotQualityProfile(fnFs[1:4]) # Forward sequences quality 
dev.off()

## 1.4 Filter & Trim ####

# Create a folder to store filtered and trimmed sequences

# Create a list of sample names
sample.names <- strsplit(basename(fnFs), "_") %>% # Creates a list of character vectors for each sample split at the underscore
  sapply(`[`, 1) # sapply uses the first element in each character vector as the sample name
sample.names # View sample.names to ensure our list is correct

# Pastes samples names and filtered specification together to make new file names, then stores them in the new "filtered" folder at our path destination.
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz")) # Forward filtered sequences
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz")) # Reverse filtered sequences

names(filtFs) <- sample.names # Applies the sample names as derived in the previous step to the filtered files
names(filtRs) <- sample.names # Applies the sample names as derived in the previous step to the filtered files

# Filter and trim
out <- filterAndTrim(fwd = fnFs, filt = filtFs, rev = fnRs, filt.rev = filtRs, truncLen = c(240,160),
                     maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE) # On Windows, set multithread = FALSE
head(out)

# Play around with truncation length and expected errors
out2 <- filterAndTrim(fwd = fnFs, filt = filtFs, rev = fnRs, filt.rev = filtRs, truncLen = c(240,155),
                      maxN = 0, maxEE = c(2,5), truncQ = 2, rm.phix = TRUE,
                      compress = TRUE, multithread = TRUE) 
head(out2)
# In this case, the least amount of reads are lost in the first output - we will therefore use it in the subsequent steps.

# After you filter and trim, it is good practice to go back and make sure the quality is still good using the plotQualityProfile function.
plotQualityProfile(filtFs) # Check quality of filtered forward sequences
plotQualityProfile(filtRs) # Check quality of filtered reverse sequences

# Compare unfiltered and filtered files
plotQualityProfile(c(fnFs[1],filtFs[1])) # Compares first forward sequence raw vs. filtered
plotQualityProfile(c(fnRs[1],filtRs[1])) # Compares first reverse sequence raw vs. filtered

# Assess the percentage of reads lost with filtering and trimming
percentage_reads <- function(output){
  as.data.frame(output) %>% mutate(percentage = 100*reads.out/reads.in)
} # Creates a function for calculating percentage reads using output objects
percentage_reads(out) # Shows you the percentage of reads maintained after filtering and trimming

## 1.5 Error Rates ####

# Learn the error rates
errF <- learnErrors(filtFs, multithread = TRUE) # Creates list with forward error rates
errR <- learnErrors(filtRs, multithread = TRUE) # Creates list with reverse error rates

# Visualize estimated error rates
plotErrors(errF, nominalQ = TRUE) # Forward error rates

## 1.6 Sample Inference ####

# Determine the number of sequence variants from the total unique sequences in each sample
dadaFs <- dada(filtFs, err = errF, multithread = TRUE) # Creates list with sample inference for each forward sequence
dadaRs <- dada(filtRs, err = errR, multithread = TRUE) # Creates list with sample inference for each reverse sequence

# Inspect the returned dada class object for the first forward and reverse sequences
dadaFs[[1]] 
dadaRs[[1]]

## 1.7 Construct Sequence Table ####

# First, we need to merge the paired reads (forward and reverse sequences) to generate denoised sequences.
# This requires at least 12 basepairs to overlap between the forward and reverse reads.
merged <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE) # Creates list of data frames for each sample

# Construct sequence table from the merged reads
# Also known as an ASV or amplicon sequence variant table
seqtab <- makeSequenceTable(merged) # Constructs sequence table (matrix)
dim(seqtab) # Check dimensions - the first number should be the number of samples you have, the second is the number of unique ASVs

table(nchar(getSequences(seqtab))) # Inspect distribution of sequence lengths

# Remove chimeras
# These are single sequences generated from multiple parent sequences and are considered an artifact of sequencing technologies.
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim) # Check dimensions

1-sum(seqtab.nochim)/sum(seqtab) # Examine percent lost to chimeric sequences 

## 1.8 Track Reads Through Pipeline ####

# The next quality control step is to track reads through the pipeline by generating a table of how many reads are lost after each step in the pipeline.

# Create function to tell you the abundance of each sequence
getN <- function(N){
  sum(getUniques(N))
} # getUniques function looks at object, gets sequence name, and tells you how abundant it is

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(merged, getN), rowSums(seqtab.nochim))
# cbind combines the information derived from all of these functions into an overview table
# When only processing one sample, remove sapply (e.g. replace sapply(dadaFs, getN) with getN(dadaFs)). sapply allows you to apply a function to every object in a list.

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim") # Naming columns of overview table
rownames(track) <- sample.names # Naming rows of overview table

track # View track to see how the number of reads change with each pipeline step

# We can see that most of our reads were maintained from step to step, apart from removing chimeras.

## 1.9 Assign Taxonomy ####

# Download the SILVA 138.1 prokaryotic SSU taxonomic training data formatted for DADA2 from https://zenodo.org/record/4587955#.ZBVL8uzMJa0
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread = TRUE) # Assigns taxonomy to species level
# This step may take a while to run

taxa.print <- taxa 
rownames(taxa.print) <- NULL # Removing sequence rownames for display only
head(taxa.print) # Examine taxonomic assignments 

# Examine how many sequences are identified at different taxonomic levels
table(taxa[,2], useNA = "always") # Tells you how many taxa were assigned to each phylum (taxonomic level 2)

# Can add decreasing = TRUE to get most abundant at the top
sort(table(taxa[,5], useNA = "always"), decreasing = TRUE) # Tells you how many taxa were assigned to each family (taxonomic level 5)

## 1.10 Create phyloseq Object ####

# Install phyloseq
#BiocManager::install("phyloseq") # Installs phyloseq package (only required once)
library(phyloseq); packageVersion("phyloseq") # Loads phyloseq package

# Load metadata
mouse_metadata <- read.csv("data/mouse_metadata.csv") 
rownames(mouse_metadata) <- mouse_metadata$Sample_ID # Change row names to reflect sample IDs

# Create phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               sample_data(mouse_metadata),
               tax_table(taxa))
ps # Examine object

# Save phyloseq object
saveRDS(ps,"data/ps.rds")

# Save the workspace
save.image("BioinfoWorkshop_DADA2.RData")
