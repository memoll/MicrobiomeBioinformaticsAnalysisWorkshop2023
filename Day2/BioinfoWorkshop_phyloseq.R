###############################################################
# Introduction to phyloseq & data exploration                 #
# Data: Mice - Antibiotic experience                          #
# By: ArrietaLab - University of Calgary                      #
# Author: Mackenzie Gutierrez                                 # 
# Dates: 3-5 April 2023                                       #
# Location: Oswaldo-Cruz Institute - Rio de Janeiro (Brazil)  #
###############################################################

# Installation and set up ####
# We will be using several packages alongside phyloseq. First we need to install these packages (if they are not already installed).
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager") # Installs BiocManager, a package that facilitates the download of packages from the BioConductor repository.

#BiocManager::install("phyloseq") # Installs the phyloseq package (McMurdie & Holmes, 2013)
install.packages("dyplr") #Installs the dyplr package for data manipulation (Wickham et al., 2014)

# Next we need to load these packages to be able to use their functions in this session. 
library(phyloseq); packageVersion("phyloseq")
library(dplyr); packageVersion("dplyr")

# Set the working directory to the appropriate folder (where the input files are stored). 
setwd("~/Downloads/results/")

# Load phyloseq object ####
ps = readRDS("ps.rds")

# Evaluating the phyloseq object
# We can evaluate the characteristics of the phyloseq object with various accessor functions. 
#Here we will focus on the functions that are relevant for the downstream analysis in this tutorial. 
#For additional functions please see phyloseq documentation. 
head(otu_table(ps)) #for both constructing and accessing the table ASV abundance 
head(tax_table(ps)) #for both constructing and accessing the table of taxonomic names
head(sample_data(ps)) #for both constructing and accessing the table of sample-level variables
head(taxa_sums(ps)) # tells the count of the ASVs summed over the samples 
head(sample_sums(ps)) # number of ASVs observed for each sample
rank_names(ps) # gives the taxonomic ranks included in the taxa table - useful to know available ranks to analyze

# There are also a few functions for editing the phyloseq object that will be used throughout the tutorial. Here is a basic description of the functions to give some context. Please see documentation for argument details and additional functions. 
??prune_taxa() # for filtering out unwanted ASVs in phyloseq object
??prune_samples() # for filtering out unwanted samples in phyloseq object
??filter_taxa() # filter ASVs based on across-sample abundance with a specified function
??subset_samples() # subset samples based on the provided conditions
??transform_sample_counts() # transform ASV abundance counts in each sample with a function of choice
??psmelt() # melt the phyloseq object to a dataframe 

# We can also extract the components from the phyloseq object if needed. 
taxatable <- as.matrix(tax_table(ps)) # extract the taxonomy table as a matrix
asvtable <- as.matrix(otu_table(ps)) # extract the ASV table as a matrix
sampledata <- data.frame(sample_data(ps)) # extract the sample metadata as a data frame

# Pre-processing ####
# First let's measure the sequencing depth for the samples and save these values to the data frame for reference. NOTE: @ for phyloseq objects is similar to $ for data frames in that it specifies the data component to be accessed.  
ps@sam_data$depth <- sample_sums(ps) # create a new column in the sample data for the total number of ASVs observed for each sample (sequencing depth)

# Now let's filter out spurious ASVs. Here is an example of commonly used filtering parameters, but these parameters may vary based on sample type and sequencing run. 
ps2 <- prune_taxa(taxa_sums(ps) > 1, ps) # Remove singleton ASVs to eliminate spurious taxa.
ps2 # Check how this changed the ps object
ps3 <- prune_samples(sample_sums(ps2) >= 1000, ps) # Remove samples with less than 1000 reads, which is a sign of poor quality sequences.
ps3 # Check how this changed the ps object
ps_clean <- filter_taxa(ps3, function(x) sum(x > 3) > (0.2*length(x)), prune = TRUE) # This returns a phyloseq object filtered to include only those taxa seen greater than 3 times in at least 20% of the samples. This threshold has been selected to remove ASVs with small means and large coefficients of variation.
ps_clean # Check how this changed the ps object. This is the final filtered ps object we will work with.

# Save phyloseq object
saveRDS(ps_clean,"ps_clean.rds")

# We can now calculate the new sequencing depth and compare this to the original to determine the reads lost from filtering. 
ps_clean@sam_data$depth_filt <- sample_sums(ps_clean) # save number of ASVs per sample after filtering to new column
sample_data <- data.frame(sample_data(ps_clean)) # save sample data frame with new and old sequencing depth columns to perform stats below

# Calculate % of reads lost from filtering
sdf_read_lost <- sample_data %>%
  mutate(read_lost = (depth - depth_filt)) %>% # new column with difference in reads pre- and post-filtering
  mutate(percent_lost = 100*(read_lost/depth)) %>% # convert reads lost to percentage for easier interpretation
  summarise(mean_lost=mean(read_lost), sd = sd(read_lost), mean_percent_lost=mean(percent_lost)) # create dataframe with stats
print(sdf_read_lost)

# Save the workspace
save.image("BioinfoWorkshop_phyloseq.RData")
