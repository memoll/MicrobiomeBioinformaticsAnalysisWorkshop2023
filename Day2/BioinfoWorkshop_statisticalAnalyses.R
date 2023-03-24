###############################################################
# Statistical analyses of 16S data                            #
# Data: Mice - Antibiotic experience                          #
# By: ArrietaLab - University of Calgary                      #
# Author: Mackenzie Gutierrez                                 # 
# Dates: 3-5 April 2023                                       #
# Location: Oswaldo-Cruz Institute - Rio de Janeiro (Brazil)  #
###############################################################

#Package installation ####
# We will be using several packages alongside phyloseq. First we need to install these packages (if they are not already installed).
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager") # Installs BiocManager, a package that facilitates the download of packages from the BioConductor repository.

#BiocManager::install("phyloseq") # Installs the phyloseq package (McMurdie & Holmes, 2013)
#install.packages("dyplr") #Installs the dyplr package for data manipulation (Wickham et al., 2014)
#BiocManager::install("DESeq2") # Installs DESeq2 for differential abundance analysis (Love et al., 2014)
#BiocManager::install("BiocGenerics") # Installs BiocGenerics for variance stabilizing transformatiom (Huber et al., 2015)
#BiocManager::install("SummarizedExperiment") # Installs SummarizedExperiment for variance stabilizing transformatiom (Morgan et al., 2022)
#install.packages("vegan") # Installs vegan for ecological statistical functions (Oksanen et al., 2020)
#install.packages("rstatix") # Installs rstatix for pipe friendly statistical analysis (Kassambara, 2021)
#install.packages("tidyverse") # Installs tidyverse for tidier coding (Wickham et al., 2019)
#install.packages("ggplot2") # Installs ggplot2 for plotting (Wickham, 2016)

# Next we need to load these packages to be able to use their functions in this session. 
library(phyloseq); packageVersion("phyloseq")
library(dplyr); packageVersion("dplyr")
library(tidyverse); packageVersion("tidyverse")
library(ggplot2); packageVersion("ggplot2")
library(vegan); packageVersion("vegan")
library(DESeq2); packageVersion("DESeq2")
library(rstatix); packageVersion("rstatix")
library(SummarizedExperiment); packageVersion("SummarizedExperiment")
library(BiocGenerics); packageVersion("BiocGenerics")

# Load filtered phyloseq object ####
setwd("~/Downloads/results/")
ps_clean = readRDS("ps_clean.rds")

# Alpha diversity ####
# Create alpha diversity variables and prepare data for plotting
richness <- estimate_richness(ps_clean, measures = c("Shannon", "Chao1")) # Here you can select from the several available diversity metrics (see documentation for more options).
colnames(richness) # Check the column names for the dataframe that was made by estimate_richness().
richness2 <- cbind(richness, ps_clean@sam_data) # Merge the results with the existing sample data.

# Plot the Shannon diversity between treatment groups
richness2$Treatment_Group <- factor(richness2$Treatment_Group, levels = c("Control", "Abx", "Abx+C.albicans")) # order categories 
fig <- ggplot(richness2, aes(x= Treatment_Group, y = Shannon, color = Treatment_Group, fill = Treatment_Group)) + 
  theme_bw() +
  geom_boxplot(color = "black", alpha = 0.5) +
  geom_jitter(aes(color = Treatment_Group), position = position_jitter(0.2),  size = 1.2) +
  labs(x = "Treatment", y = "Shannon diversity")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  scale_color_manual(name = "Treatment", values = c("#797979","#F16400", "#531B92"))+ #dot color
  scale_fill_manual(name = "Treatment", values = c("#797979","#F16400", "#531B92")) #fill color
fig

# Run the statistics
# Summary statistics (sample size, mean, standard deviation) for Shannon diversity grouped by treatment group.
richness2 %>%
  group_by(Treatment_Group) %>%
  get_summary_stats(Shannon, type = "mean_sd")

# Check for outliers (will only get results if outliers are present).
richness2 %>% 
  group_by(Treatment_Group) %>%
  select(Treatment_Group, Shannon) %>% # select variables of interest
  identify_outliers(Shannon) # outlier test

# Check if data is normally distributed (p-value < 0.05 is not normally distributed).
richness2 %>%
  group_by(Treatment_Group) %>%
  shapiro_test(Shannon) # Shapiro-Wilk normality test

# Assess homogeneity of sample variance (p-value < 0.05 indicates unequal variances).
richness2 %>% 
  levene_test(Shannon ~ Treatment_Group) # Levene's test for equality of variances

# Run the ANOVA.
richness2 %>% 
  anova_test(Shannon ~ Treatment_Group)

# test for multiple comparisons 
richness2 %>% 
  tukey_hsd(Shannon ~ Treatment_Group, conf.level=.95) %>% 
  filter(p.adj < 0.05) # filter for significant comparisons 


#Activity: reproduce this plot and statistics for Chao1

# Plot the Chao1 between treatment groups
fig <- ggplot(richness2, aes(x= Treatment_Group, y = Chao1, color = Treatment_Group, fill = Treatment_Group)) + 
  theme_bw() +
  geom_boxplot(color = "black", alpha = 0.5) +
  geom_jitter(aes(color = Treatment_Group), position = position_jitter(0.2),  size = 1.2) +
  labs(x = "Treatment", y = "Chao1")+
  theme(plot.title = element_text(face = "bold", hjust = 0.5))+
  scale_color_manual(name = "Treatment", values = c("#797979","#F16400", "#531B92"))+ #dot color
  scale_fill_manual(name = "Treatment", values = c("#797979","#F16400", "#531B92")) #fill color
fig

# Run the statistics
# Summary statistics (sample size, mean, standard deviation) for Shannon diversity grouped by treatment group.
richness2 %>%
  group_by(Treatment_Group) %>%
  get_summary_stats(Chao1, type = "mean_sd")

# Check for outliers (will only get results if outliers are present).
richness2 %>% 
  group_by(Treatment_Group) %>%
  select(Treatment_Group, Chao1) %>% # select variables of interest
  identify_outliers(Chao1) # outlier test

# Check if data is normally distributed (p-value < 0.05 is not normally distributed).
richness2 %>%
  group_by(Treatment_Group) %>%
  shapiro_test(Chao1) # Shapiro-Wilk normality test

# Assess homogeneity of sample variance (p-value < 0.05 indicates unequal variances).
richness2 %>% 
  levene_test(Chao1 ~ Treatment_Group) # Levene's test for equality of variances

# Run the ANOVA.
richness2 %>% 
  kruskal_test(Chao1~ Treatment_Group)

# test for multiple comparisons 
richness2 %>% 
  tukey_hsd(Chao1 ~ Treatment_Group, conf.level=.95) %>% 
  filter(p.adj < 0.05) # filter for significant comparisons 

# Beta diversity ####
# save sample data from ps_clean 
sdf <- as(sample_data(ps_clean), "data.frame")
sdf$Treatment_Group <- factor(sdf$Treatment_Group, levels = c("Control", "Abx", "Abx+C.albicans")) # order categories 

# Create function geo means for Variance Stabilizing Transformation (based on sample size).
gm_mean = function(x, na.rm = TRUE){exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))}

# Variance Stabilizing Transformation
ps_clean_deseq <- phyloseq_to_deseq2(ps_clean, ~Treatment_Group)

# Convert counts to integer
ps_clean_deseq = estimateSizeFactors(ps_clean_deseq, geoMeans = apply(counts(ps_clean_deseq), 1, gm_mean))
vst_blind <- DESeq2::varianceStabilizingTransformation(ps_clean_deseq, blind = TRUE)
vst_blind_mat <- SummarizedExperiment::assay(vst_blind)
vst_blind_mat <- t(vst_blind_mat) 
vst_blind_mat[which(vst_blind_mat < 0)] <- 0 
dists <- dist(t(assay(ps_clean_deseq)))

# Computing Bray-Curtis Dissimilarities and PCoA
comm_vst_blind_mat<- vegdist(vst_blind_mat, "bray")
PCoA_comm_vst_blind_mat<- capscale(comm_vst_blind_mat ~ 1, distance = "bray")
PCoA_comm_vst_blind_mat$CA$eig[1:3]/sum(PCoA_comm_vst_blind_mat$CA$eig)

PCoA_scores <- scores(PCoA_comm_vst_blind_mat)$sites

# Save scores into metadata tables
row.names(sdf) == row.names(scores(PCoA_comm_vst_blind_mat)$sites)
sdf$PCoA1 <- scores(PCoA_comm_vst_blind_mat)$sites[,1]
sdf$PCoA2 <- scores(PCoA_comm_vst_blind_mat)$sites[,2]

# Variance stabilized PCoA plot by colonization 
PCoA <- qplot(PCoA1, PCoA2,
              size = I(2), fill = Treatment_Group, color = Treatment_Group, data = (sdf))
# Customize plot with ggplot
fig <- PCoA +
  stat_ellipse(level = 0.95, geom = "polygon", alpha = 1/6, linetype = 2, size = 0.5, 
               aes(fill = Treatment_Group, color = Treatment_Group)) +
  theme_bw() + 
  labs(title = "Bray-Curtis Dissimilarity", x = "PCoA1 (48.6%)", y = "PCoA2 (31.7%)")+ #PLS DOUBLE-CHECK THE PERCENTAGES
  theme(legend.title = element_text(colour = "black", size = 9.5, face = "bold"),
        legend.text = element_text(colour = "black", size = 9.5),
        legend.position = "right",
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        axis.text = element_text(size = 9.5),
        strip.text.x = element_text(face = "bold"),
        plot.title = element_text(colour = "black", size = 12, face = "bold", hjust=0.5)) +
  scale_color_manual(name = "Treatment", values = c("#797979","#F16400", "#531B92"))+ #outer shape color
  scale_fill_manual(name = "Treatment", values = c("#797979","#F16400", "#531B92"))  #inside shape color
fig

# PERMANOVA for sample time
set.seed(111)
permanova = adonis2(comm_vst_blind_mat ~ Treatment_Group, sdf, permutations = 999)
permanova

# Activity: repeat the variance stabilization, plot and PERMANOVA but with sex as the variable of interest.
# Variance Stabilizing Transformation
ps_clean_deseq <- phyloseq_to_deseq2(ps_clean, ~Sex)

# Convert counts to integer
ps_clean_deseq = estimateSizeFactors(ps_clean_deseq, geoMeans = apply(counts(ps_clean_deseq), 1, gm_mean))
vst_blind <- DESeq2::varianceStabilizingTransformation(ps_clean_deseq, blind = TRUE)
vst_blind_mat <- SummarizedExperiment::assay(vst_blind)
vst_blind_mat <- t(vst_blind_mat) 
vst_blind_mat[which(vst_blind_mat < 0)] <- 0 
dists <- dist(t(assay(ps_clean_deseq)))

# Computing Bray-Curtis Dissimilarities and PCoA
comm_vst_blind_mat<- vegdist(vst_blind_mat, "bray")
PCoA_comm_vst_blind_mat<- capscale(comm_vst_blind_mat ~ 1, distance = "bray")
PCoA_comm_vst_blind_mat$CA$eig[1:3]/sum(PCoA_comm_vst_blind_mat$CA$eig)

PCoA_scores <- scores(PCoA_comm_vst_blind_mat)$sites

# Save scores into metadata tables
row.names(sdf) == row.names(scores(PCoA_comm_vst_blind_mat)$sites)
sdf$PCoA1 <- scores(PCoA_comm_vst_blind_mat)$sites[,1]
sdf$PCoA2 <- scores(PCoA_comm_vst_blind_mat)$sites[,2]

# Variance stabilized PCoA plot by colonization 
PCoA <- qplot(PCoA1, PCoA2,
              size = I(2), fill = Sex, color = Sex, data = (sdf))
# Customize plot with ggplot
fig <- PCoA +
  stat_ellipse(level = 0.95, geom = "polygon", alpha = 1/6, linetype = 2, size = 0.5, 
               aes(fill = Sex, color = Sex)) +
  theme_bw() + 
  labs(title = "Bray-Curtis Dissimilarity", x = "PCoA1 (48.6%)", y = "PCoA2 (31.7%)")+ #PLS CORRECT THE PERCENTAGES
  theme(legend.title = element_text(colour = "black", size = 9.5, face = "bold"),
        legend.text = element_text(colour = "black", size = 9.5),
        legend.position = "right",
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        axis.text = element_text(size = 9.5),
        strip.text.x = element_text(face = "bold"),
        plot.title = element_text(colour = "black", size = 12, face = "bold", hjust=0.5)) +
  scale_color_manual(name = "Sex", values = c("#797979","#F16400"))+ #outer shape color
  scale_fill_manual(name = "Sex", values = c("#797979","#F16400"))  #inside shape color
fig

# PERMANOVA for sample time
set.seed(112)
permanova = adonis2(comm_vst_blind_mat ~ Sex, sdf, permutations = 999)
permanova

# Relative Abundance ####
# Agglomerate at Species level and relativise
ps_genus <- ps_clean %>%
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x) x*100 / sum(x)) %>%
  psmelt()

# Plot relativea abundance
ps_genus$Treatment_Group <- factor(ps_genus$Treatment_Group, levels = c("Control", "Abx", "Abx+C.albicans")) # order categories 

fig <- ggplot(ps_genus, aes (x = Treatment_Group, y = Abundance, color = Genus, fill = Genus)) +
  geom_bar(stat = "identity", position = "fill") + 
  theme_bw() + 
  scale_y_continuous(labels = scales::percent)+
  labs(x = "Treatment", y = "Relative abundance") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12), 
        legend.text = element_text(face = "italic")) +
  scale_fill_manual(values=c("#fde725", "#d0e11c", "#a0da39", "#73d056", "#4ac16d", "#2db27d", "#1fa187", "#21918c", "#277f8e", "#2e6e8e", "#365c8d", "#3f4788", "#46327e", "#481b6d", "#440154"))+
  scale_color_manual(values = c("#fde725", "#d0e11c", "#a0da39", "#73d056", "#4ac16d", "#2db27d", "#1fa187", "#21918c", "#277f8e", "#2e6e8e", "#365c8d", "#3f4788", "#46327e", "#481b6d", "#440154"))
fig

# Activity: Plot relative abundance for the treatment groups at species level (hint - you need x colours).
ps_species <- ps_clean %>%
  tax_glom(taxrank = "Species") %>%
  transform_sample_counts(function(x) x*100 / sum(x)) %>%
  psmelt()

# Plot relativea abundance
ps_species$Treatment_Group <- factor(ps_species$Treatment_Group, levels = c("Control", "Abx", "Abx+C.albicans")) # order categories 

fig <- ggplot(ps_species, aes (x = Treatment_Group, y = Abundance, color = Species, fill = Species)) +
  geom_bar(stat = "identity", position = "fill") + 
  theme_bw() + 
  scale_y_continuous(labels = scales::percent)+
  labs(x = "Treatment", y = "Relative abundance") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12), 
        legend.text = element_text(face = "italic")) +
  scale_fill_manual(values=c("#fde725", "#d0e11c", "#a0da39", "#73d056", "#4ac16d", "#2db27d", "#1fa187", "#21918c", "#277f8e"))+
  scale_color_manual(values = c("#fde725", "#d0e11c", "#a0da39", "#73d056", "#4ac16d", "#2db27d", "#1fa187", "#21918c", "#277f8e"))
fig

# Differential abundance with DESeq ####
# Let's look at taxa that are differentially abundant in Control vs Abx
ps_clean_filt <- ps_clean %>%
  subset_samples(Treatment_Group == "Control" | Treatment_Group == "Abx")

# First we need to convert the phyloseq object to the DESeq2 format.
sample_data(ps_clean_filt)$Treatment_Group <- relevel(sample_data(ps_clean_filt)$Treatment_Group, "Control") # re-orderdering the levels of treatment group to set control as the reference in downstream tests
ps_clean_deseq = phyloseq_to_deseq2(ps_clean_filt, ~ Treatment_Group) # convert to DESeq2 format and set the independent variable as SampleType

# This dataset contains several zeros so we need to use a zero-tolerant variant of geometric mean.
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(ps_clean_deseq), 1, gm_mean)

# Estimation of size factors and dispersion and fits the model.
ps_clean_deseq = estimateSizeFactors(ps_clean_deseq, geoMeans = geoMeans)
ps_clean_deseq = DESeq(ps_clean_deseq, test="Wald", fitType="parametric") 

# Lets take a look at the results from the test above. 
res = results(ps_clean_deseq, cooksCutoff = FALSE) # extract results without applying Cook's cut off (distance threshhold) 
alpha = 0.05 # set significance threshhold
sigtab = res[which(res$padj < alpha), ] # filter for significant results with the adjusted p-values smaller than 0.05
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps_clean)[rownames(sigtab), ], "matrix")) # combine significant results with the phyloseq taxonomy table
head(sigtab %>% select(log2FoldChange, padj, Genus)) # view Log2FoldChange, adjusted p-vlaue and phylum for significant results

# Now we can plot the results. Again we will represent the results at the phylum level.
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x)) # select max log2FoldChange for each phylum
x = sort(x, TRUE) # sort from highest to lowest
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x)) # convert Phylum to factor with the levels specified as x to organize the plot
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Genus)) + geom_point(size=6) + 
  labs(title = "Differentially abundant genera in control vs Abx") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Activity: Assess differentially abundant taxa between Abx and Abx+C.albicans
ps_clean_filt <- ps_clean %>%
  subset_samples(Treatment_Group == "Abx" | Treatment_Group == "Abx+C.albicans")

# First we need to convert the phyloseq object to the DESeq2 format.
sample_data(ps_clean_filt)$Treatment_Group <- relevel(sample_data(ps_clean_filt)$Treatment_Group, "Abx") # re-orderdering the levels of treatment group to set control as the reference in downstream tests
ps_clean_deseq = phyloseq_to_deseq2(ps_clean_filt, ~ Treatment_Group) # convert to DESeq2 format and set the independent variable as SampleType

# This dataset contains several zeros so we need to use a zero-tolerant variant of geometric mean.
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(ps_clean_deseq), 1, gm_mean)

# Estimation of size factors and dispersion and fits the model.
ps_clean_deseq = estimateSizeFactors(ps_clean_deseq, geoMeans = geoMeans)
ps_clean_deseq = DESeq(ps_clean_deseq, test="Wald", fitType="parametric") 

# Lets take a look at the results from the test above. 
res = results(ps_clean_deseq, cooksCutoff = FALSE) # extract results without applying Cook's cut off (distance threshhold) 
alpha = 0.05 # set significance threshhold
sigtab = res[which(res$padj < alpha), ] # filter for significant results with the adjusted p-values smaller than 0.05
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps_clean)[rownames(sigtab), ], "matrix")) # combine significant results with the phyloseq taxonomy table
head(sigtab %>% select(log2FoldChange, padj, Genus)) # view Log2FoldChange, adjusted p-vlaue and phylum for significant results

# Now we can plot the results. Again we will represent the results at the phylum level.
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x)) # select max log2FoldChange for each phylum
x = sort(x, TRUE) # sort from highest to lowest
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x)) # convert Phylum to factor with the levels specified as x to organize the plot
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Genus)) + geom_point(size=6) + 
  labs(title = "Differentially abundant genera in Abx vs Abx+C.albicans") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

