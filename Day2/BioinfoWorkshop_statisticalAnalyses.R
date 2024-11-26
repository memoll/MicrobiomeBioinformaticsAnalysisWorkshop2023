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
#install.packages("ggpubr") # Installs ggpubr for plotting statistics from rstatix (Kassambara, 2023)
#install.packages("RColorBrewer") # Installs RColorBrewer for ggplot color palettes (Neuwirth, 2022)

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
library(ggpubr); packageVersion("ggpubr")
library(RColorBrewer); packageVersion("RColorBrewer")

# Load filtered phyloseq object ####
setwd("~/Downloads/results/")
ps_clean = readRDS("ps_clean.rds")

# Alpha diversity ####
# Create alpha diversity variables and prepare data for plotting.
richness <- estimate_richness(ps_clean, measures = c("Shannon", "Chao1")) # here you can select from the several available diversity metrics (see documentation for more options)
colnames(richness) # check the column names for the data frame that was made by estimate_richness()
richness2 <- cbind(richness, ps_clean@sam_data) # merge the results with the existing sample data

# Plot the Shannon diversity between treatment groups
richness2$Treatment_Group <- factor(richness2$Treatment_Group, levels = c("Control", "Abx", "Abx+C.albicans")) # order categories 
fig_shn_trt <- ggplot(richness2, aes(x= Treatment_Group, y = Shannon, color = Treatment_Group, fill = Treatment_Group)) + 
  geom_boxplot(color = "black", alpha = 0.5) +
  geom_jitter(aes(color = Treatment_Group), position = position_jitter(0.2),  size = 1.2) +
  labs(x = "Treatment", y = "Shannon diversity")+
  scale_color_manual(name = "Treatment", values = c("#797979","#F16400", "#531B92"))+ # dot color
  scale_fill_manual(name = "Treatment", values = c("#797979","#F16400", "#531B92"))+ # fill color
  theme_bw() + 
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = "none") # axis titles and legend aesthetics
fig_shn_trt

# Run the statistics.
# Summary statistics (sample size, mean, standard deviation) for Shannon diversity grouped by treatment group to get an overview.
richness2 %>%
  group_by(Treatment_Group) %>%
  rstatix::get_summary_stats(Shannon, type = "mean_sd")

# Check if data is normally distributed (p-value < 0.05 is not normally distributed).
richness2 %>%
  group_by(Treatment_Group) %>%
  shapiro_test(Shannon) # Shapiro-Wilk normality test

# Assess homogeneity of sample variance (p-value < 0.05 indicates variance is not equal).
richness2 %>% 
  levene_test(Shannon ~ Treatment_Group) # Levene's test for equality of variances

# Run the ANOVA.
richness2 %>% 
  anova_test(Shannon ~ Treatment_Group)

# Test for multiple comparisons. 
stat.test_shannon <- richness2 %>% 
  tukey_hsd(Shannon ~ Treatment_Group, conf.level=.95) %>% 
  filter(p.adj < 0.05) # filter for significant comparisons 
stat.test_shannon

# Add statistics to the plot - here we need to use ggpubr functions that are very similar to ggplot.
stat.test_shannon <- add_xy_position(stat.test_shannon, x="Treatment_Group") # adds the y-axis position for the significance bars
fig_shn_trt <- ggboxplot(richness2, x = "Treatment_Group", y = "Shannon", color = "Treatment_Group", fill = "Treatment_Group", alpha = 0.5)+
  theme_bw() +
  geom_jitter(aes(color = Treatment_Group), position = position_jitter(0.2),  size = 1.2) +
  labs(x = "Treatment", y = "Shannon diversity")+
  scale_color_manual(name = "Treatment", values = c("#797979","#F16400", "#531B92"))+ #dot color
  scale_fill_manual(name = "Treatment", values = c("#797979","#F16400", "#531B92"))+ #fill color
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = "none") +
  stat_pvalue_manual(stat.test_shannon, label = "p.adj.signif", tip.length = 0.01, size = 8) # add test results to the plot
fig_shn_trt

# Activity 1: reproduce this plot and statistics for Chao1.

# Beta diversity ####
# Save sample data from ps_clean. 
sdf <- as(sample_data(ps_clean), "data.frame")
sdf$Treatment_Group <- factor(sdf$Treatment_Group, levels = c("Control", "Abx", "Abx+C.albicans")) # order categories 

# Create function geo means for Variance Stabilizing Transformation (based on sample size).
gm_mean = function(x, na.rm = TRUE){exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))}

# Variance Stabilizing Transformation.
ps_clean_deseq <- phyloseq_to_deseq2(ps_clean, ~ Treatment_Group) # convert to DESeq2 format 

# Convert counts to integer.
ps_clean_deseq = estimateSizeFactors(ps_clean_deseq, geoMeans = apply(counts(ps_clean_deseq), 1, gm_mean))
vst_blind <- DESeq2::varianceStabilizingTransformation(ps_clean_deseq, blind = TRUE)
vst_blind_mat <- SummarizedExperiment::assay(vst_blind)
vst_blind_mat <- t(vst_blind_mat) 
vst_blind_mat[which(vst_blind_mat < 0)] <- 0 
#dists <- dist(t(assay(ps_clean_deseq)))

# Computing Bray-Curtis Dissimilarities and PCoA.
comm_vst_blind_mat <- vegdist(vst_blind_mat, "bray")
PCoA_comm_vst_blind_mat <- capscale(comm_vst_blind_mat ~ 1, distance = "bray")
PCoA_comm_vst_blind_mat$CA$eig[1:3]/sum(PCoA_comm_vst_blind_mat$CA$eig)
PCoA_scores <- scores(PCoA_comm_vst_blind_mat)$sites

# Save scores into metadata table.
row.names(sdf) == row.names(scores(PCoA_comm_vst_blind_mat)$sites)
sdf$PCoA1 <- scores(PCoA_comm_vst_blind_mat)$sites[,1]
sdf$PCoA2 <- scores(PCoA_comm_vst_blind_mat)$sites[,2]

# Variance stabilized PCoA plot by Treatment_Group. 
PCoA <- qplot(PCoA1, PCoA2,
              size = I(2), fill = Treatment_Group, color = Treatment_Group, data = (sdf))
# Customize plot with ggplot.
fig_pcoa_trt <- PCoA +
  stat_ellipse(level = 0.95, geom = "polygon", alpha = 1/6, linetype = 2, linewidth = 0.5, 
               aes(fill = Treatment_Group, color = Treatment_Group)) +
  theme_bw() + 
  labs(title = "Bray-Curtis Dissimilarity",
       x = paste0("PCoA1 (", round(PCoA_comm_vst_blind_mat$CA$eig[1:1]/sum(PCoA_comm_vst_blind_mat$CA$eig)*100,digits=1), "%)"), 
       y = paste0("PCoA2 (", round(PCoA_comm_vst_blind_mat$CA$eig[2:1]/sum(PCoA_comm_vst_blind_mat$CA$eig)*100,digits=1), "%)"))+ 
  theme(legend.title = element_text(colour = "black", size = 9.5, face = "bold"),
        legend.text = element_text(colour = "black", size = 9.5),
        legend.position = "right",
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        axis.text = element_text(size = 9.5),
        strip.text.x = element_text(face = "bold"),
        plot.title = element_text(colour = "black", size = 12, face = "bold", hjust=0.5)) +
  scale_color_manual(name = "Treatment", values = c("#797979","#F16400", "#531B92"))+ # outer shape color
  scale_fill_manual(name = "Treatment", values = c("#797979","#F16400", "#531B92"))  # inside shape color
fig_pcoa_trt

# PERMANOVA for Treatment_Group.
set.seed(111) # set seed for reproducible results
permanova_treatment = adonis2(comm_vst_blind_mat ~ Treatment_Group, sdf, permutations = 999)
permanova_treatment

# Activity 2: repeat the variance stabilization, plot and PERMANOVA but with sex as the variable of interest.

# Relative Abundance ####
# Prepare the taxa table for plotting. 
taxa <- data.frame(tax_table(ps_clean)) # extract taxa table as a data frame
taxa2 <- taxa %>%
  replace_na(list(Genus = "Unassigned", Species = "unassigned")) %>% # replace NAs as unassigned to ensure they are plotted 
  unite("Species_fullname", Genus:Species, remove = FALSE) # merge genus and species names and save as a new column

taxa2 <- as.matrix(taxa2) # make sure taxa table is in matrix form to fit back into ps pbject
tax_table(ps_clean) <- taxa2 # add taxa table back to the ps object
ps_clean

# Agglomerate at genus level and relativise.
ps_genus <- ps_clean %>%
  tax_glom(taxrank = "Genus") # agglomerate at the genus level

# Melt the phyloseq object.  
ps_melt = psmelt(ps_genus)

# Order based on abundance.
ps_ord = ps_melt %>%
  group_by(Genus) %>%
  summarize_at("Abundance", sum) %>% # add total abundance of each genus
  arrange(dplyr::desc(Abundance)) %>% # descending order
  mutate(rel_abund = Abundance/sum(Abundance)*100) # add relative abundance of each genus

# Filter for the top 10 genera.
genus_top10 = ps_ord$Genus[1:10] # create a list of top 10 genera names
ps_genus_top10_melt = subset_taxa(ps_genus, as.data.frame(tax_table(ps_genus))$Genus %in% genus_top10) %>% # select the top 10 genera
  psmelt() 
unique(ps_genus_top10_melt$Genus) # view the top 10 genera names 

# Plot relative abundance.
nb.cols <- 10 # specify how many colors you will need
mycolors <- colorRampPalette(brewer.pal(10, "RdBu"))(nb.cols) # for more palette options: https://r-graph-gallery.com/38-rcolorbrewers-palettes.html

ps_genus_top10_melt$Treatment_Group <- factor(ps_genus_top10_melt$Treatment_Group, levels = c("Control", "Abx", "Abx+C.albicans")) # order categories 
fig_genus_top10 <- ggplot(ps_genus_top10_melt, aes (x = Treatment_Group, y = Abundance, color = Genus, fill = Genus)) +
  geom_bar(stat = "identity", position = "fill") + 
  theme_bw() + 
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Treatment", y = "Relative abundance") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12), 
        legend.text = element_text(face = "italic"),
        axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title = element_text(size = 16, face = "bold")) +
  scale_fill_manual(values=mycolors) +
  scale_color_manual(values = mycolors)
fig_genus_top10

# Activity 3: Plot relative abundance for the treatment groups at species level (hint - use the Species_fullname column).

# Differential abundance with DESeq ####
# Let's look at taxa that are differentially abundant in Control vs Abx
ps_clean_filt <- ps_clean %>%
  subset_samples(Treatment_Group == "Control" | Treatment_Group == "Abx")

# First we need to convert the phyloseq object to the DESeq2 format.
ps_clean_deseq = phyloseq_to_deseq2(ps_clean_filt, ~ Treatment_Group) # convert to DESeq2 format

# This data set contains several zeros so we need to use a zero-tolerant variant of geometric mean.
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(ps_clean_deseq), 1, gm_mean)

# Estimation of size factors and dispersion and fits the model.
ps_clean_deseq = estimateSizeFactors(ps_clean_deseq, geoMeans = geoMeans)
ps_clean_deseq = DESeq(ps_clean_deseq, test="Wald", fitType="parametric") 

# Lets take a look at the results from the test above. 
res = results(ps_clean_deseq, cooksCutoff = FALSE) # extract results without applying Cook's cut off (distance threshold) 
alpha = 0.05 # set significance threshold
sigtab = res[which(res$padj < alpha), ] # filter for significant results with the adjusted p-values smaller than 0.05
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps_clean)[rownames(sigtab), ], "matrix")) # combine significant results with the phyloseq taxonomy table
head(sigtab %>% dplyr::select(log2FoldChange, padj, Genus)) # view Log2FoldChange, adjusted p-value and genus for significant results

# Now we can plot the results. We will represent the results at the Genus level.
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x)) # select max log2FoldChange for each genus
x = sort(x, TRUE) # sort from highest to lowest
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x)) # convert Genus to factor with the levels specified as x to organize the plot

# Plot the results.
nb.cols <- 18 # specify how many colors you will need
mycolors <- colorRampPalette(brewer.pal(18, "RdBu"))(nb.cols)
fig_CtlAbx <- ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Genus)) + geom_point(size=6) + 
  labs(title = "Differentially abundant genera in Control vs Abx") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = "italic"),
        legend.text = element_text(face = "italic"))+
  scale_colour_manual(values=mycolors)+
  scale_fill_manual(values=mycolors) 
fig_CtlAbx

# Activity 4: Assess differentially abundant taxa between Abx and Abx+C.albicans

# Save the workspace
save.image("BioinfoWorkshop_statisticalAnalyses.RData")

