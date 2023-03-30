###############################################################
# Statistical analyses of 16S data - Activities               #
# Data: Mice - Antibiotic experience                          #
# By: ArrietaLab - University of Calgary                      #
# Author: Mackenzie Gutierrez                                 # 
# Dates: 3-5 April 2023                                       #
# Location: Oswaldo-Cruz Institute - Rio de Janeiro (Brazil)  #
###############################################################

# Activity 1: Alpha diversity - Chao1 ####
# Plot the Chao1 between treatment groups.
fig_chao_trt <- ggplot(richness2, aes(x= Treatment_Group, y = Chao1, color = Treatment_Group, fill = Treatment_Group)) + 
  theme_bw() +
  geom_boxplot(color = "black", alpha = 0.5) +
  geom_jitter(aes(color = Treatment_Group), position = position_jitter(0.2),  size = 1.2) +
  labs(x = "Treatment", y = "Chao1")+
  scale_color_manual(name = "Treatment", values = c("#797979","#F16400", "#531B92"))+ #dot color
  scale_fill_manual(name = "Treatment", values = c("#797979","#F16400", "#531B92"))+ #fill color
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = "none") #axis titles and legend aesthetics
fig_chao_trt

# Run the statistics.
# Summary statistics (sample size, mean, standard deviation) for Shannon diversity grouped by treatment group.
richness2 %>%
  group_by(Treatment_Group) %>%
  rstatix::get_summary_stats(Chao1, type = "mean_sd")

# Check if data is normally distributed (p-value < 0.05 is not normally distributed).
richness2 %>%
  group_by(Treatment_Group) %>%
  shapiro_test(Chao1) # Shapiro-Wilk normality test

# Assess homogeneity of sample variance (p-value < 0.05 indicates variance is not equal).
richness2 %>% 
  levene_test(Chao1 ~ Treatment_Group) # Levene's test for equality of variances

# Run the Kruskal-Wallis test.
richness2 %>% 
  kruskal_test(Chao1~ Treatment_Group)

# Test for multiple comparisons.
stat.test_chao1 <- richness2 %>% 
  tukey_hsd(Chao1 ~ Treatment_Group, conf.level=.95) %>% 
  filter(p.adj < 0.05) # filter for significant comparisons 
stat.test_chao1

# Add statistics to the plot - here we need to use ggpubr functions that are very similar to ggplot.
stat.test_chao1 <- add_xy_position(stat.test_chao1, x="Treatment_Group") # adds the y-axis position for the significance bars
fig_chao_trt <- ggboxplot(richness2, x = "Treatment_Group", y = "Chao1", color = "Treatment_Group", fill = "Treatment_Group", alpha = 0.5)+
  theme_bw() +
  geom_jitter(aes(color = Treatment_Group), position = position_jitter(0.2),  size = 1.2) +
  labs(x = "Treatment", y = "Richness (Chao1)")+
  scale_color_manual(name = "Treatment", values = c("#797979","#F16400", "#531B92"))+ # dot color
  scale_fill_manual(name = "Treatment", values = c("#797979","#F16400", "#531B92"))+ # fill color
  theme(axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = "none") +
  stat_pvalue_manual(stat.test_chao1, label = "p.adj.signif", tip.length = 0.01, size = 8) # add test results to the plot
fig_chao_trt

# Activity 2: Beta diversity - variance stabilization, PCoA plot and PERMANOVA - sex ####
# Variance Stabilizing Transformation.
ps_clean_deseq <- phyloseq_to_deseq2(ps_clean, ~Sex) # convert to DESeq2 format and set the independent variable as Sex

# Convert counts to integer.
ps_clean_deseq = estimateSizeFactors(ps_clean_deseq, geoMeans = apply(counts(ps_clean_deseq), 1, gm_mean))
vst_blind <- DESeq2::varianceStabilizingTransformation(ps_clean_deseq, blind = TRUE)
vst_blind_mat <- SummarizedExperiment::assay(vst_blind)
vst_blind_mat <- t(vst_blind_mat) 
vst_blind_mat[which(vst_blind_mat < 0)] <- 0 
dists <- dist(t(assay(ps_clean_deseq)))

# Computing Bray-Curtis Dissimilarities and PCoA.
comm_vst_blind_mat<- vegdist(vst_blind_mat, "bray")
PCoA_comm_vst_blind_mat<- capscale(comm_vst_blind_mat ~ 1, distance = "bray")
PCoA_comm_vst_blind_mat$CA$eig[1:3]/sum(PCoA_comm_vst_blind_mat$CA$eig)
PCoA_scores <- scores(PCoA_comm_vst_blind_mat)$sites

# Save scores into metadata tables.
row.names(sdf) == row.names(scores(PCoA_comm_vst_blind_mat)$sites)
sdf$PCoA1 <- scores(PCoA_comm_vst_blind_mat)$sites[,1]
sdf$PCoA2 <- scores(PCoA_comm_vst_blind_mat)$sites[,2]

# Variance stabilized PCoA plot by Sex. 
PCoA <- qplot(PCoA1, PCoA2,
              size = I(2), fill = Sex, color = Sex, data = (sdf))
# Customize plot with ggplot.
fig_pcoa_sex <- PCoA +
  stat_ellipse(level = 0.95, geom = "polygon", alpha = 1/6, linetype = 2, size = 0.5, 
               aes(fill = Sex, color = Sex)) +
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
  scale_color_manual(name = "Sex", values = c("#797979","#F16400"))+ #outer shape color
  scale_fill_manual(name = "Sex", values = c("#797979","#F16400"))  #inside shape color
fig_pcoa_sex

# PERMANOVA for sample time.
set.seed(112) # set seed for reproducible results
permanova_sex = adonis2(comm_vst_blind_mat ~ Sex, sdf, permutations = 999)
permanova_sex

# Activity 3: Relative abundance for the treatment groups at species level ####
# Agglomerate at species level and relativise.
ps_species <- ps_clean %>%
  tax_glom(taxrank = "Species_fullname") # agglomerate at the species level

# Melt the phyloseq object.  
ps_melt_species = psmelt(ps_species)

# Order based on abundance.
ps_ord_species = ps_melt_species %>%
  group_by(Species_fullname) %>%
  summarize_at("Abundance", sum) %>% # add total abundance of each species
  arrange(dplyr::desc(Abundance)) %>% # descending order
  mutate(rel_abund = Abundance/sum(Abundance)*100) # add relative abundance of each species

# Filter for the top 10 species.
species_top10 = ps_ord_species$Species_fullname[1:10] # create a list of top 10 species names
ps_species_top10_melt = subset_taxa(ps_species, as.data.frame(tax_table(ps_species))$Species_fullname %in% species_top10) %>% # select the top 10 species
  psmelt()
unique(ps_species_top10_melt$Species_fullname) # view the top 10 species names 

# Plot relative abundance.
nb.cols <- 10 # specify how many colors you will need
mycolors <- colorRampPalette(brewer.pal(10, "RdBu"))(nb.cols) # for more palette options: https://r-graph-gallery.com/38-rcolorbrewers-palettes.html

ps_species_top10_melt$Treatment_Group <- factor(ps_species_top10_melt$Treatment_Group, levels = c("Control", "Abx", "Abx+C.albicans")) # order categories 
fig_species_trt <- ggplot(ps_species_top10_melt, aes (x = Treatment_Group, y = Abundance, color = Species_fullname, fill = Species_fullname)) +
  geom_bar(stat = "identity", position = "fill") + 
  theme_bw() + 
  scale_y_continuous(labels = scales::percent)+
  labs(x = "Treatment", y = "Relative abundance") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12), 
        legend.text = element_text(face = "italic"),
        axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title = element_text(size = 16, face = "bold")) +
  scale_fill_manual(values=mycolors)+
  scale_color_manual(values = mycolors)
fig_species_trt

# Activity 4: Differentially abundant taxa between Abx and Abx+C.albicans ####
ps_clean_filt <- ps_clean %>%
  subset_samples(Treatment_Group == "Abx" | Treatment_Group == "Abx+C.albicans")

# First we need to convert the phyloseq object to the DESeq2 format.
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
res = results(ps_clean_deseq, cooksCutoff = FALSE) # extract results without applying Cook's cut off (distance threshold) 
alpha = 0.05 # set significance threshold
sigtab = res[which(res$padj < alpha), ] # filter for significant results with the adjusted p-values smaller than 0.05
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps_clean)[rownames(sigtab), ], "matrix")) # combine significant results with the phyloseq taxonomy table
head(sigtab %>% dplyr::select(log2FoldChange, padj, Genus)) # view Log2FoldChange, adjusted p-value and genus for significant results

# Now we can plot the results. Again we will represent the results at the genus level.
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x)) # select max log2FoldChange for each genus
x = sort(x, TRUE) # sort from highest to lowest
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x)) # convert Genus to factor with the levels specified as x to organize the plot

# Plot the results.
nb.cols <- 15 # specify how many colors you will need
mycolors <- colorRampPalette(brewer.pal(15, "RdBu"))(nb.cols)
fig_AbxAbxC <- ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Genus)) + geom_point(size=6) + 
  labs(title = "Differentially abundant genera in Abx vs Abx+C.albicans") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = "italic"),
  legend.text = element_text(face = "italic"))+
  scale_colour_manual(values=mycolors)+
  scale_fill_manual(values=mycolors) 
fig_AbxAbxC


