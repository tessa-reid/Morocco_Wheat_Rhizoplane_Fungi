############################################################
#R Scripts written by Tessa Reid tessa.reid@rothamsted.ac.uk
############################################################
#R 4.2.2
#setwd()
#Plot inspiration from: https://github.com/marie-simonin/Wheat_Microbiome
#Publication: https://doi.org/10.1093/femsec/fiaa067



######## Alpha diversity analysis

#input files: 
#1. phyloseq object with ASVs classified by UNITE with rarefied ASV abundances (physeq_all_rare.rds)
#2. phyloseq object with only ASVs classified to strain level by Isolate sequences with rarefied ASV abundances (physeq_isl_rare.rds)

library(phyloseq)
physeq.all.rare <- readRDS("Data/physeq_all_rare.rds")
physeq.all.rare
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 4457 taxa and 135 samples ]
#sample_data() Sample Data:       [ 135 samples by 7 sample variables ]
#tax_table()   Taxonomy Table:    [ 4457 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 4457 tips and 4428 internal nodes ]
physeq.isl.rare <- readRDS("Data/physeq_str_rare.rds")
physeq.isl.rare
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 230 taxa and 135 samples ]
#sample_data() Sample Data:       [ 135 samples by 7 sample variables ]
#tax_table()   Taxonomy Table:    [ 230 taxa by 10 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 230 tips and 227 internal nodes ]


#Alpha diversity metrics with microbiome
library(microbiome) #v1.19.1

hmp.div.all <- microbiome::alpha(physeq.all.rare, index = c("all"))
hmp.div.isl <- microbiome::alpha(physeq.isl.rare, index = c("all"))

hmp.meta <- meta(physeq.all.rare)
hmp.meta$sam_name <- rownames(hmp.meta)

hmp.div.all$sam_name <- rownames(hmp.div.all)
hmp.div.isl$sam_name <- rownames(hmp.div.isl)

div.df.all <- merge(hmp.div.all, hmp.meta, by = "sam_name")
div.df.isl <- merge(hmp.div.isl, hmp.meta, by = "sam_name")


#Faith's Phylogenetic Diversity
library(picante) #v1.8.2

physeq.all.rare.asvtab <- as.data.frame(physeq.all.rare@otu_table)
physeq.isl.rare.asvtab <- as.data.frame(physeq.isl.rare@otu_table)

physeq.all.rare.tree <- physeq.all.rare@phy_tree
physeq.isl.rare.tree <- physeq.isl.rare@phy_tree


physeq.all.rare@phy_tree
#
#Phylogenetic tree with 4457 tips and 4428 internal nodes.
#
#Tip labels:
#  ab006d72c8e8045baef3b0834cd78ba4, 32c005d100dfd11959d64c808a83a8fe, 5a7dcdab8a5af50ae445226a1a26cf53, dd2f04a2d5835ae9b8c11b6c465b0489, d7e3aad568731ef7e2d66a428ba81091, dad15b48e976322dd2fbd9b11bf9d25e, ...
#Node labels:
#  root, 0.383, 0.913, 0.622, 0.962, 0.890, ...
#
#Rooted; includes branch lengths.
physeq.isl.rare@phy_tree
#
#Phylogenetic tree with 230 tips and 227 internal nodes.
#
#Tip labels:
#  222901f0b3a1f22c5d280636cb18f670, 7564168ed0eaac53ac910a953cb17443, c6ec2234ba5f1af2bc4b1c2d88771376, 7165563ff158b2927c79e680d7a3b71f, 752ee4adaf696313d1fb0cf5aad1d852, a6d476eed5034e046eebc8ffa3713f26, ...
#Node labels:
#  0.924, 0.459, 1.000, 0.413, 0.375, 0.441, ...
#
#Rooted; includes branch lengths.

df.pd.all <- pd(t(physeq.all.rare.asvtab), physeq.all.rare.tree, include.root=T)
df.pd.isl <- pd(t(physeq.isl.rare.asvtab), physeq.isl.rare.tree, include.root=T)

div.df.all$Phylogenetic_Diversity <- df.pd.all$PD
div.df.isl$Phylogenetic_Diversity <- df.pd.isl$PD

#write to excel
writexl::write_xlsx(div.df.all, "Data/div.df.all.xlsx")
writexl::write_xlsx(div.df.isl, "Data/div.df.isl.xlsx")




#### Figures 1A and B: Richness and Diversity

#Mean and standard deviation per Region
div.df.all.1 <- div.df.all[, c("observed", "Phylogenetic_Diversity","Region", "Zone")]
div.df.isl.1 <- div.df.isl[, c("observed", "Phylogenetic_Diversity","Region", "Zone")]

library(dplyr) #v1.1.2
grouped.all <- group_by(div.df.all.1, Region)
grouped.isl <- group_by(div.df.isl.1, Region)

richness.all <- summarise(grouped.all, mean=mean(observed), sd=sd(observed))
richness.isl <- summarise(grouped.isl, mean=mean(observed), sd=sd(observed))

diversity.all <- summarise(grouped.all, mean=mean(Phylogenetic_Diversity), sd=sd(Phylogenetic_Diversity))
diversity.isl <- summarise(grouped.isl, mean=mean(Phylogenetic_Diversity), sd=sd(Phylogenetic_Diversity))


#Plots
library(ggplot2) #v3.4.2


p1 <- ggplot() + 
  geom_point(data=richness.all, aes(x=reorder(Region, mean), y=mean), size=4, shape=16, color="navy", fill="navy") +
  geom_point(data=richness.isl, aes(x=Region, y=mean), size=4, shape=1, color="navy") +
  geom_errorbar(data=richness.all, aes(x=Region, ymin=mean-sd, ymax=mean+sd), width=.2, color="navy") +
  geom_errorbar(data=richness.isl, aes(x=Region, ymin=mean-sd, ymax=mean+sd), width=.2, color="navy") +
  theme_classic() +
  xlab("") +
  ylab("Observed ASV richness") +
  theme(axis.title = element_text(color="black", size=12),
        axis.text = element_text(color="black", size=12),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p1


p2 <- ggplot() + 
  geom_point(data=diversity.all, aes(x=reorder(Region, mean), y=mean), size=4, shape=16, color="navy", fill="navy") +
  geom_point(data=diversity.isl, aes(x=Region, y=mean), size=4, shape=1, color="navy") +
  geom_errorbar(data=diversity.all, aes(x=Region, ymin=mean-sd, ymax=mean+sd), width=.2, color="navy") +
  geom_errorbar(data=diversity.isl, aes(x=Region, ymin=mean-sd, ymax=mean+sd), width=.2, color="navy") +
  theme_classic() +
  xlab("Region") +
  ylab("Faith's phylogenetic diversity") +
  theme(axis.title = element_text(color="black", size=12),
        axis.text = element_text(color="black", size=12),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p2




##Statistical analysis

#Compute the analysis of variance
res.aov.richness.all <- aov(observed ~ Region, data = div.df.all.1)
res.aov.diversity.all <- aov(Phylogenetic_Diversity ~ Region, data = div.df.all.1)
res.aov.richness.isl <- aov(observed ~ Region, data = div.df.isl.1)
res.aov.diversity.isl <- aov(Phylogenetic_Diversity ~ Region, data = div.df.isl.1)

#Summary of the analysis
summary(res.aov.richness.all)
#             Df Sum Sq Mean Sq F value Pr(>F)    
#Region        6 209474   34912   22.79 <2e-16 ***
#Residuals   128 196117    1532                   
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(res.aov.diversity.all)
#             Df Sum Sq Mean Sq F value   Pr(>F)    
#Region        6   4329   721.5   12.74 3.03e-11 ***
#Residuals   128   7247    56.6                     
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(res.aov.richness.isl)
#             Df Sum Sq Mean Sq F value Pr(>F)    
#Region        6   2856   475.9   21.83 <2e-16 ***
#Residuals   128   2791    21.8                   
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(res.aov.diversity.isl)
#             Df Sum Sq Mean Sq F value   Pr(>F)    
#Region        6  18.13  3.0210   11.41 3.42e-10 ***
#Residuals   128  33.90  0.2648                     
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##Check assumptions
#Homogeneity of variances
plot(res.aov.richness.all, 1)
plot(res.aov.diversity.all, 1)
plot(res.aov.richness.isl, 1)
plot(res.aov.diversity.isl, 1)
library(car) #v3.1-2
leveneTest(observed ~ Region, data = div.df.all.1)
#Levene's Test for Homogeneity of Variance (center = median)
#       Df F value Pr(>F)
#group   6  0.6936 0.6552
#      128 
leveneTest(Phylogenetic_Diversity ~ Region, data = div.df.all.1)
#Levene's Test for Homogeneity of Variance (center = median)
#       Df F value Pr(>F)
#group   6   1.819 0.1003
#      128   
leveneTest(observed ~ Region, data = div.df.isl.1)
#Levene's Test for Homogeneity of Variance (center = median)
#       Df F value Pr(>F)
#group   6  3.7133 0.001948 **
#      128
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1  
leveneTest(Phylogenetic_Diversity ~ Region, data = div.df.isl.1)
#Levene's Test for Homogeneity of Variance (center = median)
#       Df F value Pr(>F)
#group   6  1.0346  0.406
#      128
       

#Normality
plot(res.aov.richness.all, 2)
plot(res.aov.diversity.all, 2)
plot(res.aov.richness.isl, 2)
plot(res.aov.diversity.isl, 2)

# Extract the residuals
aov.residuals.rich.all <- residuals(object = res.aov.richness.all)
aov.residuals.div.all <- residuals(object = res.aov.diversity.all)
aov.residuals.rich.isl <- residuals(object = res.aov.richness.isl)
aov.residuals.div.isl <- residuals(object = res.aov.diversity.isl)

# Run Shapiro-Wilk test
shapiro.test(x = aov.residuals.rich.all)
#
#Shapiro-Wilk normality test
#
#data:  aov.residuals.rich.all
#W = 0.98898, p-value = 0.3602
#
shapiro.test(x = aov.residuals.div.all)
#
#Shapiro-Wilk normality test
#
#data:  aov.residuals.div.all
#W = 0.98658, p-value = 0.2112
#
shapiro.test(x = aov.residuals.rich.isl)
#
#Shapiro-Wilk normality test
#
#data:  aov.residuals.rich.isl
#W = 0.98651, p-value = 0.2208
#
shapiro.test(x = aov.residuals.div.isl)
#
#Shapiro-Wilk normality test
#
#data:  aov.residuals.div.isl
#W = 0.96045, p-value = 0.0006008
#

##isolate richness did not pass homogeneity of variance test
##isolate diversity did not pass normality test

#While isolate richness showed heteroscedasticity and isolate diversity deviated from normality, we proceeded with ANOVA as it is robust to moderate violations of assumptions with large, balanced samples (Blanca et al., 2017; Schmider et al., 2010), and the highly significant results (p < 0.001) suggest these violations would not affect our overall conclusions.

#Tukey's Honest Significance Multiple Comparison Tests
library(agricolae) #v1.3-5
tukey.richness.all <- HSD.test(res.aov.richness.all, "Region", group = TRUE)
tukey.diversity.all <- HSD.test(res.aov.diversity.all, "Region", group = TRUE)
tukey.richness.isl <- HSD.test(res.aov.richness.isl, "Region", group = TRUE)
tukey.diversity.isl <- HSD.test(res.aov.diversity.isl, "Region", group = TRUE)

print(tukey.richness.all)
print(tukey.diversity.all)
print(tukey.richness.isl)
print(tukey.diversity.isl)


#Add grouping letters to plots
tukey.richness.means.all <- data.frame(tukey.richness.all$means)
tukey.richness.means.all <- tukey.richness.means.all[order(tukey.richness.means.all$observed),]
tukey.richness.means.all
#            observed      std  r Min Max    Q25   Q50    Q75
#Meknes        214.95 40.82114 20 157 315 187.25 202.5 242.75
#Beni Mellal   216.45 43.34012 20 114 277 203.50 222.5 244.00
#Safi          226.80 40.98292 15 153 295 209.50 223.0 252.50
#El Jadida     261.95 34.44519 20 184 315 245.75 262.0 284.25
#Settat        282.35 39.60233 20 222 371 261.50 273.0 296.50
#Sidi Kacem    306.90 42.13625 20 180 371 291.50 313.5 333.00
#Azilal        308.15 26.28843 20 267 358 291.00 303.0 327.75

tukey.diversity.means.all <- data.frame(tukey.diversity.all$means)
tukey.diversity.means.all <- tukey.diversity.means.all[order(tukey.diversity.means.all$Phylogenetic_Diversity),]
tukey.diversity.means.all
#           Phylogenetic_Diversity       std  r      Min      Max      Q25      Q50      Q75
#Beni Mellal               45.21707  7.757261 20 27.46439 57.56071 43.56857 44.75799 48.91090
#Meknes                    47.24622  8.155538 20 36.96991 64.66993 40.43116 45.59013 52.39775
#Safi                      49.51022 10.042488 15 30.68785 67.69746 43.58586 48.00714 58.03299
#El Jadida                 54.52187  7.103530 20 39.89716 63.48863 52.92544 55.05997 60.42666
#Settat                    56.80121  6.048004 20 47.46040 71.25658 52.71857 55.29098 60.20371
#Azilal                    58.56418  3.741129 20 53.34356 68.39926 56.31479 57.51505 61.33129
#Sidi Kacem                59.64888  6.772547 20 39.07121 71.11261 57.23610 60.28085 63.63472


tukey.richness.means.isl <- data.frame(tukey.richness.isl$means)
tukey.richness.means.isl <- tukey.richness.means.isl[order(tukey.richness.means.isl$observed),]
tukey.richness.means.isl
#            observed      std  r Min Max   Q25  Q50   Q75
#Safi        27.26667 2.282438 15  23  31 26.00 27.0 29.00
#Meknes      30.10000 4.700504 20  22  39 26.75 31.0 32.25
#Beni Mellal 34.75000 4.216197 20  27  42 32.50 35.5 37.25
#Settat      34.80000 6.622132 20  25  45 30.00 33.5 40.25
#El Jadida   35.35000 4.749238 20  28  43 31.00 35.5 39.25
#Sidi Kacem  37.00000 4.242641 20  29  46 35.00 36.0 38.25
#Azilal      43.35000 4.331950 20  33  52 40.75 43.0 46.25


tukey.diversity.means.isl <- data.frame(tukey.diversity.isl$means)
tukey.diversity.means.isl <- tukey.diversity.means.isl[order(tukey.diversity.means.isl$Phylogenetic_Diversity),]
tukey.diversity.means.isl
#            Phylogenetic_Diversity       std  r      Min      Max      Q25      Q50      Q75
#Meknes                    3.994503 0.5394646 20 2.940505 4.825906 3.560868 3.904097 4.485351
#Safi                      4.039385 0.5646650 15 2.835085 4.688938 3.611993 4.244582 4.491075
#Beni Mellal               4.346832 0.4654125 20 3.745587 5.094192 3.947251 4.224882 4.807986
#El Jadida                 4.629439 0.4806668 20 3.660344 5.294195 4.366629 4.695268 4.956893
#Settat                    4.650175 0.6329534 20 3.146044 5.593371 4.284745 4.930470 5.058847
#Sidi Kacem                4.710522 0.5429015 20 3.587833 5.403611 4.170170 4.898231 5.158743
#Azilal                    5.099377 0.3755491 20 4.193176 5.516296 5.036280 5.137523 5.385927


sig.groups.obs.all <- data.frame(tukey.richness.all$groups)
sig.groups.obs.all$observed <- rev(sig.groups.obs.all$observed)
sig.groups.div.all <- data.frame(tukey.diversity.all$groups)
sig.groups.div.all$Phylogenetic_Diversity <- rev(sig.groups.div.all$Phylogenetic_Diversity)
sig.groups.obs.isl <- data.frame(tukey.richness.isl$groups)
sig.groups.obs.isl$observed <- rev(sig.groups.obs.isl$observed)
sig.groups.div.isl <- data.frame(tukey.diversity.isl$groups)
sig.groups.div.isl$Phylogenetic_Diversity <- rev(sig.groups.div.isl$Phylogenetic_Diversity)

tukey.richness.means.all$groups <- sig.groups.obs.all$groups
tukey.diversity.means.all$groups <- sig.groups.div.all$groups
tukey.richness.means.isl$groups <- sig.groups.obs.isl$groups
tukey.diversity.means.isl$groups <- sig.groups.div.isl$groups

tukey.richness.means.all$Region <- row.names(tukey.richness.means.all)
tukey.diversity.means.all$Region <- row.names(tukey.diversity.means.all)
tukey.richness.means.isl$Region <- row.names(tukey.richness.means.isl)
tukey.diversity.means.isl$Region <- row.names(tukey.diversity.means.isl)



##Final plots
Fig.1A <- p1 + geom_text(data = tukey.richness.means.all, aes(x=Region, y=observed+std+11, label=groups)) + geom_text(data = tukey.richness.means.isl, aes(x=Region, y=observed+std+18, label=groups))
Fig.1A
Fig.1B <- p2 + geom_text(data = tukey.diversity.means.all, aes(x=Region, y=Phylogenetic_Diversity+std+3, label=groups)) + geom_text(data = tukey.diversity.means.isl, aes(x=Region, y=Phylogenetic_Diversity+std+4, label=groups))
Fig.1B


library(ggpubr)
Fig.1 <- ggarrange(Fig.1A, Fig.1B, nrow = 2, align = "v", labels = c("A", "B"))
Fig.1

# Save as pdf
pdf(file = "Plots/alpha_diversity_plots.pdf",   #file name
    width = 4, # The width of the plot in inches
    height = 8) # The height of the plot in inches
Fig.1
dev.off()
#RStudioGD 
#2 

############ END ###########