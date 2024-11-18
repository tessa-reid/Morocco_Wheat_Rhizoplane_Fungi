############################################################
R Scripts written by Tessa Reid tessa.reid@rothamsted.ac.uk
############################################################

######## Strain abundance in zones

#input files: 
#1. phyloseq object with ASVs classified by UNITE with abundances normalised by DeSEQ2 (physeq_all_deseq.rds)
#2. phyloseq object with ASVs classified by Isolate sequences with abundances normalised by DeSEQ2 (physeq_isl_deseq.rds)


#Import data
library(phyloseq)
physeq.all <- readRDS("Data/physeq_all_deseq.rds")
physeq.isolates <- readRDS("Data/physeq_isl_deseq.rds")
physeq.all
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 4393 taxa and 135 samples ]
#sample_data() Sample Data:       [ 135 samples by 7 sample variables ]
#tax_table()   Taxonomy Table:    [ 4393 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 4393 tips and 4364 internal nodes ]
physeq.isolates
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 4393 taxa and 135 samples ]
#sample_data() Sample Data:       [ 135 samples by 7 sample variables ]
#tax_table()   Taxonomy Table:    [ 4393 taxa by 10 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 4393 tips and 4364 internal nodes ]



#1. Kingdom - Unite vs. isolate

#Unite
#devtools::install_github("vmikk/metagMisc")
library("metagMisc") #v.0.0.4
taxa.abd.k <- phyloseq_to_df(physeq.all, addtax = T, addtot = F, addmaxrank = F,sorting = "abundance")
taxa.abd.k<-taxa.abd.k[,-c(1:2, 4:9)] #remove other taxa columns
taxa.abd.k<-setNames(data.frame(t(taxa.abd.k[ , - 1])), taxa.abd.k[ , 1]) 
#taxa.abd.k$sam_name <- rownames(taxa.abd.k)


#Isolate
taxa.abd.k.isl <- phyloseq_to_df(physeq.isolates, addtax = T, addtot = F, addmaxrank = F,sorting = "abundance")
taxa.abd.k.isl<-taxa.abd.k.isl[,-c(1:2, 4:11)]
taxa.abd.k.isl<-setNames(data.frame(t(taxa.abd.k.isl[ , - 1])), taxa.abd.k.isl[ , 1]) 
#taxa.abd.k.isl$sam_name <- rownames(taxa.abd.k.isl)


#Sample data
taxa.meta <- sample_data(physeq.all)


#Dataframes
taxa.df.k <- merge(taxa.abd.k, taxa.meta)
taxa.df.k.isl <- merge(taxa.abd.k.isl, taxa.meta)




#Prepare for plotting

#Unite
library(dplyr)
taxa.df.k <- physeq.all %>%
  tax_glom(taxrank = "Kingdom") %>%
  psmelt()

head(taxa.df.k)
str(taxa.df.k)

taxa.df.k <- taxa.df.k[, c("Niche","Code","Region", "Zone",  "Kingdom", "Abundance")]
colnames(taxa.df.k) <- c("Niche","Sample","Region","Zone","ASV", "Abundance")
taxa.df.k$ASV = factor(taxa.df.k$ASV)
taxa.df.k$Dataset <- "Unite"
head(taxa.df.k)


#Isolates
taxa.df.k.isl <- physeq.isolates %>%
  tax_glom(taxrank = "Kingdom") %>%
  psmelt()

head(taxa.df.k.isl)
str(taxa.df.k.isl)

taxa.df.k.isl <- taxa.df.k.isl[, c("Niche", "Code", "Region", "Zone", "Kingdom", "Abundance")]
colnames(taxa.df.k.isl) <- c("Niche", "Sample", "Region", "Zone", "ASV", "Abundance")
taxa.df.k.isl$ASV = factor(taxa.df.k.isl$ASV)
taxa.df.k.isl$Dataset <- "Isolate"
head(taxa.df.k.isl)


#Plot
library(ggplot2)
taxa.bar.k <- ggplot() + 
  geom_bar(data=taxa.df.k, aes(x=Dataset, y=Abundance, fill=ASV),position="fill", stat="identity") +
  geom_bar(data=taxa.df.k.isl, aes(x=Dataset, y=Abundance, fill=ASV),position="fill", stat="identity") +
  theme_classic() +
  xlab("Dataset") +
  ylab("") +
  scale_fill_manual(values=c("black", "grey")) +
  theme(axis.title = element_text(color="black", size=12)) + 
  theme(axis.text = element_text(color="black", size=8)) + 
  theme(legend.text = element_text(colour="black", size = 12)) + 
  theme(legend.title = element_text(colour="black", size=12)) +
  theme(legend.position = "bottom") +
  theme(axis.line.x = element_blank(), axis.line.y = element_blank(), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())
taxa.bar.k

# Save as pdf
pdf(file = "Plots/kingdom_abundance.pdf",   #file name
    width = 4, # The width of the plot in inches
    height = 5 # The height of the plot in inches
) 
taxa.bar.k
dev.off()
#RStudioGD 
#2 


# Calculate percentages by ASV
total.abundance <- sum(taxa.df.k$Abundance) 

#Unite
asv.percentages.k <- aggregate(Abundance ~ ASV, data=taxa.df.k, FUN=sum)
asv.percentages.k$Percentage <- (asv.percentages.k$Abundance / total.abundance) * 100
asv.percentages.k
#    ASV Abundance Percentage
#1 Fungi  188865.2        100

#Isolate
asv.percentages.k.isl <- aggregate(Abundance ~ ASV, data=taxa.df.k.isl, FUN=sum)
asv.percentages.k.isl$Percentage <- (asv.percentages.k.isl$Abundance / total.abundance) * 100
asv.percentages.k.isl
#         ASV Abundance Percentage
#1      Fungi   41265.1   21.84896
#2 Unassigned  147600.1   78.15104



#2. What isolates were identified in the dataset?

physeq.isolates.s <- readRDS("Data/physeq_str_deseq.rds")
physeq.isolates.s
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 231 taxa and 135 samples ]
#sample_data() Sample Data:       [ 135 samples by 7 sample variables ]
#tax_table()   Taxonomy Table:    [ 231 taxa by 10 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 231 tips and 228 internal nodes ]


#get strain levels
get_taxa_unique(physeq.isolates.s, "Strain")
#[1] "Actinomucor elegans strain UAMH 10931 [F1]"    
#[2] "Rhizopus arrhizus strain CMRC 585 [F12]"       
#[3] "Neodidymelliopsis ranunculi MFLU 16-1870 [F18]"
#[4] "Alternaria infectoria [F19]"                   
#[5] "Aspergillus montevidensis [F20]"               
#[6] "Aspergillus niger isolate 26Crab [F6]"         
#[7] "Aspergillus ochraceus isolate Contig_TVS [F16]"
#[8] "Talaromyces oumae-annae clone SF_862 [F2]"     
#[9] "Aspergillus niger voucher HQU AR4 [F3]"        
#[10] "Aspergillus flavus isolate LUOHE [F15]"        
#[11] "Penicillium aurantiogriseum strain YIJZ-5 [F4]"
#[12] "Aspergillus calidoustus strain K 2-12 [F17]"   
#[13] "Fusarium graminearum isolate wheat [F10]"     



#Strain abundance dataframe
taxa.abd.isl <- phyloseq_to_df(physeq.isolates.s, addtax = T, addtot = F, addmaxrank = F,sorting = "abundance")
writexl::write_xlsx(taxa.abd.isl,"Data/taxa.abd.isl.xlsx")
taxa.abd.isl<-taxa.abd.isl[,-c(1:9, 11)]
taxa.abd.isl<-setNames(data.frame(t(taxa.abd.isl[ , - 1])), taxa.abd.isl[ , 1]) 
taxa.df.isl <- cbind(taxa.abd.isl, taxa.meta)


#Prepare for plotting
taxa.df.isl.melt <- reshape2::melt(taxa.df.isl)
taxa.df.isl.melt.1 <- taxa.df.isl.melt[, c("Niche","Code","Region", "Zone_abbrv",  "variable", "value")]
colnames(taxa.df.isl.melt.1) <- c("Niche","Sample","Region","Zone","IsolateID", "Abundance")
taxa.df.isl.melt.1$IsolateID = factor(taxa.df.isl.melt.1$IsolateID)
taxa.df.isl.melt.2 <- subset(taxa.df.isl.melt.1, Abundance>0)



taxa.bar.strains <- ggplot() + 
  geom_bar(data=taxa.df.isl.melt.2, aes(x=Zone, y=Abundance, fill=IsolateID), position="fill", stat="identity") +
  theme_bw() +
  facet_grid(. ~ Region, scales = "free", space = "free") +
  xlab("Zone") +
  ylab("relative abundance") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
  scale_x_discrete(limits=rev) + theme(axis.title = element_text(color="black", size=12)) + 
  theme(axis.text = element_text(color="black", size=12)) + 
  theme(legend.position = "none") +
  theme(legend.text = element_text(colour="black", size = 10, face = "italic")) + 
  theme(legend.title = element_text(colour="black", size=12)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(strip.background = element_blank()) +
  theme(strip.text.x = element_text(colour="black", size = 12))
taxa.bar.strains

# Save as pdf
pdf(file = "Plots/strain_abundance.pdf",   #file name
    width = 8, # The width of the plot in inches
    height = 5 # The height of the plot in inches
    ) 
taxa.bar.strains
dev.off()
#RStudioGD 
#2 




#3. Identify the 231 ASVs identified as strains in the Unite dataset to see if the abundance and classification were similar


## Get vectors of numbered OTUs/ASVs
all <- rownames(otu_table(physeq.all))
strains <- rownames(otu_table(physeq.isolates.s))


## Get the intersection between treatments
shared.asvs <- intersect(all, strains)


## Subset phyloseq object to shared taxa
subset.asvs <- subset(otu_table(physeq.all), rownames(otu_table(physeq.all)) %in%  shared.asvs)
physeq.all.shared <- merge_phyloseq(subset.asvs, tax_table(physeq.all), sample_data(physeq.all))
physeq.all.shared


#get genus levels
get_taxa_unique(physeq.all.shared, "Genus")
#[1] "Actinomucor"  "Rhizopus"     NA             "Didymella"    "Alternaria"   "Aspergillus"
#[7] "Talaromyces"  "Fusarium"     "unidentified"


## Genus dataframe
taxa.abd.s <- phyloseq_to_df(physeq.all.shared, addtax = T, addtot = F, addmaxrank = F,sorting = "abundance")
#writexl::write_xlsx(taxa.abd.s,"Data/taxa.abd.s.xlsx")
taxa.abd.s<-taxa.abd.s[,-c(1:7, 9)]
taxa.abd.s<-setNames(data.frame(t(taxa.abd.s[ , - 1])), taxa.abd.s[ , 1]) 
#taxa.abd.s$sam_name <- rownames(taxa.abd.s)
taxa.df.s <- cbind(taxa.abd.s, taxa.meta)



## Prepare for plotting
taxa.df.s.melt <- reshape2::melt(taxa.df.s)
taxa.df.s.melt.1 <- taxa.df.s.melt[, c("Niche","Code","Region", "Zone_abbrv",  "variable", "value")]
colnames(taxa.df.s.melt.1) <- c("Niche","Sample","Region","Zone","Genus", "Abundance")
taxa.df.s.melt.1$Genus = factor(taxa.df.s.melt.1$Genus, levels = c("Fusarium", "Rhizopus", "Alternaria", "Aspergillus", "Actinomucor", "Talaromyces", "Didymella", "unidentified", "NA"))
taxa.df.s.melt.2 <- subset(taxa.df.s.melt.1, Abundance>0)
taxa.df.s.melt.2$Dataset <- "Unite"


## Plot
taxa.bar.all.g <- ggplot() + 
  geom_bar(data=taxa.df.s.melt.2, aes(x=Zone, y=Abundance, fill=Genus),position="fill", stat="identity") +
  theme_bw() +
  facet_grid(. ~ Region, scales = "free", space = "free") +
  xlab("Zone") +
  ylab("relative abundance") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
  scale_x_discrete(limits=rev) + theme(axis.title = element_text(color="black", size=12)) + 
  theme(axis.text = element_text(color="black", size=12)) + 
  theme(legend.position = "none") +
  theme(legend.text = element_text(colour="black", size = 10, face = "italic")) + 
  theme(legend.title = element_text(colour="black", size=12)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(strip.background = element_blank()) +
  theme(strip.text.x = element_text(colour="black", size = 12))
taxa.bar.all.g

# Save as pdf
pdf(file = "Plots/genus_abundance.pdf",   #file name
    width = 8, # The width of the plot in inches
    height = 5 # The height of the plot in inches
) 
taxa.bar.all.g
dev.off()
#RStudioGD 
#2 




#4. Total fungal community genera

#get genus levels
get_taxa_unique(physeq.all, "Genus")
#313 genera

#get ASV numbers
table(tax_table(physeq.all)[, "Genus"], exclude = NULL)

table(tax_table(physeq.isolates.s)[, "Genus"], exclude = NULL)
#
#  Actinomucor        Alternaria       Aspergillus          Fusarium Neodidymelliopsis       Penicillium 
#           10                88                41                42                44                 1 
#     Rhizopus       Talaromyces 
#            4                 1 
table(tax_table(physeq.all.shared)[, "Genus"], exclude = NULL)
#
#  Actinomucor   Alternaria  Aspergillus    Didymella     Fusarium     Rhizopus  Talaromyces unidentified 
#           10           88           42            2           27            4            1            9 
#         <NA> 
#           48 

############ END ###########