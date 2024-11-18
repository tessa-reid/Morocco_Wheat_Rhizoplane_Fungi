############################################################
#R Scripts written by Tessa Reid tessa.reid@rothamsted.ac.uk
############################################################
#R 4.2.2
#setwd()

library(phyloseq) #v1.42.0
library(rncl) #v0.8.7

dir.create("Data")
dir.create("Plots")


####### Phyloseq objects for amplicon sequences classified by Unite (physeq.all) and by isolates sequences (physeq.isl) #######

#Input files: Files produced by qiime2 (ASVs_abundance.txt, ASVs_taxonomy_unite.txt, ASVs_taxonomy_isolates.txt, tree.nwk) and sample metadata (metadata.txt).



#####Importing the files into R:

#ASV abundance table
otu <- read.table("ASVs_abundance.txt", header=TRUE, row.names=1)

#Taxonomy files
taxonomy.all <- read.csv(file="ASVs_taxonomy_unite.csv", header = TRUE, row.names = 1)
taxonomy.isl <- read.csv(file="ASVs_taxonomy_isolates.csv", header=TRUE, row.names = 1)

#Sample data
sample <- read.csv(file="metadata.csv", header = TRUE, row.names = 1)

#Phylogenetic tree
tree  <- read_newick_phylo(file="tree.nwk",simplify=FALSE)



######Creating the phyloseq object:

#otu table
otu <- as.matrix(otu)
otu.table <- otu_table(otu,taxa_are_rows=TRUE,errorIfNULL = TRUE)


#taxonomy tables

#Unite
row.names(taxonomy.all) <- taxonomy.all$ASV
taxonomy.all <- as.matrix(taxonomy.all)
taxonomy.all.table <- tax_table(taxonomy.all)

#Isolates
row.names(taxonomy.isl) <- taxonomy.isl$ASV
taxonomy.isl <- as.matrix(taxonomy.isl)
taxonomy.isl.table <- tax_table(taxonomy.isl)


#sample data
sample=sample_data(sample,errorIfNULL = TRUE)
# Define factors
sample$Niche = factor(sample$Niche, levels = c("Bulk_soil","Rhizoplane"))
sample$Region = factor(sample$Region, levels = c("Sidi Kacem", "Meknes", "Beni Mellal", "El Jadida", "Settat", "Azilal", "Safi"))
sample$Zone = factor(sample$Zone)
sample$Zone_abbrv = factor(sample$Zone_abbrv)
sample$SequencingRun = factor(sample$SequencingRun)


#1. physeq object - Unite
physeq.all <- phyloseq(otu.table, taxonomy.all.table, sample, tree)
physeq.all
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 4568 taxa and 135 samples ]
#sample_data() Sample Data:       [ 135 samples by 7 sample variables ]
#tax_table()   Taxonomy Table:    [ 4568 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 4568 tips and 4538 internal nodes ]
saveRDS(physeq.all, "Data/physeq_all.rds")

#2. physeq object - Isolates
physeq.isl <- phyloseq(otu.table, taxonomy.isl.table, sample, tree)
physeq.isl
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 4568 taxa and 135 samples ]
#sample_data() Sample Data:       [ 135 samples by 7 sample variables ]
#tax_table()   Taxonomy Table:    [ 4568 taxa by 10 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 4568 tips and 4538 internal nodes ]
saveRDS(physeq.isl, "Data/physeq_isl.rds")



#Check Data

sum(colSums(otu_table(physeq.all)))
#[1] 7675765
sum(colSums(otu_table(physeq.isl)))
#[1] 7675765

#get phylum levels
get_taxa_unique(physeq.all, "Phylum")
#[1] "Mucoromycota"       "Chytridiomycota"    "Ascomycota"         "Blastocladiomycota"
#[5] "Rozellomycota"      "Basidiomycota"      "unidentified"       "Zoopagomycota"     
#[9] "Aphelidiomycota"    "Monoblepharomycota" "Glomeromycota"      "Mortierellomycota" 
#[13] NA                   "Kickxellomycota"    "Olpidiomycota"      "Basidiobolomycota" 

get_taxa_unique(physeq.isl, "Phylum")
#[1] "Unassigned"   "Ascomycota"   "Mucoromycota"


# Create table, number of features for each kingdom
table(tax_table(physeq.all)[, "Kingdom"], exclude = NULL)
#     Fungi 
#     4568 

table(tax_table(physeq.isl)[, "Kingdom"], exclude = NULL)
#     Fungi Unassigned 
#       348       4220 



### Phyloseq object with only ASVs classifed as Strains

#Filter to only include ASVs identified as Strains

#remove unassigned ASVs
physeq.str <- subset_taxa(physeq.isl, !is.na(Strain) & !Strain %in% c("Unassigned"))


#check phlum levels
get_taxa_unique(physeq.str, "Strain")
#[1] "Actinomucor elegans strain UAMH 10931 [F1]"     "Rhizopus arrhizus strain CMRC 585 [F12]"       
#[3] "Neodidymelliopsis ranunculi MFLU 16-1870 [F18]" "Alternaria infectoria [F19]"                   
#[5] "Aspergillus montevidensis [F20]"                "Aspergillus niger isolate 26Crab [F6]"         
#[7] "Aspergillus ochraceus isolate Contig_TVS [F16]" "Talaromyces oumae-annae clone SF_862 [F2]"     
#[9] "Aspergillus niger voucher HQU AR4 [F3]"         "Aspergillus flavus isolate LUOHE [F15]"        
#[11] "Penicillium aurantiogriseum strain YIJZ-5 [F4]" "Aspergillus calidoustus strain K 2-12 [F17]"
#[13] "Fusarium graminearum isolate wheat [F10]"


#3. physeq object - Strains
physeq.str
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 232 taxa and 135 samples ]
#sample_data() Sample Data:       [ 135 samples by 7 sample variables ]
#tax_table()   Taxonomy Table:    [ 232 taxa by 10 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 232 tips and 229 internal nodes ]
saveRDS(physeq.str, "Data/physeq_str.rds")



# Normalise abundances by rarefication

min.lib.all<-min(sample_sums(physeq.all)) #24766
min.lib.str<-min(sample_sums(physeq.str)) #7873


physeq.all.rare<-rarefy_even_depth(physeq.all, min.lib.all, replace=TRUE, rngseed = 9242)
#`set.seed(9242)` was used to initialize repeatable random subsampling.
#Please record this for your records so others can reproduce.
#Try `set.seed(9242); .Random.seed` for the full vector
#...
#111OTUs were removed because they are no longer
#present in any sample after random subsampling

#...

#4. physeq object - rarefied asv abundances - unite
physeq.all.rare
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 4457 taxa and 135 samples ]
#sample_data() Sample Data:       [ 135 samples by 7 sample variables ]
#tax_table()   Taxonomy Table:    [ 4457 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 4457 tips and 4428 internal nodes ]
saveRDS(physeq.all.rare, "Data/physeq_all_rare.rds")


physeq.str.rare<-rarefy_even_depth(physeq.str, min.lib.str, replace=TRUE, rngseed = 9242)
#`set.seed(9242)` was used to initialize repeatable random subsampling.
#Please record this for your records so others can reproduce.
#Try `set.seed(9242); .Random.seed` for the full vector
#...
#2OTUs were removed because they are no longer
#present in any sample after random subsampling
#
#...

#5. physeq object - rarefied asv abundances - strains
physeq.str.rare
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 230 taxa and 135 samples ]
#sample_data() Sample Data:       [ 135 samples by 7 sample variables ]
#tax_table()   Taxonomy Table:    [ 230 taxa by 10 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 230 tips and 226 internal nodes ]
saveRDS(physeq.str.rare, "Data/physeq_str_rare.rds")



# Normalise abundances by deseq2 variance stabilising transforming

library(DESeq2)

#UNITE
dds = phyloseq_to_deseq2(physeq.all, ~  Region)

dds <- dds[ rowSums(counts(dds)) > 3, ]
cts <- counts(dds)
geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(sum(log(row[row != 0]))/length(row)))
dds <- estimateSizeFactors(dds, geoMeans=geoMeans)


NormVST<- t(assay(varianceStabilizingTransformation(dds, blind=TRUE)))
NormVST[NormVST < 0] <- 0.0


physeq.all.deseq = physeq.all
otu_table(physeq.all.deseq) <- otu_table(NormVST, FALSE)

# transpose otu table
otu.norm = as(otu_table(physeq.all.deseq), "matrix")
otu.norm = t(otu.norm)
otu.norm.table <- otu_table(otu.norm,taxa_are_rows=TRUE,errorIfNULL = TRUE)


#6. physeq object - asv abundances normalised by VST - unite
physeq.all.deseq <- phyloseq(otu.norm.table, taxonomy.all.table, sample, tree)
physeq.all.deseq
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 4393 taxa and 135 samples ]
#sample_data() Sample Data:       [ 135 samples by 7 sample variables ]
#tax_table()   Taxonomy Table:    [ 4393 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 4393 tips and 4364 internal nodes ]
saveRDS(physeq.all.deseq, "Data/physeq_all_deseq.rds")



#ISOLATES
dds = phyloseq_to_deseq2(physeq.isl, ~  Region)

dds <- dds[ rowSums(counts(dds)) > 3, ]
cts <- counts(dds)
geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(sum(log(row[row != 0]))/length(row)))
dds <- estimateSizeFactors(dds, geoMeans=geoMeans)


NormVST<- t(assay(varianceStabilizingTransformation(dds, blind=TRUE)))
NormVST[NormVST < 0] <- 0.0


physeq.isl.deseq = physeq.isl
otu_table(physeq.isl.deseq) <- otu_table(NormVST, FALSE)


# transpose otu table
otu.norm = as(otu_table(physeq.isl.deseq), "matrix")
otu.norm = t(otu.norm)
otu.norm.table <- otu_table(otu.norm,taxa_are_rows=TRUE,errorIfNULL = TRUE)


#7. physeq object - asv abundances normalised by VST - isolates
physeq.isl.deseq <- phyloseq(otu.norm.table, taxonomy.isl.table, sample, tree)
physeq.isl.deseq
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 4393 taxa and 135 samples ]
#sample_data() Sample Data:       [ 135 samples by 7 sample variables ]
#tax_table()   Taxonomy Table:    [ 4393 taxa by 10 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 4393 tips and 4364 internal nodes ]
saveRDS(physeq.isl.deseq, "Data/physeq_isl_deseq.rds")



#8. physeq object - asv abundances normalised by VST - strains
physeq.str.deseq <- subset_taxa(physeq.isl.deseq, !is.na(Strain) & !Strain %in% c("Unassigned"))
physeq.str.deseq
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 231 taxa and 135 samples ]
#sample_data() Sample Data:       [ 135 samples by 7 sample variables ]
#tax_table()   Taxonomy Table:    [ 231 taxa by 10 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 231 tips and 228 internal nodes ]
saveRDS(physeq.str.deseq, "Data/physeq_str_deseq.rds")


### ASV abundances pre- and post normalization
# Make a data frame with a column for the read counts of each sample
sample_sum_1 <- data.frame(sum = sample_sums(physeq.all))
sample_sum_2 <- data.frame(sum = sample_sums(physeq.all.rare))
sample_sum_3 <- data.frame(sum = sample_sums(physeq.all.deseq))
sample_sum_4 <- data.frame(sum = sample_sums(physeq.str))
sample_sum_5 <- data.frame(sum = sample_sums(physeq.str.rare))
sample_sum_6 <- data.frame(sum = sample_sums(physeq.str.deseq))


# Histogram of sample read counts
library(ggplot2)
theme_set(theme_bw())

p1 = ggplot(sample_sum_1, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("physeq.all") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

p2 = ggplot(sample_sum_2, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("physeq.all.rare") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

p3 = ggplot(sample_sum_3, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("physeq.all.deseq") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

p4 = ggplot(sample_sum_4, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("physeq.str") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

p5 = ggplot(sample_sum_5, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("physeq.str.rare") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

p6 = ggplot(sample_sum_6, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("physeq.str.deseq") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

#Final plot
library(ggpubr) #v0.6.0
norm.plot <- ggarrange(p1, p2, p3, p4, p5, p6)

# Save as pdf
pdf(file = "Plots/normalisation_plots.pdf",   #file name
    width = 10, # The width of the plot in inches
    height = 5) # The height of the plot in inches
norm.plot
dev.off()
#RStudioGD 
#2



#Rarefication Curves
suppressMessages(library(MicrobiotaProcess)); packageVersion("MicrobiotaProcess") #v1.8.2

#rarefaction curve
rareres.all <- get_rarecurve(obj=physeq.all, chunks=400)
#There were 50 or more warnings (use warnings() to see the first 50)
colour <- c("#4CB5F5",  "#D70026", "#EDB83D",  "navy",  "#F77604", "#B3C100", "#EC96A4")
prare1 <- ggrarecurve(obj=rareres.all,
                      factorNames="Region",
                      shadow=FALSE,
                      indexNames=c("Observe")) + 
  scale_color_manual(values = colour) +
  labs(title = "Unite", x="Number of reads", y="Observed ASV richness") +
  theme(axis.text=element_text(size=8), 
        panel.grid=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
prare1

#rarefaction curve
rareres.str <- get_rarecurve(obj=physeq.str, chunks=400)
#There were 50 or more warnings (use warnings() to see the first 50)
prare2 <- ggrarecurve(obj=rareres.str,
                      factorNames="Region",
                      shadow=FALSE,
                      indexNames=c("Observe")) + 
  scale_color_manual(values = colour) +
  labs(title = "Strains", x="Number of reads", y="Observed ASV richness") +
  theme(axis.text=element_text(size=8), 
        panel.grid=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
prare2

library(ggpubr); packageVersion("ggpubr") #v0.6.0
rrare.plot <- ggarrange(prare1, prare2, common.legend = TRUE, legend = "right", labels = c("A", "B"), font.label = list(size=12, face="plain"))
rrare.plot

# Save as pdf
pdf(file = "Plots/rareplots.pdf",   #file name
    width = 8, # The width of the plot in inches
    height = 4) # The height of the plot in inches
rrare.plot
dev.off()
#RStudioGD 
#2 


################## END ##################