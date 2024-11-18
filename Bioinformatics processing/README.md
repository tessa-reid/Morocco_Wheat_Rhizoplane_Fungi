---
title: "Benbrik et al., 2025 - Code for bioinformatic analysis"
author: "Tessa Reid"
date: "22/08/2023"
---
# Code for Bioinformatic Processing of Amplicon Sequence Data


# 1. Processing of ITS region amplicon dataset from merged reads
### In Qiime2
### Import fastq files in Qiime2 as single-end reads
```{r}
source activate qiime2-2019.7

qiime tools import \
    --type   'SampleData[SequencesWithQuality]' \
    --input-path 'Data to import' \
    --input-format CasavaOneEightSingleLanePerSampleDirFmt \
    --output-path demux-single-end.qza
```

## Visualize and verify sequence quality. Open qzv file on Qiime2View online
```{r}
qiime demux summarize \
    --i-data demux-single-end.qza \
    --o-visualization demux-single-end.qzv
```


## Sequence Denoising with DADA2
```{r}
qiime dada2 denoise-single \
    --i-demultiplexed-seqs demux-single-end.qza \
    --p-trunc-len 0 \
    --p-trim-left 2 \
    --o-representative-sequences rep-seqs.qza \
    --o-table table.qza \
    --o-denoising-stats stats.qza \
    --p-n-threads 32 \
    --verbose
```

## Visualize the outputs of DADA2 in Qiime2View online
```{r}
qiime metadata tabulate \
    --m-input-file stats.qza \
    --o-visualization stats.qzv

qiime feature-table summarize \
    --i-table table.qza \
    --o-visualization table.qzv \
    --m-sample-metadata-file metadata.tsv

qiime feature-table tabulate-seqs \
    --i-data rep-seqs.qza \
    --o-visualization rep-seqs.qzv
```

## Remove ASVs with low abundance (less than 10 total observation count) from ASV table and rep-seq file
```{r}
qiime feature-table filter-features \
    --i-table table.qza \
    --p-min-frequency 10 \
    --o-filtered-table table-filtered.qza

qiime feature-table filter-seqs \
    --i-data rep-seqs.qza \
    --i-table table-filtered.qza \
    --o-filtered-data rep-seqs-filtered.qza
```

## Visualize the outputs in Qiime2View online
```{r}
qiime feature-table summarize \
    --i-table table-filtered.qza \
    --o-visualization table-filtered.qzv \
    --m-sample-metadata-file metadata.tsv
```

## Assign taxonomy with Qiime2 using Unite 2022 database
```{r}
qiime feature-classifier classify-consensus-vsearch \
    --i-query rep-seqs-filtered.qza \
    --i-reference-reads unite_2022_ASVs.qza \
    --i-reference-taxonomy unite_2022_taxonomy.qza \
    --o-classification taxonomy.qza \
    --o-search-results blast.qza \
    --p-threads 30 \
    --verbose
```

## Visualise taxa results to check for non-fungal ASVs. Open qzv file on Qiime2View online
```{r}
qiime taxa barplot \
    --i-table table-filtered.qza \
    --i-taxonomy taxonomy.qza \
    --m-metadata-file metadata.tsv \
    --o-visualization taxa-bar-plots.qzv
```

## Remove Unassigned ASVs from ASV table and rep-seqs file
## Re-assign taxonomy with filtered rep-seqs file to produce a taxonomy table with only fungal ASVs


## Make phylogenic tree in Qiime2 from filtered rep-seqs file (unassigned, low abundance sequences removed)
```{r}
qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences rep-seqs-fungi.qza \
    --o-alignment aligned-rep-seqs.qza \
    --o-masked-alignment masked-aligned-rep-seqs.qza \
    --o-tree unrooted-tree.qza \
    --o-rooted-tree rooted-tree.qza \
    --verbose
```


## Export ASV table, taxonomy table, and tree for analysis in phyloseq

### Taxonomy Table
```{r}
qiime tools export \
  --input-path taxonomy-fungi.qza \
  --output-path phyloseq_unite
```

### Phylogenetic Tree
```{r}
qiime tools export \
  --input-path rooted-tree.qza \
  --output-path phyloseq_unite
```

### ASV Table
```{r}
qiime tools export \
  --input-path table-filtered-fungi.qza \
  --output-path phyloseq_unite
```

### Convert Biom files to tsv
```{r}
biom convert \
  -i phyloseq_unite/feature-table.biom \
  -o phyloseq_unite/otu_table_unite.txt \
  --to-tsv
```


# 2. Assigning reads to cultured fungal isolate ASVs

## Isolate database creation in Qiime2

### Convert the taxonomy file into a Qiime2 file
```{r}
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path isolate-taxonomy.txt \
  --output-path isolate-taxonomy.qza
```

### Convert the ASV file into a Qiime2 file
```{r}
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path isolate-ASVs.txt \
  --output-path isolates-ASVs.qza
```


## Assign taxonomy with Qiime2 using the isolate database
```{r}
qiime feature-classifier classify-consensus-vsearch \
    --i-query rep-seqs-fungi.qza \
    --i-reference-reads isolates-ASVs.qza \
    --i-reference-taxonomy isolate-taxonomy.qza \
    --o-classification isolate_taxonomy.qza \
    --o-search-results isolate_blast.qza \
    --p-threads 32 \
    --verbose
```


## Export taxonomy table for analysis in phyloseq

### Taxonomy Table
```{r}
qiime tools export \
  --input-path taxonomy-fungi.qza \
  --output-path phyloseq_isolates
```
