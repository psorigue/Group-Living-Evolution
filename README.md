# Molecular Evolution of Social Living in Tanganyika Cichlids
Pol Sorigue(1,2), XXXXXXXXXX, Rui Oliveira (1,2,*)

1 GIMM - Gulbenkian Institute for Molecular Medicine, Rua Quinta Grande 6, 2780-156 Oeiras, Portugal  
2 ISPA - University Institute for Psychological, Social and Life Sciences, Rua do Jardim do Tabaco 34, 1149-041 Lisbon, Portugal  
(*) Corresponding author


## Data availability
The raw sequencing data for this study have been deposited in the European Nucleotide Archive (ENA) at EMBL-EBI under accession number **PRJEB111213**. 

## Annotations Used:
- Nile tilapia: **O_niloticus_UMD_NMBU** (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001858045.2/)


## Scripts index
The scripts assume that the required input files are available in the project home directory structure used by the repository.

## 01. Trait Reconstruction at Node

### [01.1. phenotype_characterization.R](Scripts_clean/01.Trait%20Reconstruction%20at%20Node/01.1.phenotype_characterization.R)
Assigns each species to a binary phenotype (social-living vs solitary) based on the median Nmax value. It writes a phenotype table and prunes the phylogenetic tree so that only species with phenotype information remain.

### [01.2. model_selection_and_node_reconstruction_phytools.R](Scripts_clean/01.Trait%20Reconstruction%20at%20Node/01.2.model_selection_and_node_reconstruction_phytools.R)
Fits ancestral state reconstruction models using phytools, compares the ER and ARD models with AIC and likelihood-ratio tests, reconstructs ancestral states at nodes, and saves the probability tables and a plotted tree.

### [01.3. node_reconstruction_mbasr.R](Scripts_clean/01.Trait%20Reconstruction%20at%20Node/01.3.node_reconstruction_mbasr.R)
Runs MBASR for Bayesian ancestral state reconstruction. This script estimates node-state probabilities using a Bayesian framework and produces the output files and plots for trait reconstruction.

## 02. Read Mapping and Counting

### [02.1. STAR_genome_index_generation.sh](Scripts_clean/02.Read%20Mapping%20and%20Counting/02.1.STAR_genome_index_generation.sh)
Builds the STAR genome index for the reference genome before read alignment.

### [02.2. STAR_mapping.sh](Scripts_clean/02.Read%20Mapping%20and%20Counting/02.2.STAR_mapping.sh)
Maps sequencing reads to the reference genome with STAR, generates BAM files, and indexes them for downstream analyses.

### [02.3. count_table.sh](Scripts_clean/02.Read%20Mapping%20and%20Counting/02.3.count_table.sh)
Uses featureCounts to count reads overlapping annotated features and generates a sample-by-gene count matrix.

## 03. Differential Gene Expression

### [03.1. expression_PCA_and_outliers.R](Scripts_clean/03.Differential%20Gene%20Expression/03.1.expression_PCA_and_outliers.R)
Performs PCA on normalized count data and detects potential outlier samples before differential expression analysis.

### [03.2. differential_gene_expression.R](Scripts_clean/03.Differential%20Gene%20Expression/03.2.differential_gene_expression.R)
Runs DESeq2 to identify differentially expressed genes between specified species/phenotype contrasts and writes the resulting DEG table with gene annotations.

## 04. Functional Enrichment

### [04.1. GSEA_enrichment.R](Scripts_clean/04.Functional%20Enrichment/04.1.GSEA_enrichment.R)
Performs gene set enrichment analysis for GO terms and KEGG pathways using ranked gene lists derived from DESeq2 results.

### [04.2. intraspp_GO_similarity.R](Scripts_clean/04.Functional%20Enrichment/04.2.intraspp_GO_similarity.R)
Computes semantic similarity among enriched GO terms within species-level comparisons, clusters related terms, and summarizes shared clusters.

### [04.3. intraspp_KEGG_similarity.R](Scripts_clean/04.Functional%20Enrichment/04.3.intraspp_KEGG_similarity.R)
Finds KEGG pathways shared between two intraspecies enrichment analyses and reports their descriptions and enrichment statistics.

### [04.4. interspp_GO_similarity.R](Scripts_clean/04.Functional%20Enrichment/04.4.interspp_GO_similarity.R)
Compares GO terms across interspecies contrasts, groups them by similarity, and highlights patterns shared between different comparisons.



## Notes

- Some scripts depend on external R packages and annotation resources, such as AnnotationHub and species-specific annotation databases.
- Output files are written into the directories expected by the workflow, typically under the project home folder structure used in this repository.
