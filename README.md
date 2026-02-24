# Cn_MHC_Evolution

# Genome variant pipeline
# There is a python script per haplotype pair to be cleaner for students to run.
# Files: GenomeAnnot_VariantQ_020240327.py, GenomeAnnot_VariantK_020240327.py, GenomeAnnot_VariantB_020240327.py
Python scripts sort SNPEff vcf files first for quality variants in parent and both evolved strains (e.g., B1 and B2).
Then parent good quality variants are removed from both evolved stain good quality variants.
Then unique and common between the haplotype pair are identified along with the SNPEff effect (hig,h moderate, low) for each gene within the variant' annotation 
Note: Since some genes and regulatory regions overlap, there are typically multiple genes per variant.

# Transcriptome pipeline
# There is a single R script to detect outliers in the transcriptome data sets (when 3 replicates are present).
# File: outlier_minimal_20230404 (1).R
   # This code is my pared down version of the R code by Loraine et al (2015) which is AWESOME! 
   # You should cite them: Their paper is here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4387895/
R code performs clustering of gene expression counts where all genes are in column 1 and all sample IDs are in row 1.
Output is a NMD plot and a hierarchical clustering (tree-like) plot.
