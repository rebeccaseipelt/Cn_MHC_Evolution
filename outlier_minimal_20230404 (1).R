################################################################################
###    Background                  
################################################################################ 
# This is a general  file to generate a MDS plot and cluster dendrogram for RNA sequencing data 
# My name is Rebecca Seipelt-Thiemann [rebecca.seipelt@mtsu.edu], a faculty member in the Biology Department at Middle Tennessee State University
# Today's date is April 5, 2023
   # This code is my pared down version of the R code by Loraine et al (2015) which is AWESOME! 
   # You should cite them: Their paper is here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4387895/

# Input 
   # A single txt file with gene names/ID in column 1 and sample counts as the remaining columns with a header row that has a short name indicating the group (e.g., control/treatment) and replicate number.  
   # You'll want to write some notes about your file and how it was generated. 

# Questions to answer in this data analysis
   # How similar were the samples within or across groups?

# Analysis
   # To answer these questions, we'll import the edgeR package from the Bioconductor suite of tools. 
   # Also see the edgeR manual and vignettes.
   # good resource for bioconductor and packages: https://www.bioconductor.org/install/
  
################################################################################
###    Getting started - installing and loading packages you'll need                  
################################################################################ 
  
# new bioconductor - latest release
# check here for the version number and replace 3.16 with that number: https://www.bioconductor.org/install/
# let it install whatever it needs; may take a while if this is your first R code
# you may need to answer yes or y to questions in the console panel
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version = "3.16")


# install specific package edgeR ; this only needs done once (like buying the ingredients for cookies)
BiocManager::install("edgeR")

# load the library (this needs done each time you run the code - it is like taking the cookie ingredients out of the pantry)
library(edgeR)

################################################################################
###    Setting the working directory                   
################################################################################
 
#IMPORTANT STEP THAT IS NOT CODE!!!!  To set your working directly, go to the session tab at top > set working directory > source file location
# If this is not set correctly, you will get an error saying something like "can't find the file"

#check your working directory - your R file and data file should BOTH be in this folder
getwd()

################################################################################
###    Clear out all variables/storage                  
################################################################################

# to make sure you have all the space and clean variables, remove everything from your environment
rm(list=ls())  

################################################################################
###    Read from your file into R                 
################################################################################

#now read the data into a dataframe that we'll call d
#BE SURE TO replace the file name below with your file name including the file extension
# EXAMPLE 1 d=read.delim('HEK_STAR_joined_reoganized_20220216.txt')
# EXAMPLE 2 d=read.delim('HEK_FeatureCount_joined_reorganized_20220216.txt')
d=read.delim('CnRNAccombine_result_withnames.txt', header = TRUE)

# if you get errors with formatting here, open your file in excel and format the count data cells to numbers with 0 decimals

#to look at the beginning of the file to check it looks ok
head(d)

#to make sure d is the correct data type, should be list
mode(d)

################################################################################
###    remove column 1 with the gene names so there is just numerical data                
################################################################################

#to get just values; replace 37 with the column number appropriate for your data; 
# example has 36 data columns so 37
# if example has 9 data columns use 2:10
just_values<-d[,2:37]

#now look at the data again, gene name column should be gone
head(just_values)

################################################################################
###    define the experiment group names for your data columns                
################################################################################

# define a new character vector using the "c" function
# must be in the column order in the file!
# example: I have 2 replicates for 12600 and three replicates for 12602, etc...
group=c('12600','12600','12602','12602','12602','12658','12658','12658','12703','12703','12703','12823','12823','12824','12938','12938','12938','12949','12949','12949','14010','14010','14012','14012','14195','14195','14195','14203','14203','14203','14332','14332','14332','H99W','H99W','H99W')

################################################################################
###    define your dge list and associate with group names and remove rows with all zeros                
################################################################################
dge=DGEList(just_values,group=group,remove.zeros=TRUE)

################################################################################
###    calculate the normalization factor for each library so we have a fair comparison among groups                
################################################################################

dge=calcNormFactors(dge)

#look at the normalization factors for each library
dge$samples 

################################################################################
###    pick the colors for your groups - replicates wil be same color to help you see outliers          
################################################################################

# you can use lots of colors, each sample type should have its own color
# example for 14 condition experiment with some "pairs" of data so similar colors
# a nice color "cheatsheet": https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf


#b
s12600.color='blue'
s12602.color='cyan3'

#f  
s12658.color='pink'

#k
s12703.color='red'
s12938.color='firebrick'    

#BALBc  
s12823.color='black'
s12824.color='gray'
  
#Alt_q_b_q
s12949.color='darkorchid4'
#Alt_d-q
s14203.color='darkorchid1'
  
#q
s14010.color='gold'
s14012.color='orange'

#d  
s14195.color='burlywood4'
s14332.color='brown'
  
#parent
h99w.color = 'green'

################################################################################
###    making the MDS plot          
################################################################################

# define a title for the plot
main='MDS Plot for CnEvol Count Data (RLST)'
# makes y axis labels horizontal not vertical
par(las=1) 

# Assign the colors to the samples and replicates
# example for a 5 condition experiment with 3 replicates each
# colors=c(rep(m1m.color,3),rep(m1p.color,3),rep(m2m.color,3),rep(m2p.color,3),rep(m0.color,3))

colors=c(rep(s12600.color,2),rep(s12602.color,3),rep(s12658.color,3),rep(s12703.color,3),rep(s12823.color,2),rep(s12824.color,1),rep(s12938.color,3),rep(s12949.color,3),rep(s14010.color,2),rep(s14012.color,2),rep(s14195.color,3),rep(s14203.color,3),rep(s14332.color,3),rep(h99w.color,3))

#finally, make the plot (will show in window at right under Plot tab)
plotMDS(dge,main=main,labels=colnames(dge$counts),
        col=colors,las=1)

################################################################################
###    making the cluster dendogram - another way to look at it (tree-like)          
################################################################################

### Hierarchical clustering
#You can use hierarchical clustering to view the relationships between samples. This approach calculates a distance value between each sample using measurements from all the genes. Samples are then clustered according to how similar they are with respect to that distance metric.

#Note we need to use normalized "counts per million" instead of the raw count values. Also, note that we have to transpose the matrix so that columns become rows and rows become columns. In other words, to do the clustering, the entities that should be clustered need to be rows, not columns.
normalized.counts=cpm(dge)

# transposes the counts matrix
transposed=t(normalized.counts) 

# calculates distance
distance=dist(transposed) 

# does hierarchical clustering
clusters=hclust(distance) 

# plots the clusters as a dendrogram
plot(clusters) 
  
