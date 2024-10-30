#Load necessary packages
library(Seurat)
library(SeuratObject)
#Load .RDS file into R and update object (only if prompted to)
df <- readRDS('MDAMB231.rds')
df <- UpdateSeuratObject(df)
#Inspect the Seurat Object
df
#Obtain percentage of mitochondrial genes within each sample and save to meta.data
# for quality control metrics
df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^MT-")
#Generate violin plots of the nFeature_RNA(unique genes per cell), nCount_RNA(total counts) and percent.mt
VlnPlot(df, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#Subset the data if necessary (in this case, remove all cells in which percent.mt is 10>)
df_subset <- subset(df, subset = percent.mt < 10)
#Inspect the subset 
df_subset
#Normalize the subsetted data
df_subset<-NormalizeData(df_subset,normalization.method = "LogNormalize",scale.factor =  10000)
#Load a csv which contains glycolytic genes of interest
genelist<-read.csv('glycolyticgenes.csv')
#Extract gene symbols from the genelist
genes_of_interest<-genelist[1]
#Extract the average of each gene expression into a new dataframe, any genes not found in the 
#dataset will be output in the console
mean_expression <- AverageExpression(df_subset, features = c(genes_of_interest), return.seurat = FALSE)
#Write the data to a csv
write.csv(mean_expression,'MDA-MB-231_Glycogenes.csv')