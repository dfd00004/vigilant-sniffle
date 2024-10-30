#Load necessary packages
library(dplyr)
library(ggplot2)
library(ggfortify)
#Load data into R ensuring rows represent samples and columns are variables (Genes,Names,Genotype etc)
df<-read_csv("KO_vs_WT.csv")
#Subset sample information and genotype from proteomics data
protvalues<-df[,4:ncol(df)]
samples<-df['SampleID']
genotype<-df['Genotype']
#Remove NaN values from proteomics data and convert to 0 (since 0 represents no detection)
protvalues<-protvalues %>% replace(is.na(.), 0)
#Perform principal component detection on protvalues while also scaling and centering the data
pca_results<-prcomp(protvalues,scale. = TRUE)
#Generate plot of the pca_results in R using PC1 and PC2 and color by genotype of samples
prot.pca.plot<-autoplot(pca_results,data=df,colour='Genotype')
#Obtain the pca_scores for plotting
pca_scores<-pca_results$x
#Normalize the pca_scores for plotting
process_pca<-preProcess(pca_scores,method = c("range"))
pca_scores_normalized<-predict(process_pca,pca_scores)
#Combine the normalized pca scores with sample and genotype data
pca_scores_normalized_annotated<-cbind(samples,genotype,pca_scores_normalized)
#Obtain variance explained by each PC
variance_explained <- pca_results$sdev^2 / sum(pca_results$sdev^2)
# Create a scree plot in R
plot(variance_explained, type = "b", 
     xlab = "Principal Component", 
     ylab = "Proportion of Variance Explained",
     main = "Scree Plot")
#Save necessary data for later analysis
write.csv(pca_scores_normalized_annotated,'PCAscores.csv')
write.csv(variance_explained,'VarianceExplained.csv')

