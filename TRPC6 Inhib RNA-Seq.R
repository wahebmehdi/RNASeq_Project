#Setting up libraries

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(readxl)


#loading up the data
count_data = read_excel("C:\\Users\\mehdi\\Downloads\\GSE242592_COUNT.xlsx")


#Trimming the dataset
 
#Cleaning the Gene column
count_data$Gene = gsub("^ENSG", "", count_data$Gene)

#Remove the "Name" and "Type" columns
count_data = count_data[, !(names(count_data) %in% c("Name", "Type"))]

#Check for duplicates 
sum(duplicated(count_data))

#Checking for missing values 
which(is.na(count_data))



#Creation of data-frame for GeneNames
genenames = count_data$Gene


#Creation of count_data 
count_data = count_data[,2:7]
rownames(count_data) = genenames


#Changing colnames 
coldata = data.frame("condition" = as.factor(c(rep("Inhibitor",3), rep("Control", 3))), row.names = colnames(count_data))

#Checking if the rows of colData and the columns of count_data are in the same order
all(colnames(count_data) == rownames(colData))



#Creating a dds object for DESeqAnalysis
dds = DESeqDataSetFromMatrix(colData = coldata, countData = count_data, design= ~ condition)


#Pre-filtering 
#Genes with a total read count >= 10 are kept for Differential Expression Analysis 

# Check dimensions
n_rows_dds = nrow(dds)
n_rows_count_data = nrow(count_data)

if (n_rows_dds != n_rows_count_data) {
  cat("Dimensions of 'dds' and 'count_data' do not match. Adjusting...\n")
  
  if (n_rows_dds > n_rows_count_data) {
    # Trim 'dds' to match the number of rows in 'count_data'
    dds <- dds[1:n_rows_count_data, ]
  } else {
    # Trim 'count_data' to match the number of rows in 'dds'
    count_data <- count_data[1:n_rows_dds, ]
  }
}
# Now 'dds' and 'count_data' have the same number of rows
# Subsetting based on row sums
keep <- rowSums(count_data) >= 10
dds <- dds[keep, ]
dds

#Setting a factor level i.e which level(factor) to compare against- setting Control as our reference
dds$condition = relevel(dds$condition, ref = 'Control')


# Assuming 'dds' is your DESeqDataSet object
# Replace 10 and 0.5 with your chosen count and proportion thresholds
dds = dds[rowSums(counts(dds) >= 20) >= ncol(dds) * 0.5,]


#Performing DESeq Analysis on dds object 

dds = DESeq(dds)


#Save normalized read counts
normalizedcounts = counts(dds, normalized = TRUE)
write.csv(normalizedcounts, "normalized_counts.csv")


#Extracting results of DEAnalysis from dds object 
res = results(dds, alpha =0.01)
summary(res)

#Visualizations 

#PlotMA
plotMA(res)

#Dispersion Plot 
plotDispEsts(dds, main = "dispersion plot")


#PCA Plot
#Perform rlog transformation
rld = rlog(dds, blind = FALSE)

#PCA analysis
PCAA = plotPCA(rld, intgroup = "condition")

#Create PCA plot
PCAA + geom_text(aes(label = name), size = 2.5) + ggtitle('PCA Plot')


#HeatMap
sampleDists = dist(t(assay(rld)))

library("RColorBrewer")
sampleDistMatrix = as.matrix(sampleDists)

colnames(sampleDistMatrix)

colors = colorRampPalette(rev(brewer.pal(9, "Blues")) ) (255)

pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col = colors)
