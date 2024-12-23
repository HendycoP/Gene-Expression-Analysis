---
title: "Gene Expression"
author: "Hendyco Pratama"
date: "2024-11-18"
output:
  html_document: default
  pdf_document: default
  md: default
  
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Final Project

Description of data set and all

```{r}
df <- read.table("GDS1028.soft", header = T, sep = "\t", row.names = 1, skip = 51, comment.char = "!")
``` 

```{r}
head(df)
```
The identifier column will not be used, hence it can be excluded.

```{r}
expr <- df[, 2:length(df)]
head(expr)
```

```{r}
colnames(expr)
```

```{r}
head(rownames(expr))
```

Based on the annotations found here  (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE1739), the dataset contains 14 samples in which there are 4 control and 10 SARS samples. 

Renaming the column names might make it easier to analyze the data.

```{r}
colnames(expr) <- c(paste0("control_", 1:4), paste0("sars_", 1:10))
colnames(expr)
```

Check whether there are missing data in the data set or not.

```{r}
which(is.na(expr))
which(expr == 0)
```

The data does not show any NAs nor Zeroes.

transofmring the data since it is not transformed yet
```{r}
expr <- log2(expr)
head(expr)
```
```{r}
sum(expr == 0)
which(expr == 0)
```

## Normalization

Normalization will be done in couple of ways to be compared.

```{r}
library(limma)
```

```{r}
expr.qNorm <- normalizeBetweenArrays(expr, method = "quantile") #perform quantile normalization
expr.lNorm <- normalizeBetweenArrays(expr, method = "cyclicloess") #perfrom cyclic loess normalization
```

visualizing the result of the normalization using density and scatter plots

```{r, fig.height=28}
par(mfrow=c(7,2))
for (i in 1:14) {
  max.expr <- max(expr[, i], expr.qNorm[, i], expr.lNorm[, i], na.rm =TRUE)
  min.expr <- min(expr[, i], expr.qNorm[, i], expr.lNorm[, i], na.rm =TRUE)
  x.range <- c(min.expr, max.expr)

  dens.expr <- density(expr[, i], na.rm = TRUE)
  dens.qNorm <- density(expr.qNorm[, i], na.rm = TRUE)
  dens.lNorm <- density(expr.lNorm[, i], na.rm = TRUE)
  max.y <- max(dens.expr$y, dens.qNorm$y, dens.lNorm$y)

  plot(dens.expr, xlim = x.range, ylim = c(0, max.y), col = "blue", lwd = 2, main = paste("Density Plot of Sample", i), xlab = "Expression", ylab = "Density")
  lines(dens.qNorm, lwd = 2, col = "red")
  lines(dens.lNorm, lwd = 2, col = "green")

  abline(v = mean(expr[, i], na.rm = TRUE), col = "black", lty = 2, lwd = 2)

  legend("topright", legend = c("pre-normalized", "quantile", "cyclic loess", "median global"), col = c("blue", "red", "green", "purple"), lwd = 2)
}
```

```{r, fig.height = 7}
#Generate boxplot for all normalization methods
par(mfrow = c(2,2))
boxplot(expr, main = "Boxplot of pre-normalized data", outline = FALSE, las = 2)
boxplot(expr.qNorm, main = "Boxplot of normalized data (quantile)", outline = FALSE, las = 2)
boxplot(expr.lNorm, main = "Boxplot of normalized data (cyclic loess)", outline = FALSE, las = 2)
```

Based on the plots (density and boxplot), quantile normalization seems like the best normalization method. Hence, the data that has been subjected to quantile normalization will be used going forward.

```{r}
expr <- expr.qNorm
summary(expr)
```
## Filtering Noise


```{r}
# filtering the genes based on the 1st quantile value.
q.threshold <- 5.0 # 5.0 was chosen because the 1st quantile for all samples are around 5.9
n.genes <- floor(0.25 * 14) # this will be used to retain genes that are expressed in at least 25% of all the sample
expr.fil <- expr[rowSums(expr > q.threshold) >= n.genes, ]
dim(expr.fil)
```

After filtering the genes, the number of genes has dropped although not significantly - from 8793 to 8286

```{r}
expr <- expr.fil
```

## Identifying Outliers

Identifying the outliers can be done by visualizing the overall data in multiple ways. For example by using correlation plot, dendrogram, CV vs mean plot, or average correlation plot. All in all, these plots will show which samples behave differently compared to each others.

# Pearson's Correlation Plot

```{r}
library(gplots)
```

```{r, fig.height = 10}
corr <- cor(expr, use = "pairwise.complete.obs", method = "pearson")
layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 2, 2), 5, 2, byrow = T))
par(oma = c(5, 7, 1, 1))

cx <- rev(colorpanel(23, "yellow", "black", "blue"))
lg <- seq(min(corr, na.rm = T), max(corr, na.rm = T), length = 10)

image(corr, main = "Pearson's Correlation Plot", axes = F, col = cx)
axis(1, at = seq(0, 1, length = ncol(corr)), label = row.names(corr),las = 2)
axis(2, at = seq(0, 1, length = ncol(corr)), label = row.names(corr), las = 2)

image(as.matrix(lg), col = cx, las = 2, axes = F)
tmp <- round(lg, 2)
axis(1, at = seq(0, 1, length = length(lg)), labels = tmp)
```

Based on the Pearson's Correlation heatmap, the data does seem to have one outlier - sample sars_3. But to be more sure we can run another test to see.

# Clustering Dendrogram

```{r, fig.height = 5}
t.expr <- t(expr)
dist <- dist(t.expr, method = "euclidian")
clust <- hclust(dist)
plot(clust, labels = names(t), ylab = "Distance", xlab = "Sample Name", main = "Cluster Dendrogram", sub = "", cex = 0.6)
```

# CV vs Mean Plot

```{r, fig.height = 5}
expr.mean <- apply(expr, 2, mean, na.rm = TRUE)
expr.sd <- sqrt(apply(expr, 2, var, na.rm = TRUE))
expr.cv <- expr.sd/expr.mean

plot(expr.mean, expr.cv, main = "CV vs Mean Plot", xlab = "Mean", ylab = "CV")
lines(lowess(expr.mean, expr.cv, f = 0.2), col = "red")
text(expr.mean, expr.cv, labels = colnames(expr), pos = 4, offset = 1, cex = 0.75)

```

# Average correlation plot

```{r, fig.height=5}
cor.mean <- apply(corr, 1, mean)
plot(c(1, length(cor.mean)), range(cor.mean), type = "n", ylab = "Mean", xlab = "", main = "Average correlation plot", xaxt = "n")
axis(1, at = c(1:length(cor.mean)), labels = names(cor.mean), las = 2)
points(cor.mean)
grid()
```

# PCA plot

The PCA plot will use only the 1st and 2nd PC because being the first two principal compnents, they will cover a more general information.

```{r, fig.height=5}
pca <- prcomp(expr, cor = F)
plot(pca$x[, 1], pca$x[, 2], main = "PCA plot for PC1 vs PC2", xlab = "PC1", ylab = "PC2")
```

From all of the plots generated above, all of them seems to come to the same conclusion that sars_3 sample is an outlier hence it will be removed from the data.

```{r}
expr <- subset(expr, select = -sars_3)
colnames(expr)
```
## Differential testing

To determine which testing will be used, the distribution and variance of the data needs to be evaluated. 

```{r}
# evaluating the distribution of the data
hist(expr, main = "Histogram of Gene Expression", xlab = "Normalized Expression Level")
```

```{r}
# evaluating the variance between groups using boxplot
f <- c(rep("Control", 4), rep("SARS", 9)) # there are only 9 sars sample because sample sars_3 has been removed
plot.data <- data.frame(Expression = colMeans(expr), Group = f) # representing the data per sample using its respective mean
boxplot(Expression ~ f, data = plot.data, main = "Boxplot between Groups", xlab = "Sample Group", ylab = "Expression", col = c("lightblue", "lightpink"))
```

```{r}
#statistically testing the variance between groups could be done with F-test
f.test <- var.test(expr[, 1:4], expr[, 5:13])
f.test
```

The histogram of the data shows that it follow normal distribution pattern hence it can be considered parametric. Meanwhile, the boxplot shows that the SARS group has more variance compared to the control group - shown by the size of the boxes. Despite that, the F test result (F = 0.983, p-value = 0.0753) shows that there is no significant difference between the variance, hence student t-test could be used.

```{r}
# perform student t test
t.test.all.genes <- function(x, s1, s2) {
  x1 <- x[s1]
  x2 <- x[s2]
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  t.out <- t.test(x1, x2, alternative = "two.sided", var.equal = T)
  out <- as.numeric(t.out$p.value)
  return(out)
}
```

```{r}
control <- colnames(expr)[1:4]
sars <- colnames(expr)[5:13]
ttest.result <- apply(expr, 1, t.test.all.genes, s1 = control, s2 = sars)
hist(ttest.result, col = "lightblue", main = "p-value distribution (Control vs SARS)", xlab = "p-value")
abline(v = .05, col = 2, lwd = 2)
```

```{r}
n.genes.p05 <- sum(ttest.result < 0.05)
print(c("# significant genes (p-val <0.05):", n.genes.p05))
```

## Multiple Testing

Based on the Student's t test, there are 1923 significant genes - genes with p-value of less than 0.05. Hence, we are going to do multiple testing to accommodate for the increased likelihood of false positive. 

```{r}
ttest.BH <- p.adjust(ttest.result, method = "BH")
ttest.BY <- p.adjust(ttest.result, method = "BY")
ttest.bonfer <- p.adjust(ttest.result, method = "bonferroni")
```
```{r, fig.height=5}
sorted.nonAdj <- sort(ttest.result)
sorted.BH <- sort(ttest.BH)
sorted.BY <- sort(ttest.BY)
sorted.bonfer <- sort(ttest.bonfer)

plot(sorted.nonAdj, type = "l", lwd = 2, col = "blue", main = "Adjusted and Non-Adjusted p-Value for Significant Genes", ylim = c(0, max(sorted.BH, sorted.BY, sorted.bonfer, sorted.nonAdj)), xlab = "Gene Index", ylab = "Sorted p-value")
lines(sorted.BH, col = "red", lwd = 2)
lines(sorted.BY, col = "green", lwd = 2)
lines(sorted.bonfer, col = "purple", lwd = 2)
legend("bottomright", legend = c("non-adjusted", "BH", "BY", "Bonferroni"), col = c("blue", "red", "green", "purple"), lty = 1, lwd = 2)
```

Based on the graph above, BY and bonferroni method shows more conservative adjustment than BH. Hence, BH will be used since the goal is for discovery analysis meaning that more potentially significant genes will be captured.

```{r}
# calculate fold change
control.mean <- apply(expr[,1:4], 1, mean, na.rm = TRUE)
sars.mean <- apply(expr[, 5:13], 1, mean, na.rm = TRUE)
FC <- control.mean - sars.mean
head(FC)
```

```{r}
sig.genes <- which(ttest.BH < 0.05 & abs(FC) > 1)
length(sig.genes)
```

After performing differential testing and using BH adjustment and utilizing fold change, the number of differentially expressed genes drop from 1923 to 256. 

```{r}
p.trans <- -1*log10(ttest.BH)
plot(range(p.trans), range(FC), type = "n", xlab = "-1*log10(p-value)", ylab = "Fold Change", main = "Volcano Plot")
points(p.trans, FC, col = "black", pch = 21, bg = 1)
points(p.trans[(p.trans > -log10(.05) & FC > log2(2))], FC[(p.trans > -log10(.05) & FC > log2(2))], col = 1, bg = 2, pch = 21)
points(p.trans[(p.trans > -log10(.05) & FC < -log2(2))], FC[(p.trans > -log10(.05) & FC < -log2(2))], col = 1, bg = 3, pch = 21)
abline(h = log2(2), col = 2)
abline(h = -log2(2), col = 2)
abline(v = -log10(.05), col = 2)
```

```{r}
layout(matrix(c(1, 2), 1, 2, byrow = TRUE))
hist(FC[sig.genes], col = "lightblue", main = "Histogram for Fold Change \nof Retained Genes", xlab = "Log 2 Fold Change", ylab = "Frequency")
hist(ttest.BH[sig.genes], col = "lightblue", main = "Histogram for p-value \nof Retained Genes", xlab = "p-value", ylab = "Frequency")
```

## Clustering

```{r}
sig.expr <- expr[sig.genes, ] #subsetting based on the retained genes
sample.group <- c(rep("Control", 4), rep("SARS", 9))
```

```{r, fig.height=5}
par(mar = c(5, 4, 4, 8) + 0.1)
sig.pca <- prcomp(t(sig.expr))
plot(range(sig.pca$x[, 1]), range(sig.pca$x[, 2]), type = "n", xlab = "PC1", ylab = "PC2", main = "PCA plot by group \nPC1 vs PC2")
points(sig.pca$x[, 1][sample.group == "Control"], sig.pca$x[, 2][sample.group == "Control"], col = "lightblue", pch = 16)
points(sig.pca$x[, 1][sample.group == "SARS"], sig.pca$x[, 2][sample.group == "SARS"], col = "lightpink", pch = 16)
legend("topright", inset = c(-0.25, 0), legend = unique(sample.group), col = c("lightblue", "lightpink"), pch = 16, title = "Sample Group", xpd = TRUE)
text(sig.pca$x[, 1], sig.pca$x[, 2], labels = colnames(sig.expr), cex = 0.45, pos = 1, xpd = TRUE)

```

```{r}
library(pheatmap)
```

```{r, fig.height=5}
# clustering using HCA
sig.dist <- dist(t(sig.expr), method = "euclidean")
sig.clust <- hclust(sig.dist, method = "complete")
plot(sig.clust, xlab = "Samples", sub = "", main = "Hierarchical Clustering Dendrogram \n(Euclidean and Complete Linkage)", cex = 0.75)
```

```{r, fig.height=7}
sig.genes.dist <- dist(sig.expr, method = "euclidean")
sig.genes.clust <- hclust(sig.genes.dist, method = "complete")
pheatmap(sig.expr, cluster_cols = sig.clust, cluster_rows = sig.genes.clust, show_colnames = TRUE, show_rownames = TRUE, main = "Gene Expression Heatmap", fontsize_col = 10, fontsize_row = 2)
```

The scatter plot for the PCA between PC1 and PC2 successfully grouped samples based on their condition, with a clear distinction between the SARS and control groups. However, the hierarchical clustering dendrogram showed that the sample SARS_9 was clustered in the same branch as the control samples. The heatmap of gene expression further supports this finding, showing that certain genes in the SARS_9 sample have expression patterns similar to those in the control group. This observation suggests that SARS_9 may exhibit a different biological variation compared to the other SARS samples. Additionally, it is important to consider that the differential expression analysis used the BH (Benjamini-Hochberg) method, which belongs to the FDR (False Discovery Rate) family. As such, there may be some false positives included in the results, which could contribute to the unexpected grouping of SARS_9 with the control samples.

## classification
```{r}
library(MASS)
```

```{r}
classes <- c(rep("control", 4), rep("SARS", 9))
classification.data <- data.frame(classes, t(sig.expr))
train.set <- rbind(classification.data[2:4, ], classification.data[5:9, ])
test.set <- rbind(classification.data[1, ], classification.data[10:13, ])

class.label <- test.set[, 1]
test.set <- test.set[, -1]

model <- lda(train.set[, 1] ~., train.set[, 2:257])
out <- predict(model, test.set)
table(out$class, class.label)
```

```{r}
colors <- c("control" = "lightblue", "SARS" = "lightpink")
plot(out$x[, 1], rep(0, length(out$x[, 1])), col = colors[class.label], xlab = "Discriminant function", ylab = "", main = "Discriminant Function Plot", pch = 19)
legend("bottomright", legend = names(colors), col = colors, pch = 19)
```

## 10 Discriminant Genes

```{r}
lda.loadings <- model$scaling

coef <- lda.loadings[, 1]

sorted.coef <- sort(coef, decreasing = TRUE)
pos5 <- head(sorted.coef, 5)
neg5 <- tail(sorted.coef, 5)

pos5.probe <- gsub("^X", "", names(pos5))
neg5.probe <- gsub("^X", "", names(neg5))
```

```{r}
# Find the gene name from the the identifier of the original data frame
pos5.gene.name <- df[pos5.probe, "IDENTIFIER"]
pos5.gene.name
```

```{r}
# Find the gene name from the the identifier of the original data frame
neg5.gene.name <- df[neg5.probe, "IDENTIFIER"]
neg5.gene.name
```
