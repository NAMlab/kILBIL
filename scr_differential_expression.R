
# 0. Install packages ----
library(doParallel)
library(edgeR)


# 1. Load the input data ----


load(file = "/data/master_table") # table of marker-based introgression mapping integrated from Ofner et al. 2016
load(file = "/data/counts_lyc") # count table for S. lycopersicum-unique reads
load(file = "/data/counts_ind") # count table for indistinguishable reads
load(file = "/data/design_matrix") # table of marker-based introgression mapping integrated from Ofner et al. 2016

# *7.2 Ind data ----

# Generate the coldata matrix

#Select Lyc and BILs as groups
group =  design.matrix[,1]
group[group == 1] = "Lyc"
group[group == 0] = "BIL"
group = factor(group)

#Generate the data object
y <- DGEList(counts.ind, group = group, genes=rownames(counts.ind.filtered))

#Calculate the norm factor
y <- calcNormFactors(y)

#Calculate dispersions
y.ind <- estimateDisp(y, design.matrix)
plotBCV(y.ind)

# 8. DIfferential expression ----

# parallel version
library(doParallel)
registerDoParallel(cores=3)
results.list.ind = foreach(i=2:ncol(design.matrix)) %dopar% {
  design.matrix.i = design.matrix[,c(1,i)]
  design.matrix.i[,1][design.matrix.i[,2] == 0] = 1
  fit <- glmQLFit(y.ind, design.matrix.i)
  lrt <- glmLRT(fit, contrast=c(-1,1))
  result.lrt = topTags(lrt, sort.by="none", n = nrow(y.ind))$table
  return(result.lrt)
}


