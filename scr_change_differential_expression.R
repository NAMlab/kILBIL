

# 1. Load the input data ----


load(file = "/data/master_table") # table of marker-based introgression mapping integrated from Ofner et al. 2016
load(file = "/data/counts_lyc") # count table for S. lycopersicum-unique reads
load(file = "/data/counts_ind") # count table for indistinguishable reads


# *7.2 Ind data ----

# Generate the coldata matrix

#Select Lyc and BILs as groups
group =  master.table[,1]
group[group == 1] = "Lyc"
group[group == 0] = "BIL"
group = factor(group)

#Generate the data object
y <- DGEList(counts.ind, group = group, genes=rownames(counts.ind.filtered))

#Calculate the norm factor
y <- calcNormFactors(y)
head(y$samples)


#Generate the design matrix

#Calculate dispersions
y.ind <- estimateDisp(y, design.matrix.lyc)
plotBCV(y.ind)

y.cpm.ind = cpm(y.ind, normalized.lib.sizes=TRUE, log=F, prior.count=0.25)
y.cpm.ind.log = cpm(y.ind, normalized.lib.sizes=TRUE, log=T, prior.count=0.25)
write.table(y.cpm.ind, "normalized_counts_ind.txt", sep="\t")
write.table(y.cpm.ind.log, "normalized_counts_ind_log.txt", sep="\t")


# 8. DIfferential expression ----


linelist = sapply(colnames(design.matrix), function(x) rownames(design.matrix)[as.vector(design.matrix[,x]) == 1])


# *8.1 Ind data ----
# pb = txtProgressBar(min = 1, max = ncol(design.matrix), style = 3)
# results.list.ind = list()
# for(i in 2:ncol(design.matrix)){
#   setTxtProgressBar(pb, i)
#   tmp = design.matrix[,i]
#   design.matrix.i = design.matrix[,c(1,i)]
#   design.matrix.i[,1][design.matrix.i[,2] == 0] = 1
#   fit <- glmQLFit(y.ind, design.matrix.i)
#   lrt <- glmLRT(fit, contrast=c(-1,1))
#   result.lrt = topTags(lrt, sort.by="none", n = nrow(y.ind))$table
#   results.list.ind[[i]] = result.lrt
# }
# close(pb)



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


t.dataset.ind = c()
for(i in 1:length(results.list.ind)){
  t.dataset.ind  = cbind(t.dataset.ind,results.list.ind[[i]][,6])
}
dimnames(t.dataset.ind) = list(rownames(y.na$counts), colnames(design.matrix)[2:ncol(design.matrix)] )


t.dataset.ind.f= c()
for(i in 1:length(results.list.ind)){
  t.dataset.ind.f  = cbind(t.dataset.ind.f,results.list.ind[[i]][,2])
}
close(pb)
dimnames(t.dataset.ind.f) = list(rownames(y.na$counts), colnames(design.matrix)[2:ncol(design.matrix)] )


save(t.dataset.ind, file = "t.dataset.ind")
save(t.dataset.ind.f, file = "t.dataset.ind.f")
save(y.cpm.ind.log, file = "y.cpm.ind.log")
save(y.ind, file = "y.ind")




# *8.2 Lyc data ----
# pb = txtProgressBar(min = 1, max = ncol(design.matrix), style = 3)
# results.list.lyc = list()
# for(i in 2:ncol(design.matrix)){
#   setTxtProgressBar(pb, i)
#   tmp = design.matrix[,i]
#   design.matrix.i = design.matrix[,c(1,i)]
#   design.matrix.i[,1][design.matrix.i[,2] == 0] = 1
#   fit <- glmQLFit(y.lyc, design.matrix.i)
#   lrt <- glmLRT(fit, contrast=c(-1,1))
#   result.lrt = topTags(lrt, sort.by="none", n = nrow(y.lyc))$table
#   results.list.lyc[[i]] = result.lrt
# }
# close(pb)

# parallel version
library(doParallel)
registerDoParallel(cores=3)
results.list.lyc = foreach(i=2:ncol(design.matrix)) %dopar% {
  design.matrix.i = design.matrix[,c(1,i)]
  design.matrix.i[,1][design.matrix.i[,2] == 0] = 1
  fit <- glmQLFit(y.lyc, design.matrix.i)
  lrt <- glmLRT(fit, contrast=c(-1,1))
  result.lrt = topTags(lrt, sort.by="none", n = nrow(y.lyc))$table
  return(result.lrt)
}



t.dataset.lyc = c()
for(i in 1:length(results.list.lyc)){
  t.dataset.lyc  = cbind(t.dataset.lyc,results.list.lyc[[i]][,6])
}
dimnames(t.dataset.lyc) = list(rownames(y.na$counts), colnames(design.matrix)[2:ncol(design.matrix)] )


t.dataset.lyc.f= c()
for(i in 1:length(results.list.lyc)){
  t.dataset.lyc.f  = cbind(t.dataset.lyc.f,results.list.lyc[[i]][,2])
}
close(pb)
dimnames(t.dataset.lyc.f) = list(rownames(y.na$counts), colnames(design.matrix)[2:ncol(design.matrix)] )


save(t.dataset.lyc, file = "t.dataset.lyc")
save(t.dataset.lyc.f, file = "t.dataset.lyc.f")
save(y.cpm.lyc.log, file = "y.cpm.lyc.log")
save(y.lyc, file = "y.lyc")



# *8.1 Trans data ----

y.na = y.lyc
for(i in 1:length(linelist)){
  tmp = itag.class$Locus[itag.class$bin == names(linelist)[i] ]
  tmp = intersect(tmp,rownames(y.na$counts))
  if(length(tmp) > 0){
    y.na$counts[tmp, linelist[[i]]] = NA
  }
}

pc <- pca(y.na$counts, nPcs=3, method="svd")
imputed <- completeObs(pc)
imputed[imputed < 0] = 0
y.na$counts = imputed
y.cpm.trans.log = cpm(y.na, normalized.lib.sizes=TRUE, log=T, prior.count=0.25)

pb = txtProgressBar(min = 1, max = ncol(design.matrix), style = 3)
results.list.trans = list()
for(i in 2:ncol(design.matrix)){
  setTxtProgressBar(pb, i)
  design.matrix.i = design.matrix[,c(1,i)]
  design.matrix.i[,1][design.matrix.i[,2] == 0] = 1
  fit <- glmQLFit(y.na, design.matrix.i)
  lrt <- glmLRT(fit, contrast=c(-1,1))
  result.lrt = topTags(lrt, sort.by="none", n = nrow(y.na))$table
  results.list.trans[[i]] = result.lrt
}
close(pb)


t.dataset.trans = c()
for(i in 2:length(results.list.trans)){
  t.dataset.trans  = cbind(t.dataset.trans,results.list.trans[[i]][,6])
}
dimnames(t.dataset.trans) = list(rownames(y.na$counts), colnames(design.matrix)[2:ncol(design.matrix)] )


t.dataset.trans.f= c()
for(i in 2:length(results.list.trans)){
  t.dataset.trans.f  = cbind(t.dataset.trans.f,results.list.trans[[i]][,2])
}
close(pb)
dimnames(t.dataset.trans.f) = list(rownames(y.na$counts), colnames(design.matrix)[2:ncol(design.matrix)] )


save(t.dataset.trans, file = "t.dataset.trans")
save(t.dataset.trans.f, file = "t.dataset.trans.f")
save(y.cpm.trans.log, file = "y.cpm.trans.log")
save(y.na, file = "y.na")

# *8.2 Mixed data ----

# Turn the normalizet cpm values from ind reads to the lyc counts
tmp = sweep(y.cpm.ind, 2, colSums(y.cpm.ind), FUN = "/")
tmp = sweep(tmp, 2, colSums(y.lyc$counts), FUN = "*")

tmp = sweep(tmp, 1, rowMeans(tmp), FUN = "/")
tmp = sweep(tmp, 1, rowMeans(y.lyc$counts), FUN = "*")
tmp.mixed = tmp
tmp.mixed[is.na(tmp.mixed)] = 0
tmp.mixed = round(tmp.mixed,0)
# tmp = sweep(y.cpm.ind, 2, colSums(y.cpm.ind), FUN = "/")
# tmp.mixed = sweep(tmp, 2, colSums(y.lyc$counts), FUN = "*")
# tmp.mixed = round(tmp.mixed,0)


y.lyc$counts[1:10,1:10]
tmp.mixed[1:10,1:10]

# Replace the values
y.mix = y.lyc
for(i in 1:length(linelist)){
  tmp = itag.class$Locus[itag.class$bin == names(linelist)[i] ]
  tmp = intersect(tmp,rownames(y.mix$counts))
  if(length(tmp) > 0){
    y.mix$counts[tmp, linelist[[i]]] = tmp.mixed[tmp, linelist[[i]]]
  }
}

#Generate the data object

group[group == 1] = "Lyc"
group[group == 0] = "BIL"
group = factor(group)

y <- DGEList(y.mix$counts, group = group, genes=rownames(y.mix$counts))
y <- calcNormFactors(y)
y.mix <- estimateDisp(y, design.matrix.lyc)

y.cpm.mix = cpm(y.mix, normalized.lib.sizes=TRUE, log=F, prior.count=0.25)
y.cpm.mix.log = cpm(y.mix, normalized.lib.sizes=TRUE, log=T, prior.count=0.25)


y <- DGEList(y.cpm.mix, group = group, genes=rownames(y.mix$counts))
y <- calcNormFactors(y)
y.mixx <- estimateDisp(y)

y.cpm.mixx = cpm(y.mixx, normalized.lib.sizes=TRUE, log=F, prior.count=0.25)
y.cpm.mixx.log = cpm(y.mixx, normalized.lib.sizes=TRUE, log=T, prior.count=0.25)

# parallel version
library(doParallel)
registerDoParallel(cores=3)
results.list.mixx = foreach(i=2:ncol(design.matrix)) %dopar% {
  design.matrix.i = design.matrix[,c(1,i)]
  design.matrix.i[,1][design.matrix.i[,2] == 0] = 1
  fit <- glmQLFit(y.mixx, design.matrix.i)
  lrt <- glmLRT(fit, contrast=c(-1,1))
  result.lrt = topTags(lrt, sort.by="none", n = nrow(y.mixx))$table
  return(result.lrt)
}

t.dataset.mixx = c()
for(i in 1:length(results.list.mixx)){
  t.dataset.mixx  = cbind(t.dataset.mixx,results.list.mixx[[i]][,6])
}
dimnames(t.dataset.mixx) = list(rownames(y.lyc$counts), colnames(design.matrix)[2:ncol(design.matrix)] )


t.dataset.mixx.f= c()
for(i in 1:length(results.list.mixx)){
  t.dataset.mixx.f  = cbind(t.dataset.mixx.f,results.list.mixx[[i]][,2])
}
close(pb)
dimnames(t.dataset.mixx.f) = list(rownames(y.lyc$counts), colnames(design.matrix)[2:ncol(design.matrix)] )
save(t.dataset.mixx, file = "t.dataset.mixx")
save(t.dataset.mixx.f, file = "t.dataset.mixx.f")
save(y.cpm.mixx.log, file = "y.cpm.mixx.log")
save(y.mixx, file = "y.mixx")

