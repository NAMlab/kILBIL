

# 1. Load the input data ----

load(file = "/data/itag30") # ITAG30 genome annotation file
load(file = "/data/counts_lyc") # count table for S. lycopersicum-unique reads
load(file = "/data/counts_ind") # count table for indistinguishable reads

# 2. Generate the placeholder list ----

changepoint.list = list()

# 3. Run the change point analysis on the count table for S. lycopersicum-unique reads

for(i in colnames(master.table)){ # Run for each line
  
  # Generate a vector of expression deviations from the mean across all quantified genes for the S. lycopersicum-unique reads
  period1 = log10(counts.lyc[,i]/itag.lyc[,"MeanCountLyc"])
  names(period1) = rownames(itag.lyc)
  period1 = period1[is.na(period1) == F & period1 != Inf & period1 != -Inf]
  
  # Generate a vector of expression deviations from the mean across all quantified genes for the indistinguishable reads
  period2 = log10(counts.ind[,i]/counts.ind[,"MeanCountInd"])
  names(period2) = rownames(itag.lyc)
  period2 = period2[names(period1)]
  
  # Generate a vector of Chromosome annotations for all quantified genes
  chromosome = itag.lyc[names(period1),"Chromosome"]
  transbin = rep(0,length(period1))
  names(transbin) = names(period1)  
  
    n = 0
    for(j in 1:12){ # Run for each chromosome
      
      x <- period1[chromosome == j]
      x = x[is.na(x) == F & x != Inf & x != -Inf]
      resultsStream <- processStream(x, cpmType = "Mann-Whitney", ARL0 = 1000, startup = 50, lambda=NA) # change point analysis with optimized parameters
      
      # Parse the change point analysis results
      if(length(resultsStream$changePoints) == 0){
        res = matrix(c(1,length(x)),ncol=2) 
        res3 = cbind(x, rep(1, length(x)))
      }else{
        res = cbind(c(1,resultsStream$changePoints), c(resultsStream$changePoints, length(x)))
        res[2:nrow(res),1] = res[2:nrow(res),1]+1
        res2 = sapply(1:nrow(res),function(y) rbind(x[res[y,1]:(res[y,2])], rep(y, length(x[res[y,1]:res[y,2]])) ))
        res3 = as.data.frame(t(as.data.frame(res2)))
      }
      transbin[rownames(res3)] = res3[,2] + n
      n = max(res3[,2] + n) 
    }
    
    transbin = transbin[names(period1)]
    transbin = transbin+1
    tmp1 = transbin[is.na(period1) == F & period1 != Inf & period1 != -Inf]
    tmp2 = period1[is.na(period1) == F & period1 != Inf & period1 != -Inf]
    
    
    # Compare the deviation from the mean of the S. lycopersicum-unique reads between the genes within the segment of interest and all remaining genes 
    # We use one-tailed test since all introgressions are expected to result in the deviation only towards negative values
    transbin.pvalue = sapply(unique(tmp1), function(y)  t.test(tmp2[tmp1 ==y],tmp2[tmp1 !=y], alternative = "less")$p.value)
    transbin.estimate = sapply(unique(tmp1), function(y)  t.test(tmp2[tmp1 ==y],tmp2[tmp1 !=y], alternative = "less")$estimate[1])
    
    changepoint.list[[i]] = data.frame(estimate = unlist(transbin.estimate), p.value = transbin.pvalue)
    
}


