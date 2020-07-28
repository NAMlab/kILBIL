

# The script shows the metabolomic data processing for "red ripe" samples in a negative ionization mode. 
# Repeat the processing for the positive ionization mode accordingly 

# 0. Install packages ----
library(xcms)
library(CAMERA)
library(BiocParallel)
bpparam <- MulticoreParam(4)



# 1. Perform peak picking in batches ----
ionization = "negative" #
PL_PATH = "."
cdfFiles = list.files("netCDF_red_stage",full.names=T,pattern="*0[1|2].CDF")
strata = findInterval(1:length(cdfFiles), seq(1,length(cdfFiles), by = 100))
stop()
for (i in 1:max(strata)) { 
  cdfFiles.i = cdfFiles[strata==i]	
  xs.i = xcmsSet(cdfFiles.i,snthresh=5,fwhm=10,mzdiff=0.05,max=200,BPPARAM=bpparam)
  save(xs.i,file=file.path(PL_PATH,paste("xs",i,sep="")))
}

load(file.path(PL_PATH,paste("xs",1,sep="")))
xs = xs.i
rm(xs.i)
for(i in unique(strata)){
  if(i != 1){
    print(i)
    load(file.path(PL_PATH,paste("xs",i,sep="")))
    xs = c(xs, xs.i)
    rm(xs.i)
  }
}


# 2. Group peaks ----
xsg = group(xs, bw=25, mzwid = 0.05, minsamp = 4, max = 200)

# 3. Perform retention time correction ----
xs_corr = retcor(xsg,smooth="loess",plottype="mdevden",span=0.75,missing=120,extra=20)

# 4. Group corrected peaks ----
xs_corrg = group(xs_corr, bw=20, mzwid = 0.05, minsamp= 4, max = 200)

# 5. Fill peaks ----
p <- SerialParam()
xs_corrg_fill = fillPeaks(xs_corrg,BPPARAM=p)

# 6. Annotate peaks with CAMERA ----
xan = annotate(xs_corrg_fill, polarity=ionization,ppm=20,multiplier=2,quick=TRUE)
PL = getPeaklist(xan)

time.RetGroup <- system.time({ # measuring time
  resultRetcorGroup <-
    optimizeRetGroup(xset = xs, 
                     params = retcorGroupParameters, 
                     nSlaves = 8, 
                     subdir = "retCorGrpOpt")
})



