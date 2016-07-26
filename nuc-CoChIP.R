WorkDir = "~/Dropbox/CoChIP/coChIP-Analysis"
DataDir = getwd()
setwd(WorkDir)


DataDir =  "~/Data/CoChIP/HisMut-K4/HisMut"  
OrigDir = DataDir
BamDir = paste0(DataDir,"/BAM")

extendDir <- function(x) {
  a = substr(x,1,1)
  if( a == "/" || a == "~")
    return(x)
  else
    return(paste0(OrigDir,"/",x))
}

print("Initializing")
suppressMessages(source("CoChip-Functions.R"))
suppressMessages(source("NucAtlas.R"))
suppressMessages(source("DensityScatter.R"))
source("Main-Functions.R")
InputName = "-input"

GetSGDandNucs()

params <- ccParams()
params$TSSRegions = SGD$TSSRegion
params$TTSRegions = SGD$TTSRegion
params$GeneRegions = SGD$GeneRegion
params$NucRegions = NucRegions
params$MaxFragLen = 220
params$MinFragLen = 50
params$cen = FALSE
params$meta = FALSE
params$nuc = TRUE
params$DataDir = DataDir

NucFile = "nucdata"
NucFN = paste0(DataDir,"/",NucFile,".rdata")

MinReads = 200000

## get the data
nuc.data =readRDS(NucFN)
nuc.counts = nuc.data$counts
nuc.mat = nuc.data$mat

names(nuc.counts) = rownames(nuc.mat)

if( 0 ) {
  # get rid of outliers
  tot = colMeans(nuc.mat)
  t = quantile(tot, .99)
  m = mean(tot[tot <= t])
  s = sqrt(var(tot[tot <= t]))
  t = m + 3*s
  Io = tot < t
  nuc.mat = nuc.mat[,Io]
  cNucRegions = NucRegions[Io]
}

ReadCounts = rowSums(nuc.mat)

if(0) {
  #normalize
  TotalNuc = dim(nuc.mat)[[2]]
  for(i in rownames(nuc.mat)) {
    s = sum(nuc.mat[i,])/TotalNuc
    nuc.mat[i,] = nuc.mat[i,]/s
  }
}

I = nuc.counts > MinReads
nuc.smat = nuc.mat[I,]
print(dim(nuc.smat))

cnucs = nuc.mat
nucs = nuc.smat

ExpNames = rownames(nuc.mat)
SplitExpName = matrix(unlist(lapply(ExpNames, function(x) strsplit(x,"-"))), nc = 3, byrow=T)
rownames(SplitExpName) = ExpNames
colnames(SplitExpName) = c("Strain", "Ab1", "Ab2")
Ab1 = factor(SplitExpName[,"Ab1"])
Ab2 = factor(SplitExpName[,"Ab2"])
Strain = factor(SplitExpName[,"Strain"])

aH3K4me3 = as.integer(Ab1 == "H3K4me3") + as.integer(Ab2 == "H3K4me3")
aH3K18ac = as.integer(Ab1 == "H3K18ac") + as.integer(Ab2 == "H3K18ac")
aH3K36me3 = as.integer(Ab1 == "H3K36me3") + as.integer(Ab2 == "H3K36me3")
aH3K79me3 = as.integer(Ab1 == "H3K79me3") + as.integer(Ab2 == "H3K79me3")
aH3 = as.integer(Ab1 == "H3") + as.integer(Ab2 == "H3")


if( 0 ) {
#PredictCounts = glm( ReadCounts ~ Strain + Ab1 + Ab2 -1, family = poisson(link=log ) )

mm = model.matrix(ReadCounts ~ Strain + aH3+aH3K4me3 + aH3K18ac + aH3K36me3 + aH3K79me3 + 
                    aH3K4me3:aH3K18ac + aH3K4me3:aH3K36me3 + aH3K4me3:aH3K79me3 +
                    aH3K18ac:aH3K36me3 + aH3K18ac:aH3K79me3 + aH3K36me3:aH3K79me3 +                  
                    Strain:(aH3+aH3K4me3 + aH3K18ac + aH3K36me3 + aH3K79me3) )
#mm = model.matrix(ReadCounts ~ Strain + Ab1 + Ab2 + Strain:(Ab1+Ab2) + Ab1:Ab2)

fit = glmnet(mm,ReadCounts,family = "poisson")
cvfit = cv.glmnet(mm,ReadCounts,family = "poisson")
plot(cvfit)
c = as.matrix(coef(cvfit, cvfit$lambda.min))
c[c != 0]
PredictedReads = exp(predict(cvfit, newx = mm, s = "lambda.min"))

for( s in levels(Strain) )
  for( a in levels(Ab1))
    for( b in levels(Ab1) ) 
      if( !(a == b) && is.element(BuildExpName(a,b,s), ExpNames) )
#        RelativeEnrichment(a,b, pre=s)
        PlotPair(a,b,pre=s,relative=FALSE)
}    