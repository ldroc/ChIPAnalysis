WorkDir = "~/Dropbox/coChIP-scripts/"
DataDir = getwd()
setwd(WorkDir)


#DataDir =  "~/Data/CoChIP/HisMut-K4/HisMut/Data"  
DataDir =  "~/Data/CoChIP/160124/K4Sym/Data"  
#DataDir =  "~/Data/CoChIP/160124/RPD3/DataLong"
OrigDir = DataDir
BamDir = paste0(DataDir,"/BAM")
Trackdir = DataDir


print("Initializing")
suppressMessages(source("CoChip-Functions.R"))
suppressMessages(source("NucAtlas.R"))
suppressMessages(source("DensityScatter.R"))
GetSGDandNucs()

params = ccParams()
params$DataDir = DataDir
params$cov = TRUE
params$meta = TRUE
params$cen = TRUE
params$nuc = TRUE
params$TSSRegions = SGD$TSSRegion
params$TTSRegions = SGD$TTSRegion
params$GeneRegions = SGD$GeneRegion
params$PlusOneRegions = Plus1Nucs
params$NucRegions = NucRegions
params$MaxFragLen = 220

extendDir <- function(x) {
  a = substr(x,1,1)
  if( a == "/" || a == "~")
    return(x)
  else
    return(paste0(OrigDir,"/",x))
}

#cnts = read.csv(paste0(DataDir,"/HisMut-all-counts.csv"), row.names = 1)
#FileNames = rownames(cnts)
FileNames = list.files( paste0(DataDir,"/Old"), "*.rdata$", full.names = FALSE )
OrigNames = FileNames
FileNames = sub(".rdata", "", FileNames)
FileNames = sub("_.*", "", FileNames)
FileNames = sub("K4Sym-", "", FileNames)

Names = sapply(FileNames, function(n) strsplit(n,"-"), simplify = "array")

for( i in 1:length(Names) ) {
  m = Names[[i]]
  if( m[2] != m[3] && m[3] != "Input") {
    j = which(FileNames == paste0(m[1],"-",m[3],"-",m[2]))
    if( length(j)==1 && j > i ) {
      print(paste(FileNames[[i]], FileNames[[j]]))
#      dat.i = ccProcessFile(filename = FileNames[[i]], param = params)
#      dat.j = ccProcessFile(filename = FileNames[[j]], param = params)
      dat.i = readRDS(paste0(DataDir,"/Old/", OrigNames[[i]]))
      dat.j = readRDS(paste0(DataDir,"/Old/", OrigNames[[j]]))
      dat = ccCombineDat(dat.i,dat.j, paste0(m[1],"-",m[2],"-",m[3]))
      ccProcessFile(dat = dat, param = params)
    } else {
      dat = readRDS(paste0(DataDir,"/Old/", OrigNames[[i]]))
      if( m[[2]] < m[[3]]) {
        fname = paste0(m[1],"-",m[2],"-",m[3])
      } else
        fname = paste0(m[1],"-",m[3],"-",m[2])
      dat$Name = fname
      ccProcessFile(dat = dat, param = params)
    }
  } else {
    dat = readRDS(paste0(DataDir,"/Old/", OrigNames[[i]]))
    fname = paste0(m[1],"-",m[2],"-",m[3])
    dat$Name = fname
    ccProcessFile(dat = dat, param = params)
  }
}
