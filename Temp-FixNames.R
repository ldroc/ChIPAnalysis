WorkDir = "~/Dropbox/CoChIP/coChIP-Analysis"
DataDir = getwd()
setwd(WorkDir)


#DataDir =  "~/Data/CoChIP/HisMut-K4/HisMut"  
DataDir =  "~/Data/CoChIP/160124/RPD3"  
OrigDir = DataDir
BamDir = paste0(DataDir,"/BAM")
Trackdir = DataDir


print("Initializing")
suppressMessages(source("CoChip-Functions.R"))
suppressMessages(source("NucAtlas.R"))
suppressMessages(source("DensityScatter.R"))
GetSGDandNucs()

params = ccParams()
params$DataDir = paste0(DataDir,"/Data")
params$cov = FALSE
params$meta = FALSE
params$cen = FALSE
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
FileNames = list.files( paste0(DataDir,"/Data"), "^RPD3.*.rdata$", full.names = FALSE )
OrigNames = FileNames
FileNames = sub(".rdata", "", FileNames)
FileNames = sub("_.*", "", FileNames)
FileNames = sub("RPD3-", "", FileNames)



for( i in 1:length(FileNames) ) {
  print(FileNames[[i]])
  dat = readRDS(paste0(DataDir,"/Data/", OrigNames[[i]]))
  dat$Name = FileNames[[i]]
  saveRDS(dat,paste0(DataDir,"/Data/",FileNames[[i]],".rdata"))
}
