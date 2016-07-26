WorkDir = "~/Dropbox/CoChIP/coChIP-Analysis"
DataDir = getwd()
setwd(WorkDir)


DataDir =  "~/Data/CoChIP/HisMut-K4/HisMut"  
OrigDir = DataDir
BamDir = paste0(DataDir,"/BAM")
Trackdir = DataDir
FigureDir = "~/Dropbox/CoChIP/figures/"

params = ccParams()
params$DataDir = paste0(DataDir,"/Data")
params$cov = FALSE
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

print("Initializing")
suppressMessages(source("CoChip-Functions.R"))
suppressMessages(source("NucAtlas.R"))
suppressMessages(source("DensityScatter.R"))
source("Main-Functions.R")
GetSGDandNucs()

Exp = c("WT-H3K18ac-Input", "WT-H3K4me3-Input", "WT-H3K36me3-Input", "WT-H3K79me3-Input",
        "WT_merge-H3K18ac-H3K36me3", "WT_merge-H3K18ac-H3K4me3", "WT_merge-H3K36me3-H3K79me3", "WT_merge-H3K4me3-H3K79me3" 
)
for( m in Exp )
{
  print(m)
  dat = ccProcessFile(filename = m, param = params )
#  fn = paste0(FigureDir,"/Heat-expr-TSS-", m,".png")
#  PlotHeatMatrix( dat, fn, type = "TSS", order = "expr", filetype = "png"  )
#  fn = paste0(FigureDir,"/Heat-length-TSS-", m,".png")
#  PlotHeatMatrix( dat, fn, type = "TSS", order = "length", filetype = "png"  )
  fn = paste0(FigureDir,"/Heat-occ-TSS-", m,".png")
  PlotHeatMatrix( dat, fn, type = "TSS", order = "occupancy", filetype = "png"  )
}
