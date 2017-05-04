WorkDir = "~/Dropbox/coChIP-scripts/"
DataDir = getwd()
setwd(WorkDir)


#DataDir =  "~/Data/CoChIP/160124/RPD3/"
OrigDir = DataDir
BamDir = paste0(DataDir,"/BAM")
Trackdir = DataDir

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
GetSGDandNucs()

suppressMessages(library(rtracklayer))
suppressMessages(library(optparse))

option_list = list(
  make_option("--maxfraglen", type = "integer", default = 600,
              help="maximal fragment length"),
  make_option("--minfraglen", type = "integer", default = 50,
              help="minimum fragment length"),
  make_option(c("-T", "--tilewidth"), type = "integer", default = 10,
              help="tile width length"),
  make_option(c("-w", "--width"), type = "integer", default = 50,
              help="fragments trimmed to this width (0 - no trimming)"),
  make_option(c("-t", "--trackdir"), type = "character", default = NULL,
              help="track directory"),
  make_option(c("-d", "--datadir"), type = "character", default = NULL,
              help="data directory"),
  make_option("--bamdir", type="character", default=NULL,
              help="location of BAM files" )
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, positional_arguments = TRUE);
## prepere SGD/nuc structures

params <- ccParams()
params$MaxFragLen = opt$options$maxfraglen
params$MinFragLen = opt$options$minfraglen

if( !is.null(opt$options$bamdir) )
  BamDir = extendDir(opt$options$bamdir)

if( !is.null(opt$options$trackdir) )
  TrackDir = extendDir(opt$options$trackdir)

if( !is.null(opt$options$datadir)) {
  DataDir = extendDir(opt$options$datadir)
}
params$DataDir = DataDir

# Files = c("WT-H3-H3", "WT-H3-H3K4ac")
# Files = c("bar1-K18ac-Input")
Files = opt$args

# TileWidth = 20
TileWidth = opt$options$tilewidth

Tiles = tileGenome(seqinfo(SGD$Genes),tilewidth = TileWidth,cut.last.tile.in.chrom = TRUE)
tparams = params
tparams$DataDir = TrackDir

DoProcessFile <- function( x ) {
  temp = Tiles
  d = ccProcessFile(paste0(BamDir,"/",x), param = params)
  I = width(d$GR) > params$MinFragLen & width(d$GR) <= params$MaxFragLen 
  GR = unique(d$GR[I])
  
  if(opt$options$width > 0 )
    GR = resize(GR, width=opt$options$width,fix="center")
  
  score(temp) = countOverlaps(Tiles,GR)
  fname = ccBuildFN(d$Name,tparams, suff = ".bw")
  export( temp, fname, format="BigWig" )  
}

dat = lapply(Files, DoProcessFile)

