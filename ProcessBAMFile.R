#!/usr/bin/env Rscript --vanilla
WorkDir = "~/Dropbox/coChIP-scripts"
#Find where are the scripts at
initial.options <- commandArgs(trailingOnly = FALSE)
#file.arg.name <- 
#script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
#WorkDir <- dirname(script.name)


DataDir = getwd()
OrigDir = getwd()
BamDir = paste0(DataDir,"/BAM")



extendDir <- function(x) {
  a = substr(x,1,1)
  if( a == "/" || a == "~")
    return(x)
  else
    return(paste0(OrigDir,"/",x))
}

print("Initializing")
setwd(WorkDir)
suppressMessages(source("CoChip-Functions.R"))
suppressMessages(source("NucAtlas.R"))
GetSGDandNucs()
setwd(DataDir)

suppressMessages(library("optparse"))

option_list = list(
  make_option(c("-q", "--QC"), action = "store_true", type="logical", default=FALSE, 
              help="run QC"),
  make_option(c("-c", "--coverage"), action = "store_true", type="logical", default=FALSE, 
              help="compute coverage"),
  make_option(c("-C", "--centers"), action = "store_true", type="logical", default=FALSE, 
              help="compute coverage of centers"),
  make_option(c("-n", "--nucleosomes"), action = "store_true", type="logical", default=FALSE, 
              help="compute over nucleosome positions"),
  make_option(c("-x", "--genesheet"), type="character", default=NULL, 
              help="Generate genesheets"),
  make_option(c("-T", "--tracks"), action = "store_true", type="logical", default=FALSE, 
              help="generate browser tracks"),
  make_option(c("-m", "--metagene"), action = "store_true", type="logical", default=FALSE, 
              help="plot metagene"),
  make_option(c("-N", "--NucFile"), type="character", default=NULL,
              help="Collect nucleosome counts to specified RDATA file"),
  make_option(c("--Nucsummary"), type="character", default=NULL,
              help="Collect nucleosome counts to specified CSV file"),
  make_option(c("-S", "--QCsummary"), type="character", default=NULL, 
              help="Summarize QC to specificed file"),
  make_option(c("--sizes"),type="character", default="1,1000", 
              help="Comma seperated list of fragment sizes for QC summary"),
  make_option("--centerwidth", type = "integer", default = 50,
              help="width of center window"),
  make_option("--tilewidth", type = "integer", default = 10,
              help="track tile width length"),
  make_option("--nucwidth", type = "integer", default = 150,
              help="width of nucleosome window"),
  make_option("--maxfraglen", type = "integer", default = 600,
              help="maximal fragment length"),
  make_option("--minfraglen", type = "integer", default = 50,
              help="minimum fragment length"),
  make_option(c("-d", "--datadir"), type = "character", default = NULL,
              help="data directory"),
  make_option(c("-b", "--bamdir"), type="character", default=NULL,
              help="location of BAM files" ),
  make_option("--QCdir", type="character", default=NULL,
              help="location to write QC files" ),
  make_option("--trackdir", type = "character", default = NULL,
              help="track directory"),
  make_option("--metadir", type="character", default=NULL,
              help="location to write metagene files" ),
  
  make_option(c("-F", "--force"), action = "store_true", type="logical", default=FALSE, 
              help="force recomputing")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, positional_arguments = TRUE);
## prepere SGD/nuc structures


params <- ccParams()
params$TSSRegions = SGD$TSSRegion
params$TTSRegions = SGD$TTSRegion
params$GeneRegions = SGD$GeneRegion
params$PlusOneRegions = Plus1Nucs
params$NucRegions = NucRegions
params$MaxFragLen = opt$options$maxfraglen
params$MinFragLen = opt$options$minfraglen
params$cov = opt$options$coverage
params$cen = opt$options$centers
params$CenWidth = opt$options$centerwidth
params$NucRegionWidth = opt$options$nucwidth
params$meta = opt$options$metagene

TileWidth = opt$options$tilewidth
if( params$meta )
  params$cen = TRUE
params$nuc = opt$options$nucleosomes

Force = opt$options$force

if( !is.null(opt$options$datadir) )
  DataDir = extendDir(opt$options$datadir)
params$DataDir = DataDir

doQC = opt$options$QC
doMeta = opt$options$metagene
doTracks = opt$options$tracks

if( !is.null(opt$options$bamdir) )
  BamDir = extendDir(opt$options$bamdir)


MetaDir = DataDir
if( !is.null(opt$options$metadir) )
  MetaDir = extendDir(opt$options$metadir)

TrackDir = DataDir
if( !is.null(opt$options$trackdir) ) {
  TrackDir = extendDir(opt$options$trackdir)
}

QCDir = DataDir
if( !is.null(opt$options$QCdir) )
  QCDir = extendDir(opt$options$QCdir)

if( !is.null(opt$options$datadir)) {
  setwd(opt$options$datadir)
  DataDir = getwd()
  params$DataDir = DataDir
}

doQCSum = FALSE
if( !is.null(opt$options$QCsummary) ) {
  doQCSum = TRUE
  QCSumFN = opt$options$QCsummary
}

doNucs = FALSE
if( !is.null(opt$options$NucFile)) {
  doNucs = TRUE
  NucFN = opt$options$NucFile
  params$nuc = TRUE
}


doNucsCSV = FALSE
if( !is.null(opt$options$Nucsummary)) {
  doNucsCSV = TRUE
  NucCSVFN = opt$options$Nucsummary
  params$nuc = TRUE
}
doGeneSheets = FALSE
if( !is.null(opt$options$genesheet)) {
  doGeneSheets = TRUE
  GeneSheetsFN = extendDir(opt$options$genesheet)
  params$meta = TRUE
  params$cen = TRUE
}

Files = opt$args

if( doQCSum )
  QCSizeRanges = as.integer(unlist(strsplit(opt$options$sizes,",")))

if( doTracks ) {
  Tiles = tileGenome(seqinfo(SGD$Genes),tilewidth = TileWidth,cut.last.tile.in.chrom = TRUE)
  params$cen = TRUE
}

QCSummary = list()
Nucs = list()
DoProcessFile <- function( x ) {
  d = ccProcessFile(paste0(BamDir,"/",x), param = params, Force = Force)
  if( doQC ) 
    ccDoQC(d, QCDir)
  if( doQCSum ){
    print(d$Name)
    s = QCSummaryRanges(d$GR, QCSizeRanges)
    QCSummary[[d$Name]] <<- s
  }
  if( doNucs || doNucsCSV ) {
    Nucs[[d$Name]] <<- d$nuc
  }
  
  if( doMeta ) 
    ccDoMeta(d, params, MetaDir = MetaDir)
  
  if( doGeneSheets ) 
    ccExportGeneSheet( d, params, Dir = GeneSheetsFN, type = "Genes" )
  
  if( doTracks ) 
    ccExportTrack( d, params, Tiles, TrackDir)
}

## do the work...
dat = lapply(Files, DoProcessFile)

if( doQCSum ) {
  D = sapply(QCSummary, unlist, simplify = "array")
  write.csv(t(D), file = extendDir(QCSumFN))
}

if( doNucs ) {
  Ns = do.call(rbind, Nucs)
  rownames(Ns) = names(Nucs)
  saveRDS(Ns, file = extendDir(NucFN))
}

if( doNucsCSV ) {
  Ns = do.call(rbind, Nucs)
  rownames(Ns) = names(Nucs)
  write.csv(Ns, file = extendDir(NucCSVFN),quote = FALSE)
}


