WorkDir = "~/Dropbox/CoChIP/coChIP-Analysis"
DataDir = getwd()
BamDir = paste0(DataDir,"/BAM")

print("Initializing")
setwd(WorkDir)
suppressMessages(source("CoChip-Functions.R"))
#suppressMessages(source("NucAtlas.R"))
#GetSGDandNucs()
setwd(DataDir)

suppressMessages(library("optparse"))

option_list = list(
  make_option(c("-u", "--uniq"), action = "store_true", type="logical", default=FALSE, 
              help="filter uniq fragments"),
  make_option(c("-m", "--min"),  type="integer", default=0, 
              help="minimum width of fragments"),
  make_option(c("-M", "--max"),  type="integer", default=1000, 
              help="maximum width of fragments"),
  make_option(c("-d", "--datadir"), type = "character", default = NULL,
              help="data directory"),
  make_option(c("-s", "--suffix"), type = "character", default = "",
              help="suffic to add to the files"),
  make_option(c("-p", "--prefix"), type = "character", default = "",
              help="prefix to add to the files")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, positional_arguments = TRUE);
## prepere SGD/nuc structures


params <- ccParams()
params$DataDir = DataDir
params$reuseSavedData = TRUE
params$Save = FALSE
outparams = params
outparams$DataDir = DataDir
outparams$cov = TRUE
outparams$Save = TRUE
if( !is.null(opt$options$datadir) )
  outparams$DataDir = opt$options$datadir

Files = opt$args

DoProcessFile <- function( x ) {
  d = ccProcessFile(x, param = params)
  d$Name = paste0(opt$options$prefix,d$Name, opt$options$suffix)
 
  print(d$Name) 
  
  GR = d$GR
  if( opt$options$uniq )
    GR = unique(GR)
  W = width(GR)
  I = W >= opt$options$min & W <= opt$options$max
  d$GR = GR[ I ]
  d$nuc = NULL
  d$cov = NULL 
  d$cen = NULL
  d$meta = NULL 
  d$vplot = NULL
  ccProcessFile( dat = d, param = outparams )
}

## do the work...
dat = lapply(Files, DoProcessFile)

