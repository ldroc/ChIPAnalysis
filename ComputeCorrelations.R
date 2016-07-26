WorkDir = "~/Dropbox/CoChIP/coChIP-Analysis"
DataDir = getwd()
setwd(WorkDir)


#DataDir =  "~/Data/CoChIP/HisMut-K4/K4"  
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
GetSGDandNucs()

suppressMessages(library("optparse"))

option_list = list(
  make_option(c("-o", "--output"), type = "character", default = NULL, 
              help = "output file"),
  make_option(c("-p", "--plot"), type = "character", default = NULL, 
              help = "plot to file"),
  make_option(c("-b", "--baseline"), type = "character", default = NULL,
              help = "comma seperated list of baseline files"),
  make_option("--maxfraglen", type = "integer", default = 600,
              help="maximal fragment length"),
  make_option("--minfraglen", type = "integer", default = 50,
              help="minimum fragment length"),
  make_option("--minreads", type = "integer", default = 100000,
              help="minimum number of reads"),
  make_option(c("-n", "--nucfile"), type = "character", default = NULL, 
              help = "name of nucleosome coverage file"),
  make_option(c("-d", "--datadir"), type = "character", default = NULL,
              help="data directory"),
  make_option("--bamdir", type="character", default=NULL,
              help="location of BAM files" )
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, positional_arguments = TRUE);
## prepere SGD/nuc structures

MinReads = opt$options$minreads
OutputFile = opt$options$output
PlotFile = opt$options$plot
BaseLineList = opt$options$baseline

params <- ccParams()
params$TSSRegions = SGD$TSSRegion
params$TTSRegions = SGD$TTSRegion
params$GeneRegions = SGD$GeneRegion
params$NucRegions = NucRegions
params$MaxFragLen = opt$options$maxfraglen
params$MinFragLen = opt$options$minfraglen
params$cen = FALSE
params$meta = FALSE
params$nuc = TRUE

if( !is.null(opt$options$bamdir) )
  BamDir = extendDir(opt$options$bamdir)

if( !is.null(opt$options$datadir)) {
  DataDir = extendDir(opt$options$datadir)
#  setwd(DataDir)
}
params$DataDir = DataDir

#NucFile = "TestNuc"
NucFile = opt$options$nucfile

# Files = c("WT-H3-H3", "WT-H3-H3K4ac")
Files = opt$args

DoProcessFile <- function( x ) {
  d = ccProcessFile(paste0(BamDir,"/",x), param = params)
  list( Name = x, Nucs = d$nuc, UniqReads = length(d$UniqGR) )
}

if( !is.null(NucFile) )
  NucFile = paste0(DataDir,"/",NucFile,".rdata")

## get the data
if( !is.null(NucFile) && file.exists(NucFile) ) {
  nuc.data =readRDS(NucFile)
  nuc.counts = nuc.data$counts
  nuc.mat = nuc.data$mat
} else {
  dat = lapply(Files, DoProcessFile)
  nuc.mat = do.call(rbind,lapply(dat, function(x) x$Nucs ))
  rownames(nuc.mat) = lapply(dat, function(x) x$Name )
  nuc.counts = unlist(lapply(dat, function(x) x$UniqReads ))
  nuc.data = list( mat = nuc.mat, counts = nuc.counts )
  if( !is.null(NucFile) )
    saveRDS(nuc.data,NucFile)
}

names(nuc.counts) = rownames(nuc.mat)

if( 1 ) {
  # get rid of outliers
  tot = colMeans(nuc.mat)
  t = quantile(tot, .99)
  m = mean(tot[tot <= t])
  s = sqrt(var(tot[tot <= t]))
  t = m + 3*s
  Io = tot < t
  nuc.mat = nuc.mat[,Io]
}

if(1) {
  #normalize
  for(i in rownames(nuc.mat)) {
    s = sum(nuc.mat[i,])/MinReads
    nuc.mat[i,] = nuc.mat[i,]/s
  }
}

I = nuc.counts > MinReads
nuc.smat = nuc.mat[I,]
print(dim(nuc.smat))

if( !is.null(BaseLineList) ) {
  l = unlist(strsplit(BaseLineList, ","))

  ref = log2(colMeans(nuc.smat[is.element(rownames(nuc.smat),l),]) + 0.1)
  d = dim(nuc.smat)
  refm = matrix(ref, nc = d[[2]], nr = d[[1]], byrow=TRUE)
  
  nuc.smat = log2(nuc.smat + 0.1) - refm
  nuc.smat = nuc.smat[!is.element(rownames(nuc.smat),l),]
}

nuc.cor = cor(t(nuc.smat), use="complete.obs")

if( !is.null(OutputFile) ) {
  fn = paste0(extendDir(OutputFile),".csv")
  write.csv(nuc.cor, file=fn)
  fn = paste0(extendDir(OutputFile),"-counts.csv")
  write.csv(nuc.counts, file=fn)
}

if( !is.null(PlotFile) ) {
  fn = paste0(extendDir(PlotFile),".png")
  png(filename = fn, width = 1535, height = 1536)
  corrplot(nuc.cor, order="hclust", hclust.method = "complete")
  dev.off()
  fn = paste0(extendDir(PlotFile),"-corr.png")
  png(filename = fn, width = 800, height = 800)
  hist(nuc.cor, breaks = "FD")
  dev.off()
    
  fn = paste0(extendDir(PlotFile),"-corr-to-count.png")
  png(filename = fn, width = 800, height = 800)
  dd = dim(nuc.cor)
  cnt = matrix(nr = dd[[1]], nc = dd[[2]] )
  rn = rownames(nuc.cor)
  for( i in 1:dd[[1]])
    if( i < dd[[2]])
    for( j in (i+1):dd[[2]] ) {
#        print(c(rn[[i]], rn[[j]]))
#        print(nuc.counts[c(rn[[i]],rn[[j]])])
        cnt[i,j] = min(nuc.counts[c(rn[[i]],rn[[j]])])
#        print(cnt[i,j])
    }
    
#  DensityScatter(cnt,nuc.cor, alpha = 1,coordeq = F, xlabel = "Uniq Reads", ylabel = "Correlation")
  plot(cnt,nuc.cor, log="x",xlab="Uniq reads", ylab="cor", type="p", cex=.5)
  dev.off()
}
