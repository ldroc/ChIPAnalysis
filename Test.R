setwd("~/Dropbox/CoChIP/coChIP-Analysis")
source("CoChip-Functions.R")
source("NucAtlas.R")

donuc = 0
docov = 0
doQC = 0
doMeta = 0
doSize = 0
doMetaSize = 0
doVPlot = 1

MinReadNumber = 500000

DataDir = "~/Data/CoChIP/coChIP-Analysis/150909"
BamDir = paste0(DataDir,"/BAM")

## prepere SGD/nuc structures

if(file.exists("atlas")) {
  print("Loading Atlas")
  load(file="atlas")
} else {
  print("Read Atlas")
## SGD + TSS/TTS annotations
  SGD <- ReadGeneAtlas("SGD_TSS.tab")
  ## read nucs
  nucs <- ReadNucAtlas("Weiner-NucAtlas.csv", SGD )
  nucs.cov <- coverage(nucs)
  plus1nucs <- GetNucAlignedRegions(nucs, SGD)

  save(SGD, nucs, plus1nucs, file="atlas")
} 

rsfile = paste0(DataDir, "/rs.rdata")

if( file.exists(rsfile) ) {
  ## read from saved file
  print("Loading BAM data structure")
  rs = readRDS(file= rsfile)
} else {
  print("Read BAM")
  
  rs <- ReadBAMDirectory(BamDir, PairedEnd = TRUE )

  ## save the data we just loaded.
  saveRDS(rs, file=rsfile  )
}


## Do QC
if( doQC ) {
  QCdir = paste0(DataDir, "/QC")
  if( !dir.exists(QCdir) )
    dir.create(QCdir)

  if( dir.exists(QCdir) ) {
    print("Run QC")
    QCGRangesMultileFiles( rs,  path = QCdir)
  }
}

## From now on we work only with samples that are large enough

rs.lengths = lapply(rs, function(r) length(r$GR))
rs.select = rs.lengths > MinReadNumber

## Reduce fragements to centers

MaxFragLen = 600

rscovfile = paste0(DataDir,"/cov.rdata")
if( file.exists(rscovfile) ) {
  print("Loading coverage")
  load(rscovfile)
} else {
  print("Compute coverage")
  
  coveragecenters <- function( dat ) {
    I = width(dat$GR) < MaxFragLen
    
    list( BAMFile = dat$BAMFile, 
          GR = resize(dat$GR[I],width=50,fix="center")) 
  }
  
  rs.cen = lapply( rs[rs.select], coveragecenters )
  rs.files = unlist(lapply(rs[rs.select], function(r) BaseFileName(r$BAMFile) ))
  ## Coverage 
  FilterThreshold = .999
  

  docoverage <- function( dat ) {  
    c = coverage(dat$GR) 
  }

  rs.cov = lapply( rs.cen, docoverage)

  save(rs.cen, rs.cov, rs.files, file=rscovfile)
}

GeneRegions <- resize(plus1nucs, width=3000, fix="start")

if( doMeta ) {
genematfile = paste0(DataDir,"/genemat.rdata")
if( file.exists(genematfile)) {
  load(genematfile)
} else {
  print("Collect gene centered coverage")
  
  rs.TSS <- lapply(rs.cov, function (cov) CollectRegions(cov, SGD$TSSRegion) )
  rs.TTS <- lapply(rs.cov, function (cov) CollectRegions(cov, SGD$TTSRegion) )
  rs.P1S <- lapply(rs.cov, function (cov) CollectRegions(cov, plus1nucs) )
  rs.Genes <- lapply(rs.cov, function (cov) CollectRegions(cov, GeneRegions, 3000) )
  save(rs.TSS, rs.TTS, rs.P1S, rs.Genes, file=genematfile)
}


  print("Generate MetaGenes")
  
  TSS <- lapply( rs.TSS, function( mat ) AvgRegionsSubgroups(mat, SGD$TSSRegion, SGD$ExprQuan) )
  TTS <- lapply( rs.TTS, function( mat ) AvgRegionsSubgroups(mat, SGD$TTSRegion, SGD$ExprQuan) )
  P1S <- lapply( rs.P1S, function( mat ) AvgRegionsSubgroups(mat, plus1nucs, SGD$ExprQuan) )

  save(TSS, TTS, P1S, file=metafile)

  metagenedir = paste0(DataDir,"/MetaGene")
  if( !dir.exists(metagenedir) )
    dir.create(metagenedir)

  if( dir.exists(metagenedir)) {
    PlotMultiCoverage( rs.files, TSS, TTS, path =metagenedir )
  }
}

if( doSize ) {
  genematsizefile = paste0(DataDir,"/genematsize.rdata")
  if( file.exists(genematsizefile)) {
    load(genematsizefile)
  } else {
    print("Collect gene centered coverage by fragment sizes")
    
    rs.size.cov = lapply( rs[rs.select], function(r) CoverageByLength(r$GR,FragmentBoundaries))
    
    rs.size.genes <- lapply(rs.size.cov, function (cov) CollectRegionsBySize(cov, GeneRegions, 3000) )
    save(rs.size.cov, rs.size.genes, file=genematsizefile)
  }
}

if( doMetaSize ) {
  metagenedir = paste0(DataDir,"/MetaGeneSize")
  if( !dir.exists(metagenedir) )
    dir.create(metagenedir)
  
  if( dir.exists(metagenedir)) {
    AA = lapply( rs.size.genes, function(r) AvgRegionsLenBySubgroups(r,GeneRegions, SGD$ExprQuan) )
    PlotMultiSizeCoverage(rs.files,AA,path=metagenedir, pos=500, label="TSS")
  }
}

if( doVPlot ) {
  rss = rs[rs.select]
  
  vplotfile = paste0(DataDir,"/vplot.rdata")
  if( file.exists(vplotfile)) {
    load(vplotfile)
  } else {
    print("Collecting Vplots")
    tempFunc <- function(gr) {
      print(gr$BAMFile)
      cov = CreateVPlotCov(gr$GR)
      CollectGroupVPlots(cov,GeneRegions, SGD$ExprQuan, width=3000)
    }
    rs.vplots <- lapply(rss, tempFunc)
    save(rs.vplots, file=vplotfile)
  }
    
  vplotdir = paste0(DataDir,"/VPlots")
  if( !dir.exists(vplotdir) )
    dir.create(vplotdir)
  
  if( dir.exists(vplotdir)) {
    tempFunc <- function(gr, Ms) {
      print(gr$BAMFile)
      PlotMultiVPlot( rev(Ms), file = gr$BAMFile, path = vplotdir )
    }
    mapply(tempFunc, rss, rs.vplots)
  }
  for( i in c(1,3,5)) {
    cv = CreateVPlotCov(rs[[i]]$GR)
    Ms = CollectGroupVPlots(cv,GeneRegions, SGD$ExprQuan, width=3000)
    print(rs[[i]]$BAMFile)
    PlotMultiVPlot( rev(Ms), file = rs[[i]]$BAMFile, path = vplotdir )
  }
    
}
stop()

tabdir = paste0(DataDir,"/GeneTab")

# writing for TreeView 
if( 0 ) {
  plus1nucs$expr = SGD$Genes$expr[match(plus1nucs$acc,SGD$Genes$acc)]
  p1s.ord = order( plus1nucs$expr , na.last = TRUE)
  p1s.ord = p1s.ord[1:(length(p1s.ord)[1] - length(which(is.na(plus1nucs$expr))))]

  WriteExpressionTab(rs.P1S[[57]][p1s.ord,], "Flag80-K18-P1S.tab", SGD)

  WriteExpressionTab(rs.TSS[[57]], "test.tab", SGD)
}

nuccoverage <- function( dat ) { 
  c = unlist(countOverlaps(nucs, dat$GR)) 
  q = quantile(c,FilterThreshold)
#  print(q)
  c[c > q] <- NA
  return(c)
}

BAMname <- function( dat ) { x = dat$BAMFile; print(x); return(BaseFileName(dat$BAMFile)); }

if( donucs ) {
  nucfile = paste0(DataDir,"/nuc.rdata")
  if( file.exists(nucfile)) {
    load(nucfile)
  } else {
    print("Computing nucleosome coverage")
    rs.nuc = data.frame(lapply(rs.cen, nuccoverage))
    names(rs.nuc) = lapply(rs, BAMname)
    ##
    
    save(rs.nuc, file=nucfile)
  }
  
  RatioThreshold = 10
  eps = .1
  
  
  BuildRatioNucs <- function( i ) {
    j = i-1
    a = rs.nuc[[i]]
    b = rs.nuc[[j]]
    I = (b > RatioThreshold) & !is.na(a) & !is.na(b)
    ratios = vector(mode="numeric", length = length(a))
    ratios[I] = a[I] / (b[I] + eps)
    ratios[!I] = NA
    ratios
  }
  
  rs.nuc.ratios = lapply( seq(2, length(rs), by=2), BuildRatioNucs )
  names(rs.nuc.ratios) = lapply( rs[seq(2, length(rs), by=2)], BAMname)
  
  rs.nuc.cor=matrix(nrow=length(rs.nuc.ratios),ncol=length(rs.nuc.ratios))
  colnames(rs.nuc.cor) = names(rs.nuc.ratios)
  rownames(rs.nuc.cor) = names(rs.nuc.ratios)
  
  for( i in 1:length(rs.nuc.ratios) )
    for( j in 1:i ) {
      a = !is.na(rs.nuc.ratios[[i]])
      b = !is.na(rs.nuc.ratios[[j]])
      c = a & b
      I = any(c)
      r = 0
      if( I )
        r = cor(rs.nuc.ratios[[i]], rs.nuc.ratios[[j]], use="complete.obs")
      rs.nuc.cor[i,j] = r
      rs.nuc.cor[j,i] = r
    }
  
  corrplot(rs.nuc.cor)
  
}


PlotNucPair <- function (i,j) {
  xlab = paste(names(rs.nuc.ratios)[i] )
  ylab = paste(names(rs.nuc.ratios)[j] )
  title = paste(names(rs.nuc.ratios)[i], "vs", names(rs.nuc.ratios)[j])
  DensityScatter(rs.nuc.ratios[[i]],rs.nuc.ratios[[j]],threshold=1, 
                 xlabel = xlab, ylabel=ylab,alpha=0.05,
                 xlim=c(0,4.2),ylim=c(0,4.2))
}

if( 0 ) {
  
  g = arrangeGrob(PlotNucPair(7,5), PlotNucPair(7,6), PlotNucPair(8,5), PlotNucPair(8,6) )
  ggsave(file="CompareNucRatios-30-10.png", dpi=150, g)
  g = arrangeGrob(PlotNucPair(9,5), PlotNucPair(9,6), PlotNucPair(10,5), PlotNucPair(10,6) )
  ggsave(file="CompareNucRatios-100-10.png", dpi=150, g)
  g = arrangeGrob(PlotNucPair(9,7), PlotNucPair(9,8), PlotNucPair(10,7), PlotNucPair(10,8) )
  ggsave(file="CompareNucRatios-100-30.png", dpi=150, g)
  g = arrangeGrob(PlotNucPair(9,3), PlotNucPair(9,4), PlotNucPair(10,3), PlotNucPair(10,4) )
  ggsave(file="CompareNucRatios-100-03.png", dpi=150, g)
  
  g = arrangeGrob(PlotNucPair(19,17), PlotNucPair(19,18), PlotNucPair(20,17), PlotNucPair(20,18) )
  ggsave(file="CompareNucRatios-Myc-100-30.png", dpi=150, g)
  g = arrangeGrob(PlotNucPair(17,15), PlotNucPair(17,16), PlotNucPair(18,15), PlotNucPair(18,16) )
  ggsave(file="CompareNucRatios-Myc-30-10.png", dpi=150, g)
  g = arrangeGrob(PlotNucPair(17,13), PlotNucPair(17,14), PlotNucPair(18,13), PlotNucPair(18,14) )
  ggsave(file="CompareNucRatios-Myc-30-03.png", dpi=150, g)
  g = arrangeGrob(PlotNucPair(17,11), PlotNucPair(17,12), PlotNucPair(18,11), PlotNucPair(18,12) )
  ggsave(file="CompareNucRatios-Myc-30-00.png", dpi=150, g)
}

## Full bp resolution
if( docov ) {
  print("Computing coverage covariances")
  
  covs = lapply(rs.cov, function(x) as.vector(unlist(x)))

  rs.cor=matrix(nrow=length(rs.cov),ncol=length(rs.cov))
  colnames(rs.cor) = rs.files
  rownames(rs.cor) = rs.files

  for( i in 1:length(rs.cov) )
    for( j in 1:i ) {
      r = cor(covs[[i]], covs[[j]],use="complete.obs")
      rs.cor[i,j] = r
      rs.cor[j,i] = r
    }

  corrplot(rs.cor)
}

if( 0 ) {
N = 10000
DensityScatter(rs.ratios[[9]],rs.ratios[[10]])

PlotPair <- function (i,j) {
  xlab = paste("K36me/Input", names(rs.ratios)[i] )
  ylab = paste("K36me/Input", names(rs.ratios)[j] )
  title = paste( "K36me/Input", names(rs.ratios)[i], "vs", names(rs.ratios)[j])
  DensityScatter(rs.ratios[[i]],rs.ratios[[j]],threshold=1, title=title, xlabel = xlab, ylabel=ylab,alpha=0.05,xlim=c(0,2.2),ylim=c(0,2.2))
}

if( 0 ) {
  
  g = arrangeGrob(PlotPair(7,5), PlotPair(7,6), PlotPair(8,5), PlotPair(8,6) )
  ggsave(file="CompareRatios-30-10.png", dpi=150, g)
  g = arrangeGrob(PlotPair(9,5), PlotPair(9,6), PlotPair(10,5), PlotPair(10,6) )
  ggsave(file="CompareRatios-100-10.png", dpi=150, g)
  g = arrangeGrob(PlotPair(9,7), PlotPair(9,8), PlotPair(10,7), PlotPair(10,8) )
  ggsave(file="CompareRatios-100-30.png", dpi=150, g)
  
}


## Compute correlation matrix
#I = apply(rs.nuc,1, function(x) !any(is.na(x)))
rs.cor = cor(rs.nuc, use="complete.obs")

ratmat = sapply(rs.ratios, function( r ) as.vector(r$ratio))
colnames(ratmat) = lapply(rs.ratios, function(x) x$name )
rs.cor = cor(ratmat, use="complete.obs")

##
## Compute linear relations between samples
##
doslope <- function(x, y) { l = lm(y ~ -1+x); coef(l)}
n = length(rs.ratios)
rs.slope = matrix(nrow=n,ncol=n)
for( i in 1:n) 
  for( j in 1:n )
    rs.slope[i,j] = doslope(rs.ratios[[i]]$ratio,rs.ratios[[j]]$ratio)

#colnames(rs.slope) = lapply(rs, BAMname)
#rownames(rs.slope) = lapply(rs, BAMname)
colnames(rs.slope) = lapply(rs.ratios, function(x) x$name )
rownames(rs.slope) = lapply(rs.ratios, function(x) x$name )

DensityScatter(rs.nuc[1],rs.nuc[2])




## Old

BuildRatio <- function( i ) {
  j = i-1
  a = covs[[i]]
  b = covs[[j]]
  I = covs[[j]] > RatioThreshold
  ratios = vector(mode="numeric", length = length(a))
  ratios[I] = a[I] / (b[I] + eps)
  ratios[!I] = NA
  ratios
}

rs.ratios = lapply( seq(2, length(rs), by=2), BuildRatio )
names(rs.ratios) = lapply( rs[seq(2, length(rs), by=2)], BAMname)

#rs.cor = cor(sapply(rs.ratios, function(r) as.vector(r)), use="complete.obs")

plot( cov1, cov3, xlim=c(0,quantile(cov1,.999)), ylim=c(0,quantile(cov3,.999)),
       col=rgb(0,10,0,10,maxColorValue=255), 
       pch=20,ps=2)

thresh = .998
I = cov1 < quantile(cov1,thresh) & cov2 < quantile(cov2,thresh)
J = cov1 < quantile(cov1,thresh) & cov3 < quantile(cov3,thresh)

ggplot(data.frame(x=cov1[I],y=cov2[I]),aes(x=x,y=y)) + 
  geom_point(colour="blue", alpha=0.01) + 
  geom_density2d(colour="black")
#+stat_density2d(aes(fill = ..level..), geom="polygon")

ggplot(data.frame(x=cov1[J],y=cov3[J]),aes(x=x,y=y)) + 
  geom_point(colour="blue", alpha=0.01) + 
  geom_density2d(colour="black")

lines(lowess(cov1,cov2),col='red')


}