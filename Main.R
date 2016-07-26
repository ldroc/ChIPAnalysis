setwd("~/Dropbox/CoChIP/coChIP-Analysis")
source("CoChip-Functions.R")
source("coChIP-GLM.R")
source("NucAtlas.R")
source("Multiplot.R")

#DataDir = "~/Data/CoChIP/151029-SizeChIP/Short"
#DataDir = "~/Data/CoChIP/151112-iChIP"
#DataDir = "~/Data/CoChIP/HisMut-K4/HisMut/Data"  
#DataDir = "~/Google Drive/Data/160124-RPD3/Data"
DataDir = "~/Data/CoChIP/160124/RPD3"
#DataDir = "~/Data/CoChIP/160124/K4Sym"
FigureDir = "~/Dropbox/CoChIP/figures/"

## prepere SGD/nuc structures
GetSGDandNucs()
source("Main-Functions.R")


TargetMods = c("H3K18ac", "H3K4ac", "H3K36me3" , "H3K79me3","H3K4me3", "H3K56ac" )
mod.nucs = ReadHistoneModeAtlas("Weiner-HisMod.csv", NucRegions, TargetMods )
m = mcols(mod.nucs)
weiner.nucs = do.call(rbind,lapply(TargetMods, function(mod) as.vector(exp(as.matrix(m["Input"]) + as.matrix(m[mod])))))
weiner.nucs = rbind(weiner.nucs, as.vector(exp(as.matrix(m["Input"]))))
rownames(weiner.nucs) = c(TargetMods, "Input")
nuc.1 = which(mcols(NucRegions)$gene_pos == 1)
nuc.3 = which(mcols(NucRegions)$gene_pos == 3)
nuc.5 = which(mcols(NucRegions)$gene_pos == 5)

if(0) {
  p = DensityScatter(weiner.nucs["H3K56ac",],weiner.nucs["H3K18ac",],
                     xlabel = "H3K56ac", ylabel = "H3K18ac", 
                     coordeq = T, threshold = .995, 
                     logscale = T, alpha=.15, 
                     xlim =c(0.05,15),ylim =c(0.05,15))
  ggsave(filename = paste0(FigureDir,"/Weiner-K56ac-K18ac.png"), 
         p, width = 3, height=3)
  p = DensityScatter(weiner.nucs["H3K56ac",],weiner.nucs["H3K4me3",],
                     xlabel = "H3K56ac", ylabel = "H3K4me3", 
                     coordeq = T, threshold = .995, 
                     logscale = T, alpha=.15, 
                     xlim =c(0.05,15),ylim =c(0.05,15))
  ggsave(filename = paste0(FigureDir,"/Weiner-K56ac-K4me3.png"), 
         p, width = 3, height=3)
  p = DensityScatter(weiner.nucs["H3K18ac",],weiner.nucs["H3K4me3",],
                     xlabel = "H3K18ac", ylabel = "H3K4me3", 
                     coordeq = T, threshold = .995, 
                     logscale = T, alpha=.15, 
                     xlim =c(0.05,15),ylim =c(0.05,15))
  ggsave(filename = paste0(FigureDir,"/Weiner-K18ac-K4me3.png"), 
         p, width = 3, height=3)
  
  p = DensityScatter(weiner.nucs["H3K4ac",],weiner.nucs["H3K4me3",],
                     xlabel = "H3K4ac", ylabel = "H3K4me3", 
                     coordeq = T, threshold = .995, 
                     logscale = T, alpha=.15, 
                     xlim =c(0.05,15),ylim =c(0.05,15))
  ggsave(filename = paste0(FigureDir,"/Weiner-K4ac-K4me3.png"), 
         p, width = 3, height=3)
  p = DensityScatter(weiner.nucs["H3K36me3",],weiner.nucs["H3K4me3",],
                     xlabel = "H3K36me3", ylabel = "H3K4me3", 
                     coordeq = T, threshold = .995, 
                     logscale = T, alpha=.15, 
                     xlim =c(0.05,15),ylim =c(0.05,15))
  ggsave(filename = paste0(FigureDir,"/Weiner-K36me3-K4me3.png"), 
         p, width = 3, height=3)
  
}

params <- ccParams()
params$cov = TRUE
params$cen = TRUE
params$DataDir = paste0(DataDir,"/Data")
params$nuc = TRUE
params$meta = TRUE
params$TSSRegions = SGD$TSSRegion
params$TTSRegions = SGD$TTSRegion
params$NucRegions = NucRegions
params$PlusOneRegions = Plus1Nucs
params$GeneRegions = SGD$GeneRegion
params$MaxFragLen = 220
params$MinFragLen = 50

Files = list.files( paste0(DataDir,"/Data"), "rdata$", full.names = F )
#Files = Files[grep("H4", Files, invert = T)]

if(0) {
#  Files = Files[grep("K18ac.rdata", Files, invert = T)]
  ReadNucs <- function(f) {
    print(f)
    d = readRDS(paste0(DataDir,"/Data/",f))
    if( is.null(d$nuc) )
      d = ccProcessFile(d,param = params)
    d$nuc
  }
  nucs.list = lapply(Files,ReadNucs)
  nucs <- do.call(rbind,nucs.list)

#if(0)
#  dat <- lapply(dat, function(d) {d$meta = NULL; d})

#if( 0 )
#  dat <- lapply(dat, function(x) ccProcessFile(dat=x, param = params))
#nucs <- do.call(rbind,lapply( dat, function(x) x$nuc ) )
#  Names = unlist(lapply(dat, function(x) x$Name))
  Names = Files
  Names = sub("bar1", "WT", Names)
  Names = sub(".rdata", "", Names)
  Names = sub("_merge", "", Names)
  
  #  Names = unlist(lapply(Names, function(x) {z = strsplit(x,"_"); z[[1]][1]}))

  Strain <- factor(unlist(lapply(strsplit(Names,"-"), function(x) x[[1]])))
  Ab1 <- factor(unlist(lapply(strsplit(Names,"-"), function(x) x[[2]])))
  Ab2 <- factor(unlist(lapply(strsplit(Names,"-"), function(x) x[[3]])))
  Names <- mapply(function(s,x,y) paste0(s,"-",x,"-", y), Strain,Ab1,Ab2)
  rownames(nucs) <- Names
  saveRDS(nucs,paste0(DataDir,"nucs.rdata"))
}

nucs = readRDS(paste0(DataDir,"nucs.rdata"))
Names = rownames(nucs)
Strain <- factor(unlist(lapply(strsplit(Names,"-"), function(x) x[[1]])))
Ab1 <- factor(unlist(lapply(strsplit(Names,"-"), function(x) x[[2]])))
Ab2 <- factor(unlist(lapply(strsplit(Names,"-"), function(x) x[[3]])))

if(0) {
cnucs <-rbind(nucs, weiner.nucs)

Inputs = weiner.nucs["Input",]
t.high = quantile(Inputs,0.99, na.rm=TRUE)
t.low = quantile(Inputs,0.01, na.rm=TRUE)
nucs.abnormal = Inputs > t.high | Inputs < t.low
cnucs[,nucs.abnormal] = NA
Inputs[nucs.abnormal] = NA
Names = rownames(cnucs)
} else
  cnucs <- nucs

if( 0 ) {
  Sizes = rowSums(nucs)
  I = Ab2 != "K18ac" & Ab1 != "H4" & Ab2 != "H4" & Sizes > 20000
#  I = Ab1 != "K4un" & Ab2 != "K4un" & Sizes > 5000
  cnucs <- nucs[I,]
  Names <- rownames(cnucs)
}
# Compute correlation matrix
if( 0 ) {
  Cor = matrix(nc = length(Names), nr = length(Names))
  rownames(Cor) = Names
  colnames(Cor) = Names
  for( i in Names) 
    for (j in Names) 
      Cor[i,j] = cor(cnucs[i,], cnucs[j,])
  
  pdf(paste0(FigureDir, "/K4Sym-NucCorrelations.pdf"),width=25,height=25)
  corrplot(corr = Cor, order="hclust")
  dev.off()
}

if( 0 ) {
  I = Ab2 == "Input" &  Ab1 != "H4" 
  Dat = list()
  for( i in which(I) ) {
    dat = readRDS(paste0(DataDir,"/Data/", Files[[i]]))
    #   params$meta = TRUE
    #    dat$meta = NULL
    Dat[[i]] = ccProcessFile(dat = dat, param = params)
  }
}

if( 0 ) {
  params$meta = TRUE
  for( i in 1:length(Files) ) {
    print(c(i,Files[[i]]))
    dat = readRDS(paste0(DataDir,"/Data/", Files[[i]]))
    if( is.null(dat$meta[[3]])) {
      dat$meta = NULL
      ccProcessFile(dat = dat, param = params)
    }
  }
}

if( 0 ) {
  for( a in levels(Ab1)) {
    i.wt = which( Names == paste0("WT-",a,"-Input"))
    dat.wt =  readRDS(paste0(DataDir,"/Data/", Files[[i.wt]]))
    for( s in levels(Strain) )
      if( s != "WT") {
        i.s = which( Names == paste0(s,"-",a,"-Input"))
        dat.s =  readRDS(paste0(DataDir,"/Data/", Files[[i.s]]))
        PlotComparison(dat.wt,dat.s, paste0("WT-",s,"-",a,"-long"), order = "length", 
                       type="Genes", MinLength = 1800) 
      }
  }
}  

if(0) {
  NETDataDir = "~/Data/CoChIP/160124/NET-seq"
  Set2Transcripts = readRDS(paste0(NETDataDir,"/set2Trans.rdata"))
  Set2Regions = resize(resize(Set2Transcripts, width=1000, fix="start"), width=1500, fix="end" ) 
  J = GeneLengths < 1500 & GeneLengths > 500
  for( a in levels(Ab1)) {
    Ms = list()
    Ts = list()
    for( s in c("bar1", "set2", "eaf3") ) {
      print(paste0(s,"-",a))
      dat = readRDS(paste0(DataDir,"/Data/",s,"-",a,"-Input.rdata"))
      meta = CollectRegions(dat$cov, Set2Regions)
      mTSS  =dat$meta[[1]]
      mGenes = dat$meta[[3]]
      N = mean(mGenes[J,500:2000], na.rm = T)
      Ms[[s]] = colMeans(meta)/N
      Ts[[s]] = colMeans(mTSS[is.element(SGD$TSSRegion$acc, SGD$ExprQuan[[4]]),1:1500])/N
    }
    Mlen = length(Ms)
    df = data.frame( x = rep(-499:1000, 2*Mlen), y = do.call("c",c(Ms, Ts)), 
                     var = factor(c(rep(names(Ms), each = 1500),
                                    rep(paste(names(Ts),"TSS"), each = 1500))))
    p = ggplot(df) 
    p = p + xlab("Position") + ylab(a)
    p = p + geom_line(aes(x=x,y=y,colour = var),size=1.5)
    labelat = seq(-500,1000, 500)
    labelval = labelat
    labelval[[2]] = "Start"
    p = p+scale_x_continuous(name="Position", breaks = labelat, labels = labelval )
    ggsave(filename = paste0(FigureDir,"/set2-trans-meta-",a,".pdf"),
           plot=p,
           width=10,height=6, unit="in")
  }
}

if(0) {
  m = "K18ac"
  r = "K36me3"
  s = "eaf3"
  for( s in c("eaf3","set2") )
    for( m in c("K18ac", "K56ac", "K4me3"))
      for( r in  c("K18ac", "K56ac", "K4me3", "K36me3"))
        if( r != m ) {
          WTmI = paste0("WT-",m,"-Input")
          SmI = paste0(s,"-",m,"-Input")
          WTrI = paste0("WT-",r,"-Input")
          SrI = paste0(s,"-",r,"-Input")
          if( m < r ) {
            WTmr = paste0("WT-", m, "-", r)
          } else 
            WTmr = paste0("WT-", r, "-", m)
          
          WTthresh = quantile(nucs[WTmI,], 0)
          Sthresh = quantile(nucs[SmI,], 0.4)
          Sdiff = quantile(nucs[SmI,], 0.9) - Sthresh
          lo <- rlm(y~x-1, data.frame(x = nucs[WTmI,], y = nucs[SmI,]), method="MM")
          a = coef(lo)
          Modup = nucs[SmI,] > pmax(Sthresh,  a*nucs[WTmI,]+Sdiff) &
            nucs[WTmI,] > WTthresh
          #  Sample = sample(dim(nucs)[[2]], 10000 )
          Sample = 1:dim(nucs)[[2]]
          UpGroups = list(Modup[Sample])
          names(UpGroups) = list(paste(m, " up in ", s))
          p = DensityScatter(nucs[WTmI,Sample],nucs[SmI,Sample],
                             threshold=.99,coordeq = F, alpha=.2, jitter = T, 
                             Groups = UpGroups, pointsize = 1,galpha=.2,
                             title = paste(m," - WT", s), 
                             xlabel =  paste("WT",m), ylabel = paste(s,m))
          p = p+geom_abline(slope = a, interscept=0, color = "red", size = .5)
          p = p+scale_colour_brewer(palette = "Set1")
          ggsave(filename = paste0(FigureDir,"/",m,"UP-", s, "def.png"),p,width=6,height=5)
          p = DensityScatter(nucs[WTmI,Sample],nucs[WTrI,Sample],
                             threshold=.99,coordeq = F, alpha=.2, jitter = T, 
                             Groups = UpGroups, pointsize = 1,galpha=.2,
                             title = paste("WT", m, "vs",r), 
                             xlabel = paste("WT",m), ylabel = paste("WT",r),
                             GroupContour = T)
          p = p+scale_colour_brewer(palette = "Set1")
          ggsave(filename = paste0(FigureDir,"/",m,"UP-", s, "-", m,"vs",r,".png"),p,width=6,height=5)
          
          p = DensityScatter(nucs[WTrI,Sample],nucs[SrI,Sample],
                             threshold=.99,coordeq = F, alpha=.2, jitter = T, 
                             Groups = UpGroups, pointsize = 1,galpha=.2,
                             title = paste(m, "- WT vs",s), 
                             xlabel = paste("WT",r), ylabel = paste(s,r),
                             GroupContour = T)
          p = p+scale_colour_brewer(palette = "Set1")
          ggsave(filename = paste0(FigureDir,"/",m,"UP-", s, "-", r,"WT-vs-",s,".png"),p,width=6,height=5)
        }
}
if( 0 ) {
  GeneLengths = width(SGD$Genes)
  GeneExpr = SGD$Genes$expr
  I = GeneLengths > 2000 & GeneExpr > 10
  J = GeneLengths < 1500 & GeneLengths > 500
  for( a in levels(Ab1)) {
    print(paste("First IP = ", a))
    for( a2 in levels(Ab2) ){
      print(paste("Second IP = ", a2))
      Ms = list()
      for( s in levels(Strain) ) {
#        if( length(which( Names == paste0(s,"-",a,"-",a2))) != 1 )
#          print(paste0(s,"-",a,"-",a2))
        is = which( Names == paste0(s,"-",a,"-",a2))
        if( length(is) > 0 ) {
          i.s = is[[1]]
          dat = readRDS(paste0(DataDir,"/Data/", Files[[i.s]]))
          mat = dat$meta[[3]]
          m = colMeans(mat[I,1:3000], na.rm = T)
          N = mean(colMeans(mat[J,500:2000], na.rm = T))
          m = m / N      
          Ms[[s]] = m
        } else
          Ms[[s]] = rep(0,3000)
      }
      Mlen = length(Ms)
      df = data.frame( x = rep(-499:2500, Mlen), y = do.call(c,Ms), 
                       var = factor(rep(names(Ms), each = 3000)) )
      p = ggplot(df) 
      p = p + xlab("Position") + ylab(paste0(a,"-",a2))
      p = p + geom_line(aes(x=x,y=y,colour = var),size=1.5)
      labelat = seq(-500,4500, 500)
      labelval = labelat
      labelval[[2]] = "TSS"
      p = p+scale_x_continuous(name="Position", breaks = labelat, labels = labelval )
      p = p+ylim(0,NA)
      p = p +scale_color_manual(values = c("red","blue","black"))
      ggsave(filename = paste0(FigureDir,"/Long-gene-meta-",a,"-",a2,".pdf"),
             plot=p,
             width=7,height=4, unit="in")
    }
  }
}  

if( 0 ) {
  for( a in c("K18ac", "K4me3", "K36me3", "K56ac", "K79me3") )
    for( b in c("K18ac", "K4me3", "K36me3", "K56ac", "K79me3") ) 
      if( a != b ) {
        Ps = list()
        for( s in c("WT", "eaf3", "set2")) {
          X = cnucs[paste0(s,"-", a,"-Input"),]
          if( a < b ) {
            Y = cnucs[paste0(s,"-", a,"-",b),]
          } else
            Y = cnucs[paste0(s,"-", b,"-",a),]
          Ps[[s]] = DensityScatter(X,Y,threshold=.99,coordeq = F,
                                   xlabel=a,ylabel=paste(a,b), title=s,
                                   jitter = T, pointsize = 1)
        }
        png(paste0(FigureDir,"/conditional-",a,"-",b,".png"),height=600, width=1224)
        multiplot(Ps[[1]], Ps[[2]],Ps[[3]], layout = matrix(c(1,2,3),nr=1))
        dev.off()
      }
  
  for( a in c("K18ac", "K4me3", "K36me3", "K56ac", "K79me3") )
    for( b in c("K18ac", "K4me3", "K36me3", "K56ac", "K79me3") ) 
      if( a < b ) {
        Ps = list()
        for( s in c("WT", "eaf3", "set2")) {
          X = paste0(s,"-", a,"-Input")
          Y = paste0(s,"-", b,"-Input")
          
          if( a < b ) {
            Z = paste0(s,"-", a,"-",b)
          } else
            Z = paste0(s,"-", b,"-",a)
          Ps[[s]] = PlotXYZ(X,Y,Z)
        }
        png(paste0(FigureDir,"/pairwise-",a,"-",b,".png"),height=400, width=1224)
        multiplot(Ps[[1]], Ps[[2]],Ps[[3]], layout = matrix(c(1,2,3),nr=1))
        dev.off()
      }
}

if(0) {
  threshold = 0.95
  Mi = list()
  for( s in c("WT", "eaf3", "set2"))
    for( m in c("K18ac", "K4me3", "K36me3", "K56ac")) {
      ename = paste0(s,"-",m,"-Input")
      if( is.element(ename, rownames(nucs)) ) {
        N = nucs[ename,]
        Upper = quantile(N, threshold, na.rm = TRUE)
        Mi[[paste0(s,"-", m)]] = pmin(N/Upper,1)
      } else
        print(paste("Error SetupMarginal: cannot find experiment - ", ename))
    }
  Marginals = do.call(rbind,Mi)
}
abort()

if(0) {
  
  GeneCov = list()
  for( a in levels(Ab1))
  for( s in levels(Strain)) {
    i.s = which( Names == paste0(s,"-",a,"-Input"))[[1]]
    dat = readRDS(paste0(DataDir,"/Data/", Files[[i.s]]))
    if( s %in% names(GeneCov)) 
      GeneCov[[s]] = GeneCov[[s]] + countOverlaps(SGD$Genes, dat$UniqGR)
    else
      GeneCov[[s]] = countOverlaps(SGD$Genes, dat$UniqGR)
  }
  A = do.call(cbind, GeneCov)
}
if( 0 ) {
ccBuildMatrix <- function(x) {
  m = matrix(0,nr = length(levels(Ab1)), nc = length(levels(Ab2)))
  rownames(m) = levels(Ab1)
  colnames(m) = levels(Ab2)
  for( i in 1:length(x) )
    m[Ab1[i],Ab2[i]] = x[i]
  return(m)
}
}

if(0) {
focus <- GRanges(c("chrIII"), IRanges(31250,52750) ) 
#focus <- GRanges(c("chrIII"), IRanges(31250,32000) ) 
focus.nucids <- queryHits(findOverlaps(NucRegions, focus))
focus.nucs <- nucs[,focus.nucids]
focus.NucRegions <- NucRegions[focus.nucids]
}

if(0){

  for( a in levels(Ab1)) {
    i.wt = which( Names == paste0("WT-",a,"-Input"))
    for( s in levels(Strain) )
      if( s != "WT") {
        i.s = which( Names == paste0(s,"-",a,"-Input"))
        PlotComparison(Dat[[i.wt]],Dat[[i.s]], paste0("WT-",s,"-",a), order = "length", type="Plus1") 
        PlotComparison(Dat[[i.wt]],Dat[[i.s]], paste0("WT-",s,"-",a), type="Plus1",order = "occupancy-s") 
        PlotComparison(Dat[[i.wt]],Dat[[i.s]], paste0("WT-",s,"-",a), order = "length", type="Genes") 
        PlotComparison(Dat[[i.wt]],Dat[[i.s]], paste0("WT-",s,"-",a), order = "occupancy-s", type="Genes") 
      }
    }
}

if( 0 ) {
  for( s in levels(Strain) )
  for( a in levels(Ab1))
    for( b in levels(Ab1) ) 
      if( a < b )
        PlotPair(a,b, pre = s, relative = F)
  
  for( a in levels(Ab1)) {
    fn = paste0(DataDir, "/", a,"-input.png")
    png(fn,height=1024, width=600)
    p1 = DensityScatter( cnucs["Input",], 
                         cnucs[paste0(a,"-input"),],
                         xlabel = "Input", ylabel = paste0(a,"-input"), 
                         threshold = .99, coordeq = FALSE )
    p2 = DensityScatter( cnucs["Input",], 
                         cnucs[paste0(a,"-input"),],
                         xlabel = "Input", ylabel = paste0(a,"-input"), 
                         threshold = .99, coordeq = FALSE, logscale = TRUE )
    multiplot(p1,p2)
    dev.off()  
  }
}

if( 0 ) {
  for( a in levels(Ab1))
    for( b in levels(Ab1) ) 
      if( a < b) {
        PairExpectations(a,b)
        PairExpectationRatio(a,b)
      }
}

if(0)
  for( s in levels(Strain))
    for( a in levels(Ab1))
      for( b in levels(Ab1) ) 
          PairScatter(a,b,pre = s)

if(0)
for( a in levels(Ab1))
  for( b in levels(Ab1) ) 
    if( !(a == b) ){
      RelativeEnrichment(a,b)
    }

if( 0 ) 
  for( a in levels(Ab1))
    CheckBackground(a,plot=T)

if(0) {
  sacCer3 = Seqinfo(genome="sacCer3")
  for( d in dat) {
    seqinfo(d$cen) <- sacCer3
    fname = ccBuildFN(paste0("track",d$Name),params, suff = ".bw")
    export( d$cen, fname, format="BigWig" )
  }
  export( SGD$Transcripts, 
          ccBuildFN("trackTSS", params, suff = ".bed"), 
          format = "bigBed" )
}

if( 0 ) {
  for( d in dat ) 
    ccDoMeta(d, params, paste0(DataDir,"/Meta"))
}
## iChIP experiment specific stuff
if(0 ) {
  cdat = c()
  cdat[["mono-0.125"]] = ccCombineDat(dat[[1]],dat[[3]],"MNase-0-125-mono")
  cdat[["mono-0.25"]] = ccCombineDat(dat[[5]],dat[[7]],"MNase-0-25-mono")
  cdat[["mono-0.5"]] = ccCombineDat(dat[[9]],dat[[11]],"MNase-0-5-mono")
  cdat[["sub-0.125"]] = ccCombineDat(dat[[2]],dat[[4]],"MNase-0-125-sub")
  cdat[["sub-0.25"]] = ccCombineDat(dat[[6]],dat[[8]],"MNase-0-25-sub")
  cdat[["sub-0.5"]] = ccCombineDat(dat[[10]],dat[[12]],"MNase-0-5-sub")
}
if( 0 )
  Files = c( "MNase-0-125-mono",
             "MNase-0-25-mono", 
             "MNase-0-5-mono", 
             "MNase-0-125-sub",
             "MNase-0-25-sub",
             "MNase-0-5-sub" )
  cdat = lapply(Files, function(x) readRDS(paste0(DataDir,"/",x,".rdata")))
  names(cdat) =   c( "mono-0.125",
                     "mono-0.25",
                     "mono-0.5",
                     "sub-0.125",
                     "sub-0.25",
                     "sub-0.5" )
  
                     
  for( t in c("TSS", "TTS", "Plus1", "PostTTS"))
    for( o in c("expr", "length", "occupancy-m", "occupancy-s")) {
      for( m in c( "0.125", "0.25", "0.5"))
      PlotComparison( cdat[[paste0("mono-",m)]], cdat[[paste0("sub-",m)]], 
                      paste0("Mnase-", m), type=t, order = o,filetype = "pdf")
    }


if( 0 ) {
  RegionsTSS = list( "TSS-150" = 300:400, "TSS" = 450:550, "TSS+150" = 600:700)
  RegionsTTS = list( "TTS-50" = 900:1000, "TTS+50" = 1000:1100, "TTS+150" = 1100:1200)
  m = data.frame(ORF = SGD$Genes$acc, Name = SGD$Genes$name)
  for( f in names(cdat)) {
    o = CollectOccupancy(cdat[[f]], RegionsTSS, RegionsTTS )
    colnames(o) = unlist(lapply(colnames(o), function(x) paste(f,x)))
    m = cbind(m,o)
  }
  write.csv(m, paste0(DataDir,"/occupancy.csv"), row.names = F)
}

## old
abort()


Files = c( "K4me3-Input.bam",
           "K4me3-K4me3.bam", 
           "K4me3-K36me3.bam",
           "K36me3-Input.bam",
           "K36me3-K4me3.bam", 
           "K36me3-K36me3.bam" )

# meta design
#
# parameters - s4, b4, s36, b36, s4**2, , b4**2, s36**2, b36**2, s4s36, s4b36, b4s36, b4b36

basepar = c( s4 = .5, b4 = 0.001, s36 = .2, b36 = 0.005)
ExpandPar <- function( basepar ) {
  s4 = basepar[["s4"]]
  b4 = basepar[["b4"]]
  s36 = basepar[["s36"]]
  b36 = basepar[["b36"]]
  
  c( s4, b4, s36, b36, 
     s4**2, b4**2, s36**2, b36**2,
     s4*s36, s4*b36, b4*s36, b4*s36 )
}

par = ExpandPar( basepar)

# unknowns - +4+36, +4-36, -4+36, -4-36 
#
UnknownsNames = c("K4me3 K36me3", "K4me3 alone", "K36me3 alone", "Neither")

D.4.Input = matrix( c( 1, 1, 0, 0,
                     0, 0, 1, 1,
                     0, 0, 0, 0,
                     0, 0, 0, 0,
                     0, 0, 0, 0,
                     0, 0, 0, 0,
                     0, 0, 0, 0,
                     0, 0, 0, 0,
                     0, 0, 0, 0,
                     0, 0, 0, 0,
                     0, 0, 0, 0,
                     0, 0, 0, 0 ),
                  ncol=4, byrow = TRUE)
D.4.4 = matrix( c(  0, 0, 0, 0,
                    0, 0, 0, 0,
                    0, 0, 0, 0,
                    0, 0, 0, 0,
                    1, 1, 0, 0,
                    0, 0, 1, 1,
                    0, 0, 0, 0,
                    0, 0, 0, 0,
                    0, 0, 0, 0,
                    0, 0, 0, 0,
                    0, 0, 0, 0,
                    0, 0, 0, 0 ),
                  ncol=4, byrow = TRUE)
D.36.Input = matrix( c( 0, 0, 0, 0,
                      0, 0, 0, 0,
                      1, 0, 1, 0,
                      0, 1, 0, 1,
                      0, 0, 0, 0,
                      0, 0, 0, 0,
                      0, 0, 0, 0,
                      0, 0, 0, 0,
                      0, 0, 0, 0,
                      0, 0, 0, 0,
                      0, 0, 0, 0,
                      0, 0, 0, 0),
                  ncol=4, byrow = TRUE)
D.36.36 = matrix( c( 0, 0, 0, 0,
                     0, 0, 0, 0,
                     0, 0, 0, 0,
                     0, 0, 0, 0,
                     0, 0, 0, 0,
                     0, 0, 0, 0,
                     1, 0, 1, 0,
                     0, 1, 0, 1,
                     0, 0, 0, 0,
                     0, 0, 0, 0,
                     0, 0, 0, 0,
                     0, 0, 0, 0),
                   ncol=4, byrow = TRUE)
D.4.36 = matrix( c( 0, 0, 0, 0,
                  0, 0, 0, 0,
                  0, 0, 0, 0,
                  0, 0, 0, 0,
                  0, 0, 0, 0,
                  0, 0, 0, 0,
                  0, 0, 0, 0,
                  0, 0, 0, 0,
                  1, 0, 0, 0,
                  0, 1, 0, 0,
                  0, 0, 1, 0,
                  0, 0, 0, 1 ),
                  ncol=4, byrow = TRUE)
D.36.4 = D.4.36

multiDesign = list( D.4.Input, D.4.4, D.4.36, D.36.Input, D.36.4, D.36.36)

TransProb = matrix( c( 0.5, 0.0, 0.5, 0.0,  #both
               0.3, 0.7, 0.0, 0.0,  #K4
               0.0, 0.0, 0.8, 0.2,  #K36
               0.0, 0.7, 0.0, 0.3  # none
), ncol = 4, nrow = 4, byrow = TRUE )

set.seed(0)
simNus <- c(100000,100000,10000,100000,50000,20000)
sim <- ccSimulateCounts(TransProb, 1000, par, multiDesign, simNus)
simCounts = sim$Counts
simXs = sim$Xs

GLMFitOut = ccGLMFit(simCounts, multiDesign, initialPar = par, iter = 5)
fitXs = GLMFitOut[[2]]
DensityScatter(simXs,fitXs,xlabel="source", ylabel = "estimated",threshold=1,alpha = .8)
for( i in 1:4 ) DensityScatter(simXs[,i],fitXs[,i],
                               xlabel="source", ylabel = "estimated", 
                               title = UnknownsNames[i], alpha = .5)

if( 0 ) {
countsFile = ccBuildFN("TileCounts", params )
if( file.exists( countsFile )) {
  print("Loading tile counts from file")
  load(countsFile)
} else {
  dat = lapply(Files, function(x) ccProcessFile(paste0(BamDir,"/",x), param = params))
  print("Computing tile counts")
  Tiles = tileGenome(seqinfo(dat[[1]]$GR),tilewidth = 100,cut.last.tile.in.chrom = TRUE)
  tempcounts = lapply(dat, function(d) countOverlaps(Tiles,d$GR))
  counts = do.call(rbind,tempcounts)
  rownames(counts) = lapply(Files, BaseFileName)
  save(Tiles,counts,file=countsFile)
}
  
#counts = counts[,1:10000]
XsFile = ccBuildFN("GLMFit", params)
if( file.exists( XsFile) ) {
  print("Loading regressed estimates")
  load( file = XsFile )
} else {
  print("Computing regressed estimates")
  GLMFitOut = ccGLMFit(counts, multiDesign, initialPar = par)
  save(GLMFitOut, file=XsFile )
}

for( i in 1:length(UnknownsNames)) {
  temp = Tiles
  score(temp) = GLMFitOut[[2]][,i]
  fname = ccBuildFN(paste0("track",i),params, suff = ".bw")
  export( temp, fname, format="BigWig" )
}

}

if( 0 ) {
  # check gene coverage to detect which KO...
  
  KORegions = SGD$Genes
  KORegions = resize(KORegions, width = width(KORegions) - 400, fix = "center")

  countCov <- function(x) {
    g = unique(x$GR)
    g = g[width(g) < 200]
    countOverlaps(KORegions, g)
  }

  ko.cov = do.call(rbind,lapply( dat, countCov ))

  rownames(ko.cov) = Names
  colnames(ko.cov) = KORegions$name
}
