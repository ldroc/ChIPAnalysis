setwd("~/Dropbox/CoChIP/coChIP-Analysis")
source("DensityScatter.R")
source("Multiplot.R")
options(digits = 5)

Strains = c("WT", "set2", "eaf3")
Abs = c("H3", "H3K4me3", "H3K9ac", "H3K14ac", "H3K36me3", "H3K56ac", "H3K79me3")

DataDir = "~/Google Drive/CoChIPAnalysis/160331"
NucFile = "~/Google Drive/CoChIPAnalysis/160331/nuc-merge.rdata"

if( file.exists(NucFile) ) {
  NucMerge = readRDS( NucFile )
} else {
  Nuc = readRDS("~/Google Drive/CoChIPAnalysis/160331/nuc-occ.rdata")

  sNuc = list()
  for( s in Strains ) 
    for( a in Abs )
      for( b in c(Abs, "Flag") )
        sNuc[[paste0(s,"-", a, "-", b)]] = Nuc[paste0(s,"-1-", a, "-", b),] + Nuc[paste0(s,"-2-", a, "-", b),]+Nuc[paste0(s,"-3-", a, "-", b),]

  NucMerge = do.call(rbind,sNuc)
  
  NucFlag = Nuc[grep("Flag-",rownames(Nuc)),]
  rownames(NucFlag) = sub("-1-", "-", rownames(NucFlag))
  NucMerge = rbind(NucMerge, NucFlag)
  
  NucInput = readRDS("~/Google Drive/CoChIPAnalysis/160331/nuc-occ-Input.rdata")
  iNuc = list()
  for( s in Strains ) 
    for( a in Abs )
        iNuc[[paste0(s,"-", a, "-Input")]] = NucInput[paste0(s,".1-", a, "-Input"),] + NucInput[paste0(s,".2-", a, "-Input"),]+NucInput[paste0(s,".3-", a, "-Input"),]
  NucMerge = rbind( NucMerge, do.call(rbind,iNuc))
  
  saveRDS(NucMerge, NucFile )
}

FigureDir = "~/Google Drive/CoChIPAnalysis/160331/Reciprocal"

if( 0 ) {
  Counts = rowSums(NucMerge)
  hist(Counts[grep("Input",names(Counts))],breaks="FD")
  plot(sort(Counts), 1:length(Counts))
  sum(Counts < 100000)
  
  WT.Counts = Counts[grep("WT",names(Counts))]
  I = WT.Counts > 100000
  R = cor( t(NucMerge[names(WT.Counts)[I],]))
  png(paste0(DataDir,"/WTcor.png"),width=1024,height=1024)
  corrplot::corrplot(R,order="hclust")
  dev.off()
}

s = "WT"
a = "H3K9ac"
b = "H3K36me3"

Plot2D <- function( X, Y, diagonal = TRUE) {
  threshold = .99
  t1 = quantile(X,threshold, na.rm=TRUE)
  t2 = quantile(Y,threshold, na.rm=TRUE)
#  I = X > 0.025 & Y > 0.025 & X < t1 & Y < t2
  I = X < t1 & Y < t2
  p = ggplot(data.frame(x  = X[I], y= Y[I]),aes(x=x,y=y))
  p = p+geom_bin2d(binwidth=1)
  #    p = p+scale_fill_gradient(low="black",high="red")
  p = p+scale_fill_distiller(palette = "YlOrBr", direction=1, values=c(-0.1,.05,0.15,1))
  if( diagonal )
    p = p+geom_abline(slope=1,intercept=0,color="blue",size=2)
#  p = p + coord_equal(ratio=1)
  p = p+theme_classic()
  r = cor( X[I], Y[I])
  s = lm(Y~X-1)
  p = p + ggtitle(paste("r = ",sprintf("%.3f",r), "slope = ", sprintf("%.3f",coef(s))))
#  print(c(r,coef(s)))
  #  p = p + scale_x_continuous(expand=c(0,0),limits=c(0,)) 
#  p = p + scale_y_continuous(expand=c(0,0),limits=c(0,1.1)) 
  p
}

AbSlope = c( H3  = 0.069508, 
             H3K4me3 = 0.025615,
             H3K9ac = 0.24587,
             H3K14ac = 0.24875,
             H3K36me3 = 0.029901,
             H3K56ac = 0.091796, 
             H3K79me3 = 0.011131
)

if(0)
for( s in Strains )
  for( a in Abs )
    for( b in Abs )
      if( a < b ) {
        ab = paste0(s,"-",a,"-",b)
        ba = paste0(s,"-",b,"-",a)
        aa = paste0(s,"-",a,"-",a)
        bb = paste0(s,"-",b,"-",b)
        aF = paste0(s,"-",a,"-Flag")
        bF = paste0(s,"-",b,"-Flag")
        
        p1 = Plot2D(NucMerge[ab,],NucMerge[ba,]) + xlab(ab)+ylab(ba)
        p2 = Plot2D(NucMerge[ab,]-NucMerge[aF,],NucMerge[ba,]-NucMerge[bF,]) +
          xlab(paste(ab,"-",aF))+ylab(paste(ba,"-",bF))
        p3 = Plot2D(pmax(NucMerge[ab,]+NucMerge[bF,],0),
                       pmax(NucMerge[ba,]+NucMerge[aF,],0)) + 
          xlab(paste(ab,"+",bF))+ylab(paste(ba,"+",aF))
#        p4 = Plot2D(pmax(NucMerge[ab,]-AbSlope[a]*NucMerge[aa,],0),
#                                  pmax(NucMerge[ba,]-AbSlope[b]*NucMerge[bb,],0)) + xlab(paste(ab,"fix"))+ylab(paste(ba,"fix"))
        p4 = Plot2D(pmax(NucMerge[ab,]+NucMerge[aF,],0),
                    pmax(NucMerge[ba,]+NucMerge[bF,],0)) + 
          xlab(paste(ab,"+",aF))+ylab(paste(ba,"+",bF))
        png(paste0(FigureDir,"/",ab,"-fix.png"), width=600,height=950)
        multiplot(p1,p2,p3, p4,cols = 2)
        dev.off()
      }

if(0)
for( s in Strains )
  for( a in Abs ) {
    print(a)
    aF = paste0(s,"-",a,"-Flag")
    aa = paste0(s,"-",a,"-",a)
    ggsave(filename = paste0(FigureDir,"/",aF,".png"), Plot2D(NucMerge[aa,],NucMerge[aF,]) + xlab(aa)+ylab(aF) )
    
  }

{
  aF = paste0(s,"-",a,"-Flag")
  aI = paste0(s,"-",a,"-Input")
  Plot2D(NucMerge[aF,],NucMerge[aI,]) + xlab(aF)+ylab(aI)
  
}

LLfit <- function( p ) {
  L = 0
  s = "WT"
  for( a in Abs )
    for(b in Abs )
      if( a < b ) {
        ab = paste0(s,"-",a,"-",b)
        ba = paste0(s,"-",b,"-",a)
        aa = paste0(s,"-",a,"-",a)
        bb = paste0(s,"-",b,"-",b)
        aF = paste0(s,"-",a,"-Flag")
        bF = paste0(s,"-",b,"-Flag")
        wa = p[[a]]
        wb = p[[b]]
        xa = (wb*NucMerge[ab,] - NucMerge[aF,])
        xb = (wa*NucMerge[ba,] - NucMerge[bF,])
        I = xa >= 0 & xb>=0
        J = xa < 0 | xb < 0
        L = L + sum((xa[I]-xb[I])**2) +sum(xa[J]**2) + sum(xb[J]**2)
      }
  return(L)
}

if(0) {
  p0 = rep(1/3,length(Abs))
  names(p0) = Abs
  
  opt = optim(par = p0, fn = LLfit, lower = 0)
  p = opt$par
}


if(1) {
  # Build regression matrix
  
  s = "WT"
  RNames = rownames(NucMerge)
  RNames = RNames[grep(s,RNames)]
  RNames = RNames[grep("Flag-",RNames, invert = TRUE)]
  RNames = RNames[grep("Input",RNames, invert = TRUE)]
  mm = strsplit(RNames,"-")
  n = 1
  CNames = list()
  for( a in Abs )
  {
    CNames[[paste0(s,"-",a,"-Flag")]] = n
    n = n+1
  }    
  for( a in Abs )
    for( b in Abs)
      if( a <= b) {
        CNames[[paste0(s,"-",a,"-",b)]] = n
        n = n+1
      }
  
  Ws0 = c(rep(1, length(Abs)),1/3)
  names(Ws0) = c(Abs, "Flag")
  Ws = Ws0
  
  BuildsXs <- function( Ws ) {
    Xs = matrix(0,nr=length(RNames),nc=length(CNames) )
    rownames(Xs) = RNames
    colnames(Xs) = names(CNames)
    wF = Ws[["Flag"]]
    for( i in 1:length(RNames) ) {
      s = mm[[i]][[1]]
      a = mm[[i]][[2]]
      b = mm[[i]][[3]]
      
      if( b == "Flag") {
        Xs[i,CNames[[RNames[[i]]]]] = wF
      } else
      {
        wb = Ws[[b]]
        Xs[i,CNames[[paste0(s,"-",a,"-Flag")]]] = wb
        if( any(names(CNames) == paste0(s,"-",a,"-",b)) )
          Xs[i,CNames[[paste0(s,"-",a,"-",b)]]] = wb
        if( any(names(CNames) == paste0(s,"-",b,"-",a)) )
          Xs[i,CNames[[paste0(s,"-",b,"-",a)]]] = wb
      }
    }
    Xs
  }    
  
  BuildCNucs <- function(Ws, N = 5000) {

    Xs = BuildsXs(Ws)
    
    NucPoissonLL <- function(p) {
      lam = Xs %*% p
      -sum(Y*log(lam) - lam)
    }
    
    dNucPoissonLL <- function(p) {
      lam = Xs %*% p
      d = Y/lam - 1
      -colSums(Xs * as.vector(d))
    }
    
    CNucs = matrix(0,nr=length(CNames), nc=dim(NucMerge)[[2]])
    rownames(CNucs) = names(CNames)
    colnames(CNucs) = colnames(NucMerge)
    p = rep(1,length(CNames))
    
    Value = 0
    for( i in 1:N ){ 
      Y = NucMerge[RNames,i]
      opt = optim(par=p,fn=NucPoissonLL,gr=dNucPoissonLL,method="L-BFGS-B", lower=10^-8)
      CNucs[,i] = opt$par
      Value = Value + opt$value
    }
    print(paste("OptNucs:",Value))
    
    CNucs
  }
  
  BuildWs <- function(CNucs, N = 5000) {
    Zs = matrix(0, nr = N*length(RNames), nc = length(Ws) )
    colnames(Zs) = names(Ws)

    Os = seq(0,N*length(RNames)-1, by=length(RNames))
    
    for( i in 1:length(RNames) ) {
      s = mm[[i]][[1]]
      a = mm[[i]][[2]]
      b = mm[[i]][[3]]
      
      Zs[Os + i,b] = CNucs[paste0(s,"-", a, "-Flag"),1:N]
      if( b != "Flag") {
        if( any(names(CNames) == paste0(s,"-",a,"-",b)) )
          Zs[Os+i,b] = Zs[Os+i,b]+ CNucs[paste0(s,"-",a,"-",b),1:N]
        if( any(names(CNames) == paste0(s,"-",b,"-",a)) )
          Zs[Os+i,b] = Zs[Os+i,b]+CNucs[paste0(s,"-",b,"-",a),1:N]
      }
    }
    
    YY = as.vector(NucMerge[RNames,1:N])
    
    WsPoissonLL <- function(p) {
      lam = Zs %*% p 
      -sum(YY*log(lam) - lam)
    }
    
    dWsPoissonLL <- function(p) {
      lam = Zs %*% p 
      d = YY/lam - 1
      -colSums(Zs * as.vector(d))
    }
    print(paste("StartWs: ",WsPoissonLL(Ws)))
   opt = optim(par=Ws0,fn=WsPoissonLL,gr=dWsPoissonLL,method="L-BFGS-B", lower=10^-8)
    Ws = opt$par
    Ws = Ws * length(Ws) / sum(Ws)
    print(paste("OptWS:",opt$value))
    Ws
  }
  
  N = 5000
  if(1) {
    CNucs = BuildCNucs(Ws0, N)
    Plots = list()
    for( i in 1:1 ) {
      print(i)
      Ws = BuildWs(CNucs,N) 
      print(Ws)
      CNucs = BuildCNucs(Ws,N)
      
      for( t in c( "WT-H3K79me3-H3K4me3", "WT-H3K4me3-H3K79me3", "WT-H3K36me3-H3K14ac", "WT-H3K14ac-H3K36me3") )
        Plots[[t]] = c(Plots[[t]],list(DensityScatter(NucMerge[t, 1:N], Xs[t,] %*% CNucs[, 1:N], diagonal = TRUE) + labs(title=paste("Iter = ", i,t))))
    }
  }
  CNucs = BuildCNucs(Ws,dim(NucMerge)[[2]])
  
  NucFix = NucMerge * 0
  for( i in 1:length(RNames) ) {
    s = mm[[i]][[1]]
    a = mm[[i]][[2]]
    b = mm[[i]][[3]]
    NucFix[RNames[[i]],] = NucMerge[RNames[[i]],] - Ws[[b]]*CNucs[paste0(s,"-",a,"-Flag"),] 
  }
  for( s in Strains )
    for( a in Abs )
      for( b in Abs )
        if( a < b ) {
          ab = paste0(s,"-",a,"-",b)
          ba = paste0(s,"-",b,"-",a)

          p1 = Plot2D(NucMerge[ab,],NucMerge[ba,]) + xlab(ab)+ylab(ba)
          p2 = Plot2D(NucFix[ab,],NucFix[ba,]) +
            xlab(paste(ab," Fix"))+ylab(paste(ba," Fix"))
          png(paste0(FigureDir,"/",ab,"-fix.png"), width=600,height=450)
          multiplot(p1,p2,cols = 2)
          dev.off()
        }
}