if( !exists("InputName") ) 
  InputName = "Input"
if( !exists("cNucRegions") )
  cNucRegions = NucRegions

source("Multiplot.R")

BuildExpName <- function( x, y = InputName, pre = "") {
  if( pre != "" )
    paste0(pre,"-",x,"-",y)
  else
    paste0(x,"-",y)
}

PlotXYZ <- function(xl,yl,zl, alpha = 1, threshold = .995, 
#                    nucs = nucs, 
                    logscale = FALSE, relative = FALSE ) {
  x = nucs[xl,]
  y = nucs[yl,]
  z = nucs[zl,]
  if( relative ) {
    x = x/Inputs
    y = y/Inputs
    z = z/Inputs
  }
  if( relative ) {
    xl = paste0(xl,"/input")
    yl = paste0(yl,"/input") 
    zl = paste0(zl,"/input") 
  }
  
  return(Plot3D( x,y,z, xl, yl, zl, alpha, threshold, logscale ))
}

PlotPair <- function(x,y,pre="", relative = TRUE) {
  if( x != y ) {
    nx = BuildExpName(x,pre=pre)
    ny = BuildExpName(y,pre=pre)
    nxy = BuildExpName(x,y,pre=pre)
    print(c(x,y,nx,ny))
    if( length(which(nucs[nx,] > 0)) < 100 ||
        length(which(nucs[ny,] > 0)) < 100
        )
      return(FALSE)
    fn = paste0(DataDir, "/", BuildExpName(x,y,pre=pre), ".png")
    png(fn,height=1024, width=600)
    p1 = (PlotXYZ(nx, ny, nxy, relative=relative))
    p2 = (PlotXYZ(nx, ny, nxy, logscale = TRUE,relative=relative))
    multiplot(p1,p2)
    dev.off()  
  }
}


ComputeExpectations <- function(x, y,pre="") {
  o = cnucs[BuildExpName(x,pre=pre),] * cnucs[BuildExpName(y,pre=pre),] / Inputs
}

PairExpectations <- function(x,y,pre="") {
  #  Xs = ccPairIPInit(BuildPairwiseMatrix(x,y), Inputs) $ Xs
  Exp = ComputeExpectations(x,y)
  fn = paste0(DataDir, "/Obs-Exp-", BuildExpName(x,y,pre),".png")
  png(fn,height=1024, width=600)
  p1 = DensityScatter( Exp, cnucs[BuildExpName(x,y,pre),], 
                       xlabel = "Expected", ylabel = "observed", title = BuildExpName(x,y,pre),
                       coordeq = FALSE, alpha = .1, smooth = T,logscale = F)
  
  p2 = DensityScatter( Exp/Inputs, cnucs[BuildExpName(x,y,pre),]/Inputs, 
                       xlabel = "Expected/input", ylabel = "observed/input", title = paste0(x,"-",y, "/input"),
                       coordeq = FALSE, alpha = .1, smooth = T,logscale = F )
  
  
  multiplot(p1,p2)
  dev.off()  
}

PairScatter <- function(x,y,pre="") {
  xn = BuildExpName(x,pre=pre)
  if( x < y)
    xyn = BuildExpName(x,y,pre=pre)
  else
    xyn = BuildExpName(y,x,pre=pre)
  fn = paste0(DataDir, "/Scatter-", BuildExpName(x,y,pre),".png")
  png(fn,height=1024, width=600)
  p1 = DensityScatter( cnucs[xn,], cnucs[xyn,], 
                       xlabel = x, ylabel = BuildExpName(x,y,pre),
                       coordeq = FALSE, alpha = .1, logscale = F)
  
if(0 ) {
   p2 = DensityScatter( cnucs[xn,]/Inputs, cnucs[xyn,]/Inputs, 
                       xlabel = paste0(x, "/input"), ylabel = paste0(x,"-",y,"/input"),
                       coordeq = FALSE, alpha = .1, logscale = F )
} else
  p2 = p1 + scale_x_log10()+scale_y_log10()  
  
  multiplot(p1,p2)

  
  dev.off()  
}

PairExpectationRatio <- function(x,y,pre="") {
  #  Xs = ccPairIPInit(BuildPairwiseMatrix(x,y), Inputs) $ Xs
  fn = paste0(DataDir, "/Obs-Exp-Ratio-", BuildExpName(x,y,pre),".png")
  png(fn,height=1024, width=600)
  Obs = cnucs[BuildExpName(x,y,pre),]
  #  Exp = cnucs[BuildExpName(x,pre=pre),]* Xs[,1]/(Xs[,1]+Xs[,2])
  Exp = ComputeExpectations(x,y)
  ratio = log( Obs/Exp, base = 2 )
  ratio = ratio - quantile(ratio, .75, na.rm=T)
  xl = BuildExpName(x,pre=pre)
  yl = BuildExpName(y,pre=pre)
  et = quantile(Exp,.99,na.rm=T)
  rt = quantile(ratio,.99,na.rm=T)
  I = !is.na(ratio) & !is.infinite(ratio) & Exp < et & ratio < rt
  p1 = Plot3D( cnucs[xl,I], cnucs[yl,I], ratio[I], 
               xl = x, yl = y, zl = "observed/expected (log_2)",
               alpha = .1, logscale = T)
  
  p2 = Plot3D( cnucs[xl,I]/Inputs[I], cnucs[yl,I]/Inputs[I], ratio[I], 
               xl = paste0(x,"/input"), yl = paste0(y,"/input"), zl = "observed/expected (log_2)",
               alpha = .1, logscale = F)
  multiplot(p1,p2)
  dev.off()  
}

CheckBackground <- function(a, plot=F, H3 = "H3", pre="") {
  x = cnucs[BuildExpName(H3,pre=pre),]/Inputs
  t = quantile(x,0.15, na.rm=TRUE)
  y = cnucs[BuildExpName(a,pre=pre),]
  I = x<t
  bg = median(y[I]/Inputs[I],na.rm = T)
  if( plot ) {
    fn = paste0(DataDir, "/Background-",a,".png")
    png(fn,height=1024, width=1000)
    xlim = c(quantile(Inputs,0.00, na.rm=TRUE), quantile(Inputs,0.99, na.rm=TRUE))
    ylim = c(quantile(y,0.00, na.rm=TRUE), quantile(y,0.95, na.rm=TRUE))
    yilim = c(quantile(y/Inputs,0.01, na.rm=TRUE), quantile(y/Inputs,0.99, na.rm=TRUE))
    
    #    p1 = DensityScatter(Inputs[I], y[I]/Inputs[I], xlabel="input", ylabel= paste(a,"/input"), coordeq = F, alpha=.1,  
    #                 horizontal = bg, xlim=xlim, ylim=yilim )
    p1 = ggplot(data = data.frame( x = y/Inputs), aes(x = x)) + geom_histogram()+xlim(ylim)+
      geom_histogram( data = data.frame( x= y[I]/Inputs[I]), color = "blue", fill = "blue") +
      geom_vline(xintercept = bg, color="orange", size = 2 ) + labs(title=paste0(a,"/input"))
    p2 = DensityScatter(Inputs, y/Inputs, xlabel="input", paste(a,"/input"), ylabel= paste(a,"/input"),
                        coordeq = F, alpha=.1, 
                        horizontal = bg, xlim=xlim, ylim=yilim ) + geom_hline(yintercept = bg, color="orange", size=2)
    #    p3 = DensityScatter(Inputs[I], y[I], xlabel="input", ylabel= a, coordeq = F, alpha=.1,  
    #                        xlim=xlim, ylim=ylim ) + geom_abline(yintercept = 0, slope=bg, color="orange", size=2)
    p3 = ggplot(data = data.frame( x = y), aes(x = x)) + geom_histogram()+xlim(ylim)+
      geom_histogram( data = data.frame( x= y[I]), color = "blue", fill = "blue") + labs(title=a)
    p4 = DensityScatter(Inputs, y, xlabel="input", ylabel= a, coordeq = F, alpha=.1, 
                        xlim=xlim, ylim=ylim ) + geom_abline(yintercept = 0, slope=bg, color="orange", size=2)
    multiplot(p1,p2,p3,p4,layout = matrix(c(1,2,3,4), nrow=2))
    dev.off()     
  }
  return(bg)
}

RelativeEnrichmentPlot <- function( a, b, nucs = cnucs, logscale = F, pre="" ) {
  x = BuildExpName(a, pre=pre) 
  y = BuildExpName(b, pre=pre) 
  z = BuildExpName(a,b, pre=pre) 
  print(c(x,y,z))
  X = nucs[x,]
  Y = nucs[y,]
  Z = nucs[z,]  
  p = DensityScatter(Z[I]/X[I], Z[I]/Y[I], 
                 xlabel=paste(z,"/", a), ylabel=paste(z,"/", b), coordeq = F,alpha=.1,
                 logscale = logscale, threshold = 0.99)
  return(p)
}

RelativeEnrichment <- function( a, b, nucs = cnucs, prefix = "Rel-Scatter-", pre = ""  ) {
  x = BuildExpName(a, pre=pre) 
  y = BuildExpName(b, pre=pre) 
  z = BuildExpName(a,b, pre=pre) 
  X = nucs[x,]
  Y = nucs[y,]
  Z = nucs[z,]  
  tx = quantile(X,.20, na.rm = TRUE)
  ty = quantile(Y,.20, na.rm = TRUE)
  print(c(tx,ty))
  
  if( tx == 0 || ty == 0 )
    return(FALSE)
  
  I = X > 0
  
  fn = paste0(DataDir, "/", prefix, z,".png")
  png(fn,height=650, width=1224)
  
  xTs = lapply( seq(0.1,1,.1), function(t) quantile(X,t, na.rm = TRUE))
  xGroups = lapply( 1:(length(xTs) - 1), function(i) X[I]> xTs[i] & X[I] <= xTs[i+1] )
  yTs = lapply( seq(0.1,1,.1), function(t) quantile(Y,t, na.rm = TRUE))
  yGroups = lapply( 1:(length(yTs) - 1), function(i) Y[I] > yTs[i] & Y[I] <= yTs[i+1] )
  
  nucGroups = lapply(list(n1 = 1, n3= 3, n5 = 5), function(i) (mcols(cNucRegions)$gene_pos == i ) )
  J = I & Reduce( "|", nucGroups)
  nucGroups <- lapply(nucGroups, function(x) x[J])
  
  rX = 100*rank(X[I])/length(which(!is.na(X[I])))
  rY = 100*rank(Y[I])/length(which(!is.na(Y[I])))
  
  p1 = DensityScatter(Z[I]/X[I], Z[I]/Y[I], 
                      xlabel=paste(z,"/", a), ylabel=paste(z,"/", b), coordeq = F,alpha=.1,
                      logscale = F, threshold = 0.99)
  p2 = DensityScatter(Z[I]/X[I], Z[I]/Y[I], 
                      xlabel=paste(z,"/", a), ylabel=paste(z,"/", b), coordeq = F,alpha=.2,
                      Value = rX,
                      logscale = F, threshold = 0.99) + scale_colour_gradient(a,low = "darkgray",high="red")
  p3 = DensityScatter(Z[I]/X[I], Z[I]/Y[I], 
                      xlabel=paste(z,"/", a), ylabel=paste(z,"/", b), coordeq = F,alpha=.2,
                      Value = rY,
                      logscale = F, threshold = 0.99) +    scale_colour_gradient(b,low = "darkgray",high="red")
  p4 = DensityScatter(Z[J]/X[J], Z[J]/Y[J], 
                      xlabel=paste(z,"/", a), ylabel=paste(z,"/", b), coordeq = F,alpha=.2,
                      Groups = nucGroups, galpha = 1,
                      logscale = F, threshold = 0.99) +
    scale_colour_discrete("nuc", na.value = "light gray")
  
  #  p2 = DensityScatter(X,Z/X, 
  #                      xlabel=a, ylabel=paste(z,"/", a), coordeq = F,alpha=.1,
  #                      logscale = F, threshold = 0.99)
  #  p3 = DensityScatter(Y,Z/Y, 
  #                      xlabel=b, ylabel=paste(z,"/", b), coordeq = F,alpha=.1,
  #                      logscale = F, threshold = 0.99)
  #
  p5 = DensityScatter(Z[I]/X[I], Z[I]/Y[I], 
                      xlabel=paste(z,"/", a), ylabel=paste(z,"/", b), coordeq = F,alpha=.1,
                      logscale = T, threshold = 0.99)
  p6 = DensityScatter(Z[I]/X[I], Z[I]/Y[I], 
                      xlabel=paste(z,"/", a), ylabel=paste(z,"/", b), coordeq = F,alpha=.1,
                      Value = rX, galpha = .2,
                      logscale = T, threshold = 0.99) +scale_colour_gradient(a,low = "darkgray",high="red")
  p7 = DensityScatter(Z[I]/X[I], Z[I]/Y[I], 
                      xlabel=paste(z,"/", a), ylabel=paste(z,"/", b), coordeq = F,alpha=.1,
                      Value = rY, galpha = .2,
                      logscale = T, threshold = 0.99)  +scale_colour_gradient(b,low = "darkgray",high="red")
  
  p8 = DensityScatter(Z[J]/X[J], Z[J]/Y[J], 
                      xlabel=paste(z,"/", a), ylabel=paste(z,"/", b), coordeq = F,alpha=.1,
                      Groups = nucGroups,galpha = 1,
                      logscale = T, threshold = 0.99) + 
    scale_colour_discrete("nuc", na.value = "light gray")
  
  multiplot(p1,p2,p3,p4,p5,p6,p7,p8,layout = matrix(c(1,2,3,4,5,6,7,8), nrow=2, byrow=T))
  dev.off()  
}

BuildPairwiseMatrix <- function(x,y, pre="") {
  m = matrix( nc = dim(cnucs)[[2]], nr = 5 )
  m[1,] = cnucs[BuildExpName(x,pre=pre),]
  m[2,] = cnucs[BuildExpName(y,pre=pre),]
  m[3,] = cnucs[BuildExpName(x,x,pre),]
  m[4,] = cnucs[BuildExpName(y,y,pre),]
  if( x < y ) {
    m[5,] = cnucs[BuildExpName(x,y,pre),]
  } else
    m[5,] = cnucs[BuildExpName(y,x,pre),]
  return(m)
}


SimulatePair <- function( a, b, Mode = "exp", 
                          specA = 1, bgA = 0.01,
                          specB = 1, bgB = 0.01, pre= "") {
  mm = BuildPairwiseMatrix(a, b, pre)
  Xs = ccPairIPInit(mm, Inputs)$Xs
  
  pA = Xs[,1]+Xs[,2]
  pB = Xs[,1]+Xs[,3]
  pT = Xs[,1]+Xs[,2] + Xs[,3] + Xs[,4]
  
  if( Mode == "exp")
  {
    pAB = Xs[,1]
  } else {
    if( Mode == "min")
    {
      pAB = pmax(pA+pB - pT, 0)  
    } else {
    if( Mode == "max" ) {
      pAB = pmin(pA, pB)
    } else
      return(NULL)
    }
  }
  pAb = pA - pAB
  paB = pB - pAB
  pab = pT - (pAB + pAb + paB)

  ps = matrix( nr = 3, nc = length(Inputs))
  ps[1,] = specA * (pAB + pAb) + bgA * (paB + pab)
  ps[2,] = specB * (pAB + paB) + bgB * (pAb + pab)
  ps[3,] = specA * specB * pAB + specA * bgB * pAb + bgA * specB * paB + bgA * bgB * pab
  
  ps[1,] = (sum( mm[1,], na.rm=TRUE )/ sum(ps[1,],na.rm=TRUE)) * ps[1,] 
  ps[2,] = (sum( mm[2,], na.rm=TRUE )/ sum(ps[2,],na.rm=TRUE)) * ps[2,] 
  ps[3,] = (sum( mm[5,], na.rm=TRUE )/ sum(ps[3,],na.rm=TRUE)) * ps[3,] 
  rownames(ps) = c( paste0(a, InputName), paste0(b,InputName), paste0(a,"-",b)) 
  
  nucs = ps
  for( i in 1:3 ) 
    nucs[i,] = rpois(length(ps[i,]),ps[i,])
  
  return(nucs)
}

TestSimulation <- function( a, b, pre ) {
  x = BuildExpName(a, pre=pre) 
  y = BuildExpName(b, pre=pre) 
  z = BuildExpName(a,b, pre=pre) 
  p1 = PlotXYZ( x, y, z) 
  nucs.min = SimulatePair(a, b, Mode = "min", pre = pre) 
  p2 = PlotXYZ(  x, y, z, nucs = nucs.min ) + labs(title="Minimal")
  nucs.exp = SimulatePair(a, b, Mode = "exp", pre = pre)
  p3 = PlotXYZ(  x, y, z, nucs = nucs.exp ) + labs(title="Random")
  nucs.max = SimulatePair(a, b, Mode = "max", pre = pre)
  p4 = PlotXYZ(  x, y, z, nucs = nucs.max ) + labs(title="Maximal")
  fn = paste0(DataDir, "/", "Sim-Scatter", x,"-",y,".png")
  png(fn,height=600, width=1224)
  multiplot(p1,p1+ scale_x_log10()+scale_y_log10(),
            p2, p2+scale_x_log10()+scale_y_log10(),
            p3,p3 + scale_x_log10()+scale_y_log10(), 
            p4, p4 +  scale_x_log10()+scale_y_log10(), layout = matrix(c(1,2,3,4, 5, 6, 7, 8), nrow=2))
  dev.off()

  p1 = RelativeEnrichmentPlot( a, b, pre = pre) + labs(title="Observed")
  p2 = RelativeEnrichmentPlot( a, b, nucs = nucs.min, pre = pre ) + labs(title="Minimal")
  p3 = RelativeEnrichmentPlot( a, b, nucs = nucs.exp, pre = pre ) + labs(title="Random")
  p4 = RelativeEnrichmentPlot( a, b, nucs = nucs.max, pre = pre ) + labs(title="Maximal")
  
  fn = paste0(DataDir, "/", "Sim-Rel-Scatter", x,"-",y,".png")
  png(fn,height=600, width=1224)
  multiplot(p1,p1+ scale_x_log10()+scale_y_log10(),
            p2, p2+scale_x_log10()+scale_y_log10(),
            p3,p3 + scale_x_log10()+scale_y_log10(), 
            p4, p4 +  scale_x_log10()+scale_y_log10(), layout = matrix(c(1,2,3,4, 5, 6, 7, 8), nrow=2))
  dev.off()
  
}

PlotComparison <- function( dmono, dsub,fn, type = "TSS", order = "expr", filetype = "png", 
                            MinLength = 500, MaxLength = 10000 ) {
  o = NULL
  if( type == "TSS") {
    offset= 500
    mi = 1
  }
  if( type == "Genes") {
    offset= 500
    mi = 3
  }
  if( type == "TTS") {
    offset= 1000
    mi = 2
  } 
  if( type == "PostTTS") {
    offset= 1050
    mi = 2
  } 
  if( type == "Plus1") {
    offset = 500
    mi = 4
  }
  if( order == "expr") {
    v = log(SGD$Genes$expr)
  } else
    if( order == "length" ) {
      v = width(SGD$Genes)
    } else
      if( order == "occupancy-m" || order == "occupancy-s") {
        if( order == "occupancy-m")
          meta = dmono$meta[[mi]]
        else
          meta = dsub$meta[[mi]]
        v = rowMeans(meta[,offset-50:offset+50])
      }
    
  o = order(v, na.last = F)
  
  S = SGD$Genes[o]
  Ii = which(width(S) > MinLength & width(S) < MaxLength)
  I = S[Ii]$acc
  
  p1 = ccPlotAlignedMatrix(dmono,params, order = o, group = I, type = type )
  p2 = ccPlotAlignedMatrix(dsub,params, order  = o, group = I, type = type )

  p1 = p1 + geom_vline(xintercept = offset,colour="yellow")
  p2 = p2 + geom_vline(xintercept = offset,colour="yellow")
  
  p3 = ggplot(data.frame(x = v[o[Ii]], y = 1:length(o[Ii])), aes(x = x, y = y)) + geom_line() 
  p3 = p3 + labs(title = order, x = "") 
  p3 = p3 + xlim(c(min(v[o[Ii]]),quantile(v[o[Ii]],.95,na.rm=T)))
  fname = ccBuildFN(paste0(fn,"-", type, "-by-", order),params, suff = paste0(".", filetype))
  if( filetype == "png")
    png(fname, width=1024, height = 1024)
  else
    if( filetype == "pdf")
      pdf(fname,width = 8, height = 8, pointsize = 10)
  
  multiplot(p1,p2, p3, layout = matrix(c(1,1,1,2,2,2,3), nr = 1))
  dev.off()
}

PlotHeatMatrix <- function( dat, fn, type = "TSS", order = "expr", filetype = "png"  ) {
  o = NULL
  if( type == "TSS") {
    offset= 500
    mi = 1
  }
  if( type == "TTS") {
    offset= 1000
    mi = 2
  } 
  if( type == "PostTTS") {
    offset= 1050
    mi = 2
  } 
  if( type == "Plus1") {
    offset = 500
    mi = 4
  }
  if( order == "expr") {
    v = log(SGD$Genes$expr)
  } else
    if( order == "length" ) {
      v = width(SGD$Genes)
    } else
      if( order == "occupancy" ) {
          meta = dat$meta[[mi]]
          v = rowMeans(meta[,offset-50:offset+50])
      }
  
  o = order(v, na.last = F)
  
  S = SGD$Genes[o]
  Ii = which(width(S) > 500)
  I = S[Ii]$acc
  
  p1 = ccPlotAlignedMatrix(dat,params, order = o, group = I, type = type )

  p1 = p1 + geom_vline(xintercept = offset,colour="yellow")

  p3 = ggplot(data.frame(x = v[o[Ii]], y = 1:length(o[Ii])), aes(x = x, y = y)) + geom_line() 
  p3 = p3 + labs(title = order, x = "") 
  p3 = p3 + xlim(c(min(v[o[Ii]]),quantile(v[o[Ii]],.95,na.rm=T)))
  #fname = ccBuildFN(paste0(fn,"-", type, "-by-", order),params, suff = paste0(".", filetype))
  fname = fn
  if( filetype == "png")
    png(fname, width=1024, height = 1024)
  else
    if( filetype == "pdf")
      pdf(fname,width = 8, height = 8, pointsize = 10)
  
  multiplot(p1, p3, layout = matrix(c(1,1,1,1,1,2), nr = 1))
  dev.off()
}

CollectOccupancy <- function(d, RegionsTSS, RegionsTTS) {
 means.TSS = do.call("cbind", lapply(RegionsTSS, function(r) rowMeans(d$meta[[1]][,r])))
 means.TTS = do.call("cbind", lapply(RegionsTTS, function(r) rowMeans(d$meta[[2]][,r])))
 m = cbind( means.TSS, means.TTS)
 return(m)
}


CompareStrains <- function(x,y, s1,s2) {
  DensityScatter( cnucs[BuildExpName(x,y,pre=s1),], cnucs[BuildExpName(x,y,pre=s2),], 
                  xlabel = s1, ylabel =s2, title = x, threshold = .99,
#                  smooth = T,
                  coordeq = FALSE, alpha = .1, logscale = F)
}

ComputeCorrelationMatrix <- function( nucs, Indexes = NULL, Ref = NULL, Relative = FALSE) {
  if( !is.null(Indexes) )
    nucs = nucs[Indexes,]
  
  if(is.null(Ref))
    Ref = as.vector(rep(1,dim(nucs)[2]))
  
  if( Relative )
    Ref = colMeans(nucs, na.rm = TRUE)*Ref
  
  N = dim(nucs)[[1]]
  Cor = matrix(nc = N, nr = N)
  rownames(Cor) = rownames(nucs)
  colnames(Cor) = rownames(nucs)
  for( i in 1:N) 
    for (j in 1:N) 
      Cor[i,j] = cor(nucs[i,]/Ref, nucs[j,]/Ref, use="complete.obs")
  
  return(Cor)
}