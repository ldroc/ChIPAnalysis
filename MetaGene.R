library(reshape)

CollectRegions <- function( cov, regions, width = 1500 ) {
  if( is.null(regions) || is.null(cov) ) 
    return (NULL)

  mat <- matrix(nrow = length(regions), ncol = width)
  mat[,] = 0
  rownames(mat) <- regions$acc
  colnames(mat) <- 1:width
  chr <- levels(seqnames(regions))
  for( x in chr ) {
 #   print(x)
    if( any(seqnames(regions) == x)) {
      m  = length(cov[[x]])
      I = start(regions) > 0 & end(regions) < m & seqnames(regions) == x
      Ip = strand(regions) == '+' & I
      Im = strand(regions) == '-' & I
      
      v = Views(cov[[x]], ranges(regions)[Ip])
      np = length(which(Ip))
      if( np > 0 ) {
        mp = matrix(nrow=np, ncol=width )
        for( i in 1:np )
          mp[i,] = as.vector(v[[i]])
        mat[which(Ip),] <- mp
      }
      v = Views(cov[x], ranges(regions)[Im])[[1]]
      nm = length(which(Im))
      if( nm > 0 ) {
        mm = matrix(nrow=nm, ncol=width )
        for( i in 1:nm )
          ## Note the reversal
          mm[i,width:1] = as.vector(v[[i]])
        mat[which(Im),] <- mm
      }
    }
  }
  return(mat)
}

AvgFromCollection <- function( mat, I ) {
#  print(length(I))
#  print(length(which(I)))
#  print(dim(mat))
  a = colMeans(mat[I,],na.rm = TRUE)
}


AvgRegionsSubgroups <- function(mat, rs, subs){
  A = lapply(subs, function(A) AvgFromCollection(mat, is.element(rs$acc, A) ) )
  names(A) = names(subs)
  return(A)
}


PlotColors = c( rgb(0.843,0.145,0.1843), rgb(0.855,0.145,0.561), rgb(0.5922,0.301,0.601), rgb(0.333,0.321,0.6275), rgb(0.255,0.404,0.675), "darkgreen", "maroon", "purple", "orange", "green", "")


PlotCovergeGroups = function(A, pos=500, label = "TSS", main = "", ylim = NULL)
{
  cs = PlotColors[1:length(A)]
  
#  cs = rainbow((length(A)*1.5))
  if( is.null(ylim) ) {
    M = max(unlist(lapply(A,max)),na.rm=TRUE)*1.1
    m = min(0,min(unlist(lapply(A,min)),na.rm=TRUE)*1.1)
    ylim=c(m,M)
  }
  w = length(A[[1]])
  if( is.nan(w) || is.infinite(w))
    w = 1500
  plot( (1:w)-pos, A[[1]], type='l', col=cs[1], axes=FALSE, main=main, ylim = ylim, xlim=c(-pos,w-pos+500),
        xlab="Relative position", ylab="Avg occupancy", lty="solid",lwd=3)

  ipos = as.integer(pos/500)*500
  iwidth = as.integer(w/500)*500
  labelat = seq(-ipos,iwidth, 500)
  labelval = labelat
  labelval[ipos/500+1] = label
  axis(1,col="black", labels = labelval, at = labelat)
  axis(2,col="black")
  
  for( i in 2:length(A))
    lines((1:w)-pos, A[[i]], type='l', col=cs[i], lty="solid",lwd=2)
  legend("bottomright", legend = names(A), col=cs, lty="solid",lwd=2)
  abline(v=0,lty="dotted", col="gray",lwd=2)
}


ggPlotCovergeGroups = function(A, pos=500, label = "TSS", main = "", ylim = NULL)
{
  if( is.null(ylim) ) {
    M = max(0,max(sapply(A,max),na.rm=TRUE)*1.1)
    m = min(0,min(sapply(A,min), na.rm = TRUE) *1.1)
    ylim=c(m,M)
  }
  w = length(A[[1]])
  if( is.nan(w) || is.infinite(w))
    w = 1500
  
  p = ggplot() + theme(legend.position = c(1,0), legend.justification=c(1,0)) 
  p = p + scale_color_brewer( palette = "Set1",
                              guide = guide_legend(title=NULL, reverse = TRUE))
  
  # plot lines
  p = p+geom_vline(xintercept = 0, colour="gray", size=1)
  if( ylim[[1]] < 0 && ylim[[2]] > 0 )
    p = p+geom_hline(yintercept = 0, colour="gray", size=1)

  for( i in 1:length(A)) {
    p = p + geom_line( data=data.frame(x = (1:w)-pos, y = A[[i]], color = rep(names(A)[[i]], 1)), 
                       aes(x = x, y = y, colour = color),size = .5, 
                       show.legend = TRUE)
  }
  
  ipos = as.integer(pos/500)*500
  iwidth = as.integer(w/500)*500
  labelat = seq(-ipos,iwidth, 500)
  labelval = labelat
  labelval[ipos/500+1] = label
  
  p = p + ylim(ylim)
  p = p+scale_x_continuous(name="Position", breaks = labelat, labels = labelval, limits = c(-pos,w-pos+500) )
  p = p + labs(title=main,  y="Avg Occupancy") 
  
  return(p)
}

PlotMultiCoverage <- function (Files, Gr1, Gr2, path=NULL, label1="TSS", label2="TTS", pos1=500, pos2=1000, PDF = F, NameAdd = "" ) {
  tf <- function( i ) {
    f = Files[i]
    suff = ".png"
    if( PDF )
      suff = ".pdf"
    if( !is.null(path) )
      out = paste0(path,"/", BaseFileName(f),NameAdd,"-Meta",suff)
    else
      out = paste0(file_path_sans_ext(f),NameAdd,"-Meta", suff)
    if( PDF)
      pdf(out)
    else
      png(out,height = 1024, width=1024)
    parparam = par()
    par(mfrow = c(2,1),mar=c(2.1,3.1,2.1,2.1))
    PlotCovergeGroups(Gr1[[i]], main=paste(label1, "aligned -", Files[i]), label=label1, pos=pos1)
    PlotCovergeGroups(Gr2[[i]], main=paste(label2, "aligned -", Files[i]), label=label2, pos=pos2)
    dev.off()
    # par(parparm)
    suppressMessages(par( mfrow = parparam$mfrow))
    suppressMessages(par( mar = parparam$mar))
  }
 
  sapply(1:length(Gr1), tf)
}


PlotHeatMap <- function( m, xlab = "", ylab = "Y", zlab = "Z", title = "", 
                         zlim = NULL,
                         offset = 500, base = "TSS" ) {
  n = dim(m)[[1]]
  colnames(m) <- 1:dim(m)[[2]]
  breaks = seq(0,n,500) 
  labels = breaks - offset
  labels[labels == 0] = base
  ybreaks = seq(0,dim(m)[[2]], 500)
    df  = melt(m, varnames = c( "x", "y"))
  t1 = quantile(m, probs = .95, na.rm = T)
  t2 = quantile(m, probs = .05, na.rm = T)
  if( t2 > 0 )
    t2 = 0
  print(c(t2,t1))
  if( is.null(zlim) ) 
    zlim = c(t2,t1)
  df$value[df$value < zlim[[1]]] = zlim[[1]]
  df$value[df$value > zlim[[2]]] = zlim[[2]]
#  p = ggplot(df, aes(xlab, ylab, fill=value))
  p = ggplot(df, aes(x=x,y=y, fill=value))
#  p = p + geom_raster(color = "white")
  p = p + geom_raster()
  if( is.null(zlim)){
#    p = p + scale_fill_gradient(low="black", high="red", na.value="black" )
    p = p + scale_fill_gradient2(mid="black", low="green", high="red", na.value="gray")
  } else
    p = p + scale_fill_gradient2(limits=zlim, mid="black", low="green", high="red", na.value="gray")
  
  p = p + scale_x_continuous( breaks = breaks, labels = labels, limits=c(0,n),expand=c(0,0) )
  p = p + scale_y_continuous( breaks = ybreaks, labels = ybreaks, expand=c(0,0) )
  p = p + labs(x = xlab, y = ylab, fill=zlab, title = title)
  p = p + theme_minimal()
#  p = p + theme( legend.position = "none",
#                axis.ticks = element_blank()
#  )
  return(p)
}