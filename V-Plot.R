CreateVPlotCov <- function( GR, MaxLen = 600, Step = 10) {
 
  L = seq(0,MaxLen-Step,Step) 
  W = width(GR)
  
  tempFunc <- function(l) {
    I = W >= l & W <= l+Step
    if( any(which(I)) ) {
      return(coverage(GR[I]))
    } else
      return(coverage(GR[1]))
  }
  C = lapply(L, tempFunc)  
  names(C) <- L+Step/2
  C
}

CreateVPlot <- function( VCov, Chr, Start, End ) {
#  print(Chr)
#  print(Start)
#  print(End)
  L = names(VCov)
  M = matrix(0 , nrow = length(L), ncol = End - Start+1)
  rownames(M) <- L
  colnames(M) <- seq(Start, End)
  for( l in L ) {
    M[l,] = as.vector(VCov[[l]][[Chr]][Start:End])
  }
  M
}

AddVPlot <- function( M, VCov, Chr, Start, End, Strand ) 
{
  if( !all(Strand)) {
    t = Start
    Start = End
    End = t
  }
  for( l in 1:length(VCov) ) {
    M[l,] = M[l,] + as.vector(VCov[[l]][[Chr]][Start:End])
  }
  M
}

CollectVPlots <- function( VCov, regs, width = 1500 ) {
  M = matrix(0, nrow = length(names(VCov)), ncol = width)
  rownames(M) <- names(VCov)
  Chr = seqnames(regs)
  Start = start(regs)
  End = end(regs)
  Strand = strand(regs) == '+'
  if( length(regs) == 0 )
    return(M)
  n = 0
  for( i in 1:length(regs) ) 
    if( Start[i] > 0 && End[i] < length(VCov[[1]][[as.character(Chr[i])]])) {
      M = AddVPlot( M, VCov, as.character(Chr[i]), Start[i], End[i], Strand[i])
      n = n+1
      if( n %% 50 == 0)
        print(paste(n,i))
    }
  return(M/n)
}

CollectGroupVPlots <- function(cv, regs, groups, width = 1500) {
  tempFunc <- function(i) {
    I = is.element(regs$acc, groups[[i]])
    M = CollectVPlots(cv, regs[I], width=width)
  }
  Ms = lapply(1:length(groups), tempFunc )
  names(Ms) <- names(groups)
  
  return(Ms)
}

PlotVPlot <- function( M, Pos = 500, Label = "TSS", width=1500, Filter = 0, 
                       main = "", zlim=NULL, col = rev(heat.colors(255)) ) {
#  x = 1:width+1 - Pos
  x = 1:width - Pos
  L = rownames(M)
  Step = as.integer(L[[length(L)]])-as.integer(L[[length(L)-1]])
#  y = seq(as.integer(L[[1]]), as.integer(L[[length(L)]])+Step, Step)
  y = seq(as.integer(L[[1]]), as.integer(L[[length(L)]]), Step)
#  print(length(x))
#  print(length(y))
#  print(dim(M))
  if( Filter > 0 ) {
    f = makeBrush(2*(Filter%/%2)+1,shape="disc", step = FALSE)
    f = f/sum(f)
    M = filter2(M,f)
  }
  if( is.null(zlim) )
    zlim = c(0,max(M))
  
  image(x,y, t(M), useRaster=TRUE,xlab=main, ylab = "Length", zlim = zlim, col = col)
  title(main=main)
}

PlotMultiVPlot <- function( Ms, file = NULL, path = NULL, Pos = 500, Label = "TSS", Filter = 0) {
  if( !is.null(file) ) {
    if( !is.null(path) )
      out = paste0(path,"/", BaseFileName(file),"-Meta.pdf")
    else
      out = paste0(file_path_sans_ext(file),"-Meta.pdf")
    pdf(out,paper="a4",pointsize = 10)
  }
  
  ns = names(Ms)
  print(ns)
  width = dim(Ms[[1]])[2]
  par(mfrow = c(length(ns),1))
  parparam =  par(mar=c(3.1,4.1,1.5,2.1))
  
  m = max(unlist(Ms))
  Ms = lapply(Ms, function(M) 256*M/m)
  col = hcl(70,seq(0,100,100/256),seq(0,100,100/256))
  lapply(1:length(ns), function(i) PlotVPlot( Ms[[i]], width = width, Label=Label, 
                                              zlim = c(0,255), Pos=Pos, main = ns[[i]],
                                              col=col, Filter = Filter))
  par(parparam,mfrow=c(1,1))
  if( !is.null(file) )
    dev.off()
}