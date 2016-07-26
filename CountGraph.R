#DataDir = "~/Data/CoChIP/160124/RPD3"
#DataDir = "~/Data/CoChIP/160124/K4Sym"
DataDir = "~/Data/CoChIP/HisMut-K4/HisMut/"
#ExpName = "RPD3"
#ExpName = "K4Sym"
ExpName = "HisMut"
ExpOffset = 0
#ExpOffset = 1
WT = "WT"
#WT= "bar1"
FigureDir = "~/Dropbox/CoChIP/figures/"

PLotCounts <- function(M,zrange=c(1,6),colorRange = NULL, color = FALSE, text= TRUE, log=TRUE, crange=c(-6,6)) {
  Nx = dim(M)[[1]]
  Ny = dim(M)[[2]]
  if( log ) {
    Z = log10(M)
  } else
    Z = M
  Z[ Z < zrange[[1]] ] = zrange[[1]]+.1
  Z[ Z > zrange[[2]] ] = zrange[[2]]-.1
  if( is.null(colorRange)) {
    C = M
  } else {
    color = TRUE
    C = colorRange
  }
  
  C[ C < crange[[1]]] = crange[[1]]
  C[ C > crange[[2]]] = crange[[2]]
  
  df = data.frame( x = rep(1:Ny, each=Nx), 
                   y = rep(1:Nx,Ny), 
                   z = as.vector(Z),
                   c=as.vector(C),
                   t = sapply(M, function(m) format(m/1000,digits = 2)),
                   ty = rep(1:Nx,Ny) - .35
            )
  if ( color)
    p = ggplot(df, aes(x = x, y=y, size=z,colour = c))
  else
    p = ggplot(df, aes(x = x, y=y, size=z))
  
  xlabels = sub("H3K","K", rownames(M))
  ylabels = sub("H3K","K", colnames(M))
  p = p+scale_y_continuous(breaks = 1:Nx, labels = xlabels, lim=c(0.5,Nx+.5))
  p = p+scale_x_continuous(breaks = 1:Ny, labels = ylabels, lim=c(0.5,Ny+.5))
  p = p+geom_point()
  p = p+ scale_size("Reads (x 1000)", range = c(0,15), limits = zrange)
  if( text )
    p = p + geom_text(aes(label=t,y=ty,x=x), size=8)
  p = p + coord_fixed(ratio=1) 
  p = p + theme_bw()
  p = p+xlab("2nd IP")
  p = p+ylab("1st IP")
  if( color )
    p = p+scale_color_gradient2("Ratio to WT (log2)", low="green", high="red", mid="black", limits = crange)
  p
}

cnts = read.csv(paste0(DataDir,"/QCsummary.csv"),header=T)

CNames = as.character(cnts[[1]])
CNames = unlist(lapply(strsplit(CNames,"_"), function(x) x[[1]]))
CValue = cnts[[6]]
CStrain = unlist(lapply(strsplit(CNames,"-"), function(x) x[[1+ExpOffset]]))
CAb1 = unlist(lapply(strsplit(CNames,"-"), function(x) x[[2+ExpOffset]]))
CAb2 = unlist(lapply(strsplit(CNames,"-"), function(x) x[[3+ExpOffset]]))

CMs = list()
for( s in unique(CStrain)) {
  M = matrix(0, nr = length(unique(CAb1)), nc = length(unique(CAb2)))
  rownames(M) = unique(CAb1)
  colnames(M) = unique(CAb2)
  rownames(M) = c("H3", "H3K79me3", "H3K36me3","H3K18ac", "H3K4me3")
  colnames(M) = c("Input","H3K4me3", "H3K18ac","H3K36me3",  "H3K79me3", "H3" )
  for( i in which(CStrain == s) )
      M[CAb1[[i]],CAb2[[i]]] = CValue[[i]] 

#  Input = M[,"Input"]
#  M = M[,rownames(M)]
#  AbO = rownames(M)
#  M = M[AbO,c("Input", AbO)]
  
  CMs[[s]] = M
} 

for( s in unique(CStrain)) {
  if( s == WT){
    p = PLotCounts((CMs[[s]]),zrange=c(2,6.5)) + labs(title=s)
  } else
    p = PLotCounts((CMs[[s]]),zrange=c(2,6.5), colorRange =log2(CMs[[s]]/CMs[[WT]])) + labs(title=s)
  ggsave(filename=paste0(FigureDir,"/",ExpName,"-",s,"-counts.pdf"),plot = p)
}

if( 0 )
for( s in names(CMs))
  if( s != WT)
  {
    p = PLotCounts(log2(CMs[[s]]/CMs[[WT]]),zrange=c(-6,6), log=FALSE,crange = c(-3,3),TRUE) +
      labs(title=paste(s,"vs WT"))
    ggsave(filename=paste0(FigureDir,"/", ExpName,"-",s,"-counts-ratio.png"),plot = p)
  }