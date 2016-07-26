CompareToWT <- function(Strain, ip1, ip2 = "Input") {
  n1 = paste0("WT-",ip1,"-",ip2)
  n2 = paste0(Strain,"-",ip1,"-",ip2)
#  print(n1)
#  print(n2)
  x = nuc.mat[n1,]
  y = nuc.mat[n2,]
  p = DensityScatter(x,y, threshold = .99, diagonal = T, coordeq = F, 
                     linearfit = F, jitter = T,
                 xlabel = "WT", ylabel = Strain, title = paste0(ip1,"-",ip2))
  b = ComputeTopRatio(x,y, threshold = .1)
  p + geom_abline(slope = b, intercept = 0,color="red", size=2)
} 

GLMCorrelation <- function( x, y ) {
  x = x 
  y = y 
  df = data.frame( x = x, y = y)
  gl = glm( y ~ x-1, family = poisson())
  print(gl)
  Nx = sum(x)
  Ny = sum(y)
  N = length(x)
  print(paste("Poisson best slope =", Ny/(Nx+Ny+N)))
}

BuildCorMatrix <- function() {
  Names = sapply(rownames(nuc.mat), function(x) strsplit(x,"-")[[1]])
  
  Strains = unique(Names[1,])
  Exp = unique(sapply(1:length(rownames(nuc.mat)), function(i) paste0(Names[2,i], "-", Names[3,i])))
  
  cor = matrix( nr = length(Exp), nc = length(Strains), dimnames = list(Exp, Strains))
  
  for( s in Strains )
    for( e in Exp )
      cor[e,s] = cor(nuc.mat[paste0("WT-",e),],nuc.mat[paste0(s,"-", e),])
  
  cor
}