##
## Model Design:
##
## Exps:  Strain, IP1, IP2
## Strains: Strain, Interactions with each Ab
## Loading factor: Loading amount (compared to Input)
eps = 1e-6

p.GradNames = TRUE
library(abind)


# Utilities ---------------------------------------------------------------

blkdiag.wnames <- function(Ms) {
  rows = sapply(Ms, function(M) dim(M)[1])
  cols = sapply(Ms, function(M) dim(M)[2])
  A = matrix(0,nr=sum(rows), nc = sum(cols))
  rownames(A) = sapply(Ms, rownames)
  colnames(A) = sapply(Ms, colnames)
  R = 0
  C = 0
  for( M in Ms ) {
    A[R+1:dim(M)[1],C+1:dim(M)[2]] = M
    R = R + dim(M)[1]
    C = C + dim(M)[2]
  }
  A
}

blkdiag3.wnames <- function(Ms) {
  rows = sapply(Ms, function(M) dim(M)[1])
  cols = sapply(Ms, function(M) dim(M)[2])
  depth = dim(Ms[[1]])[3]
  
  A = array(0,dim = c(sum(rows),sum(cols),depth),
            dimnames = list( sapply(Ms, rownames),
                             sapply(Ms, colnames),
                             dimnames(Ms[[1]])[[3]] )
  )
  
  R = 0
  C = 0
  for( M in Ms ) {
    A[R+1:dim(M)[1],C+1:dim(M)[2],] = M
    R = R + dim(M)[1]
    C = C + dim(M)[2]
  }
  A
}
# Data structure Constructors ------------------------------------------------------

CreateDesign <- function( ExperimentList ) {
  Experiments = do.call(rbind,strsplit(ExperimentList,"-"))
  colnames(Experiments) = c("Strain","IP1", "IP2")
  rownames(Experiments) = ExperimentList
  ABs = paste0("ab-",unique(Experiments[,"IP1"]))
  Strains = unique(Experiments[,"Strain"])
  
  Mods = unique(Experiments[,"IP1"])[grep("H3", ABs, invert = TRUE)]
#  Mods = c(Mods, "K9ac")
  
  Experiments[,"IP1"] = paste0("ab-",Experiments[,"IP1"])
  Experiments[Experiments[,"IP2"] == "Input","IP2"] = NA
  Experiments[!is.na(Experiments[,"IP2"]),"IP2"] = paste0("ab-",Experiments[!is.na(Experiments[,"IP2"]),"IP2"])
  Design = list(
    # list of modifications
    Mods = Mods,
    # list of IPs
    IPs = ABs,
    # list of Strains/Samples
    Strains = Strains,
    
    Exps = data.frame( 
      Strain = Experiments[,"Strain"],
      IP1 =  Experiments[,"IP1"],
      IP2 = Experiments[,"IP2"],
      row.names = rownames(Experiments)
    ),
    
    # constrains on strain/modification combination
    Constraints = data.frame( 
      Strain = c(), 
      Mod = c(), 
      I = c() 
    )
  )
}

##
## Model Parameters: 
##
## IPs: Ip vs Mark affinity
## Base: non-specific, background
## Recovery: sequencing recovery for each second IP (compared to Input)
## Ms: MNase success for each strain
## Tags: tagging success for each first IP
##

CreateModel <- function(Design) {
  ModNumber = length(Design$Mods)
  StrainNumber = length(Design$Strains)
  IPsNumber = length(Design$IPs)
  
  Model = list(
    IPs = array( 1e-4, dim= c(IPsNumber, ModNumber),
                 dimnames = list( Design$IPs, Design$Mods )  
    ),
    Base = array( 1e-4, 
                  dim=c(IPsNumber), 
                  dimnames = list( Design$IPs )
    ),
    
    Recovery = array( 1,
                      dim=c(IPsNumber), 
                      dimnames = list( Design$IPs )
    ),
    Total = 1e7,
    Ms = array(1, dim=c(StrainNumber),  dimnames = list( Design$Strains) ),
    
    Tags = array(0.01,dim=c(IPsNumber), dimnames = list( Design$IPs ))
  )
}


##
## Nucleosome parameters
##
## Occ: Occupancy in each strain
## Marginals: value marginal for Mod in each strain
## Interactions: list of parameters for each (Mod,Mod,Strain)
##

CreateNuc <- function(Design, NucName = "") {
  ModNumber = length(Design$Mods)
  StrainNumber = length(Design$Strains)

  Nuc.initial =  list(
    Label = NucName,
    Occ = array(1,dim = c(StrainNumber), 
                dimnames = list(Design$Strains) ),
    Marginals = array( 1, dim = c(ModNumber, StrainNumber),
                       dimnames = list( Design$Mods, Design$Strains)),
    #
    # Upper diagonal part is ignored
    #
    Interactions = array( 0, dim = c(ModNumber,ModNumber,StrainNumber), 
                          dimnames = list(Design$Mods, Design$Mods, Design$Strains))
  )
}


# Basic packing/unpacking operations --------------------------------------------------------


negorder = c(1,3,2,4)  

abModelPackNucParamsStrain <- function( Nucs, Strain ) {
  A = Nucs$Interactions[,,Strain]
  c( Nucs$Occ[[Strain]], Nucs$Marginals[,Strain], A[lower.tri(A)] )
}


abModelUpperNucParamsStrain <- function( Nucs, Strain, Compile ) {
  Ni = dim(Nucs$Interactions)[1]
  c( 1, Compile$MarginUpper[[Strain]], rep(10, Ni*(Ni-1)/2))
}

abModelPackNucParams <- function( Nucs ) {
  sapply( 1:dim(Nucs$Occ)[1], function(s) abModelPackNucParamsStrain( Nucs, s ))
}

abModelLowerNucParamsStrain <- function( Nucs, Strain, Compile ) {
  Ni = dim(Nucs$Interactions)[1]
  c( 1e-5, Compile$MarginLower[[Strain]], rep(-10, Ni*(Ni-1)/2))
}

abModelUpperNucParams <- function( Nucs, Compile ) {
  sapply( 1:dim(Nucs$Occ)[[1]], function(s) abModelUpperNucParamsStrain( Nucs, s, Compile ))
}

abModelLowerNucParams <- function( Nucs, Compile ) {
  sapply( 1:dim(Nucs$Occ)[[1]], function(s) abModelLowerNucParamsStrain( Nucs, s, Compile ))
}

abModelUnPackNucParamsStrain <- function( Nucs, Strain, x ) {
  n = 1
  l = dim(Nucs$Marginals)[1]
  m = l*(l-1)/2
  Nucs$Occ[[Strain]] = x[1:n]
  Nucs$Marginals[,Strain] = x[(n+1):(n+l)]
  A = Nucs$Interactions[,,Strain]
  A[lower.tri(A)] = x[(n+l+1):(l+n+m)]
  A[upper.tri(A)] = t(A)[upper.tri(A)]
  Nucs$Interactions[,,Strain] = A
  return(Nucs)
}

abModelUnPackNucParams <- function(Nucs, x) {
  l = length( abModelPackNucParamsStrain(Nucs, 1))
  for( s in  1:length(Nucs$Occ) ) {
    Nucs = abModelUnPackNucParamsStrain(Nucs, s, x[(l*(s-1)+1):(l*s)])    
  }
  return(Nucs)
}

abModelPackModelParams <- function( Model ) {
  c( Model$Total, Model$Ms, Model$Tags, Model$Recovery, Model$Base, Model$IPs )
}

abModelUpperModelParams <- function( Model ) {
  up = c( Inf, 
          rep(2, length(Model$Ms)), 
          rep(1,length(Model$Tags)),
          rep(Inf, length(Model$Recovery)),  
          rep(1, length(Model$Base)), 
          rep(1, length(Model$IPs)))
}

abModelLowerModelParams <- function( Model ) {
  low = c( eps,
           rep(0.5, length(Model$Ms)),
           rep(eps, length(Model$Tags)),
           rep(eps, length(Model$Recovery)), 
           rep(eps, length(Model$Base)), 
           rep(eps, length(Model$IPs)))
}

abModelUnpackModelParams <- function( Model, x ) {
  l = length(Model$IPs)
  k = length(Model$Base)
  n = length(Model$Recovery)
  m = length(Model$Ms)
  t = length(Model$Tags)

  Model$Total = x[[1]]
  Model$Ms[] = x[2:(m+1)]
  Model$Tags[] = x[(m+2):(m+t+1)]
  Model$Recovery[] = x[(m+t+2):(m+t+n+1)]
  Model$Base[] = x[(m+t+n+2):(m+t+n+k+1)]
  Model$IPs[] = x[(m+t+n+k+2):(m+t+n+k+l+1)]
  return(Model)
}

# Interaction model -------------------------------------------------------


##
## Build pairwise occupancies for all pairs
##
## Assumes row/col of Interactions is the same order as Marginals
##
## Max = min(p,q)
## Min = p+q-1
## E = p*q
## if( a >= 0) X = E + 2* (Max-E)*(sigmoid(a,1/(Max-E)) - .5)
## if( a < 0) X = E + 2* (E-Min)*(sigmoid(a,1/(E-Min)) - .5)
##
## d/dp Max = (p < q) ? 1 : 0
## d/dp Min = 1
## d/dp E = q
## d/dp P = d/dp (Max - E) = (p < q) ? 1 : 0 - q
## d/dp N = d/dp (E - Min) = q - 1
## d/dz s(a,1/z) = -1/(z**2) * s(a,1/z)*(1-s(a,1/z))
## if( a >= 0 & p < q) d/dp X = q + 2*(1-q)(s(a,1/P) - 0.5) - 2/P *s(a,1/P)*(1-s(a,1/P))
## if( a >= 0 & p >= q) d/dp X = q + 2*(s(a,1/P) - 0.5) - 2/P *s(a,1/P)*(1-s(a,1/P))
## if( a < 0 ) d/dp X = q + 2* (q-1)*(sigmoid(a,1/N) - 0.5) - 2/N * s(a,1/N)*(1-s(a,1/N))

library(gtools)
library(pracma)

sigmoid2 <- function(x) {
  2*sigmoid(x) - 1
}

logit2 <- function(x) {
  x = pmin(x,1)
  x = pmax(x,-1)
  logit((x+1)/2)
}
dsigmoid2 <- function( x ) {
  (1+sigmoid2(x))*(1-sigmoid2(x))/2
}

CombineInteraction <- function( p, q, I) {
  E = p*q
  Z = array(0,dim=c(length(I)))
  Z[E <= eps] = 0
  P = pmin(p,q)-E
  P = pmax(P,eps)
  Z[I > 0 & E > eps] = (E + P*sigmoid2(I/P)) [ I > 0 & E > eps ]
  N = E - pmax(0,(p+q-1))
  N = pmax(N,eps)
  Z[I <= 0 & E>eps] = (E + N*sigmoid2(I/N))[I <= 0 & E > eps ]
#  print(paste("ComputeI = E = ", E, " I = ", I, "Z = ", Z))
  return(Z)
}

ExtractInteraction <- function( p, q, Z) {
  E = p*q
  I = array(0,dim=c(length(Z)))
  P = pmin(p,q)-E
  P = pmax(P,eps)
  Z = Z - E
  I[Z > 0 & E > eps] = (P* logit2(Z/P))[ Z > 0 & E > eps ]
  N = E - pmax(0,(p+q-1))
  N = pmax(N,eps)
  I[Z <= 0 & Z>eps] = (N*logit2(Z/N))[Z <= 0 & E > eps ]
  return(I)
}

CombineInteractionGradP <- function ( p, q, I ) {
  E = p*q
  Z[!is.nan(I)] = q
  P = pmin(p,q)-E
  P = pmax(P,eps)
  Z = array(0,dim=c(length(I)))
  Z[I > 0] = (Z + q*(I/P)*dsigmoid2(I/P) - q*sigmoid2(I/P))[I > 0]
  Z[I > 0 & p < q] = (Z -(I/P)*dsigmoid2(I/P) + sigmoid2(I/P))[I > 0 & p < q]
  
  N = E - pmax(0,(p+q-1))
  N = pmax(N,eps)
  Z[I <= 0] = (Z - q*(I/N)*dsigmoid2(I/N) + q*sigmoid2(I/N))[I <= 0]
  Z[I <= 0 & p+q > 1 ] = (Z + (I/N)*dsigmoid2(I/N) - sigmoid2(I/N))[I <= 0 & p+q>1]
  return(Z)
}

CombineInteractionGradI <- function ( p, q, I ) {
  E = p*q
  if( E < eps )
    return(0)
  Z = 0
  if( I > 0 ) {
    P = pmin(p,q)-E
    P = pmax(P,eps)
    Z = (dsigmoid2(I/P))
  } else {
    N = E - pmax(0,(p+q-1))
    N = pmax(N,eps)
    Z = (dsigmoid2(I/N) )
  }
#  print(paste("ComputeIGrad = E = ", E, " I = ", I, "Z = ", Z))
  return(Z)
  
}


EmpiricalGrad <- function(f,x) {
  (f(x+eps/2)-f(x-eps/2))/eps
}

# Nuc Predictions and Gradients -----------------------------------------------

abModelBuildXsStrain <- function( Nucs, strain ) {
  Ni = dim(Nucs$Interactions)[1]
  Z = matrix(Nucs$Marginals[,strain], nc = Ni, nr = Ni)
  # initialize with product of marginals..
  Xs = Z * t(Z)
  rownames(Xs) = dimnames(Nucs$Interactions)[[1]]
  colnames(Xs) = dimnames(Nucs$Interactions)[[2]]
  Min =  pmax(Z + t(Z) - 1, 0)
  Max = pmin(Z,t(Z))
  A = Nucs$Interactions[,,strain]
  # Using only lower triangular as parameters
  A[upper.tri(A)] = t(A)[upper.tri(A)]
  Ip = A>0 & Xs > eps
  In = A<=0 & Xs > eps
  Pslope = pmax(eps,(Max-Xs)[Ip])
  Nslope = pmax(eps,(Xs - Min)[In])
  Xs[Ip] = Xs[Ip] + sigmoid2(A[Ip]/Pslope) * Pslope
  Xs[In] = Xs[In] + sigmoid2(A[In]/Nslope) * Nslope
  diag(Xs) = Nucs$Marginals[,strain]
  return(Xs * Nucs$Occ[strain])
}

#
# convert a Nuc to an array of pairwise marginal distributions.
#
abModelBuildXs <- function( Nucs ) {
  sapply(dimnames(Nucs$Occ)[[1]], function(i) abModelBuildXsStrain(Nucs, i), simplify = "array")
}

abModelBuildXsGradStrain <- function(Nucs, Strain) {
  Xs = abModelBuildXsStrain(Nucs, Strain )
  Ni = dim(Nucs$Interactions)[1]
  OGrad = array(1,dim=c(Ni,Ni, 1))
  MGrad = array(0,dim=c(Ni, Ni, Ni ))
  IGrad = array(0,dim=c(Ni, Ni, Ni*(Ni-1)/2 ))
  if( p.GradNames ) {
    dimnames(MGrad)[[1]] = dimnames(Nucs$Interactions)[[1]]
    dimnames(MGrad)[[2]] = dimnames(Nucs$Interactions)[[1]]
    dimnames(MGrad)[[3]] = dimnames(Nucs$Interactions)[[1]]
    dimnames(IGrad)[[1]] = dimnames(Nucs$Interactions)[[1]]
    dimnames(IGrad)[[2]] = dimnames(Nucs$Interactions)[[1]]
    dimnames(IGrad)[[3]] = matrix(1:(Ni*Ni), nc = Ni, nr = Ni) [lower.tri(Xs)]
    dimnames(OGrad)[[1]] = dimnames(Nucs$Interactions)[[1]]
    dimnames(OGrad)[[2]] = dimnames(Nucs$Interactions)[[1]]
    dimnames(OGrad)[[3]] = c("Occ")
  }
  
  OGrad[,,1] = Xs/Nucs$Occ[Strain]
  Z = matrix(Nucs$Marginals[,Strain], nc = Ni, nr = Ni)
  # initialize with product of marginals..
  E = Z * t(Z)
  Min =  pmax(Z + t(Z) - 1, 0)
  Max = pmin(Z,t(Z))
  A = Nucs$Interactions[,,Strain]
  # Using only lower triangular as parameters
  A[upper.tri(A)] = t(A)[upper.tri(A)]
  #  print(A)
  Ip = A>0
  In = A<=0
  Pslope = pmax(eps,(Max-E))
  Nslope = pmax(eps,(E - Min))
  
  for( i in 1:Ni ) {
    p = Nucs$Marginals[i,Strain]
    V = rep(0,Ni)
    V[[i]] = 1
    dZ = matrix(V, nc = Ni, nr = Ni)
    Ys = dZ * t(Z) + t(dZ) * Z
    J = (Ys > 0)
    dXs = matrix(0,nr = Ni, nc = Ni)
    dXs[J & Ip] = (Ys + Z*(A/Pslope)*dsigmoid2(A/Pslope) - Z*sigmoid2(A/Pslope))[J&Ip]
    dXs[J & Ip & p < Ys ] = (dXs - (A/Pslope)*dsigmoid2(A/Pslope) + sigmoid2(A/Pslope))[J&Ip & p < Ys ]
    
    dXs[J & In] = (Ys - Z*(A/Nslope)*dsigmoid2(A/Nslope) + Z*sigmoid2(A/Nslope))[J&In]
    dXs[J & In & p+ Ys > 1 ] = (dXs + (A/Nslope)*dsigmoid2(A/Pslope) - sigmoid2(A/Nslope))[J&In& p+Ys > 1 ]
    
    #    print(list(i,dXs))
    MGrad[,,i] = dXs 
    MGrad[i,i,i] = 1
  }
  MGrad = MGrad * Nucs$Occ[Strain]
  
  iMatrix = matrix(1:(Ni**2), nr = Ni, nc = Ni)
  inds = iMatrix[lower.tri(iMatrix)]
  for( k in 1:length(inds) ) {
    j = (inds[[k]]-1)%% Ni +1
    i = as.integer((inds[[k]] - j) / Ni + 1)
    #    print( c(inds[[k]], i, j, iMatrix[j,i]))
    a = Nucs$Interactions[j,i,Strain]
    p = Nucs$Marginals[i,Strain]
    q = Nucs$Marginals[j,Strain]
    #    print(list(i,j,k,p,q,a))
    #    print(CombineInteractionGradI(p,q,a))
    IGrad[j, i,k] = CombineInteractionGradI(p,q,a) 
  }
  IGrad = IGrad * Nucs$Occ[Strain]
  #  print(IGrad)
  for( i in 1:dim(IGrad)[[3]])
    IGrad[,,i][upper.tri(Xs)] = t(IGrad[,,i])[upper.tri(Xs)]
  #  print(IGrad)
  
  return(list(val = Xs, grad = abind(OGrad, MGrad,IGrad)))
}

abModelBuildXsGrad <- function(Nucs) {
  L = lapply(1:length(Nucs$Occ), function(i) abModelBuildXsGradStrain(Nucs, i))
  val = sapply(L, function(l) l$val, simplify = "array")
  grad = sapply(L, function(l) l$grad, simplify = "array")
  list(val = val, grad = grad)
}

abModelNucStatStrain <- function(Nucs, Strain ) {
  Ni = dim(Nucs$Marginals)[1]
  Xs = abModelBuildXsStrain( Nucs, Strain )
  Stat = array(0, dim=c(1+Ni+Ni*(Ni-1)/2))
  if( p.GradNames ) {
    MNames = dimnames(Nucs$Marginals)[[1]]
    names(Stat) = paste(Strain,
                                c( "Occ", MNames,
                                   sapply(MNames,function(n) paste(n,MNames))[lower.tri(Xs)] )
    )
  }
  
  Stat[1] = Nucs$Occ[Strain]
  Stat[2:(Ni+1)] = diag(Xs)
  Stat[(Ni+2):(Ni+1+Ni*(Ni-1)/2)] = Xs[lower.tri(Xs)]
  
  return(Stat)
}

abModelNucStat <- function( Nucs ) {
  do.call(c,lapply(dimnames(Nucs$Occ)[[1]], function(i) abModelNucStatStrain(Nucs, i)))
}

abModelNucStatStrainGrad <- function( Nucs, Strain ) {
  Ni = dim(Nucs$Marginals)[1]
  Xs = abModelBuildXsStrain( Nucs, Strain )
  
  G = abModelBuildXsGradStrain( Nucs, Strain )$grad
  Grad = array(0, dim=c( 1+Ni+Ni*(Ni-1)/2, dim(G)[3]))
  if( p.GradNames ) {
    MNames = dimnames(Nucs$Marginals)[[1]]
    dimnames(Grad)[[1]] = paste(Strain,
                                c( "Occ", MNames,
                                   sapply(MNames,function(n) paste(n,MNames))[lower.tri(Xs)] )
    )
    dimnames(Grad)[[2]] = paste(Strain,dimnames(G)[[3]])
  }
  
  A = matrix(nc=Ni, nr=Ni)
  Rows = matrix( 1:Ni, nc = Ni, nr=Ni )[lower.tri(A)]
  Cols = matrix( 1:Ni, nc = Ni, nr=Ni, byrow = TRUE )[lower.tri(A)]
  
  Grad[1,1] = 1
  for( i in 1:Ni)
    Grad[i+1,] = G[i,i,]
  for( i in 1:(Ni*(Ni-1)/2))
    Grad[Ni+1+i,] = G[Rows[i], Cols[i], ]
  
  return(Grad)
}


abModelNucStatGrad <- function( Nucs ) {
  blkdiag.wnames(lapply(dimnames(Nucs$Occ)[[1]], function(i) abModelNucStatStrainGrad(Nucs, i)))
}

# Model Predictions & Gradiants -------------------------------------------

abModelCompileDesign <- function( Design ) {
  
  IpPairs = data.frame(unique(Design$Exps[c("IP1","IP2")]))
  rownames(IpPairs) = sapply(1:length(IpPairs$IP1), 
                             function(i) paste0(IpPairs[i,"IP1"], "-", IpPairs[i,"IP2"]))
  
  ExpStrainsNames = as.character(Design$Exps[,"Strain"])
  ExpStrains = sapply(ExpStrainsNames, function(n) which(Design$Strains == n) )
  ExpIP1 = sapply(as.character(Design$Exps[,"IP1"]), function(n) which(Design$IPs == n))
  
  StrainExp = lapply(Design$Strains, function(n) Design$Exp$Strain == n)
  
  ExpPairsNames = sapply(1:length(rownames(Design$Exps)), 
                         function(i) paste0(Design$Exps[i,"IP1"], "-", Design$Exps[i,"IP2"]))
  ExpPairs = sapply(seq(nrow(Design$Exps)),
                    function(i) (ExpStrains[i]-1)*nrow(IpPairs) + which(rownames(IpPairs) == ExpPairsNames[[i]]))
  Single = lapply(Design$IPs,function(m) which(as.character(IpPairs$IP1) == m & is.na(IpPairs$IP2)))
  PairIP1 = lapply(Design$IPs,function(m) which(as.character(IpPairs$IP1) == m & !is.na(IpPairs$IP2)))
  PairIP2 = lapply(Design$IPs,function(m) which(as.character(IpPairs$IP2) == m & !is.na(IpPairs$IP2)))
  
  Ni = length(Design$Mods)
  Na = length(Design$IPs)

  YsParamMask = array(F, dim=c(Ni,Ni,Ni*Na))
  SingleYsMask = lapply(1:length(Design$IPs),
                        function(j) {
                          Mask = YsParamMask
                          for( i in 1:Ni) {
                            Mask[i,i,(i-1)*Na + j] = T
                          }
                          Mask
                        } )
  
  IP2YsMask = lapply(1:length(Design$IPs),
                     function(j) {
                       Mask = YsParamMask
                       for( i in 1:Ni)
                         Mask[,i,(i-1)*Na + j] = T
                       Mask
                     } )  
  IP1YsMask = lapply(1:length(Design$IPs),
                     function(j) {
                       Mask = YsParamMask
                       for( i in 1:Ni)
                         Mask[i,,(i-1)*Na + j] = T
                       Mask
                     } )  
  
 
  ModPairNames = sapply(Design$Mods,
                        function(n) paste(n,Design$Mods), simplify = "array")
  
  ModPairNames = ModPairNames[lower.tri(ModPairNames, diag = FALSE)] 

  NucStatNames = c( "Occ", Design$Mods, ModPairNames )
  
  Ns = length(Design$Strains)
  Nss = 1+Ni+Ni*(Ni-1)/2
  
  SingleOffset = 1+1:Ni
  ModDoubleOffset = matrix(NA, nr = Ni, nc = Ni)
  rownames(ModDoubleOffset) = Design$Mods
  colnames(ModDoubleOffset) = Design$Mods
  ModDoubleOffset[lower.tri(ModDoubleOffset, diag = FALSE)] = 1+Ni + 1:(Ni*(Ni-1)/2)
  ModDoubleOffset[upper.tri(ModDoubleOffset)] = t(ModDoubleOffset)[upper.tri(ModDoubleOffset)]
  
  IndexStart = 0
  
  ModSingleIndex = 2:(Ni+1)
  ModPairIndex = (Ni+2):Nss
  ModPairSourceIndex = lower.tri(ModDoubleOffset,diag = FALSE)
  
  BaseIndex = IndexStart+1
  
  ModZeroVector = rep(0,Ni)
  names(ModZeroVector) = Design$Mods
  
  MarginLower = sapply(Design$Strains, function(s) ModZeroVector, USE.NAMES = TRUE, simplify = FALSE )
  
  ModOneVector = ModZeroVector
  ModOneVector[] = 1
  
  MarginUpper = sapply(Design$Strains, function(s) ModOneVector, USE.NAMES = TRUE, simplify = FALSE )
  for( i in 1:nrow(Design$Constraints)) 
    MarginUpper[[as.character(Design$Constraints[i,"Strain"])]][[as.character(Design$Constraints[i,"Mod"])]] = Design$Constraints[i,"I"]
  
  list(
    Exps = rownames(Design$Exps),
    Strains = Design$Strains,
    Mods= Design$Mods,
    IPs = Design$IPs,
    IpPairs = IpPairs,
    ExpStrains = ExpStrains,
    StrainExp = StrainExp,
    ExpPairs = ExpPairs,
    ExpIP1 = ExpIP1,
    Single = Single,
    PairIP1 = PairIP1,
    PairIP2 = PairIP2,
    SingleYsMask = SingleYsMask,
    IP1YsMask = IP1YsMask,
    IP2YsMask = IP2YsMask,
    NucStatNames = NucStatNames,
    BaseIndex = BaseIndex,
    ModPairSourceIndex = ModPairSourceIndex,
    ModPairIndex = ModPairIndex,
    ModSingleIndex = ModSingleIndex,
    MarginUpper = MarginUpper,
    MarginLower = MarginLower
  )
}

printVector <- function(v) {
  l = lapply(1:length(v), function(i) paste0(names(v)[[i]],": ",v[[i]]))
  paste("[", do.call(paste,l),"]")
}

abModelModelStat <- function( Model, DesignCompile) {
  Np = length(rownames(DesignCompile$IpPairs))
  Nm = length(DesignCompile$NucStatNames)
  Ys = array(0,dim = c(Np, Nm), 
             dimnames = list( rownames(DesignCompile$IpPairs), DesignCompile$NucStatNames ))

  for( i in 1:Np) {
    ip1 = as.character(DesignCompile$IpPairs[i,"IP1"])
    ip2 = as.character(DesignCompile$IpPairs[i,"IP2"])
    v1 = Model$IPs[ip1,]
    Base = Model$Base[ip1]
    if( is.na(ip2) ) {
      Ys[i, DesignCompile$ModSingleIndex] = v1 
      Recovery = 1
    } else {
      v2 = Model$IPs[ip2,]
      Ys[i,DesignCompile$ModPairIndex] = (outer(v1,v2)+outer(v2,v1))[DesignCompile$ModPairSourceIndex]
      Ys[i,DesignCompile$ModSingleIndex] = Model$Base[ip1] * v2 + Model$Base[ip2] * v1 + v1*v2
      Recovery = Model$Recovery[[ip2]]
      Base = Base*Model$Base[ip2]
    }
    Ys[i,DesignCompile$BaseIndex] = Base
    Ys[i,] = Ys[i,] * Model$Total * Recovery * Model$Tags[[ip1]]
  }
  blkdiag.wnames(lapply(1:length(DesignCompile$Strains), 
                        function(s) { 
                          YY = Model$Ms[s]*Ys 
                          rownames(YY) = paste(DesignCompile$Strains[s], rownames(Ys))
                          colnames(YY) = paste(DesignCompile$Strains[s], colnames(Ys))
                          YY} ))
}

abModelNormalizeModel <- function( Model ) {
  meanMs = mean(Model$Ms)
  Model$Ms = Model$Ms/meanMs
  Model$Total = Model$Total * meanMs
  meanT = mean(Model$Tags)/0.01
  Model$Tags = Model$Tags/meanT
  Model$Total=Model$Total * meanT
  
  Model
}

abModelModelStatGrad <- function( Model, DesignCompile ) {
  Np = length(rownames(DesignCompile$IpPairs))
  Nm = length(DesignCompile$NucStatNames)
  Ni = length(DesignCompile$Mods)
  Ns = length(DesignCompile$Strains)
  Na = length(DesignCompile$IPs)
  Nparam = 1 + (Ni+3)*Na + Ns
  
  Ys = abModelModelStat( Model, DesignCompile )
  val = Ys
  Ys = Ys[1:Np,1:Nm]/Model$Ms[[1]]
  
  G = array(0, dim=c(Np, Nm, Nparam))
  if( p.GradNames ) {
    dimnames(G)[[1]] = rownames(DesignCompile$IpPairs)
    dimnames(G)[[2]] = DesignCompile$NucStatNames
    dimnames(G)[[3]] = c( "Total", # 1
                          paste("Ms", DesignCompile$Strains), # Ns
                          paste("Tag", DesignCompile$IPs ), # Na
                          paste("Recovery", DesignCompile$IPs), #Na
                          paste("Base", DesignCompile$IPs), # Na
                          sapply(DesignCompile$IPs, function(ip) paste(ip, DesignCompile$Mods) ) # Na * Ni
    )
  }
  
  IPParamIndex <- function(ip,mod) 1+Ns+3*Na+(mod-1)*Na + ip
  
  G[,,1] = Ys[,]/Model$Total
  for( i in 1:Np) {
    ip1 = as.character(DesignCompile$IpPairs[i,"IP1"])
    ip2 = as.character(DesignCompile$IpPairs[i,"IP2"])
    ip1i = which(DesignCompile$IPs == ip1)
    v1 = Model$IPs[ip1i,]
    
    G[i,DesignCompile$BaseIndex,1+Ns+2*Na+ip1i] = Ys[i,DesignCompile$BaseIndex] / Model$Base[[ip1i]]
    if( is.na(ip2) ) {
      # Ys[i, DesignCompile$ModSingleIndex] = v1 
      for( j in 1:Ni)
        G[i,DesignCompile$ModSingleIndex[j], IPParamIndex(ip1i,j) ] = Ys[i,DesignCompile$ModSingleIndex[j]]/v1[j]
      G[i,,1+Ns+ip1i] = Ys[i,] / Model$Tags[[ip1i]]
      
    } else {
      ip2i = which(DesignCompile$IPs == ip2)
      v2 = Model$IPs[ip2i,]
      # Ys[i,DesignCompile$ModPairIndex] = (outer(v1,v2)+outer(v2,v1))[DesignCompile$ModPairSourceIndex]
      m = 1
      for( j in 2:Ni)
        for( k in 1:(j-1) ) {
          
          mpi = DesignCompile$ModPairIndex[m]
          Val = Ys[i,mpi]/(v1[j]*v2[k] + v1[k]*v2[j])
          G[i,mpi, IPParamIndex(ip1i,j) ] = Val * v2[k]
          G[i,mpi, IPParamIndex(ip1i,k) ] = G[i,mpi, IPParamIndex(ip1i,k) ]+ Val * v2[j]
          G[i,mpi, IPParamIndex(ip2i,j) ] = G[i,mpi, IPParamIndex(ip2i,j) ]+ Val * v1[k]
          G[i,mpi, IPParamIndex(ip2i,k) ] = G[i,mpi, IPParamIndex(ip2i,k) ]+ Val * v1[j]
          m = m+1
        }
      # Ys[i,DesignCompile$ModSingleIndex] = Model$Base[ip1] * v2 + Model$Base[ip2] * v1 + v1*v2
      for( j in 1:Ni) {
        msi = DesignCompile$ModSingleIndex[j]
        Val = Ys[i,msi]/(Model$Base[ip1i]*v2[j] + Model$Base[ip2i]*v1[j] + v1[j]*v2[j])
        G[i,msi, IPParamIndex(ip1i,j)] = Val*(Model$Base[ip2i] + v2[j])
        G[i,msi, IPParamIndex(ip2i,j)] = G[i,msi, IPParamIndex(ip2i,j)]+Val*(Model$Base[ip1i] + v1[j])
        G[i,msi,1+Ns+2*Na+ip1i] = G[i,msi,1+Ns+2*Na+ip1i] + Val*v2[j]
        G[i,msi,1+Ns+2*Na+ip2i] = G[i,msi,1+Ns+2*Na+ip2i] + Val*v1[j]
      }
      G[i,DesignCompile$BaseIndex,1+Ns+2*Na+ip2i] = 
        G[i,DesignCompile$BaseIndex,1+Ns+2*Na+ip2i] +
        Ys[i,DesignCompile$BaseIndex] / Model$Base[[ip2i]]
      
      G[i,,1+Ns+Na+ip2i] = Ys[i,]/Model$Recovery[[ip2i]]
    }
    G[i,,1+Ns+ip1i] = Ys[i,] / Model$Tags[[ip1i]]
  }
  
  Gs = blkdiag3.wnames(lapply(1:Ns, 
                         function(s) {
                           Gs = G * Model$Ms[[s]]; 
                           Gs[,,1+s] = Ys;
                           dimnames(Gs)[[1]] = paste(DesignCompile$Strains[s], dimnames(G)[[1]])
                           dimnames(Gs)[[2]] = paste(DesignCompile$Strains[s], dimnames(G)[[2]])
                           Gs})
  )
  list(Ys = val, Gs = Gs)
}

abModelPredictInternal <- function( NucStat, ModelStat, DesignCompile ) {
  Spec = ModelStat[DesignCompile$ExpPairs,] %*% NucStat 
}

abModelPredict <- function( Model, Nuc, DesignCompile ) {
  ModelStat = abModelModelStat(Model, DesignCompile )
  
  abModelPredictInternal(abModelNucStat(Nuc), ModelStat, DesignCompile)
}

abModelPredictNucList <- function( Model, Nucs, DesignCompile ) {
  ModelStat = abModelModelStat(Model, DesignCompile )

  sapply(Nucs, function(n) abModelPredictInternal(abModelNucStat(n), ModelStat, DesignCompile))
}

abModelPredictNucGrad <- function( Model, Nucs, DesignCompile, MStat = NULL ){
  if( is.null(MStat) )
    MStat = abModelModelStat(Model, DesignCompile )
  
  NStat = abModelNucStat(Nucs)
  NStatG = abModelNucStatGrad(Nucs)
  
  Y = MStat[DesignCompile$ExpPairs,]
  Grad = Y %*% NStatG
  
  val =  Y %*% NStat
  list( val = val, grad = Grad)
}

abModelPredictModelGrad <- function( NStat, Model, DesignCompile ) {
  
#  print(NStat)
  tmp = abModelModelStatGrad(Model, DesignCompile)  
  MStat = tmp$Ys
  MStatG = tmp$Gs
  
  Grad = do.call(abind,c(lapply(seq(dim(MStatG)[3]), 
                              function(k) MStatG[DesignCompile$ExpPairs,,k] %*% NStat),
                 list(along = 3)))
  #              USE.NAMES = TRUE))
  if( p.GradNames) {
    dimnames(Grad)[[3]] = dimnames(MStatG)[[3]]
    dimnames(Grad)[[1]] = DesignCompile$Exps
  }
#  if( dim(Grad)[2] == 1)
#    Grad = Grad[,1,]
  val = MStat[DesignCompile$ExpPairs,] %*% NStat
  list( val = val, grad = Grad)
}

abModelPoissonLoss <- function( N, Y ) {
  if(any(is.nan(log(Y))))
    abort()
  return(-sum( N*log(Y) - Y ) )
}

abModelPoissonLossGrad <- function( N, Y, G ) {
#  print(list(dim(as.array(N)), dim(as.array(Y)), dim(G)))
#  print(G)
  T = N/Y-1
  T = T*G
#  print(list(N = dim(N), Y = dim(Y), G = dim(G), NY=dim(N/Y-1), T = dim(T)))
  -apply(T, c(2), sum)
}

abModelPoissonLossMultiGrad <- function( N, Y, G ) {
#  print(list(N = dim(N), Y = dim(Y), G = dim(G)))
#  -apply(rep((N/Y-1), dim(G)[[2]])*G, c(2), sum)
  -apply(rep((N/Y-1), dim(G)[3])*G, c(3), sum)
}

# Optimization ------------------------------------------------------------


abModelSelectRandomNucs <- function( Nucs, n = 1000) {
  # select a random subset...
  N = length(Nucs)
  if( n < N ) {
    I = sample(N, n)
  } else
    I = 1:N
  return(I)
}

abModelOptimizeModelMultiNucs <- function( Ns, Model, Nucs, 
                                           DesignCompile,
                                           n = 1000,
                                           Ids = abModelSelectRandomNucs(Nucs, n),
                                           Mask = NULL) {
  x = abModelPackModelParams(Model)
  Ni = length(x)
  Nn = length(Ids)
  if( is.null(Mask) )
    Mask = rep(1,Ni)
  Xlist = sapply(Ids, function(i) abModelNucStat(Nucs[[i]]))
  
  loss <- function(x) {
#    print(paste("Loss(",do.call(paste,as.list(x)),")"))
#    print("loss")
    M = abModelUnpackModelParams(Model, x)
    Am = abModelModelStat( M, DesignCompile )
    Ps = abModelPredictInternal(Xlist,Am,DesignCompile)
    abModelPoissonLoss(Ns[,Ids], Ps )
  }
 
  grad <- function(x) {
#    print(paste("Grad(",do.call(paste,as.list(x)),")"))
    M = abModelUnpackModelParams(Model,x)
    G = abModelPredictModelGrad(Xlist, M, DesignCompile )
    Grad = abModelPoissonLossMultiGrad(Ns[,Ids], G$val, G$grad)
    Grad * Mask
  }
  
  up = abModelUpperModelParams(Model)
  lower = abModelLowerModelParams(Model)
  opt = optim(x, fn = loss, gr = grad, lower = lower, upper = up, method = "L-BFGS-B")
  print(opt)
#  print(opt$value)
  x = opt$par
  M = abModelUnpackModelParams( Model, x)
  return(M)
}

#
# Optimize a single Nuc
#
abModelOptimizeNuc <- function( N, Nuc, DesignCompile, MStat,
                                 MaxIter, Mask) {
 

  loss <- function( x ) {
    NStat = abModelNucStat(abModelUnPackNucParams( Nuc, x))
    Y = abModelPredictInternal(NStat, MStat, DesignCompile )
    if(any(Y < 0)) {
      print(NStat)
      print(MStat)
      abort()
    }
    abModelPoissonLoss( N, Y)
  }
  grad <- function( x ) {
    
    G = abModelPredictNucGrad(Model, abModelUnPackNucParams( Nuc, x), DesignCompile, MStat )
    abModelPoissonLossGrad( N, G$val[,1], G$grad)*Mask
  }
  
  up = abModelUpperNucParams(Nuc, DesignCompile)
  lower = abModelLowerNucParams(Nuc, DesignCompile)
  x = abModelPackNucParams( Nuc )
  opt = optim(x, fn = loss, gr = grad, 
              lower = lower, upper = up, 
              method = "L-BFGS-B",
              control = list(maxit = MaxIter)#, factr = 1e10)
  )
#  print(opt)
#  print(opt$value)
  x = opt$par
  Nuc = abModelUnPackNucParams( Nuc, x)
  return(Nuc)
}

#
# Optimize a list of Nucs
#
abModelOptimizeNucs <- function( N, Model, Nucs, DesignCompile,
                                MaxIter = 100,
                                Mask = NULL) {
#  print("Opt Nucs")
  MStat = abModelModelStat(Model, DesignCompile)
  x = abModelPackNucParams( Nucs[[1]] )
  if( is.null(Mask) )
    Mask = rep(1,length(x))
  lapply(seq(length(Nucs)), function(n) abModelOptimizeNuc(N[,n], Nucs[[n]], DesignCompile, MStat, MaxIter, Mask))
}

abModelOptimizeStepWise <- function( N, Model, Nucs, DesignCompile, 
                                     NC = NULL,
                                     NucFirst = TRUE, MaxIter = 10, Delta = 1e-3 ) {
  l = abModelPoissonLoss( N, abModelPredict( Model, Nucs, DesignCompile ))
  print(l)
  Oldl = Inf
  Iter = 1
  while( Iter <= MaxIter && l+Delta < Oldl ) {
    print(paste("Stepwise iter = ", Iter, "l = ",l, "[", Oldl, "]"))
    Oldl = l
    if( Iter > 1 || NucFirst )
      Nucs = abModelOptimizeNucs(N, Model, Nucs, DesignCompile )
    Model =  abModelOptimizeModel(N, Model, Nucs, DesignCompile )
    l = abModelPoissonLoss( N, abModelPredict( Model, Nucs, DesignCompile))
    Iter = Iter+1
  }
  return(list(Model = Model, Nucs = Nucs))
}

# Setup functions ---------------------------------------------------------


FindBestModel <- function( N, Design, Model, Nucs, TotalRange = NULL,
                           MaxIter = 5, Delta = 1e-3 ) {
  DesignCompile = abModelCompileDesign( Design, Model, Nucs)
  
  BestModel = NULL
  BestLoss = Inf
  
  if( is.null(TotalRange) ) {
    q.99 = quantile(N, .99)
    q.50 = quantile(N, 80)
    maxval = q.50 + (q.99-q50)*100
    minval = max(10,q.50)
    TotalRange = seq(log10(minval), log10(maxval), length.out = 5)
  }
  
  for( p in TotalRange ) {
    tModel = Model
    tNucs = Nucs
    tModel$Total = 10 ** p
    print(paste("Power =",p))
    Iter = 1
    OldL = Inf
    l = 0
    tryCatch( {
      while( Iter <= MaxIter && l + Delta < OldL) {
        print(paste("FindBestModel Iter = ", Iter))
        opt = abModelOptimizeStepWise( N, tModel, tNucs, DesignCompile, Delta = Delta )
        tNucs = opt$Nucs
        tModel = opt$Model
        opt = abModelOptimize(N, tModel, tNucs, DesignCompile, NC)
        tNucs = opt$Nucs
        tModel = opt$Model
        OldL = l
        l = abModelPoissonLoss( N, abModelPredict( tModel, tNucs, DesignCompile ))
        print(paste("FindBestModel Iter ", Iter, " l = ", l, "[", OldL, "]"))
        Iter = Iter + 1
      }
    }, error = function(e) print(e) )
    
    Loss = abModelPoissonLoss( N, abModelPredict(tModel, tNucs, DesignCompile ))
    print(paste("Loss = ", Loss))
    if( Loss < BestLoss ) {
      BestModel = list( Design = Design, Model = tModel, Nucs = tNucs)
      BestLoss = Loss
      BestPow = p
    }
    #    DensityScatter(N,abModelPredict(Design, opt$Model, opt$Nucs), alpha=1, logscale = T, 
    #                   title = paste("Model",p, "loss", l) )
  }
  print(paste("BestLoss = ",BestLoss, "p =", BestPow))
  return(BestModel)   
}

abModelLossNucList <- function( Ns, Model, NucList, DesignCompile ) {
  
  Nn = length(NucList)
  A = abModelBuildYMod( Model, DesignCompile )
  Ps = sapply( NucList, function(x) 
    abModelPredict(Model, x, DesignCompile, A), simplify = "array")
  abModelPoissonLoss(Ns, Ps)
}

FindBestModelNucList <- function( Ns, Design, Model, NucList, TotalRange = NULL ) {
  
  DesignCompile = abModelCompileDesign( Design, Model, NucList[[1]])

  if( is.null(TotalRange) ) {
    q.99 = quantile(Ns, .99)
    q.50 = quantile(Ns, .5)
    maxval = q.50 + (q.99-q.50)*100
    minval = max(10,q.50)
    TotalRange = seq(log10(minval), log10(maxval), length.out = 4)
  }
  
  BestModel = NULL
  BestLoss = Inf
  
  for( p in TotalRange ) {
    tModel = Model
    tModel$Total = 10 ** p
    print(paste("Power =",p))
    tryCatch( {
      tModel = abModelOptimizeModelMultiNucs(Ns, tModel, NucList, DesignCompile )
      }, error = function(e) print(e) )
    
    Loss = abModelLossNucList( Ns, tModel, NucList, DesignCompile )
    print(paste("Loss = ", Loss))
    if( Loss < BestLoss ) {
      BestModel = tModel
      BestLoss = Loss
      BestPow = p
    }
  }
  print(paste("BestLoss = ",BestLoss, "p =", BestPow))
  return(BestModel)   
}

OptimizeNucList <- function( mat, mod, NucList, NucParamFN = NULL, MaxIter = 100 ) {
  
  DesignCompile = abModelCompileDesign( mod$Design, mod$Model, mod$Nucs )
  
  Nnuc = length(NucList)
  A = abModelBuildYMod(mod$Model, DesignCompile)
  
  for( i in 1:Nnuc) {
    N = mat[,i]
    Nucs = NucList[[i]]
    NucList[[i]] = abModelOptimizeNucs(N, mod$Model, Nucs, DesignCompile, MaxIter, A = A )
    if( i %% 500 == 0 )
    {
      print(i)
      if( !is.null(NucParamFN))
        saveRDS(NucList[1:i],NucParamFN)
    }
  }
  if( !is.null(NucParamFN))
    saveRDS(NucList,NucParamFN)
#  print(length(NucList))
  return(NucList)
}


InitializeNucList <- function(Marginals, Pairwise, mod ) {
  NucList = list()
  mList = sapply( mod$Design$Strains, function( strain) paste0(strain, "-", mod$Design$Mods))
  iList = sapply( mod$Design$Strains, function( strain) paste0(strain, "-", sapply(mod$Design$Mods, function(m) paste0(m,"-",mod$Design$Mods))))
  
  nS = length(mod$Design$Strains)
  nM = length(mod$Design$Mods)
  pw = array(0,dim=c(nS, nM, nM))
  
  rowIndex = rep( rep(1:nM, nM), nS)
  colIndex = rep( rep(1:nM, each = nM))
  
  for( i in dim(Marginals)[[2]]) {
    Nucs = mod$Nucs
    Nucs$lable = i
    Nucs$Occ[] = 1
    Nucs$Marginals = Marginals[mList, i]
    Nucs$Marginals[Nucs$Marginals>1] = 1
    Nucs$Marginals[Nucs$Marginals<0] = 0
  
    
    pw = Pairwise[iList, i]
    Nucs$Interactions[] = ExtractInteraction(Nucs$Marginals[rowIndex],Nucs$Marginals[colIndex], pw)
    
    NucList[[i]] = Nucs
  }
  
  NucList
}


NucListIterativelyOptimizeModel <- function( mat, mod, NucList, n = 1000, m = 5000,
                                             MaxIter = 5, Delta = 1e-3,
                                             modFN = NULL ) {

  DesignCompile = abModelCompileDesign( mod$Design, mod$Model, mod$Nucs )

  Iter = 1
  Model = mod$Model
  OldL = 0
  L = - 2*Delta
  
  Test = abModelSelectRandomNucs( mat, NucList, m )
  BestLoss = abModelLossNucList( Test$Ns, mod$Design, Model, Test$NucList, DesignCompile)
  BestModel = Model
  print(paste("Iter ", 0, " TestL = ", BestLoss))
  
  while( Iter <= MaxIter && (OldL - L) > Delta ) {
    Model = BestModel
    mod$Model = Model
    
    Select = abModelSelectRandomNucs( mat, NucList, n )
    OldL = abModelLossNucList( Select$Ns, mod$Design, Model, Select$NucList, DesignCompile)
    print(paste("Iter ", Iter, " OldL = ", OldL))

    # Try recalibrating model
    RecalModel = RecalibrateModelLoadingFactors( Select$Ns, mod, Select$NucList)
    RecalL = abModelLossNucList( Select$Ns, mod$Design, RecalModel, Select$NucList, DesignCompile )
    if( RecalL < OldL ) {
      print(paste("After recalibration L = ", RecalL ))        
      mod$Model = RecalModel
    }
    Model = mod$Model
    
    print("optimize Nucs")
    Select$NucList = OptimizeNucList(Select$Ns, mod, Select$NucList, MaxIter = 10)
#    OldL = abModelLossNucList( Select$Ns, mod$Design, Model, Select$NucList, DesignCompile )
#    print(paste("Iter ", Iter, " OldL = ", OldL))
    L = abModelLossNucList( Select$Ns, mod$Design, mod$Model, Select$NucList, DesignCompile)
    print(paste("After optimizing Nucs L = ", L))
    # optimize loading parameters
    Model = OptimizeModelLoadingFactors(Select$Ns, mod, Select$NucList)
    L = abModelLossNucList( Select$Ns, mod$Design, Model, Select$NucList, DesignCompile )
    mod$Model = Model
    print(paste("After optimizing loading parameters L = ", L))

    print("optimize Model")
    tryCatch( {
      Model = abModelOptimizeModelMultiNucs(Select$Ns,mod$Design, Model, Select$NucList)
      mod$Model = Model
      L = abModelLossNucList( Select$Ns, mod$Design, Model, Select$NucList, DesignCompile )
    }, error = function(e) print(e) ) 
    print("improving test nucs")
    Test$NucList = OptimizeNucList(Test$Ns, mod, Test$NucList, MaxIter = 10)
    
    TestL = abModelLossNucList( Test$Ns, mod$Design, Model, Test$NucList, DesignCompile )
    print(paste("Iter ", Iter, " L = ", L, " TestL = ", TestL))
    
    Iter = Iter + 1
    if( TestL < BestLoss ) {
      BestModel = Model
      BestLoss = TestL
      print(paste("New Best Loss = ", BestLoss))
      if( !is.null(modFN) ) 
        saveRDS(mod, modFN)
    }
  }
  return(BestModel)
}

# Test model --------------------------------------------------------------

Test = 0
if(Test) {
  
  Test.Experiments = c( "WT-K4me3-Input", 
                        "WT-K18ac-Input", 
                        "WT-K12ac-Input",
                        "WT-H3-K12ac", 
                        "WT-K12ac-K18ac", 
                        "WT-K4me3-K18ac",
                        "WT-K4me3-K12ac",
                        "WT-K4me3-K4me3",
                        "K18R-K18ac-Input",
                        "K18R-K4me3-Input",
                        "K18R-K4me3-K12ac",
                        "K18R-K12ac-Input")
  
  Test.Design = CreateDesign( Test.Experiments)
  Test.Design$Constraints =  data.frame( 
    Strain = c("K18R"), 
    Mod = c("K18ac"), 
    I = c(0) 
  )
  
  
  Test.Compile = abModelCompileDesign(Test.Design)
  
  Test.Model = CreateModel(Test.Design)
  
  Test.Model$IPs[,] = eps
  Test.Model$IPs[1,1] = .15
  Test.Model$IPs[2,2] = .1
  Test.Model$IPs[3,3] = .05
  Test.Model$Base[] = 1e-3
  Test.Model$Recovery[] = .1
  
  upperMarginals = unlist(Test.Compile$MarginUpper)
  lowerMarginals = unlist(Test.Compile$MarginLower)
  Test.NucNum = 1000
  Test.Nucs = lapply(1:Test.NucNum, function(i) { 
    temp = CreateNuc(Test.Design,i)
    temp$Occ[] = rbeta(length(temp$Occ), 5, 1)
    tempM = runif(length(Test.Design$Mods))
    temp$Marginals[] = pmin( rep(tempM, length(Test.Design$Strains)), upperMarginals )
    temp$Interactions[,,1] = rnorm(length(temp$Interactions[,,1]))
    temp$Interactions[,,1][upper.tri(temp$Interactions[,,1])] = t(temp$Interactions[,,1])[upper.tri(temp$Interactions[,,1])]
    
    for( s in Test.Compile$Strains)
      temp$Interactions[,,s] = temp$Interactions[,,1]
    
    temp
  })
  names(Test.Nucs) = lapply(Test.Nucs, function(n) n$Label)
  
  Test.mod = list( Design = Test.Design, Model = Test.Model, Nucs = Test.Nucs )
  
  Test.N1 = c(9526, 92629,  1030,    82,  1018) 
  Test.N2 = c(61739, 51109, 1450,   320,  1031) 
  Test.Noise=0.25
  Test.N = do.call(cbind,lapply(Test.Nucs, function(n) rpois(length(Test.Experiments), 
                                                             rgamma(length(Test.Experiments),1/(Test.Noise**2),scale=Test.Noise**2)*
                                                               abModelPredict(Test.Model, n, Test.Compile))))
  colnames(Test.N) = names(Test.Nucs)
  
  Test.ModelFromTrue = abModelOptimizeModelMultiNucs(Test.N, Test.Model, Test.Nucs, Test.Compile)
  Test.ModelFromScratch = abModelOptimizeModelMultiNucs(Test.N, CreateModel(Test.Design), Test.Nucs, Test.Compile)
  Test.PredictLearned = do.call(cbind,lapply(Test.Nucs, function(n) abModelPredict(Test.ModelFromScratch, n, Test.Compile)))
  Test.Predict = do.call(cbind,lapply(Test.Nucs, function(n) abModelPredict(Test.Model, n, Test.Compile)))
  DensityScatter(Test.N, Test.PredictLearned,logscale=T)
  DensityScatter(Test.Predict, Test.PredictLearned,logscale=T)
}