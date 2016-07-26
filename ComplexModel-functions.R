setwd("~/Dropbox/coChIP-scripts/")
library(abind)
library(compiler)
library(pracma)

source("CoChip-Functions.R")
source("NucAtlas.R")
source("DensityScatter.R")

DataDir = "~/Google Drive/CoChIPAnalysis"
FigureDir = "~/Google Drive/CoChIPAnalysis/ModelFigures"
TrackDir = "~/Google Drive/CoChIPAnalysis/ModelTracks"
source("ComplexModel.R")
source("PairwiseModel.R")

## 
## set to TRUE for debuging
##
p.GradNames = TRUE
p.Run = 0


LoadWeinerNucs <- function( TargetMods = c("H3K18ac", "H3K4ac", "H3K36me3" , "H3K79me3","H3K4me3", "H3K56ac" ) ) {
  mod.nucs = ReadHistoneModeAtlas("Weiner-HisMod.csv", NucRegions, TargetMods )
  m = mcols(mod.nucs)
  weiner.nucs <<- do.call(rbind,lapply(TargetMods, function(mod) as.vector(exp(as.matrix(m["Input"]) + as.matrix(m[mod])))))
}

BuildExpFromFilenames <- function( FileNames, Mods = NULL, baseMod = Test.mod ) {
  # split experiments to Strain-IP1-{IP2,Input}
  #
  mm = matrix(unlist(strsplit(FileNames,"-")), nc = 3, byrow = T )

  Strains = unique(mm[,1])
  IPs = unique(c(mm[,2], mm[,3]))
  IPs = IPs[IPs != "Input"]
  
  if( is.null( Mods ) )
  {
    Mods = IPs
  }
  
  Exps = data.frame( Strain = mm[,1], IP1 = mm[,2], IP2 = mm[,3], row.names = FileNames )
  Exps[Exps$IP2 == "Input","IP2"] = NA
  
  Design = baseMod$Design
  Design$Mods = Mods
  Design$IPs = IPs
  Design$Strains = Strains
  Design$Exps = Exps
  
  #Note that we leave the constraints to be manually constructed by the user..
  
  MIPs = matrix(eps, nr = length(IPs), nc = length(Mods)+1)
  rownames(MIPs) = IPs
  colnames(MIPs) = c(Base, Mods)
  for( r in rownames(MIPs)) 
    if( r %in% rownames(baseMod$Model$IPs) )
      for( c in colnames(MIPs) )
        if( c %in% colnames(baseMod$Model$IPs))
          MIPs[r,c] = baseMod$Model$IPs[r,c]
  
  Recovery = matrix(1, nc = 1, nr = length(IPs), dimnames = list(IPs, c("l")))
  for( r in IPs )
    if( r %in% rownames(baseMod$Model$Recovery) )
      Recovery[r,"l"] = baseMod$Model$Recovery[r,"l"] 
    
  Mm = unique(mm[,1:2])
  Ms = data.frame( Strain = Mm[,1], IP1 = Mm[,2], M = rep(1,length(Mm[,1])))
  
  Model = baseMod$Model
  Model$IPs = MIPs
  Model$Recovery = Recovery
  Model$Ms = Ms
  
  Occ = rep(1,length(Strains))
  names(Occ) = Strains
  
  Marginals = matrix(0.5, nc = length(Strains), nr = length(Mods))
  rownames(Marginals) = Mods
  colnames(Marginals) = Strains
  
  Interactions = array(1, dim = c(length(Mods), length(Mods), length(Strains)))
  dimnames(Interactions)[[1]] = Mods
  dimnames(Interactions)[[2]] = Mods
  dimnames(Interactions)[[3]] = Strains

  Nucs = mod$Nucs
  Nucs$Occ = Occ
  Nucs$Marginals = Marginals
  Nucs$Interactions = Interactions
  
  return( list( Design = Design, Model = Model, Nucs = Nucs ) )  
}


SimulateNucData <- function( Marginals, mod, Oalpha = 5, Obeta = 2, Imu = 0, Isigma = 1) {

  N = dim(nuc.mat)[[2]]
  M = dim(nuc.mat)[[1]]
  Occ = rbeta(N, Oalpha, Obeta)
  Marg = matrix(nr=M, nc=N)
  for( i in 1:M )
    Marg[i,] = Occ * Marginals[i,] / max(Marginals[i,])
  
  P = M*(M+1)/2
  Indexes = matrix(0, nc = M, nr = M )
  Indexes[lower.tri(Indexes,diag = TRUE)] = 1:P
  
  Pairwise = matrix(nr = P, nc = N)
  for( p in 1:P ) {
    i = which(sapply(1:M, function(r) any(Indexes[r,] == p)) > 0 )
    j = which(sapply(1:M, function(c) r*any(Indexes[,c] == p)) > 0 )
    print(c(i,j,p))
    Int = rnorm(N, Imu,Isigma)
    if( i == j )
      Pairwise[p,] = Marg[i,] 
    else
      Pairwise[p,] = CombineInteraction(Marg[i,], Marg[j,], Int)
  }
  L = length(mod$Design$Exps$Strain)
  for( l in 1:L) {
    IP1 = mod$Design$Exps[l,"IP1"]
    IP2 = mod$Design$Exps[l,"IP2"]
  }
}

# Testing Gradients -------------------------------------------------------

TestGradVec <- function(x,fn,grad) {
  
  l = do.call(fn,list(x))
  g = do.call(grad,list(x))
  g1 = g
  #  print(l)
  for( i in seq_along(x))
  {
    xt = x 
    xt[[i]] = x[[i]]+eps
    lt = do.call(fn,list(xt))
    g1[,i] = (lt -l)/eps
    #    print(list(i,lt))
  }
  Diff = g -g1
  NormDiff = abs(Diff)/pmax(abs(g),abs(g1))
  NormDiff[is.na(NormDiff)] = 0
  return(list(Computed = g,Empirical = g1, Diff = Diff, NormDiff = NormDiff))
}

TestGrad <- function(x,fn,grad) {
  
  l = do.call(fn,list(x))
  g = do.call(grad,list(x))
  g1 = g
  for( i in seq_along(x))
  {
    xt = x 
    #    xt[[i]] = x[[i]]+eps
    #    lt = do.call(fn,list(xt))
    #    g1[[i]] = (lt -l)/eps
    
    xt[[i]] = x[[i]]-eps
    lt = do.call(fn,list(xt))
    g1[[i]] = (l -lt)/eps
    
  }
  Diff = g -g1
  NormDiff = abs(Diff)/pmax(abs(g),abs(g1))
  NormDiff[is.na(NormDiff)] = 0
  return(list(Computed = g,Empirical = g1, Diff = Diff, NormDiff = NormDiff))
}

abModelOptimizeNucsTestGrad <- function( N, Model, Nucs, DesignCompile ) {
  x = abModelPackNucParams( Nucs )
  MStat = abModelModelStat(Model, DesignCompile)
  
  loss <- function( x ) {
    Nucs = abModelUnPackNucParams( Nucs, x)
    Y = abModelPredict(Model, Nucs, DesignCompile )
    if( any(is.na(Y)) ) print(list(Y,x,Model,Nucs))
    abModelPoissonLoss( N, Y)
  }
  grad <- function( x ) {
    Nucs = abModelUnPackNucParams( Nucs, x)
    G = abModelPredictNucGrad(Model, Nucs, DesignCompile, MStat )
    abModelPoissonLossMultiGrad( N, G$val, G$grad)
  }
  
  lossV <- function( x ) {
    Nucs = abModelUnPackNucParams( Nucs, x)
    Y = abModelPredict(Model, Nucs, DesignCompile )
  }
  gradV <- function( x ) {
    Nucs = abModelUnPackNucParams( Nucs, x)
    G = abModelPredictNucGrad(Model, Nucs, DesignCompile )
    G$grad
  }
  
  lossX <- function( x ) {
    Nucs = abModelUnPackNucParams( Nucs, x)
    abModelNucStat(Nucs)
  }
  gradX <- function( x ) {
    Nucs = abModelUnPackNucParams( Nucs, x)
    abModelNucStatGrad(Nucs)
  }
  
  Nucs = abModelUnPackNucParams( Nucs, x)
  #  print(Nucs)
  #  print(abModelBuildXs(Nucs))
  print("Xs Grad")
  M = TestGradVec(x,lossX,gradX)
#  print(M)
  C = which(colSums(M$NormDiff) > 0.001)
  R = which(rowSums(M$NormDiff) > 0.001)
  print(list(Computed = M$Computed[R,C], Empirical = M$Empirical[R,C], Diff = M$Diff[R,C], NormDiff = M$NormDiff[R,C]))
  print(list(max(abs(M$Diff)),max(M$NormDiff)))
  print("Ys Grad")
  M = TestGradVec(x,lossV,gradV)
  #  print(M)
  C = which(colSums(M$NormDiff) > 0.001)
  R = which(rowSums(M$NormDiff) > 0.001)
  print(list(Rows = rownames(M$Computed)[R], Cols = colnames(M$Comptued)[C]))
  print(list(Computed = M$Computed[R,C], Empirical = M$Empirical[R,C], Diff = M$Diff[R,C], NormDiff = M$NormDiff[R,C]))
  print(list(max(abs(M$Diff)),max(M$NormDiff)))
  
  M = do.call(cbind,TestGrad(x,loss,grad))
  print(M)
}

abModelOptimizeModelTestGrad <- function( N, Model, Nucs, DesignCompile, Threshold = 0.001 ) {
  doMod <- function(x) {
    Mod <- abModelUnpackModelParams(Model, x)
    Y = abModelModelStat(Mod, DesignCompile)
    as.vector(Y)
  }
  
  doModGrad <- function(x) {
    Mod <- abModelUnpackModelParams(Model, x)
    ll = abModelModelStatGrad(Mod, DesignCompile)
    G = sapply( seq(dim(ll$Gs)[3]), function(k) as.vector(ll$Gs[,,k]))
    dimnames(G)[[2]] = dimnames(ll$Gs)[[3]]
    dimnames(G)[[1]] = sapply(dimnames(ll$Gs)[[2]], function(x) paste(x,dimnames(ll$Gs)[[1]]))
    G
  }
  
  NStat = abModelNucStat( Nucs )
  loss <- function( x ) {
    Mod = abModelUnpackModelParams( Model, x)
    Y = abModelPredict(Mod, Nucs, DesignCompile)
    abModelPoissonLoss( N, Y)
  }
  grad <- function( x ) {
    Mod = abModelUnpackModelParams( Model, x)
    G = abModelPredictModelGrad(NStat, Mod, DesignCompile )
    abModelPoissonLossMultiGrad( N, G$val, G$grad)
  }
  
  x = abModelPackModelParams(Model)
  M = TestGradVec(x, doMod, doModGrad)
#  print(M)
  C = which(colSums(M$NormDiff) > 0.001)
  R = which(rowSums(M$NormDiff) > 0.001)
  print(list(Computed = M$Computed[R,C], Empirical = M$Empirical[R,C], Diff = M$Diff[R,C], NormDiff = M$NormDiff[R,C]))
#  print(list(max(abs(M$Diff)),max(M$NormDiff)))
#  return(0)
#  S = t(do.call(rbind,lapply(M,as.vector.wnames)))
#  I = S[,"Diff"] != 0 & S[,"NormDiff"] > Threshold
#  if( any(I)) {
#    tS = S[I,]
#    #    rownames(tS) = names(as.vector.wnames(M$Computed))[I]
#    #    print(names(as.vector.wnames(M$Computed))[I])
#    print(tS)
#  }
  TestG = TestGrad(x,loss,grad)
  M = do.call(rbind,TestG)
  print(t(M))
}

# Scripts -----------------------------------------------------------------



if( 0 ) {
  # deal with count file
  cnts = read.csv(paste0(DataDir,"/HisMut-all-counts.csv"), row.names = 1)
  mod = BuildExpFromFilenames( rownames( cnts ))
  Design = mod$Design
  Model = mod$Model
  Nucs = mod$Nucs
  N = cnts[,"x"]
  
  mod = FindBestModel( N, Design, Model, Nucs, seq(5.5,8.7,by=.15))  
}

FNext = "v2-recal"
modFN = paste0(DataDir,"/NucModel-", FNext, ".rdata")
NucParamFN = paste0(DataDir,"/NucParam-", FNext, ".rdata")
NucParamInitFN = paste0(DataDir,"/NucInitParam.rdata")
DistFN = paste0(DataDir,"/NucDistribution-", FNext, ".rdata")
TrackFN = paste0(DataDir,"HisMut/Tracks/DistTrack-")



if( !exists("SGD") ) {
  print("loading atlas")
  GetSGDandNucs()
}
#source("Main-Functions.R")
#InputName = "-input"

if(0) {
  
  NucFile = "nucdata"
  NucFN = paste0(DataDir,"/",NucFile,".rdata")
  
  print("reading Nuc counts")
  ## get the data
  nuc.data =readRDS(NucFN)
  MaxE = length(nuc.data$counts)
  MaxN = dim(nuc.data$mat)[[2]]
  #  MaxE = 20
  #  MaxN = 2000
  nuc.counts = nuc.data$counts[1:MaxE]
  nuc.mat = nuc.data$mat[1:MaxE,1:MaxN]
  Nnuc = dim(nuc.mat)[[2]]
  
  # get rid of outliers
  tot = colMeans(nuc.mat)
  t = quantile(tot, .95)
  m = mean(tot[tot <= t])
  s = sqrt(var(tot[tot <= t]))
  t.high = m + 3*s
  t.low = m - 2*s
  nuc.I = tot < t.high & tot > t.low  
  
  names(nuc.counts) = rownames(nuc.mat)
}

if( 0 ) {
  mod = readRDS(modFN)
  mod$Model$IPs[,"l"] = mod$Design$IP2LoadingFactor[,rownames(mod$Model$IPs)]
  mod$DesignCompile = abModelCompileDesign( mod$Design, mod$Model, mod$Nucs)
  saveRDS(mod, modFN)
}

if( 0 ) {
  mod = BuildExpFromFilenames( rownames( nuc.mat ) )
  DesignCompile = abModelCompileDesign( mod$Design, mod$Model, mod$Nucs )
  mod = FindBestModel( nuc.counts, mod$Design, mod$Model, mod$Nucs, seq(6,7,by=1))
  N = rep(1,MaxE)
  abModelOptimizeNucsTestGrad(N, mod$Design, mod$Model, mod$Nucs, DesignCompile )
  abModelOptimizeModelTestGrad(N, mod$Design, mod$Model, mod$Nucs, DesignCompile, Threshold = 0.01 )
  ## end test code
}

if( p.Run > 0 ) {
  TargetMods = c("H3K18ac", "H3K4me3", "H3K36me3" , "H3K79me3" )
  if( !exists( "WeinerInputs" ) ) {
    mod.nucs = ReadHistoneModeAtlas("Weiner-HisMod.csv", NucRegions, TargetMods )
    m = mcols(mod.nucs)
    weiner.nucs = do.call(rbind,lapply(TargetMods, function(mod) as.vector(exp(as.matrix(m["Input"]) + as.matrix(m[mod])))))
    rownames(weiner.nucs) = TargetMods
    WeinerInputs =  as.vector(exp(as.matrix(m["Input"])))
    t = quantile(WeinerInputs,0.95,na.rm=T)
    WeinerInputs = WeinerInputs/t
  }
  
  
  print("Initializing")
  mod = BuildExpFromFilenames( rownames( nuc.mat ) )
  
  if( file.exists(NucParamInitFN) ) {
    print("Reading initial nucleosome model")
    NucList = readRDS( NucParamInitFN )
    NucList = NucList[1:MaxN]
  } else {
    print("Building initial nucleosome model")
    NucList = InitNucList( nuc.mat, mod, Inputs = WeinerInputs )
    saveRDS(NucList,NucParamInitFN)
  }
}

if( 0 ) {
  weiner.99 = quantile(WeinerInputs, .99, na.rm = T)
  weiner.95 = quantile(WeinerInputs, .95, na.rm = T)
  GoodNucs = WeinerInputs[1:MaxN] < (weiner.95 + (weiner.99 - weiner.95)*10)
  NucList = NucList[GoodNucs]
}

if( p.Run > 0 ) 
  NucIds = as.integer(sapply(NucList, function(x) x$Label))

if(p.Run == 1) { 
  if( file.exists(modFN) ) {
    print("Reading initial model")
    mod = readRDS(modFN) 
  }else {
    print("Optimizing initial model")
    Select = abModelSelectRandomNucs( nuc.mat, NucList, n = 500 )
    tModel = FindBestModelNucList(Select$Ns,mod$Design, mod$Model, Select$NucList)
    mod$Model = tModel
    saveRDS(mod, modFN)
  }    
  
  if( file.exists(NucParamFN) ) {
    print("Read nucleosome model")
    NucList = readRDS(NucParamFN)
  }
  
  if( 1 ) {
    print("Optimizing model")
    tModel = NucListIterativelyOptimizeModel( nuc.mat, mod, NucList, n = 1000, m = 5000, MaxIter = 10, modFN = modFN)
    mod$Model = tModel
    saveRDS(mod, modFN)
  }
  
}

if( p.Run == 2 ) {
  mod = readRDS(modFN) 
  mod$DesignCompile = abModelCompileDesign(mod$Design, mod$Model, mod$Nucs)
  
  if( file.exists(NucParamFN)) {
    print("reading Nuc parameters")
    NucList = readRDS(NucParamFN)
  } else {
    print("Optimizing nucleosomes")
    NucList = OptimizeNucList(nuc.mat, mod, NucList, NucParamFN = NucParamFN )
    saveRDS(NucList, NucParamFN)
  }
  
  
  
  
  if( file.exists(DistFN) ) {
    print("Reading nucleosome distribution")
    Dist = readRDS(DistFN)
  } else {
    print("Computing nucleosome distribution")
    Dist = BuildNucDist( NucList )
    saveRDS(Dist, DistFN)
  }
  
  
  # compute errors
  print("Computing prediction errors")
  Preds = ComputePredictions( mod, NucList )
  
  m = "WT-H3K18ac-Input"
  DensityScatter(nuc.mat[m,NucIds],Preds[m,],  alpha = .2, coordeq = F, threshold = .95, xlabel="Obs", ylabel="Pred", title = m, diagonal = T)
  
  Interactions = sapply(NucList,function(n) n$Interactions$a, simplify = "array")
  colnames(Interactions) = NucIds
  rownames(Interactions) = sapply(1:length(rownames(mod$Nucs$Interactions)), function(i) paste0(mod$Nucs$Interactions[i,"IP1"],"-",mod$Nucs$Interactions[i,"IP2"]))
  
  ComputeRatio <- function(m) {
    B = nuc.mat[m,NucIds] > 10
    if( any(B) ) {
      data = data.frame(x = nuc.mat[m,NucIds][B], y = Preds[m, ][B])
      lm = lm(y~x-1, data = data)
      return(coef(lm))
    } else
      return(NA)
  }
  
  Ratios = sapply(rownames(nuc.mat), ComputeRatio )
  RatioN = sapply(rownames(nuc.mat), function(m) length(which(nuc.mat[m,NucIds] > 10)))
  names(Ratios) = rownames(nuc.mat)
  RRatios = Ratios[RatioN > 1000]
  
  NPred = dim(Preds)[[1]]*dim(Preds)[[2]]
  PredError = nuc.mat[,NucIds] - Preds
  S = sample(1:NPred,1e5)
  DensityScatter( nuc.mat[,NucIds][S], Preds[S], coordeq = F, logscale = T, diagonal = T, ylim = c(1,100), smooth = T)
  hist(PredError, breaks="FD")
  print(list(mean(PredError[S]), sqrt(var(PredError[S]))))
  hist(colSums(abs(PredError)), breaks="FD")
  hist(apply(abs(PredError),c(2), median), breaks="FD")
  hist(apply(abs(PredError),c(1), median), breaks="FD")
  hist(apply(PredError,c(2), median), breaks="FD")
  hist(apply(PredError,c(2), mean), breaks="FD")
  hist(apply(PredError,c(1), mean), breaks="FD")
  
  
  # Plot marginals
  plots = list()
  plots[["Occ"]] = DensityScatter(WeinerInputs[NucIds], DistGetMarginal(Dist, "Occ"), title = "Occupancy", 
                                  xlabel = "observed", ylabel = "estimated", logscale = F,  alpha = .2, coordeq = F, threshold = .95)
  
  for( m in names(mod$Nucs$Marginals))
    plots[[m]] = DensityScatter(nuc.mat[paste0("WT-", m, "-Input"),NucIds], DistGetMarginal(Dist, m), title = m,
                                xlabel = "observed (WT)", ylabel = "estimated", logscale = F,  alpha = .2, coordeq = F)
  
  for( m in rownames(mod$Nucs$Interactions) )
    plots[[m]] = DensityScatter(nuc.mat[paste0("WT-", m),NucIds], Dist[,m], title = m,
                                xlabel = "observed (WT)", ylabel = "estimated", logscale = F,  alpha = .2, coordeq = F)
  
}

if( p.Run == "wtPreds" ) {
  ## optimize to WT results
  wt.sub = SubsetModel(nuc.mat[,NucIds], mod, "WT")
  wt.mod = wt.sub$mod
  wt.mat = wt.sub$mat
  wtPredFN = paste0(DataDir,"/wtPreds.rdata")
  if( file.exists(wtPredFN) ) {
    print("Reading WT distribution")
    wt.Data = readRDS(wtPredFN)
    wt.NucList = wt.Data$NucList
    wt.Preds = wt.Data$Preds
    wt.Dist = wt.Data$Dist
  } else {
    print("Optimizing WT nucs")
    wt.NucList = OptimizeNucList( wt.mat, wt.mod, NucList)
    print("Computing WT predictions")
    wt.Preds =  ComputePredictions( wt.mod, wt.NucList )
    print("Computing WT distribution")
    wt.Dist = BuildNucDist( wt.NucList )  
    print("Saving WT distribution")
    saveRDS(list(NucList = wt.NucList, Preds= wt.Preds, Dist = wt.Dist), wtPredFN)
  }
  
}

if( p.Run == "altPreds" ) {
  mod = readRDS(modFN) 
  mod$DesignCompile = abModelCompileDesign(mod$Design, mod$Model, mod$Nucs)
  
  if( file.exists(NucParamFN)) {
    print("reading Nuc parameters")
    NucList = readRDS(NucParamFN)
  }
  print("orig")
  Preds.orig =  ComputePredictions( mod,NucList )
  print("min")
  Preds.min =  ComputePredictions( mod, lapply(NucList, function(n) { n$Interactions[,"a"] = 0; n} ) )
  print("exp")
  Preds.exp =  ComputePredictions( mod, lapply(NucList, function(n) { n$Interactions[,"a"] = 1; n} ) )
  print("max")
  Preds.max =  ComputePredictions( mod, lapply(NucList, function(n) { n$Interactions[,"a"] = 1e6; n} ) )
  altPreds = list(orig = Preds.orig, min = Preds.min, exp = Preds.exp, max = Preds.max) 
  print("saving data")
  saveRDS(altPreds, paste0(DataDir,"/altPreds.rdata"))
  
  if(0) 
    altPreds = readRDS(paste0(DataDir,"/altPreds.rdata"))
  NucIds = as.integer(sapply(NucList, function(x) x$Label))
  
  print("Building Matrix")
  altPreds.mat = abind(nuc.mat[,NucIds], altPreds$orig, altPreds$min, altPreds$exp, altPreds$max, along=3)
  dimnames(altPreds.mat)[[3]] = c("obs", "orig", "min", "exp", "max")
  
  print("saving tracks")
  
  library(rtracklayer)
  for( m in c("WT-H3K18ac-H3K4me3", "WT-H3K36me3-H3K4me3", "WT-H3K36me3-H3K79me3", "WT-H3K4me3-H3K79me3") )
    for( c in names(altPreds) )
    {
      temp = resize(NucRegions[NucIds], fix="center", width=100)
      score(temp) = altPreds[[c]][m,]
      FN = paste0(TrackFN,m,"-",c,".bw")
      export( temp, FN, format="BigWig" )  
    }
  for( m in c("WT-H3K18ac-H3K4me3", "WT-H3K36me3-H3K4me3", "WT-H3K36me3-H3K79me3") )
  {
    temp = resize(NucRegions[NucIds], fix="center", width=100)
    score(temp) = nuc.mat[m,NucIds]
    FN = paste0(TrackFN,m,"-obs",".bw")
    export( temp, FN, format="BigWig" )  
  }    
  source("Multiplot.R")
  for( m in c("WT-H3K18ac-H3K4me3", "WT-H3K36me3-H3K4me3", "WT-H3K36me3-H3K79me3", "WT-H3K4me3-H3K79me3") ) {
    k1 = strsplit(m,"-")[[1]][[2]]
    k2 = strsplit(m,"-")[[1]][[3]]
    x = altPreds.mat[paste0("WT-",k1,"-Input"),,"obs"]
    y = altPreds.mat[paste0("WT-",k2,"-Input"),,"obs"]
    for( C in c("min", "max", "exp")) {
      print(C)
      z = log2(altPreds.mat[m,,"orig"]/altPreds.mat[m,,C])
      png(paste0(FigureDir,"/3dscatter-",m,"-",C,".png"), height = 400, width = 400)
      p = Plot3DSmooth(x,y,z,k1,k2,paste0("log obs/",C), N = 5000)
      multiplot(p)
      dev.off()
    }
  }
}


if( p.Run == "ModelTracks" ) {
  ## output tracks from the model
  library(rtracklayer)
  for( m in colnames(Dist) )
  {
    temp = NucRegions
    s = array(NA, dim=c(length(temp)))
    s[NucIds] = Dist[,m]
    score(temp) = s
    FN = paste0(TrackFN,m,".bw")
    temp = temp[!is.na(as.vector(s))]
    export( resize(temp,fix="center",width = 100), FN, format="BigWig" )  
  }
}

if(0) {
  ## output nuc counts to file
  mat = nuc.mat[,GoodNucs]
  NucM = mcols(NucRegions[GoodNucs])
  Nuc.Names = as.character(NucM[,"acc"])
  I = !is.na(Nuc.Names)
  Nuc.Names[I] = AccToGeneName(Nuc.Names[I],SGD)
  Nuc.Names[!I] = " "
  Nuc.Pos = as.character(NucM[,"gene_pos"])
  Nuc.Pos[is.na(Nuc.Pos)] = ""
  Nuc.Info = data.frame(ID = paste0("N", which(GoodNucs)),
                        Name = Nuc.Names,
                        Pos = Nuc.Pos )
  
  Mat = cbind(Nuc.Info, data.frame(t(mat)))
  write.table(Mat,paste0(DataDir,"/NucCounts-orig.tab"),sep = "\t", quote = F, col.names = T, row.names = F)
}

if(0) {
  ## Generate figures for GM Jan 2016
  
  
  source("Multiplot.R")
  for(k in rownames(weiner.nucs)) {
    png(paste0(FigureDir,"/Weiner-scatter-",k,".png"), height = 400, width = 400)
    p = DensityScatter(nuc.mat[paste0("WT-",k,"-Input"),NucIds], weiner.nucs[k,NucIds], coordeq = F, threshold=.95,
                       xlabel="coChIP input (reads)", ylabel="Weiner et al (AU)", title = k, smooth = T)
    multiplot(p)
    dev.off()
  }
  
  for(k in rownames(weiner.nucs)) {
    png(paste0(FigureDir,"/Weiner-scatter-",k,"-log.png"), height = 400, width = 400)
    p = DensityScatter(nuc.mat[paste0("WT-",k,"-Input"),NucIds], weiner.nucs[k,NucIds], coordeq = F, threshold=.99,logscale=T,
                       xlabel="coChIP input (reads)", ylabel="Weiner et al (AU)", title = k, smooth = T)
    multiplot(p)
    dev.off()
  }
  for( k1 in rownames(weiner.nucs))
    for( k2 in rownames(weiner.nucs))
      if( k1 < k2 ) {
        png(paste0(FigureDir,"/SymScatter-",k1,"-",k2,".png"), height = 400, width = 400)
        p = DensityScatter(nuc.mat[paste0("WT-",k1,"-", k2),NucIds], nuc.mat[paste0("WT-",k2,"-", k1),NucIds], coordeq = F, threshold=.99,
                           xlabel=paste0(k1,"-", k2), ylabel=paste0(k2,"-", k1), title = paste0(k1," ", k2), smooth = T)
        multiplot(p)
        dev.off()
        png(paste0(FigureDir,"/SymScatter-",k1,"-",k2,"-log.png"), height = 400, width = 400)
        p = DensityScatter(nuc.mat[paste0("WT-",k1,"-", k2),NucIds], nuc.mat[paste0("WT-",k2,"-", k1),NucIds], coordeq = F, threshold=.99,
                           xlabel=paste0(k1,"-", k2), ylabel=paste0(k2,"-", k1), title = paste0(k1," ", k2), smooth = T, logscale = T)
        multiplot(p)
        dev.off()
        
      }
}

