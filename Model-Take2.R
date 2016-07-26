setwd("~/Dropbox/coChIP-scripts")
source("CoChip-Functions.R")
source("NucAtlas.R")

DataDir = "~/Google Drive/CoChIPAnalysis"
FigureDir = "~/Google Drive/CoChIPAnalysis/ModelII"


source("PairwiseModel.R")
source("ComplexModel.R")

#
# Step 1 - Build initial model from Pairwise estimates
#

ChosenStrains = c("1.WT", "1.K18R", "1.K4R", "1.K36R", "1.K79R")
ChosenExperiments = do.call(c,lapply(ChosenStrains, function(s) rownames(Nucs)[grep(s, rownames(Nucs))]))
Experiments = do.call(rbind,strsplit(ChosenExperiments,"-"))
colnames(Experiments) = c("Strain","IP1", "IP2")
rownames(Experiments) = ChosenExperiments
ChosenABs = unique(Experiments[,"IP1"])

ChosenMods = ChosenABs[grep("H3", ChosenABs, invert = TRUE)]

Experiments[Experiments[,"IP2"] == "Input","IP2"] = NA


ModNumber = length(ChosenMods)
StrainNumber = length(ChosenStrains)
ABsNumber = length(ChosenABs)

NucNumber = dim(Nucs)[2]
NucNames = colnames(Nucs)

InitialFN = paste0(DataDir,"/modelII-initial.rdata")

if(file.exists(InitialFN)) {
  InitialModel = readRDS(InitialFN)
  Design = InitialModel$Design
  Model.dummy = InitialModel$Model
  Nuc.initial = InitialModel$Nuc
} else {
  Design = CreateDesign(ChosenExperiments)
  Design$Constraints = data.frame( 
    Strain = c("1K4R", "1.K18R","1.K36R"), 
    Mod = c("K4me3", "K18ac", "K36me3"), 
    I = c(0,0,0) 
  )
  
  Model.dummy = CreateModel(Design)
  Nuc.Dummy = CreateNuc(Design)
  
  Model.Compile = abModelCompileDesign(Design)    
  
  TempMarginals = array( 1, dim = c(ModNumber, StrainNumber, NucNumber),
         dimnames = list( ChosenMods, ChosenStrains, NucNames))
  TempInteractions = array( 0, dim = c(ModNumber,ModNumber,StrainNumber,NucNumber), 
                        dimnames = list(ChosenMods, ChosenMods, ChosenStrains, NucNames))
  for( s in ChosenStrains ) {
    
    upperMarginals = Model.Compile$MarginUpper[[s]]
    for( a in ChosenMods ) 
      TempMarginals[a,s,] = pmin(ModelDistribution(s, a)[NucNames], upperMarginals[[a]])
    
    for( a in ChosenMods ) 
      for( b in ChosenMods )
        TempInteractions[a,b,s,] = ExtractInteraction(TempMarginals[a,s,],TempMarginals[b,s,],  
                                                      ModelDistribution(s, a, b)[NucNames])
  }
  TempInteractions[TempInteractions == Inf] = 10
  TempInteractions[TempInteractions == -Inf] = -10
  
  Nuc.initial =  list()
  for( i in 1:length(NucNames) ) {
    TempNuc = Nuc.Dummy
    TempNuc$Label = NucNames[[i]]
    TempNuc$Occ[] = Inputs[[i]]
    TempNuc$Marginals = TempMarginals[,,i] 
    TempNuc$Interactions = TempInteractions[,,,i]
    Nuc.initial[[i]] = TempNuc
  }
    
  
  saveRDS(list(Design=Design, Model=Model.dummy, Nuc = Nuc.initial), InitialFN)
}  
  
Model.Compile = abModelCompileDesign(Design)

if( 0 ) {
Model.start = Model.dummy
for( a in Design$Mods )
  Model.start$IPs[paste0("ab-",a),a] = 0.1
Model.start$IPs["ab-H3",] = 0.05

Model.learned = abModelOptimizeModelMultiNucs(Nucs[ChosenExperiments,], Model.start, Nuc.initial, Model.Compile)
Model.learned.1 = abModelNormalizeModel(Model.learned)
Model.learned.2 = abModelOptimizeModelMultiNucs(Nucs[ChosenExperiments,], Model.learned.1, Nuc.initial, Model.Compile, n = 2500)
Model.learned.3 = abModelNormalizeModel(Model.learned.2)
Model.learned.4 = abModelOptimizeModelMultiNucs(Nucs[ChosenExperiments,], Model.learned.3, Nuc.initial, Model.Compile, n = 5000)
Model.learned.5 = abModelNormalizeModel(Model.learned.4)

Model.Predict = do.call(cbind,lapply(Nuc.initial[1:5000], function(n) abModelPredict(Model.learned.5, n, Model.Compile)))

MStat = abModelModelStat(Model.learned.5, Model.Compile)
NewNucs = abModelOptimizeNucs(Nucs[ChosenExperiments,1:1000], Model.learned.5, Nuc.initial[1:1000], Model.Compile )

I = sample(65000,5000)
NewNucs = abModelOptimizeNucs(Nucs[ChosenExperiments,I], Model.learned.5, Nuc.initial[I], Model.Compile )

Model.Predict.1 = do.call(cbind,lapply(NewNucs, function(n) abModelPredict(Model.learned.5, n, Model.Compile)))
Model.learned.6 = abModelOptimizeModelMultiNucs(Nucs[ChosenExperiments,], Model.learned.5, NewNucs, Model.Compile, n = 5000)
Model.learned.7 = abModelNormalizeModel(Model.learned.6)

I.1 = sample(65000,5000)
NewNucs.1 = abModelOptimizeNucs(Nucs[ChosenExperiments,I.1], Model.learned.7, Nuc.initial[I.1], Model.Compile )
Model.learned.7.5 = abModelOptimizeModelMultiNucs(Nucs[ChosenExperiments,], Model.start, NewNucs.1, Model.Compile, n = 5000)

Model.learned.8 = abModelOptimizeModelMultiNucs(Nucs[ChosenExperiments,], Model.learned.7, NewNucs.1, Model.Compile, n = 5000)
Model.learned.9 = abModelNormalizeModel(Model.learned.8)

Model.Final = Model.learned.9
saveRDS(Model.Final,paste0(DataDir,"/modelII-learned.rdata"))
BatchSize = 500
#NewNucs = list()
NewNucs =   readRDS(paste0(DataDir,"/modelII-nucs.rdata"))
for( i in seq(length(NewNucs)+1,length(Nuc.initial)-1, by=BatchSize)) {
  j = min( i + BatchSize-1, length(Nuc.initial) )
  TempNucs = abModelOptimizeNucs(Nucs[ChosenExperiments,i:j], 
                                 Model.Final, Nuc.initial[i:j], Model.Compile )
  NewNucs = c(NewNucs, TempNucs)
  saveRDS(NewNucs,paste0(DataDir,"/modelII-nucs.rdata"))
  print(j)
}
Model.Final.Predict = do.call(cbind,lapply(NewNucs, function(n) abModelPredict(Model.Final, n, Model.Compile)))

plt = list()
for( a in c(Design$Mods, paste0("-", Design$Mods,"-"), Design$Strains, "Input"))
  plt[[a]] = DensityScatter(Nucs[ChosenExperiments[grep(a, ChosenExperiments)],], Model.Final.Predict[grep(a, ChosenExperiments),], title = a)
}

NewNucs = readRDS(paste0(DataDir,"/modelII-nucs.rdata"))

InitialXs = sapply(Nuc.initial, abModelNucStat)
NewXs = sapply(NewNucs, abModelNucStat)
