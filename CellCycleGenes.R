##
## CellCycle Analysis
##
## Data downloaded from cyclebase on 31/12/2015
##
source("NucAtlas.R")
source("DensityScatter.R")
source("Multiplot.R")
GetSGDandNucs()

TargetMods = c("H3K18ac", "H3K4ac", "H3K36me3" , "H3K79me3", "H3K4me3", "H2AS129ph", "H3S10ph", "Htz1", "H3K56ac", "H4K12ac", "H4K16ac", "H4K20me" )

mod.nucs = ReadHistoneModeAtlas("Weiner-HisMod.csv", NucRegions, TargetMods )
m = mcols(mod.nucs)
weiner.nucs = do.call(rbind,lapply(TargetMods, function(mod) as.vector(exp(as.matrix(m["Input"]) + as.matrix(m[mod])))))
weiner.nucs = rbind(weiner.nucs, as.vector(exp(as.matrix(m["Input"]))))
rownames(weiner.nucs) = c(TargetMods, "Input")
Inputs = weiner.nucs["Input",]
t.high = quantile(Inputs,0.99, na.rm=TRUE)
t.low = quantile(Inputs,0.01, na.rm=TRUE)
nucs.abnormal = Inputs > t.high | Inputs < t.low

nuc.1 = which(mcols(NucRegions)$gene_pos == 1 & !nucs.abnormal)
names(nuc.1) = mcols(NucRegions)$acc[nuc.1]
nuc.3 = which(mcols(NucRegions)$gene_pos == 3& !nucs.abnormal)
names(nuc.3) = mcols(NucRegions)$acc[nuc.3]
nuc.5 = which(mcols(NucRegions)$gene_pos == 5& !nucs.abnormal)
names(nuc.5) = mcols(NucRegions)$acc[nuc.5]

nuc.1.expr = unlist(lapply(names(nuc.1),
                           function(a) SGD$Genes$expr[SGD$Genes$acc == a] ))

nuc.1.90 = unlist(lapply(TargetMods, function(x) quantile(weiner.nucs[x,nuc.1]/weiner.nucs["Input",nuc.1],.90, na.rm = T)))
names(nuc.1.90) = TargetMods
nuc.1.50 = unlist(lapply(TargetMods, function(x) quantile(weiner.nucs[x,nuc.1]/weiner.nucs["Input",nuc.1],.50, na.rm = T)))
names(nuc.1.50) = TargetMods
nuc.1.25 = unlist(lapply(TargetMods, function(x) quantile(weiner.nucs[x,nuc.1]/weiner.nucs["Input",nuc.1],.25, na.rm = T)))
names(nuc.1.25) = TargetMods

nuc.1.expr.10 = quantile(nuc.1.expr, .10, na.rm = T)
nuc.1.weird = nuc.1.expr > nuc.1.expr.10 & 
  weiner.nucs["H3K36me3",nuc.1]/weiner.nucs["Input",nuc.1] > nuc.1.90["H3K36me3"] &
  weiner.nucs["H3K4me3",nuc.1]/weiner.nucs["Input",nuc.1] > nuc.1.25["H3K4me3"] &
  weiner.nucs["H3K18ac",nuc.1]/weiner.nucs["Input",nuc.1] > nuc.1.25["H3K4ac"]
  

cctab = read.delim("cerevisiae_periodic.tsv")

I = cctab$periodicity_pvalue < 0.0001 & cctab$regulation_pvalue < 0.001
cctab = cctab[I,]
# hist(cctab$peaktime, breaks = "FD")

cc.nuc.1 = names(nuc.1) %in% cctab$gene
cc.nuc.5 = names(nuc.5) %in% cctab$gene

cc.groups = lapply(seq(0,89,by=10), function(x)  cctab$gene[which(cctab$peaktime > x & cctab$peaktime <= x+10)])
cc.groups.1 = lapply(cc.groups, function(x) names(nuc.1) %in% x)
cc.groups.5 = lapply(cc.groups, function(x) names(nuc.5) %in% x)

cc.ks = -log10(unlist(lapply(TargetMods, function(x) ks.test(weiner.nucs[x, nuc.1[cc.nuc.1]],weiner.nucs[x, nuc.1[!cc.nuc.1]])$p.value )))
names(cc.ks) = TargetMods
nuc.1.weird.ks = -log10(unlist(lapply(TargetMods, function(x) ks.test(weiner.nucs[x, nuc.1[nuc.1.weird]],weiner.nucs[x, nuc.1[!nuc.1.weird]])$p.value )))
names(nuc.1.weird.ks) = TargetMods

DataDir = "~/Dropbox/CoChIP/151217-HisMut-K4/weiner-pairwise-plots"
if( 0 )
for( x in TargetMods )
  for( y in TargetMods )
    if( x < y ) {
      if(0) {
      png(paste0(DataDir,"/nuc1-",x,"-",y,".png"), height = 800, width = 800)
      p = DensityScatter(weiner.nucs[x,nuc.1], weiner.nucs[y,nuc.1],
                       xlabel = x, ylabel = y,
                       coordeq = F,  alpha = .25, 
                       Groups = cc.groups.1, logscale = T, galpha = 1 )
      multiplot(p)
      dev.off()
      png(paste0(DataDir,"/nuc5-",x,"-",y,".png"), height = 800, width = 800)
      p = DensityScatter(weiner.nucs[x,nuc.5], weiner.nucs[y,nuc.5],
                         xlabel = x, ylabel = y,
                         coordeq = F,  alpha = .25, 
                         Groups = cc.groups.5, logscale = T, galpha = 1 )
      multiplot(p)
      dev.off()
      }
      png(paste0(DataDir,"/nuc1-",x,"-",y,".png"), height = 800, width = 800)
      p = DensityScatter(weiner.nucs[x,nuc.1]/weiner.nucs["Input",nuc.1], weiner.nucs[y,nuc.1]/weiner.nucs["Input",nuc.1],
                         xlabel = paste(x,"/Input"), ylabel = paste(y,"/Input"),
                         coordeq = F,  alpha = .75, 
                         Groups = list(cc.nuc.1), logscale = T, galpha = 1, GroupContour = T )
      multiplot(p)
      dev.off()
      png(paste0(DataDir,"/nuc5-",x,"-",y,".png"), height = 800, width = 800)
      p = DensityScatter(weiner.nucs[x,nuc.5], weiner.nucs[y,nuc.5],
                         xlabel = x, ylabel = y,
                         coordeq = F,  alpha = .75, 
                         Groups = list(cc.nuc.5), logscale = T, galpha = 1,
                         GroupContour = T )
      multiplot(p)
      dev.off()
    }
if(0)
for( x in TargetMods ) {
  png(paste0(DataDir,"/nuc1-",x,"-exr",".png"), height = 800, width = 800)
  p = DensityScatter(weiner.nucs[x,nuc.1]/weiner.nucs["Input",nuc.1], nuc.1.expr,
                     xlabel = paste(x,"/Input"), ylabel = "expr",
                     coordeq = F,  alpha = .75, 
                     Groups = list(cc.nuc.1), logscale = T, galpha = 1, GroupContour = T )
  multiplot(p)
  dev.off()
}