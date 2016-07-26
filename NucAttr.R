source("NET-Seq.R")

##
## Assumes - SGD, Inputs, 
BuildNucProperties <- function( NucRegions, NucInputs ) {
  if(length(NETWigs) < 1)
    ReadNETData()
  #
  # Position within a gene
  #
  N = length(NucRegions)
  MaxPos = 10
  MinPos = -3
  PosRange = c(MinPos:-1,1:MaxPos)
  GenePos = matrix( nr = length(PosRange), nc = N)
  rownames(GenePos) = paste("Nuc ",PosRange)
  for( i in 1:length(PosRange)) {
    p = PosRange[[i]]
    GenePos[i,] = NucRegions$gene_pos == p
  }
  
  GeneRegion = matrix(nr=3,nc =N)
  rownames(GeneRegion) = c("Promoter", "5' end", "Body")
  GeneRegion[1,] = NucRegions$gene_pos %in% c(-3,-2,-1)
  GeneRegion[2,] = NucRegions$gene_pos %in% c(1, 2, 3)
  GeneRegion[3,] = NucRegions$gene_pos > 3
  
  #
  # Occupancy
  # 
  OccBin = 5
  OccT = quantile(NucInputs,seq(0,1,length.out = OccBin+1),na.rm=TRUE)
  NucOcc = matrix(nr = OccBin, nc = N)
  rownames(NucOcc) = paste0("Occ ", 100*(1:OccBin -1)/OccBin,"-",100*(1:OccBin)/OccBin,"%")
  for( i in 1:OccBin)
    NucOcc[i,] = NucInputs > OccT[[i]] & NucInputs <= OccT[[i+1]]
  
  #
  # Gene expression
  # 
  Expr = matrix(nr = length(SGD$ExprQuan), nc = N)
  rownames(Expr) = paste("ExprQ", names(SGD$ExprQuan))
  for( i in 1:length(SGD$ExprQuan)) {
    Expr[i,] = NucRegions$acc %in% SGD$ExprQuan[[i]]
    Expr[i,is.na(NucRegions$acc)] = NA
  }
  
  # Gene Expression X Region
  ExprRegion = matrix(nr = length(SGD$ExprQuan)*3, nc = N)
  rownames(ExprRegion) = sapply(rownames(GeneRegion), function(r) paste(rownames(Expr),r))
  for( j in 1:3) 
    for( i in 1:length(SGD$ExprQuan)) 
      ExprRegion[(j-1)*length(SGD$ExprQuan) + i, ] = Expr[i,] & GeneRegion[j,]
      
  #
  # Build netseq
  #
  NET.Plus = NetSeqSignal(NETWigs$WT_plus, NucRegions)
  NET.Minus = NetSeqSignal(NETWigs$WT_minus, NucRegions)
  GeneStrand = as.vector(strand(SGD$Genes))
  names(GeneStrand) = SGD$Genes$acc
  NucStrand = GeneStrand[as.character(NucRegions$acc)]
  NucStrand[is.na(NucStrand)] = "*"
  NET.Sense = vector(length=N)
  NET.Sense[] = NA
  NET.Antisense = vector(length=N)
  NET.Antisense[] = NA
  NET.Sense[NucStrand == "+"] = NET.Plus[NucStrand == "+"]
  NET.Sense[NucStrand == "-"] = NET.Minus[NucStrand == "-"]
  NET.Antisense[NucStrand == "+"] = NET.Minus[NucStrand == "+"]
  NET.Antisense[NucStrand == "-"] = NET.Plus[NucStrand == "-"]
  
  NET.Levels = 5
  NET = matrix(nr=7*NET.Levels,nc = N)
  SenseT = quantile(NET.Sense[NET.Sense > 0],seq(0,1,length.out = NET.Levels+1),na.rm = TRUE)
  SenseT[[1]] = 0
  SenseT[[NET.Levels+1]] = Inf
  AntisenseT = quantile(NET.Antisense[NET.Antisense>0],seq(0,1,length.out = NET.Levels+1),na.rm = TRUE)
  AntisenseT[[1]] = 0
  AntisenseT[[NET.Levels+1]] = Inf
  NET.ASvsS = NET.Antisense/NET.Sense 
  ASvsST = quantile(NET.ASvsS,seq(0,1,length.out = NET.Levels+1),na.rm = TRUE)
  ASvsST[[1]] = -Inf
  ASvsST[[NET.Levels+1]] = Inf
  
  rownames(NET) = c(paste0("NET S ",100*(1:NET.Levels -1)/NET.Levels,"-",100*(1:NET.Levels)/NET.Levels,"%"),
                    paste0("NET AS ",100*(1:NET.Levels -1)/NET.Levels,"-",100*(1:NET.Levels)/NET.Levels,"%"),
                    paste0("NET AS/S ",100*(1:NET.Levels -1)/NET.Levels,"-",100*(1:NET.Levels)/NET.Levels,"%"),
                    paste(paste0("NET S ",100*(1:NET.Levels -1)/NET.Levels,"-",100*(1:NET.Levels)/NET.Levels,"%"),"5' end"),
                    paste(paste0("NET S ",100*(1:NET.Levels -1)/NET.Levels,"-",100*(1:NET.Levels)/NET.Levels,"%"),"Body"),
                    paste(paste0("NET AS ",100*(1:NET.Levels -1)/NET.Levels,"-",100*(1:NET.Levels)/NET.Levels,"%"),"5' end"),
                    paste(paste0("NET AS ",100*(1:NET.Levels -1)/NET.Levels,"-",100*(1:NET.Levels)/NET.Levels,"%"),"Body") )
  for( i in 1:NET.Levels) {
    NET[i,] = NET.Sense >= SenseT[[i]] & NET.Sense < SenseT[[i+1]]
    NET[NET.Levels+i,] = NET.Antisense >= AntisenseT[[i]] & NET.Antisense < AntisenseT[[i+1]]
    NET[2*NET.Levels+i,] = NET.ASvsS >= ASvsST[[i]] & NET.ASvsS < ASvsST[[i+1]]
    NET[3*NET.Levels+i,] = NET[i,] & GeneRegion["5' end",]
    NET[4*NET.Levels+i,] = NET[i,] & GeneRegion["Body",]
    NET[5*NET.Levels+i,] = NET[NET.Levels+i,] & GeneRegion["5' end",]
    NET[6*NET.Levels+i,] = NET[NET.Levels+i,] & GeneRegion["Body",]
    
  }

  #
  # Compute 5'/3' ratio and Sense/As ratio
  #
  GeneNucdf = data.frame(nuc = NucRegions$gene_pos)
  GeneNucLength = aggregate(GeneNucdf,by=list(NucRegions$acc),max)
  GeneNucLength = GeneNucLength[GeneNucLength$nuc >6,]
  Nuc.5 = sapply(1:length(GeneNucLength$nuc), function(i) which(NucRegions$acc == GeneNucLength$Group.1[[i]] & 
                                                              NucRegions$gene_pos > 0 & 
                                                              NucRegions$gene_pos <= min(5,GeneNucLength$nuc[[i]]/3)))
  Nuc.3 = sapply(1:length(GeneNucLength$nuc), function(i) which(NucRegions$acc == GeneNucLength$Group.1[[i]] & 
                                                                  NucRegions$gene_pos >= max(GeneNucLength$nuc[[i]]-5,2*GeneNucLength$nuc[[i]]/3)))
  NETGenes.5 = sapply(Nuc.5, function(l) mean(NET.Sense[l]))
  names(NETGenes.5) = GeneNucLength$Group.1
  NETGenes.3 = sapply(Nuc.3, function(l) mean(NET.Sense[l]))
  names(NETGenes.3) = GeneNucLength$Group.1
  NETGenes.Total = sapply(GeneNucLength$Group.1, function(a) mean(NET.Sense[which(NucRegions$acc == a & NucRegions$gene_pos > 0)]))
  names(NETGenes.Total) = GeneNucLength$Group.1
  Nuc.I = pmin(sapply(Nuc.3,length),sapply(Nuc.5,length)) > 1 & NETGenes.Total > median(NETGenes.Total,na.rm=TRUE) - 2*mad(NETGenes.Total,na.rm=TRUE)
  NETGenes.5to3 = log2(NETGenes.5/NETGenes.3)[Nuc.I]
  NETGenes.mu = median(NETGenes.5to3,na.rm=TRUE)
  NETGenes.sd = mad(NETGenes.5to3,na.rm=TRUE)
  Nuc.2thirds = sapply(1:length(GeneNucLength$nuc), function(i) which(NucRegions$acc == GeneNucLength$Group.1[[i]] & 
                                                                      NucRegions$gene_pos > GeneNucLength$nuc[[i]]/3))
  Nuc.1third = sapply(1:length(GeneNucLength$nuc), function(i) which(NucRegions$acc == GeneNucLength$Group.1[[i]] & 
                                                                        NucRegions$gene_pos <= GeneNucLength$nuc[[i]]/3))
  

  UnBalancedNET = matrix(nr=4,nc=N)
  rownames(UnBalancedNET) = c("5'/3' high, body","5'/3' low, body","5'/3' high, all","5'/3' low, all")
  UnBalancedNET[,unlist(Nuc.2thirds[Nuc.I])] = FALSE
  UnBalancedNET[c(1,3),unlist(Nuc.2thirds[Nuc.I][NETGenes.5to3 > NETGenes.mu + NETGenes.sd])] = TRUE
  UnBalancedNET[c(2,4),unlist(Nuc.2thirds[Nuc.I][NETGenes.5to3 < NETGenes.mu - NETGenes.sd])] = TRUE
  UnBalancedNET[c(3,4),unlist(Nuc.1third[Nuc.I])] = FALSE
  UnBalancedNET[3,unlist(Nuc.1third[Nuc.I][NETGenes.5to3 > NETGenes.mu + NETGenes.sd])] = TRUE
  UnBalancedNET[4,unlist(Nuc.1third[Nuc.I][NETGenes.5to3 < NETGenes.mu - NETGenes.sd])] = TRUE
  
  Attr = rbind(GenePos,GeneRegion,NucOcc,Expr,ExprRegion, NET, UnBalancedNET)
}

MultiEnrichmentTest <- function(Test,Attr) {
  M = dim(Attr)[[1]]
  L = length(Test)
  Counts = data.frame(Test = rep(names(Test), M), Attr = rep(rownames(Attr), each = L))
  
  for( t in names(Test) )
    for( a in rownames(Attr))
    {
      Ref = !is.na(Test[[t]]) & !is.na(Attr[a,])
      N = length(which(Ref))
      n = length(which(Ref & Attr[a,]))
      K = length(which(Ref & Test[[t]]))
      k = length(which(Ref & Test[[t]] & Attr[a,]))
      pval = -phyper(k-1,n,N-n,K,log.p=TRUE,lower.tail = FALSE)/log(10)
      i = which(Counts$Test == t & Counts$Attr == a)
      if(length(i) != 1)
        print(paste(t,a,length(i)))
      Counts[i,"N"] = N
      Counts[i,"K"] = K
      Counts[i,"n"] = n
      Counts[i,"k"] = k
      Counts[i,"pval"] = pval
      Counts[i,"%"] = (k/K)*100
      Counts[i,"Exp %"] = (n/N)*100
      Counts[i,"Enrichment"] = (k/K)/(n/N)
    }
  Counts
}

Bonferroni <- function(Counts) {
  log10(length(which(!is.na(Counts$pval)))) - log10(0.05)
}
