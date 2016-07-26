setwd("~/Dropbox/coChIP-scripts/")
source("CoChip-Functions.R")
source("NucAtlas.R")

Turnover = read.csv("~/Work/Papers/ChIP_Analysis/Files/Agilent_lambdas_vs_probes.tab",
                    sep="\t", stringsAsFactors = FALSE)

XX = strsplit(Turnover[,"Ag_id"]," ")

lchrnames = c("chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI", "chrVII", 
              "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII", 
              "chrXIV", "chrXV", "chrXVI" )

Chr = sapply(XX, function(l) lchrnames[as.integer(substr(l[[1]],4,6))])
Pos = sapply(XX, function(l) as.integer(l[[2]]))

Probe2Nuc = rep(NA, length(Chr))

for( l in c( 1, 5,  10, 25, 50, 75,  100) ) {
  Probes = which( is.na(Probe2Nuc))
  print(paste(l,length(Probes)))
  if( length(Probes) > 0 ) {
    rs <- GRanges( seqnames = Chr[Probes],
                   ranges = IRanges( Pos[Probes] - l, Pos[Probes] + l),
                   seqinfo = Seqinfo(genome="sacCer3"))
    Probe2Nuc[Probes] = findOverlaps(rs, NucRegions, select="arbitrary")
  }
}

NucTurnoverRate = rep(NA, length(NucRegions))
names(NucTurnoverRate) = 1:length(NucRegions)
NucTurnoverZ = NucTurnoverRate

Probes = which(!is.na(Probe2Nuc))
NucTurnoverRate[Probe2Nuc[Probes]] = Turnover[Probes,"AvgOfLAMBDA"]
NucTurnoverZ[Probe2Nuc[Probes]] = Turnover[Probes,"AvgOfZscore"]

saveRDS(list(Rate = NucTurnoverRate, Z = NucTurnoverZ), "~/Google Drive/CoChIPAnalysis/NucTurnover.rdata")

