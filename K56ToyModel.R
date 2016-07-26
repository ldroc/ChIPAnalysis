setwd("~/Dropbox/CoChIP/coChIP-Analysis")
library(deSolve)
library(ggplot2)
source("Multiplot.R")

cellcycle.time = 100

Parameters <- c( Rto = log(2)/40, Rto.k36 = 0, Rme = log(2)/30, Rdeac.bg = log(2)/100, Rdeac.Mphase = log(2)/25 )


# First bit K36me3, second bit K56ac

State <- c( S00 = .25, S01 = 0.25, S10 = 0.25, S11 = 0.25)

K56K36 <- function( t, State, Parameters ) {
  with(as.list(c(State,Parameters)), {
    r.deac = Rdeac.bg
    if( t > 70 && t <= 90 )
      r.deac = Rdeac.Mphase
    
    dS00 = -S00 * (Rto + Rme) + S01*r.deac
    dS01 = -S01 * (Rme + r.deac) + S00*Rto + (S10+S11)*Rto.k36
    dS10 = -S10 * Rto.k36 + S00 * Rme + S11 * r.deac
    dS11 = -S11 * (Rto.k36 + r.deac) + S01 * Rme
    list(c(dS00, dS01, dS10, dS11))
  })
}

K56K36.replicate <- function( State ) {
  s = State/2
  s["S01"] = s["S01"] + 1/2
  s
}

SimulateLoci <- function( Parameters, rep.time = 50, N.repeat = 10 ) {

  Times1 = seq(0,rep.time,by = 2)
  Times2 = seq(rep.time, cellcycle.time, by = 2)

  s = State
  for( i in 1:N.repeat) {
    out1 = ode( y = s, times = Times1, func = K56K36, parms = Parameters )
    s = out1[length(Times1), names(State)]
    s = K56K36.replicate(s)
    out2 = ode( y = s, times = Times2, func = K56K36, parms = Parameters )
    s = out2[length(Times2), names(State)]
  }
  rbind(out1,out2)
}


PlotSimResults <- function( ODE.out ) {
  n = dim(ODE.out)[[2]] -1
  m = dim(ODE.out)[[1]]
  p = ggplot()
  p = p + scale_color_brewer( palette = "Set1",
#                             labels=c("K56ac", "K36me3", "both"),
                             guide = "legend")
  
  p = p + geom_line(data=data.frame(x = ODE.out[,1],
                                    y = ODE.out[,"S01"]+ODE.out[,"S11"], color = rep("K56ac",m)),
                      aes(x=x, y=y, color=color), 
                      show.legend = TRUE)
  
  p = p + geom_line(data=data.frame(x = ODE.out[,1],
                                    y = ODE.out[,"S10"]+ODE.out[,"S11"], color = rep("K36me3",m)),
                    aes(x=x, y=y, color=color), 
                    show.legend = TRUE)
  p = p + geom_line(data=data.frame(x = ODE.out[,1],
                                    y = ODE.out[,"S11"], color = rep("both",m)),
                    aes(x=x, y=y, color=color), 
                    show.legend = TRUE)
  p = p + labs(x=colnames(ODE.out)[[1]], y="frequency")
  p = p+theme(legend.position="right")
  p = p+ylim(c(0,1))
}


TestLoci <- function(delta = 5, Parameters = Parameters ) {
  RepT = seq( 40, 70, by = delta)
  K36T = seq( 10,80, by = delta)
  
  M = matrix( nr = length(K36T), nc = length(RepT))
  rownames(M) = as.character(K36T)
  colnames(M) = as.character(RepT)
  M.K36 = M
  M.K56 = M
  M.K56.K36 = M
  
  for( k in K36T )
    for( r in RepT ) {
      p = Parameters
      p["Rme"] = log(2)/k
      o = SimulateLoci(p, r)
      s = colMeans(o)
      kn = as.character(k)
      rn = as.character(r)
      M.K56[kn,rn] = s["S01"]+s["S11"]
      M.K36[kn,rn] = s["S10"]+s["S11"]
      M.K56.K36[kn,rn] = s["S11"]
    }
  
  return( list(K56 = M.K56, K36 = M.K36, K56.K36 = M.K56.K36))
}

PlotTestMatrix <- function(M, name = "K36me3") {
  m = dim(M)[[1]]
  n = dim(M)[[2]]
  X = rep(1:m, n)
  Y = rep(1:n, each=m)
  df = data.frame(x=X, y=Y, z = as.vector(M) )
  p = ggplot(df, aes(x=x,y=y,fill=z)
             )
  p = p + geom_raster(interpolate = TRUE)
  p = p + scale_fill_distiller(name, guide = "colorbar", palette = "Reds", direction = 1)
  p = p + scale_x_continuous(breaks = 1:m, labels = rownames(M))
  p = p + scale_y_continuous(breaks = 1:n, labels = colnames(M))
  p = p + labs(x = "K36 deposit half time", y = "Replication time")
  p = p + coord_fixed(ratio=1) 
  p  
}

do <- function() {
  turnover = c( glacial = 200, slow = 70, medium = 40, fast = 5 )
  K36TurnoverRatio = c( zero = 0, full = 1)
  for( k in 1:length(K36TurnoverRatio))
    for( i in 1:length(turnover)) {
      Parameters["Rto"] = log(2) / turnover[i]  
      Parameters["Rto.k36"] = K36TurnoverRatio[k] * Parameters["Rto"]
      iName = names(turnover)[[i]]
      kName = names(K36TurnoverRatio)[[k]]
      Ms = TestLoci(2, Parameters)
      p36 = PlotTestMatrix(Ms$K36, "K36me3") + ggtitle(paste(iName, "turnover (t1/2 = ", turnover[[i]],")", "K36 turnover - ", kName))
      p56 = PlotTestMatrix(Ms$K56, "K56ac")
      p5636 = PlotTestMatrix(Ms$K56.K36, "K56ac+K36me3")
      png(paste0("~/Google Drive/CoChIPAnalysis/K56Figures/K56ac-K36me-Simulation-",iName,"-",kName,".png"),width=700, height=1000)
      multiplot(p36,p56,p5636)
      dev.off()
      
    }
  times = c(early =40, mid = 55, late= 70)

  for( k in 1:length(K36TurnoverRatio))
    for( t in 1:length(times))
      for( i in 1:length(turnover))
        for( j in 1:length(turnover) ) {
          Parameters["Rto"] = log(2) / turnover[i]
          Parameters["Rto.k36"] = K36TurnoverRatio[k] * Parameters["Rto"]
          Parameters["Rme"] = log(2) / turnover[j]  
          iName = names(turnover)[[i]]
          jName = names(turnover)[[j]]
          tName = names(times)[[t]]
          kName = names(K36TurnoverRatio)[[k]]
          p = PlotSimResults(SimulateLoci(Parameters, times[t]))+
            ggtitle(paste("Turnover: ", iName ," K36me3: ", jName, " Rep time = ",tName, " K36 methylation = ", kName))
          ggsave( p,
                  filename = paste0("~/Google Drive/CoChIPAnalysis/K56Figures/K56ac-K36me-Simulation-",iName,"-",jName,"-", tName,"-",kName,".png"),width=5, height=4)
          
        }
  
}