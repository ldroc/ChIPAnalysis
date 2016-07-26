library(gdata)
library(ggplot2)
library(grid)
setwd("~/Dropbox/CoChIP/coChIP-Analysis")

source("CoChip-Functions.R")
source("NucAtlas.R")
source("Multiplot.R")

rc = read.xls("~/Data/CoChIP/coChIP-Analysis/KLac/SeqCounts.xlsx")
## remove bad entries
Bad = (rc$X2nd.Ab == "K36me3" & rc$Hi.Low == "Low") 
#| 
#  (rc$X1st.Ab == "FlagA" & rc$X2nd.Ab == "K36me3" & rc$X.SC == 100 & rc$Hi.Low == "Low") |
#  (rc$X1st.Ab == "FlagA" & rc$X2nd.Ab == "K18ac" & rc$X.SC == 40 & rc$Hi.Low == "Low")

rc = rc[!Bad,]
rc$Flag = (rc$X1st.Ab == "FlagA" | rc$X1st.Ab == "FlagB")
rc$K18 = (rc$X2nd.Ab == "K18ac")
rc$K36 = (rc$X2nd.Ab == "K36me3")
rc$Hi = (rc$Hi.Low == "Hi")
rc$logSC = log(rc$X.SC+.0001)


t = vector(length=length(rc$X1st.Ab))
for(y in levels( rc$Hi.Low) ) 
  for( x in levels(rc$X2nd.Ab) ) 
    t[rc$Flag & rc$X2nd.Ab == x&rc$Hi.Low==y] = paste("Flag", x,y)
for(y in levels( rc$Hi.Low) ) 
  for( x in levels(rc$X2nd.Ab) ) 
    t[!rc$Flag & rc$X2nd.Ab == x&rc$Hi.Low==y] = paste("Myc", x,y)

rc$Title = factor(t)

rc.glm = glm(SC.Only ~ logSC+Flag+Hi+K36+K18, family = poisson(link="log"), rc)
rc$fitted = fitted(rc.glm)

PlotRc <- function( rc ) {
  p = ggplot( rc, aes(x = X.SC, y = SC.Only, color=Title))
  for( t in levels(factor(rc$Title)) ) {
    df = data.frame( x = rc[rc$Title == t,]$X.SC, 
                     y = rc[rc$Title == t,]$fitted,
                     Title = rc[rc$Title == t,]$Title )
#    print(df)
    p = p + geom_line( data=df,aes(x = x, y=y) )
  }
  p = p+geom_point(alpha=1,size=4)
  p
}

p1 = PlotRc( rc[rc$Flag & rc$Hi.Low == "Hi",])
p2 = PlotRc( rc[rc$Flag & rc$Hi.Low != "Hi",])
p3= PlotRc( rc[!rc$Flag & rc$Hi.Low == "Hi",])
p4= PlotRc( rc[!rc$Flag & rc$Hi.Low != "Hi",])
multiplot(p1,p2,p3,p4,cols=2)


if(0) {


ll <- function( par  ) {
  aFlagMyc = par[1]
  aHighLow = par[2]
  aK36H3 = par[3]
  aK18H3 = par[4]
  b = par[5]
  res = 0
  for(l in levels( rc$Hi.Low) ) 
    for( a1 in levels(factor(rc$Flag)) ) 
      for( a2 in levels(rc$X2nd.Ab) ) {
#        print(paste(l,a1,a2))
        a = b
        if( l == "Hi") a = a + aHighLow
        if( a1 == TRUE ) a = a + aFlagMyc
        if( a2 == "K36me3") a = a+aK36H3
        if( a2 == "K18ac") a = a+aK18H3
        rcTemp = rc[rc$Hi.Low == l & rc$Flag == a1 & rc$X2nd.Ab == a2,] 
#        print(exp(a))
        pred = exp(a)*rcTemp$X.SC+1
        obs = rcTemp$SC.Only
#        print((rcTemp$SC.Only - pred)**2)
        res = res + sum(((obs - pred )/pred)**2)
      }
  return(res)
}

dll <- function( par  ) {
  aFlagMyc = par[1]
  aHighLow = par[2]
  aK36H3 = par[3]
  aK18H3 = par[4]
  b = par[5]
  res = c(0,0,0,0,0)
  for(l in levels( rc$Hi.Low) ) 
    for( a1 in levels(factor(rc$Flag)) ) 
      for( a2 in levels(rc$X2nd.Ab) ) {
        #        print(paste(l,a1,a2))
        a = b
        mask = c(0,0,0,0,1)
        if( l == "Hi") { a = a + aHighLow; mask[2] = 1}
         if( a1 == TRUE ) {a = a + aFlagMyc; mask[1] = 1}
        if( a2 == "K36me3") {a = a+aK36H3; mask[3] = 1}
        if( a2 == "K18ac") {a = a+aK18H3; mask[4] = 1}
        rcTemp = rc[rc$Hi.Low == l & rc$Flag == a1 & rc$X2nd.Ab == a2,] 
        pred = exp(a)*rcTemp$X.SC+1
        obs = rcTemp$SC.Only
        z = sum((obs-pred)*obs/(pred**2))
#        print(z)
        res = res - 2*z*mask
      }
  return(res)
}

optim.control = list(fnscale=10**6, maxit = 1000, reltol=10**-10, abstol=10**-10, trace=10)
init.par = c(log(8), log(1), log(1), log(10), log(16))
est = optim(init.par, ll, dll, method="BFGS", control=optim.control)
dll(est$par)
}