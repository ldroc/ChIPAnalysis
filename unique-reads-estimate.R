uniq.reads <- function(x,n) n*(1-ppois(0,x/n))

x = 10^seq(4,8,by=.025)
p = ggplot()
for( n in seq(5,10,by=.5))
  p = p + geom_line(data=data.frame(x=x,y=uniq.reads(x,10^n)), 
                    mapping=aes(x=x,y=y),
                    size=1.5)

p = p + scale_x_log10("# reads") + scale_y_log10("# unique reads")
p = p + geom_abline(slope=1,color="blue",linetype=2)
