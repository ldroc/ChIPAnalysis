# 
# Remove duplicate reads from sorted bam record
#
# Assumes pos, mpos, rname appear in the file
RemoveDuplicateFragments <- function( b, Is ) {
  
  dat <- data.frame( rname = b$rname[Is], pos = b$pos[Is], mpos = b$mpos[Is])
  
  # filter applies linear transformation to successive entries
  f <- (filter(dat,c(-1,1))!= 0)
  
  # choode the ones where there is a differnece in at least one column
  r <- c(TRUE,apply(f,1,any))
  
  # choose only unique reads
  newdat = dat[r[1:length(r)-1],]
  
  # compute difference between unique occurences to get repeat numbers
  rep <- filter( c(which(r == TRUE),length(r)+1), c(1,-1) )
  rep <- rep[1:length(rep)-1]
  newdat <- cbind( newdat , rep )
}
