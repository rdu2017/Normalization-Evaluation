
library('MASS')

## Function for correct rounding
round2 = function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}

simu <-function(p,
                spl.size,
                p.0,
                ovdisp,
                suffspl.size,
                seed1=F,
                seed2=F)
{
  if (seed1!=FALSE) set.seed(seed1)
  spls <- round2(quantile(spl.size, runif(length(spl.size))), 0) #empirical distribution of the data
  spls <- as.integer(sort(spls))
  SimuExp<-round(matrix(p, nrow=length(p),ncol=1)%*%spls)
  rownames(SimuExp)<-names(p)
  
  # NB distribution generated counts
  if (seed2!=FALSE) set.seed(seed2)
  SimuCount<-matrix(0,nrow=nrow(SimuExp),ncol=ncol(SimuExp))
  rownames(SimuCount)<-rownames(SimuExp)
  for (i in 1:nrow(SimuCount)) 
    for (j in 1:ncol(SimuCount))
    {
      if (any(p.0$expij==SimuExp[i,j]))
        SimuCount[i,j] <- (1-rbinom(1, size=1, prob=p.0$p.0_mass[p.0$expij==SimuExp[i,j]]))*rnegbin(1,SimuExp[i,j],ovdisp) else
          SimuCount[i,j] <- rnegbin(1,SimuExp[i,j],ovdisp)
    }
  
  # under-sampling simulation        
  spl.size <- colSums(SimuCount)
  for (k in 1:nrow(suffspl.size))
  {
    if (any(spl.size < suffspl.size$size[k]))
      SimuCount[rownames(SimuCount)==suffspl.size$feature[k], spl.size < suffspl.size$size[k]] <- 0
  }
  
  ############# out: SimuExp, SimuCount
  list(SimuExp = SimuExp,
       SimuCount = SimuCount)
}
