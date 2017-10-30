
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

dataLearning <- function(dat)
{
  dat <- dat[!apply(dat, 1, function(x) all(x==0)), ] # remove a feature with all zeros
  
  dat <- dat[order(-rowSums(dat)),] # order features from most to least
  
  dat <- dat[,order(colSums(dat))] # order sample from least to most 
  
  p <- rowSums(dat)/sum(rowSums(dat)) # est of relative abundance
  
  mu <- colSums(dat) # est of sample scale
  
  muij <- t(matrix(mu, ncol=1) %*% matrix(p, nrow=1)) # matrix of expected counts 
  
  muij <- round2(muij, 0) # round the expected counts to integers (1)
  
  #################### (2)
  # if cumsum > 3, the feature is thought to appear, otherwise not
  suffspl.inx <- apply(t(apply(muij, 1, cumsum)), 1, function(x) which(x>3)[1]) 
  
  suffspl.size <- data.frame(feature = rownames(dat),
                             size = as.numeric(sapply(suffspl.inx, 
                                                      function(x) if(!is.na(x)) mu[x] else NA)))
  
  #################### (3)
  ovdisp.fit <- NULL
  for (expij in unique(sort(as.vector(muij)))[-1]) # not to include 0
  {
    locij <- which(muij==expij, arr.ind = T)
    if (nrow(locij)>50) # only take the expij with the number of a expij greater than 50
    {
      obsij <- dat[locij]
      
      fitij <- fitdistr(obsij, 'negative binomial')
      
      ovdisp.fit <- rbind(ovdisp.fit,
                          data.frame(expij = expij,
                                     mu.fitted = fitij$estimate['mu'],
                                     size.fitted = fitij$estimate['size'],
                                     stringsAsFactors = F))
    }
  }
  rownames(ovdisp.fit) <- NULL
  ovdisp <- mean(ovdisp.fit$size.fitted)
  
  
  #################### (4)
  p.0 <- NULL
  for (expij in unique(sort(as.vector(muij)))[-1])
  {
    locij <- which(muij==expij, arr.ind = T)
    
    if (nrow(locij)>50)
    {
      #print(expij)
      obsij <- dat[locij]
      
      p.0_nb <- dnbinom(0, size=ovdisp, mu=expij) # probability of 0 estimated from nb distribution
      p.0_obs <- length(which(obsij==0))/length(obsij) # probability of 0 observed
      
      if (p.0_obs > p.0_nb) p.0_mass <- (p.0_obs - p.0_nb) else mass0 <- NA
      
      p.0 <- rbind(p.0,
                   data.frame(expij = expij,
                              p.0_mass = p.0_mass,
                              stringsAsFactors = F))
    }
  }
  
  ############# out: mu, p, ovdisp, suffspl.size, p.0
  list(mu = mu,
       p = p,
       ovdisp = ovdisp,
       suffspl.size = suffspl.size[!is.na(suffspl.size$size),],
       p.0 = p.0)
}
