## Percentile-based NRI: bin data based on evenly-spaced quantiles of risk
library(boot)
library(survMisc)

## Function for computing net reclassification improvement, taking censoring into account
## Due to Pencina (2011)
## risks.1 and risks.2 are predicted risks for the models being compared
## T = time to event, C = censored/uncensored (vectors of length equal to risks.1 and risks.2)
## t = length of time period under consideration (single numeric quantity)
compute.cNRI <- function(risks.1,risks.2,T,C,t) {
  n <- length(T)
  
  up.class <- (risks.2>risks.1)
  down.class <- (risks.1>risks.2)
  n.u <- sum(up.class) ## Number up-classified
  n.d <- sum(down.class) ## Number down-classified
  
  KM <- survfit(Surv(T,C)~1) ## overall
  p.KM <- 1 - KM$surv[max(which(KM$time<=t))] ## P(event)
  
  if(n.u==0) {
    p.KM.u <- 1
  } else {
    KM.u <- survfit(Surv(T,C)~1,subset=up.class) ## up-classified
    p.KM.u <- 1 - KM.u$surv[max(which(KM.u$time<=t))] ## P(event|up)
  }
  if(n.d==0) {
    p.KM.d <- 1
  } else {
    KM.d <- survfit(Surv(T,C)~1,subset=down.class) ## down-classified
    p.KM.d <- 1 - KM.d$surv[max(which(KM.d$time<=t))] ## P(event|down)
  }
  
  nri.e <- (n.u*p.KM.u - n.d*p.KM.d)/(n*p.KM)
  nri.ne <- (n.d*(1-p.KM.d) - n.u*(1-p.KM.u))/(n*(1-p.KM))
  
  c(cNRI.events=nri.e,cNRI.nonevents=nri.ne,cNRI=nri.e+nri.ne)
}

## Function for computing net reclassification improvement specific for bootstrap sampling
## Used in the following function for bootstrapping confidence intervals for the NRI
compute.cNRI.boot <- function(data,ind) {
  ## Use only data indexed by bootstrap sampling
  risks.1 <- data[ind,]$risks.1
  risks.2 <- data[ind,]$risks.2
  T <- data[ind,]$T
  C <- data[ind,]$C
  t <- data[1,]$t
  ## Compute censored NRI as above 
  compute.cNRI(risks.1, risks.2, T,C,t)
}
  
## Function to return bootstrap confidence intervals for event, non-event, and overall NRI
ci.nri <- function(risks.1, risks.2,T,C,t){
  ## Set number of cores accordingly for parallel computing
  boot <- boot(data = data.frame(risks.1,risks.2,T=T,C=C,t), 
               statistic = compute.cNRI.boot, R = 1000,
               parallel="multicore",ncpus=20)
  ## Compute confidence intervals for event, non-event, and overall NRI
  ci.event <- boot.ci(boot.out = boot, conf=0.95, type = "perc")
  ci.nonevent <- boot.ci(boot.out = boot, conf=0.95, type = "perc", index=2)
  ci.overall <- boot.ci(boot.out = boot, conf=0.95, type = "perc", index=3)
  ## Output lower and upper bounds 
  ci <- cbind(ci.event$perc[4],ci.event$perc[5],ci.nonevent$perc[4],ci.nonevent$perc[5],
              ci.overall$perc[4],ci.overall$perc[5])
  colnames(ci) <- c("Event NRI LB","Event NRI UB","Non-event NRI LB",
                    "Non-event NRI UB","Overall NRI UB","Overall NRI LB")
  ci
}

## Function to compute percentile NRI for estimates from model 2 vs. estimates from model 1
## Bin predicted risks based on evenly-spaced percentiles of risk prior to use in the previously
## stated censored-NRI formula
pct.nri <- function(risks.1,risks.2,T,C,t){
  # Create empty data frame with space for number of categories,
  # event, non-event, and overall NRI with corresponding 95% confidence intervals
  results <- data.frame("num.cat"=as.numeric(),"event"=as.numeric(),"nonevent" = as.numeric(), 
                      "overall"=as.numeric(), "e.lb"=as.numeric(), "e.ub"=as.numeric(),
                      "ne.lb"=as.numeric(), "ne.ub"=as.numeric(),
                      "ov.lb"=as.numeric(),"ov.ub"=as.numeric())
  ## Apply percentile-based NRI to a varied number of categories chosen by the user 
  ## 5, 10, and 15 categories chosen here for illustration 
  num.cat <- c(5,10,15)
  for (i in 1:length(num.cat)){
    # Sort estimates into categories based on quantiles of risk
    risks1.bin <- cut(risks.1, breaks=quantile(risks1,seq(0,1,length.out=num.cat[i]+1)),include.lowest = T,labels=F)
    risks2.bin <- cut(risks.2, breaks=quantile(risks2,seq(0,1,length.out=num.cat[i]+1)),include.lowest = T,labels=F)
  
    ## Calculate NRI
    nri <- compute.cNRI(risks1.bin, risks2.bin,T,C,t) 
    
    ## Calculate 95% bootstrap confidence intervals
    ci <- ci.nri(risks1.bin, risks2.bin,T,C,t)
      
    # Store results
    results[i,] <- c(num.cat[i],nri,ci)
  }
  results
}
