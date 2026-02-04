
CFO.simu <- function(design, target, p.true, init.level=1, ncohort, cohortsize,
                     prior.para=list(alp.prior=target, bet.prior=1-target),
                     cutoff.eli=0.95, early.stop=0.95, seed=NULL){
  ###############################################################################
  ###############define the functions used for main function#####################
  ###############################################################################
  # posterior probability of pj >= phi given data
  post.prob.fn <- function(phi, y, n, alp.prior=0.1, bet.prior=0.9){
    if(n != 0){
      alp <- alp.prior + y 
      bet <- bet.prior + n - y
      res <- 1 - pbeta(phi, alp, bet)
    }else{
      res <- NA
    }
    return(res)
  }
  
  overdose.fn <- function(phi, threshold, prior.para=list()){
    y <- prior.para$y
    n <- prior.para$n
    alp.prior <- prior.para$alp.prior
    bet.prior <- prior.para$bet.prior
    pp <- post.prob.fn(phi, y, n, alp.prior, bet.prior)
    # print(data.frame("prob of overdose" = pp))
    if ((pp >= threshold) & (prior.para$n>=3)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }
  
  MTD.level <- function(phi, p.true){
    if (p.true[1]>phi+0.1){
      MTD <- 99
      return(MTD)
    }
    MTD <- which.min(abs(phi - p.true))
    return(MTD)
  }
  ###############################################################################
  ############################MAIN DUNCTION######################################
  ###############################################################################
  if (is.null(prior.para$alp.prior)){
    prior.para <- c(prior.para, list(alp.prior=target, bet.prior=1-target))
  }
  alp.prior <- prior.para$alp.prior
  bet.prior <- prior.para$bet.prior
  
  set.seed(seed)
  p_est <- 0
  p_overdose <- 0
  earlystop <- 0
  ndose <- length(p.true)
  doselist <- rep(0, ncohort)
  currdose <- init.level
  
  ays <- rep(0, ndose) # number of responses for different doses.
  ans <- rep(0, ndose) # number of subject for different doses.
  tover.doses <- rep(0, ndose) # Whether each dose is overdosed or not, 1 yes
  DLTlist <- c()
  for (i in 1:ncohort){
    pc <- p.true[currdose]
    doselist[i] <- currdose
    # set.seed(seed+i)
    # sample from current dose
    cres <- rbinom(cohortsize, 1, pc)
    DLTlist <- c(DLTlist, cres)
    # update results
    ays[currdose] <- ays[currdose] + sum(cres)
    ans[currdose] <- ans[currdose] + cohortsize
    
    cy <- ays[currdose]
    cn <- ans[currdose]
    
    prior.para <- c(list(y=cy, n=cn), list(alp.prior=alp.prior, bet.prior=bet.prior))
    
    if (overdose.fn(target, cutoff.eli, prior.para)){
      tover.doses[currdose:ndose] <- 1
    }
    
    if (currdose == 1){
      if (cutoff.eli != early.stop) {
        cy <- ays[1]
        cn <- ans[1]
        prior.para <- c(list(y=cy, n=cn), list(alp.prior=alp.prior, bet.prior=bet.prior))
        if (overdose.fn(target, early.stop, prior.para)){
          tover.doses[1:ndose] <- 1
        }
      }
    }
    
    
    if (tover.doses[1] == 1){
      earlystop <- 1
      break()
    }
    
    prior.para <- c(list(alp.prior=alp.prior, bet.prior=bet.prior))
    if (design == 'CFO'|| design == 'rCFO' || design == 'pCFO'){
      # the results for current 3 dose levels
      if (currdose!=1){
        cys <- ays[(currdose-1):(currdose+1)]
        cns <- ans[(currdose-1):(currdose+1)]
      }else{
        cys <- c(NA, ays[1:(currdose+1)])
        cns <- c(NA, ans[1:(currdose+1)])
      }
      if (design == 'CFO'){
        currdose <- CFO.next(target, cys, cns, currdose, prior.para, cutoff.eli, early.stop)$nextdose
      } else if (design == 'rCFO'){
        if (!is.null(seed)){seed <- seed+i}
        currdose <- rCFO.next(target, cys, cns, currdose, prior.para, cutoff.eli, early.stop, seed)$nextdose
      } else if (design == 'pCFO'){
        currdose <- pCFO.next(target, cys, cns, currdose, prior.para, cutoff.eli, early.stop)$nextdose
      }
    }else if (design == 'aCFO'){
      currdose <- aCFO.next(target, ays, ans, currdose, prior.para, cutoff.eli, early.stop)$nextdose
    }else{
      stop("The input design is invalid; it can only be set as 'CFO', 'aCFO', 'rCFO' or 'pCFO'.")
    }
    
    if (i == ncohort){
      result <- CFO.selectmtd(target, ans, ays, prior.para, cutoff.eli, early.stop, verbose = TRUE)
      p_est <- result$p_est
      p_overdose <- result$p_overdose
    }
  }
  
  if (earlystop==0){
    MTD <- CFO.selectmtd(target, ans, ays, prior.para, cutoff.eli, early.stop, verbose=FALSE)$MTD
  }else{
    MTD <- 99
  }
  
  tmtd <- MTD.level(target, p.true)
  correct <- 0
  if (MTD == tmtd){
    correct <- 1
  }
  
  ptoxic <- sum(ans[which(p.true > target)])/(ncohort*cohortsize)
  
  out<-list(target=target, MTD=MTD, correct=correct, npatients=ans, ntox=ays,
            over.doses=tover.doses, cohortdose=doselist, ptoxic=ptoxic, patientDLT=DLTlist,
            sumDLT=sum(DLTlist), earlystop=earlystop, p_est = p_est, p_overdose = p_overdose)
  class(out) <- c("cfo_trial", "cfo")
  return(out)
}
