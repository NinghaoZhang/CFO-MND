#' @examples 
#' target <- 0.30; mineff <- 0.30; cohortsize = 3; ncohort = 20; init.level = 1
#' prior.para = list(alp.prior = target, bet.prior = 1 - target, 
#'                   alp.prior.eff = 0.5, bet.prior.eff = 0.5)
#' p.true=c(0.05, 0.07, 0.1, 0.12, 0.16)
#' pE.true=c(0.35, 0.45, 0.5, 0.55, 0.75)
#' result <- CFOeff.simu(target, p.true, pE.true, ncohort, init.level, cohortsize,
#'                        prior.para, mineff = mineff, seed = 1)
#' summary(result)
#' plot(result)
#' \donttest{### overly toxic
#' target <- 0.30; mineff <- 0.30; cohortsize = 3; ncohort = 20; init.level = 1
#' prior.para = list(alp.prior = target, bet.prior = 1 - target, 
#'                   alp.prior.eff = 0.5, bet.prior.eff = 0.5)
#' p.true=c(0.55, 0.57, 0.61, 0.62, 0.66)
#' pE.true=c(0.35, 0.45, 0.5, 0.55, 0.75)
#' result <- CFOeff.simu(target, p.true, pE.true, ncohort, init.level, cohortsize,
#'                        prior.para, mineff = mineff, seed = 1)
#' summary(result)
#' plot(result)
#' }
#' \donttest{### low efficacy
#' target <- 0.30; mineff <- 0.30; cohortsize = 3; ncohort = 20; init.level = 1
#' prior.para = list(alp.prior = target, bet.prior = 1 - target, 
#'                   alp.prior.eff = 0.5, bet.prior.eff = 0.5)
#' p.true=c(0.05, 0.07, 0.1, 0.12, 0.16)
#' pE.true=c(0.001, 0.003, 0.004, 0.005, 0.006)
#' result <- CFOeff.simu(target, p.true, pE.true, ncohort, init.level, cohortsize,
#'                        prior.para, mineff = mineff, seed = 1)
#' summary(result)
#' plot(result)
#' }
CFOeff.simu <- function(target, p.true, pE.true, ncohort=10, init.level=1, cohortsize=3,
                        prior.para = list(alp.prior = target, bet.prior = 1 - target, alp.prior.eff = 0.5, 
                                          bet.prior.eff = 0.5), cutoff.eli=0.95, early.stop=0.95, 
                        effearly.stop = 0.9, mineff, seed = NULL) {
  
  ###############################################################################
  ###############define the functions used for main function#####################
  ###############################################################################
  post.prob.fn <- function(target, y, n, alp.prior=0.1, bet.prior=0.1){
    if(n != 0){
      alp <- alp.prior + y 
      bet <- bet.prior + n - y
      res <- 1 - pbeta(target, alp, bet)
    }else{
      res <- NA
    }
    return(res)
  }
  
  
  under.eff.fn <- function(mineff, effearly.stop,prior.para=list())
  {
    args <- c(list(target = mineff), prior.para)
    x <- prior.para$x
    n <- prior.para$n
    alp.prior <- prior.para$alp.prior.eff
    bet.prior <- prior.para$bet.prior.eff
    ppE <- 1 - post.prob.fn(mineff, x, n, alp.prior, bet.prior)
    if ((ppE >= effearly.stop) & (n >= 3)) {
      return(TRUE)
    }else{
      return(FALSE)
    }
  }
  
  OBD.level <- function(phi, mineff, p.true, pE.true) {
    if (p.true[1] > phi + 0.1) {
      OBD <- 99
      return(OBD)
    }
    OBD <- which.min(abs(phi - p.true))
    eff.idxs <- mineff > pE.true[1:OBD]
    if (sum(eff.idxs) == OBD) {
      OBD <- 99
      return(OBD)
    }
    
    OBD <- which.max(pE.true[1:OBD])
    return(OBD)
  }
  
  # compute the marginal prob when lower < phiL/phiC/phiR < upper
  # i.e., Pr(Y=y|lower<target<upper)
  overdose.fn <- function(target, threshold, prior.para=list()){
    y <- prior.para$y
    n <- prior.para$n
    alp.prior <- prior.para$alp.prior
    bet.prior <- prior.para$bet.prior
    pp <- post.prob.fn(target, y, n, alp.prior, bet.prior)
    if ((pp >= threshold) & (prior.para$n >= 3)) {
      return(TRUE)
    }else{
      return(FALSE)
    }
  }
  
  moveprobs <- function(ad.xs, ad.ns, alp.prior, bet.prior){
    alps <- ad.xs + alp.prior
    bets <- ad.ns - ad.xs + bet.prior
    nd <- length(ad.xs)
    
    Nsps <- 10000
    sps.list <- list() 
    for (i in 1:nd){
      sps.list[[i]] <- rbeta(Nsps, alps[i], bets[i])
    }
    
    spss <- do.call(rbind, sps.list)
    argMaxs <- apply(spss, 2, which.max)
    probs <- as.vector(table(argMaxs))/Nsps
    
    probs
  }
  
  ###############################################################################
  ############################MAIN DUNCTION######################################
  ###############################################################################  
  
  
  set.seed(seed)
  tadd.args <- prior.para
  earlystop <- 0
  stopreason <- NULL
  ndose <- length(p.true)
  doselist <- rep(0, ncohort)
  cidx <- init.level
  OBDprob <- NULL
  
  tys <- rep(0, ndose) # number of DLT responses for different doses.
  txs <- rep(0, ndose) # number of efficacy responses for different doses.
  tns <- rep(0, ndose) # number of subject for different doses.
  tover.doses <- rep(0, ndose) # Whether each dose is too toxic or not, 1 yes.
  tunder.effs <- rep(0, ndose) # Whether the dose is not efficacious or not, 1 yes
  DLTlist <- c()
  EFFlist <- c()
  # if a dose is not efficacious enough or it is too toxic, it is would be eliminated from the admissible set.
  
  for (i in 1:ncohort) {
    pc <- p.true[cidx] 
    pEc <- pE.true[cidx] 
    doselist[i] <- cidx
    
    # sample from current dose
    cres <- rbinom(cohortsize, 1, pc)
    cEres <- rbinom(cohortsize, 1, pEc)
    DLTlist <- c(DLTlist, cres)
    EFFlist <- c(EFFlist, cEres)
    
    # update results
    tys[cidx] <- tys[cidx] + sum(cres)#The number of observed DLT
    txs[cidx] <- txs[cidx] + sum(cEres)#The number of efficacy outcomes
    tns[cidx] <- tns[cidx] + cohortsize#The number of patient at dose level k
    
    cy <- tys[cidx]
    cx <- txs[cidx]
    cn <- tns[cidx]
    
    prior.para <- c(list(y = cy, n = cn, x = cx, tys = tys, txs = txs, tns = tns, cidx = cidx), tadd.args)
    
    if (overdose.fn(target, cutoff.eli, prior.para)) {
      tover.doses[cidx:ndose] <- 1
    }
    
    if (under.eff.fn(mineff, effearly.stop, prior.para)) {
      tunder.effs[cidx] <- 1
    }else{
      tunder.effs[cidx] <- 0
    }
    
    
    if (tover.doses[1] == 1) {
      stopreason <- "overly_toxic"
      earlystop <- 1
      break()
    }
    if (sum(tunder.effs[tover.doses == 0]) == sum(tover.doses == 0)){
      stopreason <- "low_efficacy"
      earlystop <- 1
      break()
      
    }
    
    nextinfo <- CFOeff.next(target, axs = txs, ays = tys, ans = tns, cidx, prior.para, cutoff.eli, early.stop, effearly.stop, mineff)
    cidx <- nextinfo$nextdose
    if (cidx == 99){
      if (nextinfo$decision == "stop_for_tox"){
        stopreason <- "overly_toxic"
      }else{
        stopreason <- "low_efficacy"
      }
      earlystop <- 1
      break()
    }
    if (i == ncohort){
      para = list(alp.prior.eff = prior.para$alp.prior.eff, bet.prior.eff = prior.para$bet.prior.eff)
      result <- CFOeff.selectobd(target, txs, tys, tns, para, mineff, effearly.stop)
      OBDprob <- result$OBD.probs
    }
    
  }
  
  
  if (earlystop == 0) {
    MTD <- CFO.selectmtd(target, tns, tys)$MTD
    if ( (MTD == 99) | (sum(tunder.effs[1:MTD]) == MTD)) {
      OBD <- 99
    }else{
      OBD.probs <- moveprobs(txs[1:MTD], tns[1:MTD], prior.para$alp.prior.eff, prior.para$bet.prior.eff)
      OBD <- which.max(OBD.probs)
    }
    
  }else{
    OBD <- 99
  }
  
  tobd = OBD.level(target, mineff, p.true, pE.true)
  correct <- 0
  if (OBD == tobd){
    correct <- 1
  }
  
  
  
  ptoxic <- sum(tns[which(p.true > target)])/(ncohort*cohortsize)
  
  out <- list(OBD = OBD, MTD = MTD, target = target, npatients = tns, neff = txs, ntox = tys, pE.true = pE.true, p.true = p.true, 
              cohortdose = doselist, ptoxic = ptoxic, patientDLT = DLTlist, patienteff = EFFlist, 
              over.doses = tover.doses, under.eff = tunder.effs, correct = correct, OBDprob = OBDprob, 
              sumDLT = sum(DLTlist), sumeff = sum(EFFlist), earlystop = earlystop, stopreason = stopreason, 
              class = "phaseI/II")
  class(out) <- c("cfo_eff_trial", "cfo")
  return(out)
  
}