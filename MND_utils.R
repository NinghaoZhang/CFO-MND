source("CFOcode/CFO.selectmtd.R")
source("CFOcode/CFO.simu.R")
source("CFOcode/CFO.next.R")
source("selectMND.R")
library(ggplot2)
library(ggpubr)
OBD.level <- function(phi, mineff, p.true, pE.true) {
  #The lowest dose level is overly toxic
  if (p.true[1] > phi + 0.1) {
    OBD <- 99
    return(OBD)
  }
  #select MTD
  MTD <- which.min(abs(phi - p.true))
  
  #Rule out the situation that all dose levels have low efficacy rate
  eff.idxs <- mineff > pE.true[1:MTD]
  if (sum(eff.idxs) == MTD) {
    OBD <- 99
    return(OBD)
  }
  #Select OBD which has the highest efficacy rate
  OBD <- which.max(pE.true[1:MTD])
  return(OBD)
}

MTD.level <- function(phi, p.true){
  if (p.true[1]>phi+0.1){
    MTD <- 99
    return(MTD)
  }
  MTD <- which.min(abs(phi - p.true))
  return(MTD)
}


MND.level <- function(dose.val, q.true, gamma, p.true = NULL, target = NULL, mineff = NULL, obd = NULL, placebo = FALSE, estimate = FALSE){
  OBD <- obd
  ndose <- length(p.true)
  
  if(is.null(OBD)){
    OBD <- OBD.level(target, mineff, p.true, q.true)
    if (placebo & OBD == 1){
      MND = 99
      return(MND)
    }else if (OBD == 99){
      MND = OBD
      return(MND)
    }
  }
  
  q.tilde <- q.true/(dose.val)^gamma
  q.lb <- 0.9*q.true[OBD]
  idx <- which(q.true >= q.lb & q.true <= q.true[OBD])
  MND <- idx[which.max(q.tilde[idx])]
  return(MND)
}

pos.betasample <- function(n, cx, cn, alp.prior, bet.prior){
  alp.pos <- alp.prior + cx
  bet.pos <- bet.prior + cn - cx
  pos.sample <- rbeta(n, alp.pos, bet.pos)
  pos.sample
}

efficacy.pos.mean <- function(axs, ans, alp.prior.eff, bet.prior.eff){
  return((axs + alp.prior.eff)/(ans + alp.prior.eff + bet.prior.eff))
} 

posterior.variance <- function(axs, ans, alp.prior.eff, bet.prior.eff){
  alp.pos <- axs + alp.prior.eff
  beta.pos <- ans - axs + bet.prior.eff
  var <- (alp.pos*beta.pos)/((alp.pos+beta.pos)^2*(alp.pos+beta.pos+1))
  return (var)
}

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

overdose.fn <- function(phi, threshold, cy, cn, prior.para=list()){
  y <- cy
  n <- cn
  alp.prior <- prior.para$alp.prior
  bet.prior <- prior.para$bet.prior
  pp <- post.prob.fn(phi, y, n, alp.prior, bet.prior)
  if ((pp >= threshold) & (n>=3)){
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
  argMaxs <- factor(argMaxs, levels = 1:nd)
  probs <- as.vector(table(argMaxs))/Nsps
  return(probs)
}

under.eff.fn <- function(mineff, cx, cn, effearly.stop, prior.para=list())
{
  x <- cx
  n <- cn
  alp.prior <- prior.para$alp.prior.eff
  bet.prior <- prior.para$bet.prior.eff
  ppE <- 1 - post.prob.fn(mineff, x, n, alp.prior, bet.prior)
  if ((ppE >= effearly.stop)) {
    return(TRUE)
  }else{
    return(FALSE)
  }
}




MND.simu <- function(p.true, q.true, dose.val, cohortsize.stage1, cohortsize.stage2,
                     ncohort.stage1, ncohort.stage2, target, gamma, mineff, init.level, placebo = FALSE, realtrial = FALSE,
                     prior.para = list(alp.prior = target, bet.prior = 1 - target, alp.prior.eff = 0.5, bet.prior.eff = 0.5)){
  
  count <- 0
  earlystop.count <- 0
  ########Stage 1 Forward (Go through a round)########
  earlystop <- 0
  earlystop.eff <- 0
  prop.diff <- NULL
  decision <- NULL
  currdose <- 1
  decision.tox <- NULL
  currdose <- 1
  dose.list <- NULL
  ndose <- length(p.true)
  ays <- rep(0, ndose)#tox
  axs <- rep(0, ndose)#eff
  ans <- rep(0, ndose)#number of patients
  tover.doses <- rep(0, ndose)
  DLTlist <- c()
  EFFlist <- c()
  q.hat <- NA
  MTD <- NA
  OBD <- NA
  MND <- NA
  canset <- NA
  canans <- NA
  canaxs <- NA
  admset <- NA
  dose.val.adm <- NA
  if(realtrial){
    sink("results/output")
  }
  for (i in 1:ncohort.stage1){
    if(realtrial){
      statement <- paste0("------------------stage1_cohort: ",i,"------------------------")
      print(statement)
    }
    dose.list <- c(dose.list, currdose)
    p <- p.true[currdose]
    q <- q.true[currdose]
    currtox <- rbinom(cohortsize.stage1, 1, p)
    curreff <- rbinom(cohortsize.stage1, 1, q)
    
    
    DLTlist <- c(DLTlist, currtox)
    EFFlist <- c(EFFlist, curreff)
    
    
    ays[currdose] <- ays[currdose] + sum(currtox)
    axs[currdose] <- axs[currdose] + sum(curreff)
    ans[currdose] <- ans[currdose] + cohortsize.stage1
    
    #Cutoff and early stop
    cy <- ays[currdose]
    cn <- ans[currdose]
    
    if (placebo){
      if (overdose.fn(target, cutoff.eli, cy, cn, prior.para)){
        tover.doses[currdose:ndose] <- 1
      }
      if (currdose == 1){
        if (cutoff.eli != early.stop) {
          cy <- ays[1]
          cn <- ans[1]
          if (overdose.fn(target, early.stop, cy, cn, prior.para)){
            tover.doses[1:ndose] <- 1
          }
        }
      }
      if (tover.doses[2] == 1){
        earlystop <- 1
        break
      }
    }else{
      if (overdose.fn(target, cutoff.eli, cy, cn, prior.para)){
        tover.doses[currdose:ndose] <- 1
      }
      if (currdose == 1){
        if (cutoff.eli != early.stop) {
          cy <- ays[1]
          cn <- ans[1]
          if (overdose.fn(target, early.stop, cy, cn, prior.para)){
            tover.doses[1:ndose] <- 1
          }
        }
      }
      if (tover.doses[1] == 1){
        earlystop <- 1
        break
      }
    }
    
    #CFO algorithm
    if (currdose!=1){
      cys <- ays[(currdose-1):(currdose+1)]
      cns <- ans[(currdose-1):(currdose+1)]
    }else{
      cys <- c(NA, ays[1:(currdose+1)])
      cns <- c(NA, ans[1:(currdose+1)])
    }
    cfores <- CFO.next(target, cys, cns, currdose, prior.para, cutoff.eli = 0.95, early.stop = 0.95)
    currdose <- cfores$nextdose
    if(realtrial){
      print("The CFO rule odds ratio:\n")
      print(cfores$odds.ratio)
      print("The decision:\n")
      print(cfores$decision)
      print("Next dose:\n")
      print(currdose)
    }
    if (earlystop == 0){
      MTD <- CFO.selectmtd(target, ans, ays, prior.para, cutoff.eli, early.stop, verbose = FALSE)$MTD
    }else{
      MTD <- 99
    }
  }
  
  
  
  
  
  
  ########Stage2 Backward##########Randomization
  if (earlystop == 0) {
    #Remove the over toxicity dose based on the results given by CFO. Using admissible set to control the dose level tested each round.
    #Select the dose level that have the biggest probability to be OBD
    effprobs <- moveprobs(axs[1:MTD], ans[1:MTD], alp.prior.eff, bet.prior.eff)
    prior.para <- c(list(alp.prior.eff = alp.prior.eff, bet.prior.eff = bet.prior.eff), prior.para)
    admset <- 1:MTD
    OBD.hat <- which.max(effprobs)
    currdose <- OBD.hat
    ndose <- MTD
    alp.prior <- 0.5
    bet.prior <- 0.5
    probs <- NA
    if(realtrial){
      print("MTD是")
      print(MTD)
      print("OBD是")
      print(OBD.hat)
      print(effprobs)
    }
    if (under.eff.fn(mineff = 0.15, axs[OBD.hat], ans[OBD.hat], 0.9, prior.para)){
      earlystop <- 1
    }else{
      #Stage2
      for (i in 1:ncohort.stage2){
        if(realtrial){
          statement <- paste0("------------------stage2_cohort: ",i,"------------------------")
          print(statement)
        }
        dose.list <- c(dose.list, currdose)
        overdose <- 0
        #Generate the simulated trial outcome based on the given scenerio 
        pt <- p.true[currdose]
        q <- q.true[currdose]
        currtox <- rbinom(cohortsize.stage2, 1, pt)
        curreff <- rbinom(cohortsize.stage2, 1, q)
        
        
        DLTlist <- c(DLTlist, currtox)
        EFFlist <- c(EFFlist, curreff)
        
        
        ays[currdose] <- ays[currdose] + sum(currtox)
        axs[currdose] <- axs[currdose] + sum(curreff)
        ans[currdose] <- ans[currdose] + cohortsize.stage2
        
        #Safety rule & efficacy
        cy <- ays[currdose]
        cn <- ans[currdose]
        cx <- axs[currdose]
        
        if (overdose.fn(target, cutoff.eli, cy, cn, prior.para)) {
          overdose <- 1
        }
        
        if(placebo){
          if (currdose == 2 & overdose == 1){
            earlystop <- 1
            break()
          }else if(overdose == 1){
            admset <- admset[admset < currdose]
          }
        }else{
          if (currdose == 1 & overdose == 1){
            earlystop <- 1
            break()
          }else if(overdose == 1){
            admset <- admset[admset < currdose]
          }
        }
        
        
        # if (under.eff.fn(mineff, effthreshold, cx, cn, prior.para)) {
        #   admset <- admset[admset > currdose]
        # }
        
        axs.adm <- axs[admset]
        ans.adm <- ans[admset]
        dose.val.adm <- dose.val[admset]
        
        # Randomization by using q.tilde
        q.hat <- efficacy.pos.mean(axs.adm, ans.adm, alp.prior.eff, bet.prior.eff)
        q.tilde <- q.hat/(dose.val.adm)^gamma
        prob.sample <- q.tilde/sum(q.tilde)
        # print(admset)
        
        currdose <- sample(admset, size = 1, prob = prob.sample)
        if(realtrial){
          print("The admissible set:\n")
          print(admset)
          print("Sample probability:\n")
          print(prob.sample)
          print("The next dose level:\n")
          print(currdose)
        }
        
      }
    }
    ######Select MND after the trial finshed#######
    
    
    
    probs <- moveprobs(axs[admset], ans[admset], alp.prior.eff, bet.prior.eff)
    #tMND <- MND.level(dose.val, q.true, gamma, p.true, target, mineff)
    OBD <- admset[which.max(probs)]
    canset <- admset[admset <= OBD]
    canaxs <- axs[canset]
    canans <- ans[canset]
    dose.val.adm <- dose.val[canset]
    prior.para <- c(prior.para, list(alp.prior = alp.prior.eff, bet.prior = bet.prior.eff))
    
    res <- selectMND(canset, canans, canaxs, prior.para)
    q.hat <- res$phat
    if(realtrial){
      print("q.hat for MND")
      print(q.hat)
    }
  }
  if(realtrial){
    sink()
  }
  ATE <- NA
  if (earlystop == 0) {
    MND <- MND.level(dose.val.adm, q.hat, gamma, mineff = 0.15, obd = OBD, placebo = placebo, estimate = TRUE)
    if(placebo){
      MND.mean <- efficacy.pos.mean(axs[MND],ans[MND], alp.prior.eff, bet.prior.eff)
      placebo.mean <- efficacy.pos.mean(axs[1],ans[1], alp.prior.eff, bet.prior.eff)
      ATE <- MND.mean - placebo.mean
    }
  }
  else{
    MND <- 99
  }
  if (placebo){
    res <- list(MND = MND, p.true = p.true, q.true = q.true, ate = ATE, canset = canset,
                earlystop = earlystop, npatients = ans, ntox = ays, neff = axs, dose.list = dose.list, DLTlist = DLTlist, EFFlist = EFFlist)
  }else{
    res <- list(MND = MND, p.true = p.true, q.true = q.true, dose = dose.val, canset = canset,
                earlystop = earlystop, npatients = ans, ntox = ays, neff = axs, dose.list = dose.list, DLTlist = DLTlist, EFFlist = EFFlist)
  }
  return(res)
}



###########generate random scenario(beta uniform)#############
raw.gen.tox <- function(ndose, phi){
  gen.p.true <- rep(0, ndose)
  mtd <- sample.int(ndose, 1)
  gen.p.true[mtd] <- rbeta(1, phi, (1-phi))
  
  if (mtd > 1){
    raw.data <- runif(mtd - 1, 0, gen.p.true[mtd])
    gen.p.true[1:(mtd - 1)] <- sort(raw.data)
  }
  
  if (mtd < ndose){
    raw.data <- runif(ndose - mtd ,gen.p.true[mtd], 1)
    gen.p.true[(mtd + 1): ndose] <- sort(raw.data)
  }
  res <- list(p.true = gen.p.true, mtd = mtd)
  return(res)
}


gen.rand.doses <- function(ndose, target, mu1=0.55, mu2=0.55, inv.fn=qnorm, fn=pnorm, MTD=NA){
  sigma0 <- 0.05
  sigma1 <- 0.35
  sigma2 <- 0.35
  raw.p.trues <- rep(0, ndose)
  if (is.na(MTD)){
    mtd.level <- sample.int(ndose, 1)
  }else{
    mtd.level <- MTD
  }
  
  raw.p.trues[mtd.level] <- rnorm(1, inv.fn(target), sigma0)
  
  if (mtd.level != 1){
    eps.minus <- rnorm(1, mu1, sigma1)
    if (raw.p.trues[mtd.level] > inv.fn(target)){
      dlt <- raw.p.trues[mtd.level] - inv.fn(2*target-fn(raw.p.trues[mtd.level]))
    }else{
      dlt <- 0 
    }
    raw.p.trues[mtd.level-1] <- raw.p.trues[mtd.level] - (dlt +  eps.minus**2)
  }
  
  if (mtd.level != ndose){
    eps.plus <- rnorm(1, mu2, sigma2)
    if (raw.p.trues[mtd.level] < inv.fn(target)){
      dlt <- inv.fn(2*target-fn(raw.p.trues[mtd.level])) - raw.p.trues[mtd.level]
    }else{
      dlt <- 0 
    }
    raw.p.trues[mtd.level+1] <- raw.p.trues[mtd.level] + (dlt +  eps.plus**2)
  }
  
  if ((mtd.level-2)>0){
    for (i in (mtd.level-2):1){
      eps.minus <- rnorm(1, mu1, sigma1)
      raw.p.trues[i]  <- raw.p.trues[i+1] - eps.minus**2
    }
  }
  if ((mtd.level+2)<=ndose){
    for (j in (mtd.level+2):ndose){
      eps.plus <- rnorm(1, mu2, sigma2)
      raw.p.trues[j] <- raw.p.trues[j-1] + eps.plus**2
    }
  }
  p.trues <- fn(raw.p.trues)
  list(p.trues=p.trues, mtd.level=mtd.level)
}

gen.rand.doses.plateau <- function(ndose, phi, psi, psi.U, mu1, mu2){
  # args:
  # ndose: number of dose levels
  # phi: DLT target rate
  # psi: the minimal acceptable efficacy rate 
  # psi.U: the upper bound of the efficacy rate
  # mu1, mu2: parameters for phase I trial
  k.OBD <- sample.int(ndose, 1)
  k.MTD <- sample(k.OBD:ndose, 1)
  
  ps <- gen.rand.doses(ndose, phi, mu1, mu2, MTD=k.MTD)
  
  q.OBD <- runif(1, psi, psi.U)
  if (k.OBD == 1){
    qs <- rep(q.OBD, ndose)
  }else if (k.OBD == ndose){
    qs.l <- sort(runif(ndose-1, 0, q.OBD))
    qs <- c(qs.l, q.OBD)
  }else {
    qs.l <- sort(runif(k.OBD-1, 0, q.OBD))
    qs.u <- rep(q.OBD, ndose-k.OBD)
    qs <- c(qs.l, q.OBD, qs.u)
  }
  
  res <- list(
    qs=qs,
    ps=ps$p.trues, 
    k.MTD=k.MTD, 
    k.OBD=k.OBD
  )
  res
}

gen.rand.doses.umbrella <- function(ndose, phi, psi, psi.U, mu1, mu2){
  # args:
  # ndose: number of dose levels
  # phi: DLT target rate
  # psi: the minimal acceptable efficacy rate 
  # psi.U: the upper bound of the efficacy rate
  # mu1, mu2: parameters for phase I trial
  
  k.OBD <- sample.int(ndose, 1)
  k.MTD <- sample(k.OBD:ndose, 1)
  
  ps <- gen.rand.doses(ndose, phi, mu1, mu2, MTD=k.MTD)
  
  q.OBD <- runif(1, psi, psi.U)
  if (k.OBD == 1){
    qs.u <- sort(runif(ndose-1, 0, q.OBD), decreasing=TRUE)
    qs <- c(q.OBD, qs.u)
  }else if (k.OBD == ndose){
    qs.l <- sort(runif(ndose-1, 0, q.OBD))
    qs <- c(qs.l, q.OBD)
  }else {
    qs.l <- sort(runif(k.OBD-1, 0, q.OBD))
    qs.u <- sort(runif(ndose-k.OBD, 0, q.OBD), decreasing=TRUE)
    qs <- c(qs.l, q.OBD, qs.u)
  }
  
  res <- list(
    qs=qs,
    ps=ps$p.trues, 
    k.MTD=k.MTD, 
    k.OBD=k.OBD
  )
  res
}

gen.rand.doses.plateau <- function(ndose, phi, psi, psi.U, mu1, mu2){
  # args:
  # ndose: number of dose levels
  # phi: DLT target rate
  # psi: the minimal acceptable efficacy rate 
  # psi.U: the upper bound of the efficacy rate
  # mu1, mu2: parameters for phase I trial
  k.OBD <- sample.int(ndose, 1)
  k.MTD <- sample(k.OBD:ndose, 1)
  
  ps <- gen.rand.doses(ndose, phi, mu1, mu2, MTD=k.MTD)
  
  q.OBD <- runif(1, psi, psi.U)
  if (k.OBD == 1){
    qs <- rep(q.OBD, ndose)
  }else if (k.OBD == ndose){
    qs.l <- sort(runif(ndose-1, 0, q.OBD))
    qs <- c(qs.l, q.OBD)
  }else {
    qs.l <- sort(runif(k.OBD-1, 0, q.OBD))
    qs.u <- rep(q.OBD, ndose-k.OBD)
    qs <- c(qs.l, q.OBD, qs.u)
  }
  
  res <- list(
    qs=qs,
    ps=ps$p.trues, 
    k.MTD=k.MTD, 
    k.OBD=k.OBD
  )
  res
}

gen.rand.doses.increase <- function(ndose, phi, psi, psi.U, mu1, mu2){
  # args:
  # ndose: number of dose levels
  # phi: DLT target rate
  # psi: the minimal acceptable efficacy rate 
  # psi.U: the upper bound of the efficacy rate
  # mu1, mu2: parameters for phase I trial
  k.OBD <- sample.int(ndose, 1)
  k.MTD <- k.OBD
  
  ps <- gen.rand.doses(ndose, phi, mu1, mu2, MTD=k.MTD)
  qs <- sort(runif(ndose,psi,psi.U))
  res <- list(
    qs=qs,
    ps=ps$p.trues, 
    k.MTD=k.MTD, 
    k.OBD=k.OBD
  )
  res
}
# raw.gen.eff <- function(ndose){
#   shape <- sample.int(3,1)
#   q.true <- NA
#   if (shape == 1){
#     q.true <- gen.rand.doses.plateau
#   }else if(shape == 2){
#     q.true <- gen.rand.doses.umbrella
#   }else{
#     q.true <- sort(runif(ndose,))
#   }
# }


prob.diff.fn <- function(res, target){
  p.trues <- res$p.true
  mtd <- res$mtd
  ndose <- length(p.trues)
  
  diffs <- c()
  if (mtd!=1){
    diffs <- c(abs(p.trues[mtd-1]-target), diffs)
  }
  
  if (mtd!=ndose){
    diffs <- c(abs(p.trues[mtd+1]-target), diffs)
  }
  min(diffs)
}

scs.gen <- function(ndose, dose.val, ndata, gamma, phi, psi, psi.U, mu1, mu2){
  scs <- list()
  i = 1
  while(i<=ndata){
    shape <- sample.int(3,1)
    if (shape == 1){
      res <- gen.rand.doses.plateau(ndose, phi, psi, psi.U, mu1, mu2)
    }else if(shape == 2){
      res <- gen.rand.doses.umbrella(ndose, phi, psi, psi.U, mu1, mu2)
    }else{
      res <- gen.rand.doses.increase(ndose, phi, psi, psi.U, mu1, mu2)
    }
    MND <- MND.level(dose.val, res$qs, gamma, obd = res$k.OBD)
    res$MND <- MND
    scs[[i]] <- res
    i = i+1
  }
  return(scs)
}

scs.gen.placebo <- function(ndose, dose.val, ndata, gamma, phi, psi, psi.U, mu1, mu2){
  scs <- list()
  i = 1
  while(i<=ndata){
    shape <- sample.int(3,1)
    if (shape == 1){
      res <- gen.rand.doses.plateau(ndose, phi, psi, psi.U, mu1, mu2)
    }else if(shape == 2){
      res <- gen.rand.doses.umbrella(ndose, phi, psi, psi.U, mu1, mu2)
    }else{
      res <- gen.rand.doses.increase(ndose, phi, psi, psi.U, mu1, mu2)
    }
    placebo.dose.val <- 0.1*dose.val[1]
    placebo.q <- 0.1*res$qs[1]
    placebo.p <- 0
    dose.val <- c(placebo.dose.val, dose.val)
    res$qs <- c(placebo.q, res$qs)
    res$ps <- c(placebo.p, res$ps)
    res$MTD <- res$MTD + 1
    res$OBD <- res$OBD + 1
    MND <- MND.level(dose.val, res$qs, gamma, obd = res$k.OBD)
    res$MND <- MND
    scs[[i]] <- res
    i = i + 1
  }
  return(scs)
}

##############Results analysis###############
calculate_ci_z <- function(data_list) {
  mean_values <- sapply(data_list, mean)
  sd_values <- sapply(data_list, sd)
  n_values <- sapply(data_list, length)
  
  z <- qnorm(0.975) 
  ci_list <- lapply(1:length(data_list), function(i) {
    mean_value <- mean_values[i]
    sd_value <- sd_values[i]
    n_value <- n_values[i]
    ci_lower <- mean_value - z * (sd_value / sqrt(n_value))
    ci_upper <- mean_value + z * (sd_value / sqrt(n_value))
    return(c(ci_lower, ci_upper))
  })
  
  return(ci_list)
}

calculate_coverage_probability <- function(data_list, conf_level = 0.95) {
  # 参数检查
  if (!is.list(data_list)) {
    stop("data_list 必须是一个列表")
  }
  
  # 计算每个元素的均值、标准差和样本量
  mean_values <- sapply(data_list, mean)  # 均值作为真实值
  sd_values <- sapply(data_list, sd)      # 标准差
  n_values <- sapply(data_list, length)   # 样本量
  
  # 计算置信区间
  z <- qnorm(1 - (1 - conf_level) / 2)  # 置信水平对应的 z 值
  ci_list <- lapply(1:length(data_list), function(i) {
    mean_value <- mean_values[i]
    sd_value <- sd_values[i]
    n_value <- n_values[i]
    ci_lower <- mean_value - z * (sd_value )
    ci_upper <- mean_value + z * (sd_value )
    return(c(ci_lower, ci_upper))
  })
  
  # 检查每个置信区间是否覆盖真实值（均值），并跳过第一个剂量水平
  coverage_by_dose <- data.frame(
    dose_level = paste0("Dose ", 2:length(data_list)),  # 剂量水平标签
    ci_lower = sapply(ci_list[2:length(data_list)], function(x) x[1]),  # 置信区间下限
    ci_upper = sapply(ci_list[2:length(data_list)], function(x) x[2]),  # 置信区间上限
    coverage = sapply(2:length(data_list), function(i) {  # 是否覆盖真实值
      ci_lower <- ci_list[[i]][1]
      ci_upper <- ci_list[[i]][2]
      true_value <- mean_values[i]  # 真实值是均值
      return(true_value >= ci_lower && true_value <= ci_upper)
    })
  )
  
  # 计算每个剂量水平的覆盖概率
  coverage_by_dose$coverage_probability <- sapply(2:length(data_list), function(i) {
    ci_lower <- ci_list[[i]][1]
    ci_upper <- ci_list[[i]][2]
    true_value <- mean_values[i]  # 真实值是均值
    return(as.numeric(true_value >= ci_lower && true_value <= ci_upper))
  })
  
  # 返回结果
  result <- list(
    
    coverage_by_dose = coverage_by_dose
    
  )
  
  return(result)
}

random_scs_analysis <- function(res, ndata, mu,cohortsize.stage1, cohortsize.stage2, ncohort.stage1, ncohort.stage2){
  ndose <- 5
  n <- length(res)
  N <- cohortsize.stage1*ncohort.stage1+cohortsize.stage2*ncohort.stage2
  MND.sel <- 0
  overdose.allo <- 0
  overdose.sel <- 0
  MND.allo <- 0
  for (i in 1:n){
    result <- res[[i]]
    MND <- result$MND
    tMND <- result$tMND
    MTD <- result$tMTD
    npt <- result$npatients
    if(MND == tMND){
      MND.sel = MND.sel + 1
    }
    
    if (MND > MTD){
      overdose.sel = overdose.sel + 1
    }
    
    if (MTD < ndose){
      overdose.allo = overdose.allo + sum(npt[(MTD+1):ndose])
    }
    MND.allo = MND.allo + npt[tMND]
    
  }
  numpts <- N*ndata
  MND.sel.perc <-  MND.sel/ndata
  overdose.sel.perc <- overdose.sel/ndata
  MND.all.perc <- MND.allo/numpts
  overdose.allo.perc <-  overdose.allo/numpts
  MND.allo.perc <-  MND.allo/numpts
  out <- list(MND.sel.perc = MND.sel.perc, overdose.allo.perc = overdose.allo.perc, MND.allo.perc = MND.allo.perc, overdose.sel.perc = overdose.sel.perc)
  return(out)
}
plot_ate_scatter_with_CI <- function(ate.list, ci) {
  # Flatten the ate.list data and prepare the dose levels
  dose_levels <- unlist(lapply(1:length(ate.list), function(i) rep(paste0("Dose ", i+1), length(ate.list[[i]]))))
  ate_values <- unlist(ate.list)  # Flatten the ATE values
  
  # Extract confidence intervals from calculate_coverage_probability result
  coverage_data <- ci$coverage_by_dose
  ci_lower <- rep(coverage_data$ci_lower, sapply(ate.list, length))  # Repeat ci_lower for each dose level
  ci_upper <- rep(coverage_data$ci_upper, sapply(ate.list, length))  # Repeat ci_upper for each dose level
  
  # Create a data frame
  ate_data <- data.frame(
    dose_level = dose_levels,
    ate = ate_values,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )
  
  # Create the scatter plot with confidence intervals
  p <- ggplot(ate_data, aes(x = dose_level, y = ate)) +
    geom_point(aes(color = dose_level), size = 3) +  # ATE points for each dose level
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +  # Confidence intervals
    labs(
      title = "Scatter Plot of ATE by Dose Level with Confidence Intervals", 
      x = "Dose Level", 
      y = "ATE (Individual Values)",
      color = "Dose Level"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels
  
  # Print the plot
  print(p)
}

cover_prob <- function(sd.data, data, ate_real){
  ndose <- length(data)
  res <- rep(0,ndose)
  for (i in 2:ndose){
    count_num <- 0
    num <- 0
    ate <- na.omit(data[[i]])
    sd.list <- sd.data[[i]]
    #print(ate)
    n <- length(ate)
    ate_true <- ate_real[i]
    for (j in 1:n){
      x <- ate[j]
      sd <- sd.list[j]
      ci_lower <- x - 1.96 * sd 
      ci_upper <- x + 1.96 * sd 
      if (length(x) > 0 && length(ci_lower) > 0 && length(ci_upper) > 0 &&
          !is.na(x) && !is.na(ci_lower) && !is.na(ci_upper)){
        num = num + 1
        #print(c("upper_bpund:",ci_upper))
        if (ate_true >= ci_lower && ate_true <= ci_upper){
          count_num = count_num + 1
          #print(c("count_num:",count_num))
        }
      }
    }
    res[i] = count_num/num
  }
  return(res)
}


fix_scs_analysis <- function(res, phi, dose.val, cohortsize.stage1, cohortsize.stage2, ncohort.stage1, ncohort.stage2, placebo){
  ndose <- length(res[[1]]$p.true)
  n <- length(res)
  N <- cohortsize.stage1*ncohort.stage1+cohortsize.stage2*ncohort.stage2
  MND.sel <- rep(0,ndose)
  overdose.allo <- 0
  MND.allo <- 0
  ATE <- 0
  earlystop <- 0
  p.true = res[[1]]$p.true
  q.true = res[[1]]$q.true
  npts <- rep(0, ndose)
  ntoxs <- rep(0, ndose)
  neffs <- rep(0,ndose)
  nstop = 0
  tMND <- MND.level(dose.val, q.true,0.5,p.true, phi)
  ATE_dose_each <- rep(0, ndose)
  none.sel <- 0
  ate.list <- vector("list", ndose)
  var.pos.ate <- vector("list", ndose)
  for (i in 1:n){
    result <- res[[i]]
    MTD <- MTD.level(phi, p.true)
    MND <- result$MND
    npt <- result$npatients
    neff <- result$neff
    npts <- npts+npt
    ntoxs <- ntoxs+result$ntox
    neffs <- neffs+result$neff 
    if(MND != 99){
      MND.sel[MND] <- MND.sel[MND] + 1
      if (placebo){
        ATE <- ATE + result$ate
        ATE_dose_each[MND] <- ATE_dose_each[MND] + result$ate
        ate.list[[MND]] <- c(ate.list[[MND]], result$ate)
        var.pos.ate[[MND]] <- c(var.pos.ate[[MND]], sqrt(posterior.variance(neff[MND],npt[MND],0.5,0.5)+posterior.variance(neff[1],npt[1],0.5,0.5)))
      }
    }else{
      none.sel = none.sel + 1
    }
    if (MTD < ndose){
      overdose.allo = overdose.allo + sum(npt[(MTD+1):ndose])
    }
    if (earlystop == 0){
      nstop = nstop + 1
    }
  }
  numpts <- N*n
  MND.sel.perc <-  MND.sel/n
  if (placebo){
    #ATE <- ATE/nstop
    real_ate <- q.true - q.true[1]
    bias_ate <- ATE_dose_each/MND.sel - real_ate
  }
  overdose.allo.perc <-  overdose.allo/numpts
  MND.allo.perc <-  MND.allo/numpts
  earlystop.per <- earlystop/n
  if (placebo) {
    ate.list.sd <- lapply(ate.list, sd)
    var.pos.ate.mean <- lapply(var.pos.ate, mean)
    ci <- calculate_coverage_probability(ate.list)
    
    cp <- cover_prob(var.pos.ate, ate.list, real_ate)
    #plot_ate_scatter_with_CI(ate.list[2:length(ate.list)], ci)
    data <- ate.list
    out <- list(tMND = tMND, p.true = p.true, q.true = q.true, npt = npts/n, neff = neffs/n, ntoxs = ntoxs/n, 
                MND.sel.perc = MND.sel.perc, overdose.allo.perc = overdose.allo.perc, MND.allo.perc = MND.allo.perc, 
                ATE = ATE, none.sel = none.sel/n, ATE.dose = ATE_dose_each/MND.sel,bias_ate = bias_ate, standard.error = ate.list.sd, pos.var = var.pos.ate.mean, CI = ci, cp = cp)
    
    res <- list(out = out, data = data)
    return (res)
  }else{
    out <- list(tMND = tMND, p.true = p.true, q.true = q.true, npt = npts/n, neff = neffs/n, ntoxs = ntoxs/n, 
                MND.sel.perc = MND.sel.perc, overdose.allo.perc = overdose.allo.perc, MND.allo.perc = MND.allo.perc, 
                none.sel = none.sel/n)
    
    return(out)
  }
}

plot_ate_scatter_with_CI <- function(ate.list, ci) {
  # Flatten the ate.list data and prepare the dose levels
  dose_levels <- unlist(lapply(1:length(ate.list), function(i) rep(paste0(i), length(ate.list[[i]]))))
  ate_values <- unlist(ate.list)  # Flatten the ATE values
  
  # Extract confidence intervals from calculate_coverage_probability result
  coverage_data <- ci$coverage_by_dose
  ci_lower <- rep(coverage_data$ci_lower, sapply(ate.list, length))  # Repeat ci_lower for each dose level
  ci_upper <- rep(coverage_data$ci_upper, sapply(ate.list, length))  # Repeat ci_upper for each dose level
  
  # Create a data frame
  ate_data <- data.frame(
    dose_level = factor(dose_levels, levels = paste0(1:5)),
    ate = ate_values,
    ci_lower = ci_lower,
    ci_upper = ci_upper
  )
  adjusted_rdpu <- c("#FBB4B9", "#F768A1", "#C51B8A", "#7A0177", "#49006A")
  # Create the scatter plot with confidence intervals
  p <- ggplot(ate_data, aes(x = dose_level, y = ate)) +
    geom_point(aes(color = dose_level), size = 0.8, alpha = 0.3) +  # ATE points for each dose level
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1, size = 0.3, alpha = 0.5, color = "black") +  # Confidence intervals
    labs(
      title = "Scatter Plot of ATE by Dose Level with Confidence Intervals", 
      x = "Dose Level", 
      y = "ATE",
      color = "Dose Level"
    ) +
    scale_color_manual(values = adjusted_rdpu)+
    #scale_color_brewer(palette = adjusted_rdpu) +
    #scale_color_viridis(discrete = TRUE, option = "RdPu") +
    theme_classic() +
    theme(axis.text.x = element_text(hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 12, hjust=0.5,face = "bold"),   
          axis.title.x = element_text(size = 10),  
          axis.title.y = element_text(size = 10),
          legend.position = "bottom",
          aspect.ratio = 0.4,
          axis.line = element_blank(),
          panel.border = element_rect(fill = NA, color = "black", size = 1.5),)  # Rotate x-axis labels
  
  
}

CFOMND_simu <- function(scs, nsim, dose.val, cutoff.eli, early.stop, effthreshold, 
                        cohortsize.stage1, cohortsize.stage2, ncohort.stage1, ncohort.stage2, target,
                        mineff, init.level, alp.prior.eff, bet.prior.eff, alp.prior, bet.prior, gamma, dir_path, with_placebo = FALSE){
  #check the path                  
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  reslist <- vector("list", length(scs))
  #simulation
  for(i in 1:length(scs)){
    file.name <- paste0(dir_path,"/fix_",i,"_5000.RData")
    if(with_placebo){
      p.true <- scs[[i]]$p.true
      q.true <- scs[[i]]$q.true
      placebo.q <- 0.1*scs[[i]]$q.true[1]
      placebo.p <- 0.1*scs[[i]]$p.true[1]
      placebo.dose.val <- 0.1*dose.val[1]
      
      p.true <- c(placebo.p, p.true)
      q.true <- c(placebo.q, q.true)
      dose.val <- c(placebo.dose.val, dose.val)
    }else{
      p.true <- scs[[i]]$p.true
      q.true <- scs[[i]]$q.true
    }
    npts <- cohortsize.stage1*ncohort.stage1+cohortsize.stage2*ncohort.stage2
    plot_list <- list()
    run.fn.MND <- function(k){
      res <- MND.simu(p.true, q.true, dose.val, cohortsize.stage1, cohortsize.stage2,
                      ncohort.stage1, ncohort.stage2, target, gamma, mineff, init.level, with_placebo,
                      prior.para = list(alp.prior = target, bet.prior = 1 - target, alp.prior.eff = 0.5, bet.prior.eff = 0.5))
    }
    ress.MND <-  pbmclapply(1:nsim, run.fn.MND, mc.cores = 4)
    reslist[[i]] <- ress.MND
    save(ress.MND, file = file.name)
  }
  #analyze
  for(i in 1:length(scs)){
    cat(paste0("------------------------------Scenario_",i,"---------------------------------"))
    cat("\n")
    res <- reslist[[i]]
    if (with_placebo){
      ana <- fix_scs_analysis(res, 0.3, dose.val, cohortsize.stage1, cohortsize.stage2, ncohort.stage1, ncohort.stage2, TRUE)
      print(ana$out)
      scenario_name <- paste("Scenario", i)  
      p <- plot_ate_scatter_with_CI(ana$data[2:length(ana$data)], ana$out$CI)  
      p <- p + ggtitle(scenario_name) 
      plot_list[[i]] <- p  
    }else{
      ana <- fix_scs_analysis(res, 0.3, dose.val, cohortsize.stage1, cohortsize.stage2, ncohort.stage1, ncohort.stage2, FALSE)
      print(ana)
    }
    
  }
  if(with_placebo){
    final_plot <- ggarrange(
      plotlist = plot_list,  
      ncol = 2, nrow = ceiling(length(plot_list) / 2),    
      common.legend = TRUE,  
      legend = "bottom",     
      align = "hv"           
    )
    print(final_plot)
  }
}


CFOMND_simu <- function(scs, nsim, dose.val, cutoff.eli, early.stop, effthreshold, 
                        cohortsize.stage1, cohortsize.stage2, ncohort.stage1, ncohort.stage2, target,
                        mineff, init.level, alp.prior.eff, bet.prior.eff, alp.prior, bet.prior, gamma, dir_path, with_placebo = FALSE){
  #check the path                  
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  reslist <- vector("list", length(scs))
  #simulation
  for(i in 1:length(scs)){
    file.name <- paste0(dir_path,"/fix_",i,"_5000.RData")
    if(with_placebo){
      p.true <- scs[[i]]$p.true
      q.true <- scs[[i]]$q.true
      placebo.q <- 0.1*scs[[i]]$q.true[1]
      placebo.p <- 0.1*scs[[i]]$p.true[1]
      placebo.dose.val <- 0.1*dose.val[1]
      
      p.true <- c(placebo.p, p.true)
      q.true <- c(placebo.q, q.true)
      dose.val <- c(placebo.dose.val, dose.val)
    }else{
      p.true <- scs[[i]]$p.true
      q.true <- scs[[i]]$q.true
    }
    npts <- cohortsize.stage1*ncohort.stage1+cohortsize.stage2*ncohort.stage2
    plot_list <- list()
    run.fn.MND <- function(k){
      res <- MND.simu(p.true, q.true, dose.val, cohortsize.stage1, cohortsize.stage2,
                      ncohort.stage1, ncohort.stage2, target, gamma, mineff, init.level, with_placebo,
                      prior.para = list(alp.prior = target, bet.prior = 1 - target, alp.prior.eff = 0.5, bet.prior.eff = 0.5))
    }
    ress.MND <-  pbmclapply(1:nsim, run.fn.MND, mc.cores = 4)
    reslist[[i]] <- ress.MND
    save(ress.MND, file = file.name)
  }
  #analyze
  for(i in 1:length(scs)){
    cat(paste0("------------------------------Scenario_",i,"---------------------------------"))
    cat("\n")
    res <- reslist[[i]]
    if (with_placebo){
      ana <- fix_scs_analysis(res, 0.3, dose.val, cohortsize.stage1, cohortsize.stage2, ncohort.stage1, ncohort.stage2, TRUE)
      print(ana$out)
      scenario_name <- paste("Scenario", i)  
      p <- plot_ate_scatter_with_CI(ana$data[2:length(ana$data)], ana$out$CI)  
      p <- p + ggtitle(scenario_name) 
      plot_list[[i]] <- p  
    }else{
      ana <- fix_scs_analysis(res, 0.3, dose.val, cohortsize.stage1, cohortsize.stage2, ncohort.stage1, ncohort.stage2, FALSE)
      print(ana)
    }
    
  }
  if(with_placebo){
    final_plot <- ggarrange(
      plotlist = plot_list,  
      ncol = 2, nrow = ceiling(length(plot_list) / 2),    
      common.legend = TRUE,  
      legend = "bottom",     
      align = "hv"           
    )
    print(final_plot)
  }
}

CFOMND_simu_rand <- function(mu, ndata, ndose, dose.val, cutoff.eli, early.stop, effthreshold, 
                             cohortsize.stage1, cohortsize.stage2, ncohort.stage1, ncohort.stage2, target,
                             mineff, init.level, alp.prior.eff, bet.prior.eff, alp.prior, bet.prior, gamma, dir_path){
  nsim = 1
  if (!dir.exists(paste0(dir_path,"/scs/"))) {
    dir.create(paste0(dir_path,"/scs/"), recursive = TRUE)
  }
  
  if (!dir.exists(paste0(dir_path,"res/"))) {
    dir.create(paste0(dir_path,"res/"), recursive = TRUE)
  }
  phi <- 0.3
  psi <- 0.1
  psi.U <- 0.8
  npts <- cohortsize.stage1*ncohort.stage1+cohortsize.stage2*ncohort.stage2
  mu1 = mu2 = mu
  #generate random scenario
  scs <- scs.gen(ndose, dose.val, ndata, gamma, phi, psi, psi.U, mu1, mu2)
  save.file.name <- paste0(dir_path, "scs/random_scs",ndata,"gamma_",gamma,
                           "mu_",mu1,"target_",phi,"eff[",psi,",",psi.U,"]",".RData")
  save(scs, file = save.file.name)
  #simulation
  run.fn.MND.random <- function(k){
    scs.test <- scs[[k]]
    p.true <- scs.test$ps
    q.true <- scs.test$qs
    res <- MND.simu(p.true, q.true, dose.val, cohortsize.stage1, cohortsize.stage2,
                    ncohort.stage1, ncohort.stage2, target, gamma, mineff, init.level,
                    prior.para = list(alp.prior = target, bet.prior = 1 - target, alp.prior.eff = 0.5, bet.prior.eff = 0.5))
    res$tMND <- scs.test$MND
    res$tMTD <- scs.test$k.MTD
    return(res)
  }
  
  
  #Random scenario
  save.file.name <- paste0(dir_path, "res/",ndata,"gamma_",gamma,"npts_",npts,
                           "mu_",mu1,"target_",phi,"eff[",psi,",",psi.U,"]",".RData")
  ress.MND <- pbmclapply(1:ndata, run.fn.MND.random, mc.cores=4)
  save(ress.MND, file=save.file.name)
  ana <- random_scs_analysis(ress.MND, ndata, mu, cohortsize.stage1, cohortsize.stage2, ncohort.stage1, ncohort.stage2)
  print(ana)
}

MND.next <- function(ays, axs, ans, currdose, dose.val, cohortsize.stage1, cohortsize.stage2,
                     ncohort.stage1, ncohort.stage2, target, gamma, mineff,
                     placebo = FALSE, prior.para = list(alp.prior = target, bet.prior = 1 - target,
                                                        alp.prior.eff = 0.5, bet.prior.eff = 0.5),
                     cutoff.eli = 0.95, early.stop = 0.95, effearly.stop = 0.9,
                     admset = NULL, MTD = NULL){
  ndose <- length(ans)
  if (length(ays) != ndose || length(axs) != ndose || length(dose.val) != ndose){
    stop("ays/axs/ans/dose.val length mismatch")
  }
  if (currdose < 1 || currdose > ndose){
    stop("currdose out of range")
  }
  if (is.null(prior.para$alp.prior)){
    prior.para <- c(prior.para, list(alp.prior = target, bet.prior = 1 - target))
  }
  if (is.null(prior.para$alp.prior.eff)){
    prior.para <- c(prior.para, list(alp.prior.eff = 0.5, bet.prior.eff = 0.5))
  }
  alp.prior.eff <- prior.para$alp.prior.eff
  bet.prior.eff <- prior.para$bet.prior.eff
  
  total.stage1.n <- ncohort.stage1 * cohortsize.stage1
  stage <- if (sum(ans) < total.stage1.n) "stage1" else "stage2"
  
  neigh.vec <- function(vec, idx){
    left <- if (idx > 1) vec[idx - 1] else NA
    right <- if (idx < length(vec)) vec[idx + 1] else NA
    c(left, vec[idx], right)
  }
  
  if (stage == "stage1"){
    overtox.curr <- overdose.fn(target, cutoff.eli, ays[currdose], ans[currdose], prior.para)
    overtox.early <- FALSE
    if (currdose == 1 && cutoff.eli != early.stop){
      overtox.early <- overdose.fn(target, early.stop, ays[1], ans[1], prior.para)
    }
    
    if (placebo){
      earlystop <- (currdose == 1 && (overtox.curr || overtox.early)) || (currdose >= 2 && overtox.curr)
    }else{
      earlystop <- (currdose == 1 && (overtox.curr || overtox.early))
    }
    
    if (earlystop){
      return(list(nextdose = 99, stage = stage, earlystop = 1, MTD = 99))
    }
    
    cys <- neigh.vec(ays, currdose)
    cns <- neigh.vec(ans, currdose)
    cfores <- CFO.next(target, cys, cns, currdose, prior.para, cutoff.eli = cutoff.eli, early.stop = early.stop)
    MTD <- CFO.selectmtd(target, ans, ays, prior.para, cutoff.eli, early.stop, verbose = FALSE)$MTD
    return(list(nextdose = cfores$nextdose, stage = stage, earlystop = 0, MTD = MTD,
                decision = cfores$decision, odds.ratio = cfores$odds.ratio, toxprob = cfores$toxprob))
  }
  
  if (is.null(MTD)){
    MTD <- CFO.selectmtd(target, ans, ays, prior.para, cutoff.eli, early.stop, verbose = FALSE)$MTD
  }
  if (is.null(admset)){
    admset <- 1:MTD
  }else{
    admset <- admset[admset <= MTD]
  }
  if (length(admset) == 0 || MTD == 99){
    return(list(nextdose = 99, stage = stage, earlystop = 1, MTD = MTD, admset = admset))
  }
  
  if (sum(ans) == total.stage1.n){
    effprobs <- moveprobs(axs[1:MTD], ans[1:MTD], alp.prior.eff, bet.prior.eff)
    OBD.hat <- which.max(effprobs)
    if (under.eff.fn(mineff, axs[OBD.hat], ans[OBD.hat], effearly.stop, prior.para)){
      return(list(nextdose = 99, stage = stage, earlystop = 1, MTD = MTD, admset = admset,
                  effprobs = effprobs, OBD.hat = OBD.hat))
    }
    return(list(nextdose = OBD.hat, stage = stage, earlystop = 0, MTD = MTD, admset = admset,
                effprobs = effprobs, OBD.hat = OBD.hat))
  }
  
  overtox.curr <- overdose.fn(target, cutoff.eli, ays[currdose], ans[currdose], prior.para)
  if (placebo){
    if (currdose == 2 && overtox.curr){
      return(list(nextdose = 99, stage = stage, earlystop = 1, MTD = MTD, admset = admset))
    }
  }else{
    if (currdose == 1 && overtox.curr){
      return(list(nextdose = 99, stage = stage, earlystop = 1, MTD = MTD, admset = admset))
    }
  }
  if (overtox.curr){
    admset <- admset[admset < currdose]
    if (length(admset) == 0){
      return(list(nextdose = 99, stage = stage, earlystop = 1, MTD = MTD, admset = admset))
    }
  }
  
  q.hat <- efficacy.pos.mean(axs[admset], ans[admset], alp.prior.eff, bet.prior.eff)
  q.tilde <- q.hat/(dose.val[admset])^gamma
  prob.sample <- q.tilde/sum(q.tilde)
  nextdose <- sample(admset, size = 1, prob = prob.sample)
  
  list(nextdose = nextdose, stage = stage, earlystop = 0, MTD = MTD, admset = admset,
       q.hat = q.hat, q.tilde = q.tilde, prob.sample = prob.sample)
}


MND.select <- function(ays, axs, ans, dose.val, prior.para, gamma,
                       target, cutoff.eli = 0.95, early.stop = 0.95,
                       mineff = 0.15, placebo = FALSE){
  MTD <- CFO.selectmtd(target, ans, ays, prior.para, cutoff.eli, early.stop, verbose = FALSE)$MTD
  admset <- 1:MTD
  if (length(admset) == 0 || isTRUE(MTD == 99)) {
    return(list(OBD = 99, canset = integer(0), canaxs = numeric(0), canans = numeric(0),
                dose.val.adm = numeric(0), q.hat = NA, probs = NULL, selectMND = NULL,
                earlystop = 1, MTD = MTD, admset = admset))
  }
  alp.prior.eff <- prior.para$alp.prior.eff
  bet.prior.eff <- prior.para$bet.prior.eff
  probs <- moveprobs(axs[admset], ans[admset], alp.prior.eff, bet.prior.eff)
  OBD <- admset[which.max(probs)]
  canset <- admset[admset <= OBD]
  canaxs <- axs[canset]
  canans <- ans[canset]
  dose.val.adm <- dose.val[canset]
  prior.para.sel <- c(prior.para, list(alp.prior = alp.prior.eff, bet.prior = bet.prior.eff))
  res <- selectMND(canset, canans, canaxs, prior.para.sel)
  q.hat <- res$phat
  MND <- MND.level(dose.val.adm, q.hat, gamma, mineff = mineff, obd = OBD, placebo = placebo, estimate = TRUE)
  ATE <- NA
  if (placebo && !is.na(MND) && MND != 99) {
    MND.mean <- efficacy.pos.mean(axs[MND], ans[MND], alp.prior.eff, bet.prior.eff)
    placebo.mean <- efficacy.pos.mean(axs[1], ans[1], alp.prior.eff, bet.prior.eff)
    ATE <- MND.mean - placebo.mean
  }
  list(OBD = OBD, MND = MND, ATE = ATE, canset = canset, canaxs = canaxs, canans = canans,
       dose.val.adm = dose.val.adm, q.hat = q.hat, probs = probs, selectMND = res)
}
