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

MND.level <- function(dose.val, q.true, gamma, p.true = NULL, target = NULL, mineff = NULL, obd = NULL, PED = FALSE, estimate = FALSE){
  OBD <- obd
  ndose <- length(p.true)

  if(is.null(OBD)){
    OBD <- OBD.level(target, mineff, p.true, q.true)
    if (PED & OBD == 1){
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
                     ncohort.stage1, ncohort.stage2, target, gamma, mineff, init.level, PED = FALSE,
                     cutoff.eli = 0.95, early.stop = 0.95,
                     alp.prior.eff = 0.5, bet.prior.eff = 0.5,
                     prior.para = list(alp.prior = target, bet.prior = 1 - target, alp.prior.eff = alp.prior.eff, bet.prior.eff = bet.prior.eff)){

  count <- 0
  earlystop.count <- 0
  if (is.null(prior.para$alp.prior.eff)) prior.para$alp.prior.eff <- alp.prior.eff
  if (is.null(prior.para$bet.prior.eff)) prior.para$bet.prior.eff <- bet.prior.eff
  alp.prior.eff <- prior.para$alp.prior.eff
  bet.prior.eff <- prior.para$bet.prior.eff
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
  for (i in 1:ncohort.stage1){
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

    if (PED){
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
    if (under.eff.fn(mineff = 0.15, axs[OBD.hat], ans[OBD.hat], 0.9, prior.para)){
      earlystop <- 1
    }else{
      #Stage2
      for (i in 1:ncohort.stage2){
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

        if(PED){
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
  }
  ATE <- NA
  if (earlystop == 0) {
    MND <- MND.level(dose.val.adm, q.hat, gamma, mineff = 0.15, obd = OBD, PED = PED, estimate = TRUE)
    if(PED){
      MND.mean <- efficacy.pos.mean(axs[MND],ans[MND], alp.prior.eff, bet.prior.eff)
      PED.mean <- efficacy.pos.mean(axs[1],ans[1], alp.prior.eff, bet.prior.eff)
      ATE <- MND.mean - PED.mean
    }
  }
  else{
    MND <- 99
  }
  if (PED){
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

scs.gen.PED <- function(ndose, dose.val, ndata, gamma, phi, psi, psi.U, mu1, mu2){
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
    PED.dose.val <- 0.1*dose.val[1]
    PED.q <- 0.1*res$qs[1]
    PED.p <- 0
    dose.val <- c(PED.dose.val, dose.val)
    res$qs <- c(PED.q, res$qs)
    res$ps <- c(PED.p, res$ps)
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
 if (!is.list(data_list)) {
    stop("data_list must be a list")
  }


  mean_values <- sapply(data_list, mean)
  sd_values <- sapply(data_list, sd)
  n_values <- sapply(data_list, length)
  z <- qnorm(1 - (1 - conf_level) / 2)
  ci_list <- lapply(1:length(data_list), function(i) {
    mean_value <- mean_values[i]
    sd_value <- sd_values[i]
    n_value <- n_values[i]
    ci_lower <- mean_value - z * (sd_value )
    ci_upper <- mean_value + z * (sd_value )
    return(c(ci_lower, ci_upper))
  })

  coverage_by_dose <- data.frame(
    dose_level = paste0("Dose ", 2:length(data_list)),
    ci_lower = sapply(ci_list[2:length(data_list)], function(x) x[1]),
    ci_upper = sapply(ci_list[2:length(data_list)], function(x) x[2]),
    coverage = sapply(2:length(data_list), function(i) {
    ci_lower <- ci_list[[i]][1]
      ci_upper <- ci_list[[i]][2]
      true_value <- mean_values[i]
      return(true_value >= ci_lower && true_value <= ci_upper)
    })
  )

  coverage_by_dose$coverage_probability <- sapply(2:length(data_list), function(i) {
    ci_lower <- ci_list[[i]][1]
    ci_upper <- ci_list[[i]][2]
    true_value <- mean_values[i]
    return(as.numeric(true_value >= ci_lower && true_value <= ci_upper))
  })

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

fix_scs_analysis <- function(res, phi, dose.val, cohortsize.stage1, cohortsize.stage2, ncohort.stage1, ncohort.stage2, PED){
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
      if (PED){
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
  if (PED){
    #ATE <- ATE/nstop
    real_ate <- q.true - q.true[1]
    bias_ate <- ATE_dose_each/MND.sel - real_ate
  }
  overdose.allo.perc <-  overdose.allo/numpts
  MND.allo.perc <-  MND.allo/numpts
  earlystop.per <- earlystop/n
  if (PED) {
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
  p <- ggplot(ate_data, aes(x = .data$dose_level, y = .data$ate)) +
    geom_point(aes(color = .data$dose_level), size = 0.8, alpha = 0.3) +  # ATE points for each dose level
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

selectMND <- function(canset, npts, neff, prior.para = list(alp.prior = 0.5, bet.prior = 0.5)) {
  pava <- function(x, wt = rep(1, length(x))) {
    n <- length(x)
    if (n <= 1) return(x)
    if (any(is.na(x)) || any(is.na(wt))) {
      stop("Missing values in 'x' or 'wt' not allowed")
    }
    lvlsets <- (1:n)
    repeat {
      viol <- (as.vector(diff(x)) < 0)
      if (!(any(viol))) break
      i <- min((1:(n - 1))[viol])
      lvl1 <- lvlsets[i]
      lvl2 <- lvlsets[i + 1]
      ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
      x[ilvl] <- sum(x[ilvl] * wt[ilvl]) / sum(wt[ilvl])
      lvlsets[ilvl] <- lvl1
    }
    x
  }

  if (is.null(prior.para$alp.prior) || is.null(prior.para$bet.prior)) {
    prior.para <- c(prior.para, list(alp.prior = 0.5, bet.prior = 0.5))
  }
  alp.prior <- prior.para$alp.prior
  bet.prior <- prior.para$bet.prior

  y.adm <- neff
  n.adm <- npts
  phat <- (y.adm + alp.prior) / (n.adm + alp.prior + bet.prior)
  phat.var <- (y.adm + alp.prior) * (n.adm - y.adm + bet.prior) /
    ((n.adm + alp.prior + bet.prior)^2 * (n.adm + alp.prior + bet.prior + 1))
  phat <- pava(phat, wt = 1 / phat.var)
  phat <- phat + (1:length(phat)) * 1e-10
  list(phat = phat)
}



CFO.next <- function(target, cys, cns, currdose, prior.para=list(alp.prior=target, bet.prior=1-target),
                     cutoff.eli=0.95, early.stop=0.95){
  ###############################################################################
  ###############define the functions used for main function#####################
  ###############################################################################
  # posterior probability of pj >= phi given data
  post.prob.fn <- function(phi, y, n, alp.prior=0.1, bet.prior=0.1){
    alp <- alp.prior + y
    bet <- bet.prior + n - y
    1 - pbeta(phi, alp, bet)
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

  prob.int <- function(phi, y1, n1, y2, n2, alp.prior, bet.prior){
    alp1 <- alp.prior + y1
    alp2 <- alp.prior + y2
    bet1 <- bet.prior + n1 - y1
    bet2 <- bet.prior + n2 - y2
    fn.min <- function(x){
      dbeta(x, alp1, bet1)*(1-pbeta(x, alp2, bet2))
    }
    fn.max <- function(x){
      pbeta(x, alp1, bet1)*dbeta(x, alp2, bet2)
    }
    const.min <- integrate(fn.min, lower=0, upper=0.99, subdivisions=1000, rel.tol = 1e-10)$value
    const.max <- integrate(fn.max, lower=0, upper=1, rel.tol = 1e-10)$value
    p1 <- integrate(fn.min, lower=0, upper=phi)$value/const.min
    p2 <- integrate(fn.max, lower=0, upper=phi)$value/const.max

    list(p1=p1, p2=p2)
  }

  OR.values <- function(phi, y1, n1, y2, n2, alp.prior, bet.prior, type){
    ps <- prob.int(phi, y1, n1, y2, n2, alp.prior, bet.prior)
    if (type=="L"){
      pC <- 1 - ps$p2
      pL <- 1 - ps$p1
      oddsC <- pC/(1-pC)
      oddsL <- pL/(1-pL)
      OR <- oddsC*oddsL

    }else if (type=="R"){
      pC <- 1 - ps$p1
      pR <- 1 - ps$p2
      oddsC <- pC/(1-pC)
      oddsR <- pR/(1-pR)
      OR <- (1/oddsC)/oddsR
    }
    return(OR)
  }

  All.OR.table <- function(phi, n1, n2, type, alp.prior, bet.prior){
    ret.mat <- matrix(rep(0, (n1+1)*(n2+1)), nrow=n1+1)
    for (y1cur in 0:n1){
      for (y2cur in 0:n2){
        ret.mat[y1cur+1, y2cur+1] <- OR.values(phi, y1cur, n1, y2cur, n2, alp.prior, bet.prior, type)
      }
    }
    ret.mat
  }

  # compute the marginal prob when lower < phiL/phiC/phiR < upper
  # i.e., Pr(Y=y|lower<phi<upper)
  margin.phi <- function(y, n, lower, upper){
    C <- 1/(upper-lower)
    fn <- function(phi) {
      dbinom(y, n, phi)*C
    }
    integrate(fn, lower=lower, upper=upper)$value
  }

  # Obtain the table of marginal distribution of (y1, y2)
  # after intergrate out (phi1, phi2)
  # under H0 and H1
  # H0: phi1=phi, phi < phi2 < 2phi
  # H1: phi2=phi, 0   < phi1 < phi
  margin.ys.table <- function(n1, n2, phi, hyperthesis){
    if (hyperthesis=="H0"){
      p.y1s <- dbinom(0:n1, n1, phi)
      p.y2s <- sapply(0:n2, margin.phi, n=n2, lower=phi, upper=2*phi)
    }else if (hyperthesis=="H1"){
      p.y1s <- sapply(0:n1, margin.phi, n=n1, lower=0, upper=phi)
      p.y2s <- dbinom(0:n2, n2, phi)
    }
    p.y1s.mat <- matrix(rep(p.y1s, n2+1), nrow=n1+1)
    p.y2s.mat <- matrix(rep(p.y2s, n1+1), nrow=n1+1, byrow=TRUE)
    margin.ys <- p.y1s.mat * p.y2s.mat
    margin.ys
  }

  # Obtain the optimal gamma for the hypothesis test
  optim.gamma.fn <- function(n1, n2, phi, type, alp.prior, bet.prior){
    OR.table <- All.OR.table(phi, n1, n2, type, alp.prior, bet.prior)
    ys.table.H0 <- margin.ys.table(n1, n2, phi, "H0")
    ys.table.H1 <- margin.ys.table(n1, n2, phi, "H1")

    argidx <- order(OR.table)
    sort.OR.table <- OR.table[argidx]
    sort.ys.table.H0 <- ys.table.H0[argidx]
    sort.ys.table.H1 <- ys.table.H1[argidx]
    n.tol <- length(sort.OR.table)

    if (type=="L"){
      errs <- rep(0, n.tol-1)
      for (i in 1:(n.tol-1)){
        err1 <- sum(sort.ys.table.H0[1:i])
        err2 <- sum(sort.ys.table.H1[(i+1):n.tol])
        err <- err1 + err2
        errs[i] <- err
      }
      min.err <- min(errs)
      if (min.err > 1){
        gam <- 0
        min.err <- 1
      }else {
        minidx <- which.min(errs)
        gam <- sort.OR.table[minidx]
      }
    }else if (type=='R'){
      errs <- rep(0, n.tol-1)
      for (i in 1:(n.tol-1)){
        err1 <- sum(sort.ys.table.H1[1:i])
        err2 <- sum(sort.ys.table.H0[(i+1):n.tol])
        err <- err1 + err2
        errs[i] <- err
      }
      min.err <- min(errs)
      if (min.err > 1){
        gam <- 0
        min.err <- 1
      }else {
        minidx <- which.min(errs)
        gam <- sort.OR.table[minidx]
      }

    }
    list(gamma=gam, min.err=min.err)
  }

  ###############################################################################
  ############################MAIN DUNCTION######################################
  ###############################################################################
  if (is.null(prior.para$alp.prior)){
    prior.para <- c(prior.para, list(alp.prior=target, bet.prior=1-target))
  }
  alp.prior <- prior.para$alp.prior
  bet.prior <- prior.para$bet.prior

  cover.doses <- c(0,0,0)

  for (i in 1:3){
    cy <- cys[i]
    cn <- cns[i]
    if (is.na(cn)){
      cover.doses[i] <- NA
    }else{
      prior.para <- c(list(y=cy, n=cn),list(alp.prior=alp.prior, bet.prior=bet.prior))
      if (overdose.fn(target, cutoff.eli, prior.para)){
        cover.doses[i:3] <- 1
        break()
      }
    }
  }

  cover.prob <- c(0,0,0)
  for (i in 1:3){
    cy <- cys[i]
    cn <- cns[i]
    if (is.na(cn)){
      cover.prob[i] <- NA
    }else{
      cover.prob[i] <- post.prob.fn(target, cy, cn, alp.prior, bet.prior)
    }
  }

  if (cutoff.eli != early.stop) {
    cy <- cys[1]
    cn <- cns[1]
    if (is.na(cn)){
      cover.doses[1] <- NA
    }else{
      prior.para <- c(list(y=cy, n=cn),list(alp.prior=alp.prior, bet.prior=bet.prior))
      if (overdose.fn(target, early.stop, prior.para)){
        cover.doses[1:3] <- 1
      }
    }
  }

  cover.doses <- ifelse(is.na(cys), NA, cover.doses)
  gam1 <- NA
  gam2 <- NA
  OR.v1 <- NA
  OR.v2 <- NA
  position <- which(cover.doses == 1)[1]
  overtox <- c(-1, 0, 1)[position] + currdose
  prior.para <- c(list(alp.prior=alp.prior, bet.prior=bet.prior))
  if ((cover.doses[2] == 1)&(currdose == 1)){
    index <- NA
    decision <- "stop"
  } else {
    if (cover.doses[2] == 1){
      index <- -1
      decision <- "de-escalation"
    }
    else{
      if (is.na(cys[1]) & (cover.doses[3]==1)){
        index <- 0
        decision <- "stay"
      }
      else  if (is.na(cys[1]) & (!(cover.doses[3]==1))){
        gam2 <- optim.gamma.fn(cns[2], cns[3], target, "R", alp.prior, bet.prior)$gamma
        OR.v2 <- OR.values(target, cys[2], cns[2], cys[3], cns[3], alp.prior, bet.prior, type="R")
        if (OR.v2>gam2){
          index <- 1
          decision <- "escalation"
        }else{
          index <- 0
          decision <- "stay"
        }
      }
      else  if (is.na(cys[3]) | (cover.doses[3]==1)){
        gam1 <- optim.gamma.fn(cns[1], cns[2], target, "L", alp.prior, bet.prior)$gamma
        OR.v1 <- OR.values(target, cys[1], cns[1], cys[2], cns[2], alp.prior, bet.prior, type="L")
        if (OR.v1>gam1){
          index <- -1
          decision <- "de-escalation"
        }else{
          index <- 0
          decision <- "stay"
        }
      }
      else  if (!(is.na(cys[1]) | is.na(cys[3]) | cover.doses[3]==1)){
        gam1 <- optim.gamma.fn(cns[1], cns[2], target, "L", alp.prior, bet.prior)$gamma
        gam2 <- optim.gamma.fn(cns[2], cns[3], target, "R", alp.prior, bet.prior)$gamma
        OR.v1 <- OR.values(target, cys[1], cns[1], cys[2], cns[2], alp.prior, bet.prior, type="L")
        OR.v2 <- OR.values(target, cys[2], cns[2], cys[3], cns[3], alp.prior, bet.prior, type="R")
        v1 <- OR.v1 > gam1
        v2 <- OR.v2 > gam2
        if (v1 & !v2){
          index <- -1
          decision <- "de-escalation"
        }else if (!v1 & v2){
          index <- 1
          decision <- "escalation"
        }else{
          index <- 0
          decision <- "stay"
        }
      }
    }
  }

  if (decision=='stop'){
    nextdose <- 99
  }else{
    nextdose <- currdose+index
  }
  odds.ratio <- list(gamL = gam1, gamR = gam2, OR.L = OR.v1, OR.R = OR.v2)
  out <- list(target=target, cys=cys, cns=cns, decision=decision, currdose = currdose,
              nextdose=nextdose, overtox=overtox, toxprob=cover.prob, odds.ratio = odds.ratio)
  class(out) <- c("cfo_decision", "cfo")
  return(out)
}



CFO.selectmtd <- function (target, npts, ntox, prior.para=list(alp.prior=target, bet.prior=1-target),
                           cutoff.eli = 0.95, early.stop = 0.95, verbose = TRUE)
{
  pava <- function(x, wt = rep(1, length(x))) {
    n <- length(x)
    if (n <= 1)
      return(x)
    if (any(is.na(x)) || any(is.na(wt))) {
      stop("Missing values in 'x' or 'wt' not allowed")
    }
    lvlsets <- (1:n)
    repeat {
      viol <- (as.vector(diff(x)) < 0)
      if (!(any(viol)))
        break
      i <- min((1:(n - 1))[viol])
      lvl1 <- lvlsets[i]
      lvl2 <- lvlsets[i + 1]
      ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
      x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
      lvlsets[ilvl] <- lvl1
    }
    x
  }

  phat <- NA

  if (is.null(prior.para$alp.prior)){
    prior.para <- c(prior.para, list(alp.prior=target, bet.prior=1-target))
  }
  alp.prior <- prior.para$alp.prior
  bet.prior <- prior.para$bet.prior

  y = ntox
  n = npts
  ndose = length(n)
  elimi = rep(0, ndose)
  for (i in 1:ndose) {
    if (n[i] >= 3) {
      if (1 - pbeta(target, y[i] + alp.prior, n[i] - y[i] + bet.prior) >
          cutoff.eli) {
        elimi[i:ndose] = 1
        break
      }
    }
  }
  if (cutoff.eli != early.stop) {
    if (n[1] >= 3) {
      if (1 - pbeta(target, y[1] + alp.prior, n[1] - y[1] + bet.prior) >
          early.stop) {
        elimi[1:ndose] = 1
      }
    }
  }
  if (elimi[1] == 1 || sum(n[elimi == 0]) == 0){
    selectdose = 99
  }else {
    adm.set = (n != 0) & (elimi == 0)
    adm.index = which(adm.set == T)
    y.adm = y[adm.set]
    n.adm = n[adm.set]
    phat = (y.adm + alp.prior)/(n.adm + alp.prior + bet.prior)
    phat.var = (y.adm + alp.prior) * (n.adm - y.adm + bet.prior)/((n.adm +
                                                                     alp.prior + bet.prior)^2 * (n.adm + alp.prior + bet.prior + 1))
    phat = pava(phat, wt = 1/phat.var)
    phat = phat + (1:length(phat)) * 1e-10
    selectd = sort(abs(phat - target), index.return = T)$ix[1]
    selectdose = adm.index[selectd]
  }


  if (selectdose == 99) {
    print(npts)
    print(ntox)
    message("All tested doses are overly toxic. No MTD is selected! \n")
  }
  out = list(target = target, MTD = selectdose, p = phat)
  class(out)<-c("cfo_sel","cfo")
  return(out)
}



