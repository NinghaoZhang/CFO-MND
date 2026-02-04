
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

