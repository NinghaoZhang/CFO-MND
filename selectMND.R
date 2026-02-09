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
