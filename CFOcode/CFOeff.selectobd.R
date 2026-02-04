#' Select the optimal biological dose (OBD) for the real single-drug trials
#'
#' Select the optimal biological dose (OBD) when the real single-drug trials is completed
#' @usage CFOeff.selectobd(target, txs, tys, tns, prior.para, mineff, effearly.stop)
#' 
#' 
#' @param target the target DLT rate.
#' @param txs the cumulative counts of efficacy outcomes at all dose levels.
#' @param tys the cumulative counts of DLTs observed at all dose levels.
#' @param tns the cumulative counts of patients treated at all dose levels.
#' @param prior.para the prior parameters for two beta distributions, where set as \code{list(alp.prior.eff = 0.5, bet.prior.eff = 0.5)} by default. 
#'                   \code{alp.eff.prior} and \code{bet.eff.prior}represent the parameters of the Jeffreys' prior distribution for the efficacy probability 
#'                   at any dose level.This prior distribution is specified as Beta(\code{alpha.eff.prior}, \code{beta.eff.prior}). 
#' @param mineff the lowest acceptable efficacy rate.
#' @param effearly.stop the threshold value for early stopping due to low efficacy. The trial would be terminated
#'                      early if \eqn{Pr(q_k<\psi |y_k,m_k \ge 3)} is smaller than the value of \code{effearly.stop} where \eqn{q_k, y_k} and \eqn{m_k}
#'                      are the efficacy probability, the number of efficacy outcomes and the number of patients at dose level \eqn{k}. 
#'                      \eqn{\psi} is the the lowest acceptable efficacy rate which is set by \code{mineff} here. 
#'                      By default, \code{effearly.stop} is set as \code{0.9}. 
#'
#' @return The \code{CFOeff.selectobd()} function returns a list object comprising the following elements:
#' \itemize{
#'   \item OBD: the selected OBD. \code{OBD = 99} indicates that all tested doses are overly toxic or having low efficacy.
#'   \item MTD: MTD here is get by using function \code{CFO.selectmtd}. MTD is used as the upper bound of the admissible set.
#'   \item OBD.probs: the probability that each dose level would be selected as OBD. The probability indicates that \eqn{q_k} corresponds 
#'   to dose level \eqn{k} being the highest in the admissible set. \eqn{q_k} is efficacy probability correspond to dose level k here.
#' } 
#' @export
#' @author Jialu Fang, Ninghao Zhang, Wenliang Wang, and Guosheng Yin
#' 
#' @references Jin H, Yin G (2022). CFO: Calibration-free odds design for phase I/II clinical trials.
#'             \emph{Statistical Methods in Medical Research}, 31(6), 1051-1066. \cr
#' @examples
#' target <- 0.3; mineff<- 0.3
#' txs <- c(3, 1, 7, 11, 26); tys <- c(0, 0, 0, 0, 6); tns <- c(6, 3, 12, 17, 36)
#' prior.para = list(alp.prior.eff = 0.5, bet.prior.eff = 0.5)
#' effearly.stop <- 0.95
#' result <- CFOeff.selectobd(target, txs, tys, tns, prior.para, mineff, effearly.stop)
#' summary(result)
#' \donttest{
#' ##Low efficacy
#' target <- 0.3; mineff<- 0.3
#' txs = c(0, 0, 0, 0, 0); tys = c(2, 1, 1, 1, 6); tns = c(36, 23, 22, 27, 36)
#' prior.para = list(alp.prior.eff = 0.5, bet.prior.eff = 0.5)
#' effearly.stop <- 0.95
#' result <- CFOeff.selectobd(target, txs, tys, tns, prior.para, mineff, effearly.stop)
#' summary(result)
#' }
#' \donttest{
#' ##High toxicity
#' target <- 0.3; mineff<- 0.3
#' txs = c(3, 1, 7, 11, 26); tys = c(36, 23, 22, 27, 36); tns = c(36, 23, 22, 27, 36)
#' prior.para = list(alp.prior.eff = 0.5, bet.prior.eff = 0.5)
#' effearly.stop <- 0.95
#' result <- CFOeff.selectobd(target, txs, tys, tns, prior.para, mineff, effearly.stop)
#' summary(result)
#' }
#' 
CFOeff.selectobd <- function(target, txs, tys, tns, prior.para = list(alp.prior.eff = 0.5,
                                                              bet.prior.eff = 0.5), mineff, effearly.stop){
  ###############################################################################
  ###############define the functions used for main function#####################
  ###############################################################################
  post.prob.fn <- function(target, y, n, alp.prior=0.1, bet.prior=0.1){
    alp <- alp.prior + y 
    bet <- bet.prior + n - y
    1 - pbeta(target, alp, bet)
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
  tunder.effs <- rep(0, length(txs))
  for(dose in 1:length(txs)){
    tadd.args <- prior.para
    cy <- tys[dose]
    cx <- txs[dose]
    cn <- tns[dose]
    prior.parause <- c(list(y = cy, n = cn, x = cx, tys = tys, txs = txs, tns = tns, cidx = dose), tadd.args)
    if (under.eff.fn(mineff, effearly.stop, prior.parause)) {
      tunder.effs[dose] <- 1
    }else{
      tunder.effs[dose] <- 0
    }
  }                                                            
  MTD <- CFO.selectmtd(target, tns, tys)$MTD
  if ( (MTD == 99) || (sum(tunder.effs[1:MTD]) == MTD)) {
    OBD <- 99
    OBD.probs <- NA
    if (MTD == 99){
      message("All tested doses are overly toxic. No OBD is selected! \n")
    }
    
    if ((sum(tunder.effs[1:min(length(txs), MTD)])) == min(length(txs), MTD)){
      message("All tested dose levels show low efficacy. No OBD is selected! \n")
    }
  }else{
    OBD.probs <- moveprobs(txs[1:MTD], tns[1:MTD], prior.para$alp.prior.eff, prior.para$bet.prior.eff)
    OBD <- which.max(OBD.probs)
  }
  admiset = c(1:MTD)
  out <- list(OBD = OBD, MTD = MTD, OBD.probs = OBD.probs, admiset = admiset)
  class(out) <- c("cfo_decision", "cfo")
  return(out)
}