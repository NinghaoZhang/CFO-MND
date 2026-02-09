#' Select the next dose level by using CFO-MND method
#'
#' @param ays Number of toxicities observed at each dose level.
#' @param axs Number of efficacy responses observed at each dose level.
#' @param ans Number of patients treated at each dose level.
#' @param currdose Current dose level index.
#' @param dose.val Dose values for each level.
#' @param cohortsize.stage1 Number of patients per cohort in stage 1.
#' @param cohortsize.stage2 Number of patients per cohort in stage 2.
#' @param ncohort.stage1 Number of cohorts in stage 1.
#' @param ncohort.stage2 Number of cohorts in stage 2.
#' @param target Target DLT rate.
#' @param gamma Tuning parameter for dose-adjusted efficacy.
#' @param mineff Minimum acceptable efficacy rate.
#' @param PED Whether PED is included (TRUE or FALSE). If TRUE, ATE is estimated for each dose.
#' @param prior.para List of prior parameters for toxicity and efficacy beta models.
#' @param cutoff.eli Cutoff to eliminate overly toxic doses.
#' @param early.stop Threshold for early stopping.
#' @param effearly.stop Threshold for early stopping due to low efficacy.
#' @param admset Admissible dose set. If NULL, it is set internally to 1:MTD in stage 2.
#' @param MTD Maximum tolerated dose estimate.
#'
#' @return A list describing the next-dose decision and related quantities. Possible components include:
#' \describe{
#'   \item{nextdose}{Next dose level to assign (or 99 for stop).}
#'   \item{stage}{Current stage ("stage1" or "stage2").}
#'   \item{earlystop}{Indicator for early stopping (1 = stop, 0 = continue).}
#'   \item{MTD}{Estimated maximum tolerated dose (or 99 if no MTD is chosen).}
#'   \item{admset}{Admissible dose set in stage 2.}
#'   \item{decision}{Decision label from `CFO.next` in stage 1.}
#'   \item{odds.ratio}{Odds ratio from `CFO.next` in stage 1.}
#'   \item{toxprob}{Posterior overdose probability from `CFO.next` in stage 1.}
#'   \item{effprobs}{Posterior probabilities for efficacy maximization in stage 2.}
#'   \item{OBD.hat}{Estimated optimal biological dose in stage 2.}
#'   \item{q.hat}{Posterior mean efficacy for admissible doses.}
#'   \item{q.tilde}{Dose-adjusted efficacy for admissible doses.}
#'   \item{prob.sample}{Sampling probabilities over admissible doses.}
#' }
#' @export
#'
#' @examples
#' \donttest{
#' ays <- c(0, 1, 2, 0, 0)
#' axs <- c(1, 2, 3, 1, 0)
#' ans <- c(3, 6, 9, 3, 0)
#' currdose <- 3
#' dose.val <- c(1, 2, 3, 4, 5)
#' cohortsize.stage1 <- 3
#' cohortsize.stage2 <- 3
#' ncohort.stage1 <- 4
#' ncohort.stage2 <- 4
#' target <- 0.3
#' gamma <- 0.5
#' mineff <- 0.2
#' PED <- FALSE
#' result <- MND.next(ays, axs, ans, currdose, dose.val, cohortsize.stage1, cohortsize.stage2,
#'                    ncohort.stage1, ncohort.stage2, target, gamma, mineff, PED)
#' print(result)
#' }
MND.next <- function(ays, axs, ans, currdose, dose.val, cohortsize.stage1, cohortsize.stage2,
                     ncohort.stage1, ncohort.stage2, target, gamma, mineff,
                     PED = FALSE, prior.para = list(alp.prior = target, bet.prior = 1 - target,
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

    if (PED){
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
  if (PED){
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

  return(list(nextdose = nextdose, stage = stage, earlystop = 0, MTD = MTD, admset = admset,
       q.hat = q.hat, q.tilde = q.tilde, prob.sample = prob.sample))
}
