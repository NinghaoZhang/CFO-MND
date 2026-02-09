#' Select the MND by using CFO-MND method
#'
#' @param ays Number of toxicities observed at each dose level.
#' @param axs Number of efficacy responses observed at each dose level.
#' @param ans Number of patients treated at each dose level.
#' @param dose.val Dose values for each level.
#' @param prior.para List of prior parameters for toxicity and efficacy beta models.
#' @param gamma Tuning parameter for dose-adjusted efficacy.
#' @param target Target DLT rate.
#' @param cutoff.eli Cutoff to eliminate overly toxic doses.
#' @param early.stop Threshold for early stopping.
#' @param mineff Minimum acceptable efficacy rate.
#' @param PED Whether PED is included (TRUE or FALSE). If TRUE, ATE is estimated for each dose.
#'
#' @return A list describing the MND selection and related quantities. Possible components include:
#' \describe{
#'   \item{OBD}{Estimated optimal biological dose.}
#'   \item{MND}{Estimated maximum tolerated dose.}
#'   \item{ATE}{Estimated average treatment effect.}
#'   \item{canset}{Admissible dose set.}
#'   \item{canaxs}{Number of efficacy responses at admissible doses.}
#'   \item{canans}{Number of patients treated at admissible doses.}
#'   \item{dose.val.adm}{Dose values for admissible doses.}
#'   \item{q.hat}{Posterior mean efficacy for admissible doses.}
#'   \item{probs}{Probabilities for efficacy maximization in stage 2.}
#'   \item{selectMND}{Result from `selectMND` function.}
#' }
#' @export
#'
#' @examples
#' \donttest{
#' ays <- c(0, 1, 1, 2)
#' axs <- c(0, 1, 2, 2)
#' ans <- c(3, 3, 3, 3)
#' dose.val <- c(1, 2, 4, 6)
#' prior.para <- list(alp.prior = 0.3, bet.prior = 0.7, alp.prior.eff = 0.5, bet.prior.eff = 0.5)
#' res <- MND.select(ays, axs, ans, dose.val, prior.para, gamma = 0.5,
#'                   target = 0.3, mineff = 0.15, PED = FALSE)
#' res
#' #' }

MND.select <- function(ays, axs, ans, dose.val, prior.para, gamma,
                       target, cutoff.eli = 0.95, early.stop = 0.95,
                       mineff = 0.15, PED = FALSE){
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
  MND <- MND.level(dose.val.adm, q.hat, gamma, mineff = mineff, obd = OBD, PED = PED, estimate = TRUE)
  ATE <- NA
  if (PED && !is.na(MND) && MND != 99) {
    MND.mean <- efficacy.pos.mean(axs[MND], ans[MND], alp.prior.eff, bet.prior.eff)
    PED.mean <- efficacy.pos.mean(axs[1], ans[1], alp.prior.eff, bet.prior.eff)
    ATE <- MND.mean - PED.mean
  }
  list(OBD = OBD, MND = MND, ATE = ATE, canset = canset, canaxs = canaxs, canans = canans,
       dose.val.adm = dose.val.adm, q.hat = q.hat, probs = probs, selectMND = res)
}
