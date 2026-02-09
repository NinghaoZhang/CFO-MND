#' CFOMND Simulation with random scenarios
#'
#' @param mu The parameter controls the target values for the average probability differences.
#' @param ndata The number of random scenarios to be generated.
#' @param ndose The number of dose levels in each random scenario.
#' @param dose.val The dosage corresponding to each dose level.
#' @param cutoff.eli The cutoff to eliminate overly toxic doses for safety.
#' @param early.stop The threshold value for early stopping.
#' @param effthreshold The threshold value for early stopping due to low efficacy.
#' @param cohortsize.stage1 The number of patients per cohort in the first stage.
#' @param cohortsize.stage2 The number of patients per cohort in the second stage.
#' @param ncohort.stage1 The number of cohorts in the first stage.
#' @param ncohort.stage2 The number of cohorts in the second stage.
#' @param target The target DLT rate.
#' @param mineff The lowest acceptable efficacy rate.
#' @param init.level The dose level assigned to the first cohort.
#' @param alp.prior.eff The prior alpha parameter for the efficacy beta distribution.
#' @param bet.prior.eff The prior beta parameter for the efficacy beta distribution.
#' @param alp.prior The prior alpha parameter for the DLT beta distribution.
#' @param bet.prior The prior beta parameter for the DLT beta distribution.
#' @param gamma Tuning parameter that modulates the trade-off between efficacy and dose, used to calculate the dose-adjusted efficacy measure.
#' @param dir_path The path for storing the results. If NULL, results are not saved to disk.
#'
#' @return A list with components:
#' \describe{
#'   \item{scs}{Generated random scenarios.}
#'   \item{results}{Simulation results for each scenario.}
#'   \item{analysis}{Summary from `random_scs_analysis`, including MND selection percentage (MND.sel.perc), overdose allocation percentage (overdose.allo.perc), MND allocation percentage (MND.allo.perc), and overdose selection percentage (overdose.sel.perc).}
#' }
#' @export
#'
#' @examples
#' \donttest{
#' dose.val <- c(1, 2, 5, 10, 20)
#' cutoff.eli <- 0.95
#' early.stop <- 0.95
#' effthreshold <- 0.9
#' cohortsize.stage1 <- 3
#' cohortsize.stage2 <- 1
#' ncohort.stage1 <- 10
#' ncohort.stage2 <- 30
#' target <- 0.3
#' mineff <- 0.2
#' init.level <- 1
#' alp.prior.eff <- 0.5
#' bet.prior.eff <- 0.5
#' alp.prior <- target
#' bet.prior <- 1 - target
#' gamma <- 0.5
#'
#' ndose <- 5
#' ndata <- 3
#' mu.list <- c(0.23, 0.38, 0.53, 0.71)
#' for (mu in mu.list) {
#'   CFOMND_simu_rand(mu, ndata, ndose, dose.val, cutoff.eli, early.stop, effthreshold,
#'                    cohortsize.stage1, cohortsize.stage2, ncohort.stage1, ncohort.stage2, target,
#'                    mineff, init.level, alp.prior.eff, bet.prior.eff, alp.prior, bet.prior,
#'                    gamma, dir_path = NULL)
#' }
#' }

CFOMND_simu_rand <- function(mu, ndata, ndose, dose.val, cutoff.eli, early.stop, effthreshold,
                             cohortsize.stage1, cohortsize.stage2, ncohort.stage1, ncohort.stage2, target,
                             mineff, init.level, alp.prior.eff, bet.prior.eff, alp.prior, bet.prior, gamma, dir_path = NULL){
  nsim = 1
  if (!is.null(dir_path)) {
    if (!dir.exists(paste0(dir_path,"/scs/"))) {
      dir.create(paste0(dir_path,"/scs/"), recursive = TRUE)
    }
    if (!dir.exists(paste0(dir_path,"res/"))) {
      dir.create(paste0(dir_path,"res/"), recursive = TRUE)
    }
  }
  phi <- 0.3
  psi <- 0.1
  psi.U <- 0.8
  npts <- cohortsize.stage1*ncohort.stage1+cohortsize.stage2*ncohort.stage2
  mu1 = mu2 = mu
  #generate random scenario
  scs <- scs.gen(ndose, dose.val, ndata, gamma, phi, psi, psi.U, mu1, mu2)
  if (!is.null(dir_path)) {
    save.file.name <- paste0(dir_path, "scs/random_scs",ndata,"gamma_",gamma,
                             "mu_",mu1,"target_",phi,"eff[",psi,",",psi.U,"]",".RData")
    save(scs, file = save.file.name)
  }
  #simulation
  run.fn.MND.random <- function(k){
    scs.test <- scs[[k]]
    p.true <- scs.test$ps
    q.true <- scs.test$qs
    res <- MND.simu(p.true, q.true, dose.val, cohortsize.stage1, cohortsize.stage2,
                    ncohort.stage1, ncohort.stage2, target, gamma, mineff, init.level,
                    cutoff.eli = cutoff.eli, early.stop = early.stop,
                    alp.prior.eff = alp.prior.eff, bet.prior.eff = bet.prior.eff,
                    prior.para = list(alp.prior = target, bet.prior = 1 - target, alp.prior.eff = alp.prior.eff, bet.prior.eff = bet.prior.eff))
    res$tMND <- scs.test$MND
    res$tMTD <- scs.test$k.MTD
    return(res)
  }


  #Random scenario
  save.file.name <- NULL
  if (!is.null(dir_path)) {
    save.file.name <- paste0(dir_path, "res/",ndata,"gamma_",gamma,"npts_",npts,
                             "mu_",mu1,"target_",phi,"eff[",psi,",",psi.U,"]",".RData")
  }
  ress.MND <- lapply(1:ndata, run.fn.MND.random)
  if (!is.null(dir_path)) {
    save(ress.MND, file=save.file.name)
  }
  ana <- random_scs_analysis(ress.MND, ndata, mu, cohortsize.stage1, cohortsize.stage2, ncohort.stage1, ncohort.stage2)
  print(ana)
  return(list(scs = scs, results = ress.MND, analysis = ana))
}
