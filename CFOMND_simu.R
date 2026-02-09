#' CFOMND Simulation with fixed scenarios
#'
#' @param scs A list containing the true DLT rates and true efficacy rates for each scenario.
#' @param nsim The total number of trials to be simulated.
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
#' @param alp.prior.eff The prior alpha parameter for the efficacy rate (Beta prior).
#' @param bet.prior.eff The prior beta parameter for the efficacy rate (Beta prior).
#' @param alp.prior The prior alpha parameter for the DLT rate (Beta prior).
#' @param bet.prior The prior beta parameter for the DLT rate (Beta prior).
#' @param gamma Tuning parameter that modulates the trade-off between efficacy and dose.
#' @param dir_path The path for storing the results.
#' @param with_PED Whether to introduce PED for estimating the ATE. If TRUE, PED is introduced
#'   and the ATE is estimated; if FALSE, PED is not introduced.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{analysis}{A list of per-scenario analysis outputs}
#'   \item{plot}{A combined \code{ggarrange} plot of scenario ATE scatter plots
#'     when \code{with_PED = TRUE}; otherwise \code{NULL}.}
#' }
#'
#' @export
#' @import ggplot2 ggpubr 
#' @importFrom stats dbeta dbinom integrate na.omit pbeta pnorm qnorm rbeta rbinom rnorm runif sd
#' @examples
#' \donttest{
#' nsim <- 2
#' scs <- list()
#' scs[[1]] <- list(p.true = c(0.15, 0.25, 0.3, 0.35, 0.4),
#'                 q.true = c(0.2, 0.5, 0.5, 0.5, 0.5))
#' dose.val <- c(0.1, 1, 2, 5, 10, 20)
#' cutoff.eli <- 0.95
#' early.stop <- 0.95
#' effthreshold <- 0.9
#' cohortsize.stage1 <- 3
#' cohortsize.stage2 <- 3
#' ncohort.stage1 <- 6
#' ncohort.stage2 <- 14
#' target <- 0.3
#' mineff <- 0.2
#' init.level <- 1
#' alp.prior.eff <- 0.5
#' bet.prior.eff <- 0.5
#' alp.prior <- target
#' bet.prior <- 1 - target
#' prior.para <- list(alp.prior = target, bet.prior = 1 - target)
#' gamma <- 0.5
#'
#' CFOMND_simu(scs, nsim, dose.val, cutoff.eli, early.stop, effthreshold,
#'             cohortsize.stage1, cohortsize.stage2, ncohort.stage1, ncohort.stage2,
#'             target, mineff, init.level, alp.prior.eff, bet.prior.eff, alp.prior,
#'             bet.prior, gamma)
#' }
CFOMND_simu <- function(scs, nsim, dose.val, cutoff.eli, early.stop, effthreshold,
                        cohortsize.stage1, cohortsize.stage2, ncohort.stage1, ncohort.stage2, target,
                        mineff, init.level, alp.prior.eff, bet.prior.eff, alp.prior, bet.prior, gamma, dir_path = NULL, with_PED = FALSE){
  #check the path
  if (!is.null(dir_path)) {
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
  }
  reslist <- vector("list", length(scs))
  #simulation
  for(i in 1:length(scs)){
    file.name <- if (!is.null(dir_path)) paste0(dir_path,"/fix_",i,"_5000.RData") else NULL
    if(with_PED){
      p.true <- scs[[i]]$p.true
      q.true <- scs[[i]]$q.true
      PED.q <- 0.1*scs[[i]]$q.true[1]
      PED.p <- 0.1*scs[[i]]$p.true[1]
      PED.dose.val <- 0.1*dose.val[1]

      p.true <- c(PED.p, p.true)
      q.true <- c(PED.q, q.true)
      dose.val <- c(PED.dose.val, dose.val)
    }else{
      p.true <- scs[[i]]$p.true
      q.true <- scs[[i]]$q.true
    }
    npts <- cohortsize.stage1*ncohort.stage1+cohortsize.stage2*ncohort.stage2
    plot_list <- list()
    run.fn.MND <- function(k){
      res <- MND.simu(p.true, q.true, dose.val, cohortsize.stage1, cohortsize.stage2,
                      ncohort.stage1, ncohort.stage2, target, gamma, mineff, init.level, with_PED,
                      prior.para = list(alp.prior = target, bet.prior = 1 - target, alp.prior.eff = 0.5, bet.prior.eff = 0.5))
    }
    ress.MND <- lapply(1:nsim, run.fn.MND)
    reslist[[i]] <- ress.MND
    if (!is.null(file.name)) {
      save(ress.MND, file = file.name)
    }
  }
  analysis_list <- vector("list", length(scs))
  #analyze
  for(i in 1:length(scs)){
    res <- reslist[[i]]
    if (with_PED){
      ana <- fix_scs_analysis(res, 0.3, dose.val, cohortsize.stage1, cohortsize.stage2, ncohort.stage1, ncohort.stage2, TRUE)
      scenario_name <- paste("Scenario", i)
      p <- plot_ate_scatter_with_CI(ana$data[2:length(ana$data)], ana$out$CI)
      p <- p + ggtitle(scenario_name)
      plot_list[[i]] <- p
    }else{
      ana <- fix_scs_analysis(res, 0.3, dose.val, cohortsize.stage1, cohortsize.stage2, ncohort.stage1, ncohort.stage2, FALSE)
    }
    analysis_list[[i]] <- ana

  }
  final_plot <- NULL
  if(with_PED){
    final_plot <- ggarrange(
      plotlist = plot_list,
      ncol = 2, nrow = ceiling(length(plot_list) / 2),
      common.legend = TRUE,
      legend = "bottom",
      align = "hv"
    )
  }
  return(list(analysis = analysis_list, plot = final_plot))
}
