path <- '/Users/zhangninghao/Desktop/CFO-MND code'
setwd(path)
source("MND_utils.R")
library(pbapply)
library(pbmcapply)


#Fix Scenario example
nsim <- 200
scs <- list()
scs[[1]] <- list(p.true=c(0.15, 0.25, 0.3, 0.35, 0.4), q.true=c(0.2, 0.5, 0.5, 0.5, 0.5))
scs[[2]] <- list(p.true=c(0.05, 0.1, 0.3, 0.5, 0.6), q.true=c(0.2, 0.3, 0.5, 0.5, 0.5))
scs[[3]] <- list(p.true=c(0.05, 0.07, 0.12, 0.16, 0.18), q.true=c(0.05, 0.15, 0.39, 0.70, 0.72))
scs[[4]] <- list(p.true=c(0.05, 0.07, 0.1, 0.12, 0.16), q.true=c(0.35, 0.45, 0.5, 0.55, 0.75))
scs[[5]] <- list(p.true=c(0.05, 0.1, 0.3, 0.5, 0.6), q.true=c(0.15, 0.35, 0.5, 0.6, 0.7))
scs[[6]] <- list(p.true=c(0.05, 0.1, 0.3, 0.5, 0.6), q.true=c(0.15, 0.35, 0.6, 0.3, 0.2))
scs[[7]] <- list(p.true=c(0.05, 0.15, 0.25, 0.40, 0.45), q.true=c(0.08, 0.17, 0.45, 0.30, 0.25))
scs[[8]] <- list(p.true=c(0.10, 0.22, 0.25, 0.3, 0.4), q.true=c(0.3, 0.6, 0.55, 0.35, 0.2))
scs[[9]] <- list(p.true=c(0.40, 0.5, 0.55, 0.6, 0.7), q.true=c(0.15, 0.25, 0.5, 0.5, 0.5))
scs[[10]] <- list(p.true=c(0.35, 0.42, 0.58, 0.62, 0.7), q.true=c(0.15, 0.35, 0.6, 0.3, 0.2))



#Simulation
dose.val <- c(0.1, 1, 2, 5, 10, 20)
cutoff.eli <- 0.95
early.stop <- 0.95
effthreshold <- 0.9
cohortsize.stage1 <- 3
cohortsize.stage2 <- 3
ncohort.stage1 <- 6
ncohort.stage2 <- 14
target <- 0.3
mineff <- 0.2
init.level <- 1
alp.prior.eff <- 0.5
bet.prior.eff <- 0.5
alp.prior <-  target
bet.prior <-  1 - target
prior.para = list(alp.prior = target, bet.prior = 1 - target)
gamma <- 0.5




#Example with_out placebo
dir_path <- paste0(path, "results/test1")
CFOMND_simu(scs, nsim, dose.val, cutoff.eli, early.stop, effthreshold,
            cohortsize.stage1, cohortsize.stage2, ncohort.stage1, ncohort.stage2, target,
            mineff, init.level, alp.prior.eff, bet.prior.eff, alp.prior, bet.prior,
            gamma, dir_path)

#Example with placebo
dir_path <- paste0(path, "results/test2")
CFOMND_simu(scs, nsim, dose.val, cutoff.eli, early.stop, effthreshold,
            cohortsize.stage1, cohortsize.stage2, ncohort.stage1, ncohort.stage2, target,
            mineff, init.level, alp.prior.eff, bet.prior.eff, alp.prior, bet.prior,
            gamma, dir_path,TRUE)

#Real Trial
nsim <- 1000

dose.val <- c(1,3,6)
cutoff.eli <- 0.95
early.stop <- 0.95
effthreshold <- 0.9
cohortsize.stage1 <- 1
cohortsize.stage2 <- 1
ncohort.stage1 <- 6
ncohort.stage2 <- 12
target <- 0.2
mineff <- 0.15
init.level <- 1
alp.prior.eff <- 0.5
bet.prior.eff <- 0.5
alp.prior <-  target
bet.prior <-  1 - target
prior.para = list(alp.prior = target, bet.prior = 1 - target)
gamma <- 0.5


scs <- list()
scs[[1]] <- list(p.true=c(0.03, 0.04, 0.33), q.true=c(0.85, 0.93, 0.63))
dir_path <- paste0(path, "results/test3")
CFOMND_simu(scs, nsim, dose.val, cutoff.eli, early.stop, effthreshold,
            cohortsize.stage1, cohortsize.stage2, ncohort.stage1, ncohort.stage2, target,
            mineff, init.level, alp.prior.eff, bet.prior.eff, alp.prior, bet.prior,
            gamma, dir_path)

#Random Scenario
dose.val <- c(1, 2, 5, 10, 20)
cutoff.eli <- 0.95
early.stop <- 0.95
effthreshold <- 0.9
cohortsize.stage1 <- 3
cohortsize.stage2 <- 1
ncohort.stage1 <- 10
ncohort.stage2 <- 30
target <- 0.3
mineff <- 0.2
init.level <- 1
alp.prior.eff <- 0.5
bet.prior.eff <- 0.5
alp.prior <-  target
bet.prior <-  1 - target
prior.para = list(alp.prior = target, bet.prior = 1 - target)
gamma <- 0.5

ndose <- 5
dose.val <- dose.val <- c(1, 2, 5, 10, 20)
ndata <- 50
gamma <- 0.5
mu.list <- c(0.23, 0.38, 0.53, 0.71)

dir_path <- paste0(path, "results/test4/")
for (mu in mu.list){
  CFOMND_simu_rand(mu, ndata, ndose, dose.val, cutoff.eli, early.stop, effthreshold, 
                   cohortsize.stage1, cohortsize.stage2, ncohort.stage1, ncohort.stage2, target,
                   mineff, init.level, alp.prior.eff, bet.prior.eff, alp.prior, bet.prior, gamma, dir_path)
}


##############choose next dose
ays <- c(0, 1, 2, 0, 0)
axs <- c(1, 2, 3, 1, 0)
ans <- c(3, 6, 9, 3, 0)
currdose <- 3
dose.val <- c(1, 2, 3, 4, 5)
cohortsize.stage1 <- 3
cohortsize.stage2 <- 3
ncohort.stage1 <- 4
ncohort.stage2 <- 4
target <- 0.3
gamma <- 0.5
mineff <- 0.2
placebo <- FALSE
# result <- MND.next(ays, axs, ans, currdose, dose.val, cohortsize.stage1, cohortsize.stage2,
#                    ncohort.stage1, ncohort.stage2, target, gamma, mineff, placebo)
# print(result)