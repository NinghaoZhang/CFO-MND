# CFO-MND

CFO-MND is a two-stage dose-finding workflow to identify the **MND** by CFO and efficacy-driven allocation.

## Repository Structure

- `MND_utils.R`: core utility functions (`MND.simu`, scenario generation/analysis utilities, `selectMND`, `CFO.next`, `CFO.selectmtd`).
- `CFOMND_simu.R`: fixed-scenario simulation wrapper.
- `CFOMND_simu_rand.R`: random-scenario simulation wrapper.
- `MND.next.R`: next-dose decision function.
- `MND.select.R`: final MND selection function.
- `Example.R`: runnable examples.
- `selectMND.R`: standalone helper mirror (kept for compatibility).

## Requirements

```r
install.packages(c("ggplot2", "ggpubr", "pbapply", "pbmcapply"))
```

## Quick Start

```r
source("MND_utils.R")
source("CFOMND_simu.R")
source("CFOMND_simu_rand.R")
source("MND.next.R")
source("MND.select.R")
```

### Next-dose recommendation

```r
ays <- c(0, 1, 2, 0, 0)
axs <- c(1, 2, 3, 1, 0)
ans <- c(3, 6, 9, 3, 0)

res <- MND.next(
  ays = ays,
  axs = axs,
  ans = ans,
  currdose = 3,
  dose.val = c(1, 2, 3, 4, 5),
  cohortsize.stage1 = 3,
  cohortsize.stage2 = 3,
  ncohort.stage1 = 4,
  ncohort.stage2 = 4,
  target = 0.3,
  gamma = 0.5,
  mineff = 0.2,
  PED = FALSE
)

print(res$nextdose)
```

### Fixed-scenario simulation

```r
scs <- list(list(p.true = c(0.15, 0.25, 0.30), q.true = c(0.20, 0.45, 0.50)))

res <- CFOMND_simu(
  scs = scs,
  nsim = 1,
  dose.val = c(1, 2, 3),
  cutoff.eli = 0.95,
  early.stop = 0.95,
  effthreshold = 0.9,
  cohortsize.stage1 = 1,
  cohortsize.stage2 = 1,
  ncohort.stage1 = 1,
  ncohort.stage2 = 1,
  target = 0.3,
  mineff = 0.2,
  init.level = 1,
  alp.prior.eff = 0.5,
  bet.prior.eff = 0.5,
  alp.prior = 0.3,
  bet.prior = 0.7,
  gamma = 0.5,
  dir_path = "results/readme-check-v2",
  with_PED = FALSE
)
```

Expected output file pattern:
- `results/<dir>/fix_<scenario_id>_5000.RData`

## Main APIs

- `CFOMND_simu(...)`: run fixed scenarios.
- `CFOMND_simu_rand(...)`: generate and run random scenarios.
- `MND.next(...)`: recommend next dose during an ongoing trial.
- `MND.select(...)`: select final MND after trial completion.

## Full Example Script

```r
source("Example.R")
```
