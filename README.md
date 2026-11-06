# When can a joint species distribution model capture biotic interactions?

Code of the simulation study of Tous & Chiquet, *When can a joint species
distribution model capture biotic interactions? A simulation study across
interaction mechanisms and levels of model integration* (submitted to PCI
Ecology).

Communities are simulated under three interaction mechanisms (the lottery of
VirtualCom, a multispecies Gompertz and a multispecies Ricker), each world also
simulated without interactions as a paired control. Four estimators of
increasing integration (log-GGM, PLNPCA + glasso, PLN + glasso, PLNnetwork)
rank the pairs of species, and are scored by AUC and top-k precision, as
increments over the paired control.

## Credits

Jeanne Tous and Julien Chiquet. The species pool and the lottery are adapted
from VirtualCom (Münkemüller T. & Gallien L. 2015, VirtualCom: a simulation
model for eco-evolutionary community assembly and invasion, *Methods in
Ecology and Evolution* 6: 735–743; code at
<https://sourceforge.net/projects/virtualcom2014/>).

Licence: GNU General Public License, version 3 or later (file `LICENSE`).

## Use of generative AI

The first version of this code was written entirely by Jeanne Tous and Julien
Chiquet, and built its simulations on VirtualCom (Münkemüller & Gallien 2015).
For the second version, prepared for the resubmission, Claude Code (Anthropic)
was used to:

- reformat the code and reorganise the scripts, starting from the code of the
  first version and from that of VirtualCom;
- fix bugs in the code of VirtualCom and make it faster;
- set up the new simulation campaign, with the two new interaction mechanisms
  (multispecies Gompertz, a first-order autoregressive model on the log scale,
  and multispecies Ricker, a discrete Lotka-Volterra model);
- write this README, which it generated entirely.

Throughout, the authors read and checked the code as it was generated, and
amended it where needed to make it their own; they checked that it reproduces
the results of the paper, and take full responsibility for it.

## Requirements

R 4.6.1 and the packages PLNmodels 1.3.2, gllvm 2.0.13, geiger, ape, picante,
plyr, R.utils, RhpcBLASctl, jsonlite, ggplot2, Hmisc; Hmsc 3.3.7 and coda for
`hmsc_convergence.R`. The runs fork processes (`parallel`), hence Linux or
macOS. Pin BLAS to one thread before starting R:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
```

## Layout

```
R/                      functions, loaded by R/load.R
  simulate_species.R    niche optima (phylogeny, Brownian motion; VirtualCom)
  simulate_networks.R   interaction networks, random reference networks
  simulate_mechanisms.R lottery, Gompertz, Ricker
  simulate_data.R       one data set: world, assembly, sampling effort, covariate
  fit_models.R          PLN / PLNnetwork, PLNPCA, gllvm, log-GGM
  scoring.R             rankings of pairs, AUC, top-k precision
  replicate.R           abiotic bases, one replicate (both arms), a grid of settings
  results_io.R          result files
run_simulations.R       the campaigns of the paper
describe_simulations.R  descriptions of the simulated data (no model fitted)
hmsc_convergence.R      convergence of Hmsc
analysis/               numbers and figures of the paper
numbers/                tables of numbers cited in the paper (paper_numbers.md)
```

## Reproducing the paper

All commands run from the root of the repository. Results are written to
`results/` (or to the folder given by `VC_RESULTS`).

1. Simulations (one command per campaign; days on a 30-core server):

   ```bash
   Rscript run_simulations.R world         # main design and supplementary experiments
   Rscript run_simulations.R exploration   # abundance level K under Gompertz
   Rscript run_simulations.R emission_pois5
   Rscript run_simulations.R emission_nb5
   Rscript run_simulations.R emission_nb10
   Rscript describe_simulations.R          # a few minutes
   Rscript hmsc_convergence.R              # about 40 minutes on 16 cores
   ```

   A single grid can be run alone (`Rscript run_simulations.R world
   world_mech_smooth_lottery`); `VC_NSIMU=2` gives a quick trial.

2. Numbers and figures:

   ```bash
   Rscript analysis/paper_numbers.R        # numbers/T*.csv, numbers/paper_numbers.md
   Rscript analysis/figures_ladder.R       # figures/ladder_*.png
   Rscript analysis/figures_bases.R        # figures/smooth_{2,3,7,8}_*.png
   Rscript analysis/figures_compo.R        # figures/compo_*.png
   Rscript analysis/figures_supplement.R   # figures/smooth_{4,4b,5a,5b}_*.png, figures/annex_*.png
   ```

The result files of the paper (189 MB, zipped 53 MB) are deposited on Zenodo
(<https://doi.org/10.5281/zenodo.24001750>). Unzipped, they recompute every number
and figure without rerunning the simulations: set `VC_RESULTS` to their
`results/` folder before running the scripts of step 2.

### Exactness

Every replicate draws its world from its own seed, a function of the setting
and of the replicate: the simulated data, PLN, PLNnetwork, PLNPCA, the log-GGM
and the reference rankings are reproduced exactly with the package versions
above. PLNnetwork depends on the version of PLNmodels: the development
version 1.3.2.9500 already gives different penalty paths. The starting values of gllvm are random and drawn in forked processes,
whose random streams follow the order in which processes are started: gllvm
fits are reproduced exactly with the number of cores of the paper (the default
of `run_simulations.R`; 30 for `world`, 14 for `exploration`, 4 and 14 for
`emission_*`), and otherwise differ within the variability of gllvm's starts.
`hmsc_convergence.R` draws the same 30 data sets as for the paper; 34 of
its 60 groups of four chains are reproduced exactly, and the table of the
Supplementary material is reproduced to the precision it is printed with
(one proportion moves from 67 % to 68 %). It takes about 40 minutes on 16
cores.

## Result files

`results/<campaign>/grid_res_<grid>.csv`: one row per setting, replicate
(`simu`), arm (`interactions`, `no.interaction`) and method. Methods are named
`<family>.<basis>[.<reading>]`:

- families: `PLN` (PLNnetwork; `.marginal`: the unpenalised PLN fit, read as a
  thresholded correlation or, in the `_alt` columns, through the graphical
  lasso: PLN + glasso), `gllvm`, `empirical` (log-GGM), `plnpca.<basis>.q<rank>`;
  `informed.random`, `fully.random`, `niche.proximity`: reference rankings
  that use no counts;
- bases: `no.covar` (none), `cat.2` (categorical, 10 classes), `quadratic`,
  `poly.3`, `spline.2` (B-spline, 10 columns), `both.quadratic` (quadratic in
  the observed and the hidden gradient);
- effort corrections: `.off` (offset log effort), `.tss` (offset log total, or
  counts scaled to the median total), `.div` (Y / effort), `.clr`.

Metrics: `AUC_competition`, `AUC_facilitation`, `precision_k_competition`,
`precision_k_facilitation`; with suffix `_alt`, the graphical-lasso reading of
a covariance. In the analysis, `quick` is the log-GGM.
