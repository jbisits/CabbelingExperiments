# Cabbeling Experiments

Direct Numerical Simulation (DNS) experiments run as part of the study: *Cabbeling as a catalyst and driver of turbluent mixing*, submitted to the Journal of Fluid Mechanics.

## Reproducing the experiments

The DNS experiments are in [dns_runs](https://github.com/jbisits/CabbelingExperiments/tree/main/dns_runs) as julia scripts
with accompanying shell scripts to run them on a HPC.
They are currently configured to run on a GPU but can be configured to run on CPU by setting `architecture = CPU()` in any of the experiment scripts

## Analysis and reproducing figures

Saved copies of the analysis files can be accessed through figshare, DOI: [10.6084/m9.figshare.28012616.v2](https://doi.org/10.6084/m9.figshare.28012616.v2) (they can also be computed by running the analysis scripts in located in the experiment folders).
Figures from the manuscript can then be reproduced by running `notebooks/paper_plots_v3.jl` provided paths to the analysis files have been set correctly.
