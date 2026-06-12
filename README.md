# CellTransport

This repository contains an anonymized diffusion-advection simulation project for modeling transport dynamics in a generic system.

## Project overview

The current implementation generates synthetic time series data, fits a simple diffusion-advection model, and builds visualizations from the synthetic outputs. The code is organized into three main modules:

1. `DataGeneration/` - Creates anonymized synthetic signals that mimic transport-like dynamics.
2. `Model/` - Contains a simplified diffusion-advection model and a parameter search routine to fit the generated mean signal.
3. `Plots/` - Produces publication-style plots for the generated data and model fit.

The repository no longer depends on specific protein labels, experimental file paths, or sensitive local data references.

## How to run

From `Diffusion Advection Numerical Solution/`:

```powershell
python run_simulation.py
```

This script generates:

- `Results/synthetic_signals.csv` - anonymized synthetic time series dataset
- `Results/simulation_curves.png` - multiple synthetic series with the mean and standard deviation band
- `Results/signal_distribution.png` - distribution of synthetic signal values
- `Results/mean_fit_comparison.png` - mean synthetic signal compared to fitted model output

## Output and plots

The `Results/` directory contains the main visual outputs:

- a combined time-series plot of all synthetic signals
- a distribution plot for signal variability
- a model fit comparison plot with inferred diffusion, velocity, and decay parameters

These figures illustrate how the anonymized model behaves and how the fitted model tracks the generated data.

Inline preview of generated plots (stored under `Diffusion Advection Numerical Solution/Results/`):

![Simulation curves](Diffusion%20Advection%20Numerical%20Solution/Results/simulation_curves.png)

*Figure: Multiple anonymized synthetic series with mean and ±1 std band.*

![Signal distribution](Diffusion%20Advection%20Numerical%20Solution/Results/signal_distribution.png)

*Figure: Distribution of anonymized synthetic signal values.*

![Mean vs fit](Diffusion%20Advection%20Numerical%20Solution/Results/mean_fit_comparison.png)

*Figure: Mean synthetic signal compared to the fitted diffusion-advection model.*

## Notes

- The project has been sanitized to remove sensitive file paths and protein-specific names.
- The code now emphasizes reproducible synthetic simulation and plotting.
- The original biological reference is retained for context, but the current repository is structured around generic transport modeling.

## Reference

Dmitrieff S, Rao M, Sens P. Quantitative analysis of intra-Golgi transport shows intercisternal exchange for all cargo. Proc Natl Acad Sci U S A. 2013 Sep 24;110(39):15692-7. doi:10.1073/pnas.1303358110. PMCID: PMC3785783.
