# sparseSOMNiBUS_simu

This repository contains simulation code for the manuscript:  
**“A novel high-dimensional model for identifying regional DNA methylation QTLs.”**

## 📁 Folder Structure

- **`ToLoadFast/`** — Contains the main function `sparseSmoothFitCV()` and its supporting scripts:
  - `sparseSmoothFitCV.R` — Wrapper function for model fitting with CV
  - `sparseSmoothGrid.R` — Penalty parameter grid setup
  - `sparseSmoothPred.R` — Prediction and post-fit functions
  - `fitProxGradCpp.cpp` — C++ implementation of the proximal gradient descent algorithm
  - `sparseOmegaCr.cpp` — Computes structured penalty matrices (`$\Omega^{(1)}$` and `$\Omega^{(2)}$`) used for sparsity-smoothness regularization
  - `utils.R` — Helper utilities

- **`Experiment-smooth/`**  
  - `Scripts/` — Simulation runners (e.g., `run_100snps_samp_50.R`)
  - `Summary_scripts/` — Code for summarizing simulation output

- **`ToLoad/Simu.R`** — Code for generating simulated data, used across experiments

## 🔧 Usage

1. Source all necessary functions in `ToLoadFast/`.
2. Use `Simu.R` to generate data under specified simulation settings.
3. Run the simulation using a script from `Experiment-*/Scripts/`.
4. Summarize results using the scripts in `Summary_scripts/`.

## 📄 Reference

If you use this code, please cite:

> Zhao K, et al. (2025). *A novel high-dimensional model for identifying regional DNA methylation QTLs.* [Manuscript under review].

