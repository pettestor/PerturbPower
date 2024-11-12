
# PerturbPower: A Package for Power Simulations in Perturb-scRNAseq Screens

## Overview

The `PerturbPower` is an R package that provides tools to simulate and calculate power in CRISPR-based perturbation screens using single-cell data. This package allows users to simulate multiple samples, in silico transfect them with an gRNA library, and to analyze sensitivity and power across different experimental conditions.

## Installation

You can install the `PerturbPower` package by following these steps:

```r
devtools::install_github("https://github.com/pettestor/PerturbPower/")
```


## Example Usage

Here is a complete example of how to use the package for power simulations.

### Step 1: Load the Package

```r
library(PerturbPower)
library(foreach)
```

### Step 2: Set Parameters for Simulation

```r
# Parameters for the simulation
n_samples_vals <- c(2:8)  # Number of samples (typically 10X runs) to simulate
n_cells <- 5000           # Number of cells per sample
nGenes_vals <- c(8)       # Number of genes to modify
fraction_grnas_vals <- c(0.4, 0.6, 0.8, 1.0)  # Fraction of gRNAs to modify
fraction_cell_proportion_change_vals <- seq(1, 2, by = 0.2)  # Effect size
n_sims <- 5               # Number of simulations per condition, very low value to speed up example run.
```

### Step 3: Setup Parallel Processing

```r
# Setup parallel processing to speed up the computation
cl <- setup_parallel()
```

### Step 4: Load Data

Here we will use some example data from SeuratData and generate a synthetic gRNA library with five guides per gene for 150 genes. This data is provided within this repository.

```r

gene_universe = read.csv("https://raw.githubusercontent.com/pettestor/PerturbPower/refs/heads/main/data/gene_symbols.csv")
gRNA_counts = generate_grna_data(n_grna = 5,n_genes = 150,gene_universe = gene_universe$SYMBOL)
seurat_obj.metadata = read.csv("https://raw.githubusercontent.com/pettestor/PerturbPower/refs/heads/main/data/seurat_obj.metadata.ifnb.csv")

```

### Step 5: Run Power Simulations

The following loop runs power simulations for each number of samples and collects the results.

```r
final_results <- data.frame()

for (n_samples in n_samples_vals) {
  cat("Running power simulations for", n_samples, "samples...")

  start_time <- Sys.time()

  # Run the power simulation
  results <- power_simulation(n_samples = n_samples, 
                              n_cells = n_cells, 
                              nGenes_vals = nGenes_vals,
                              seurat_obj_metadata = seurat_obj.metadata,
                              n_sims = n_sims,
                              fraction_grnas_vals = fraction_grnas_vals, 
                              fraction_cell_proportion_change_vals = fraction_cell_proportion_change_vals,
                              gRNA_count_data = gRNA_counts)

  end_time <- Sys.time()

  # Calculate elapsed time
  elapsed_time <- difftime(end_time, start_time, units = "secs")
  cat("Time taken for", n_samples, "samples:", round(elapsed_time, 2), "seconds.")

  # Append the results
  results$n_samples <- n_samples
  final_results <- rbind(final_results, results)
}
```

### Step 6: Plot Results

```r
# Plot the power simulation results
plot_power_simulation(final_results)
```

This will generate a plot showing the sensitivity vs. the fraction of cell proportion changes across different sample sizes.

### Step 7: Stop Parallel Processing

After the simulations are completed, stop the parallel processing cluster:

```r
# Stop parallel processing
stop_parallel(cl)
```

## Additional Information

For detailed documentation on individual functions, please refer to the help files or vignettes provided in the package.

