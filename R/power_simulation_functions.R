#' Generate gRNA Data with Normally Distributed Counts
#'
#' This function generates a data frame simulating gRNA count data for a given number of genes and gRNAs.
#' It uses gene names from the `pbmc3k` dataset and assigns normally distributed counts to each gRNA.
#'
#' @param n_grna Integer. The number of gRNAs to generate for each gene.
#' @param n_genes Integer. The number of unique genes to sample.
#' @param gene_universe A character vector of gene names to sample from (default: row names of the ifnb dataset).
#' @param mean_count Numeric. The mean of the normal distribution for generating count values.
#' @param sd_count Numeric. The standard deviation of the normal distribution for generating count values.
#'
#' @return A data frame with columns:
#' \item{count}{Normally distributed count values for each gRNA}
#' \item{ID}{Unique gRNA identifier in the format gene_symbol_X}
#' \item{target_gene_symbol}{The gene symbol associated with the gRNA}
#'
#' @examples
#' # Generate gRNA data for 3 genes, each with 4 gRNAs, normally distributed counts
#' generate_grna_data(n_grna = 4, n_genes = 3, mean_count = 100000, sd_count = 20000)
#'
#' @export
generate_grna_data <- function(n_grna, n_genes, gene_universe = rownames(ifnb@assays$RNA@counts),
                               mean_count = 100000, sd_count = 20000) {

  # Sample the number of genes as required
  sampled_genes <- sample(gene_universe, n_genes, replace = FALSE)

  # Create the data frame
  df <- data.frame()
  counter <- 1

  # Loop through each gene and generate gRNAs
  for (gene in sampled_genes) {
    for (i in 1:n_grna) {
      # Create the ID in the format gene_symbol_X
      id <- paste0(gene, "_", i)
      # Generate a random count value from a normal distribution
      count <- round(rnorm(1, mean = mean_count, sd = sd_count))
      # Append to the data frame
      df <- rbind(df, data.frame(count = count, ID = id, target_gene_symbol = gene))
    }
  }

  # Return the generated data frame
  return(df)
}



#' Modify Gene Effect in a Data Frame
#'
#' Modifies the effect of selected genes by adjusting cell type frequencies in a given data frame.
#' This function simulates the effect of a gene on a specific cell type, altering the frequency of
#' the cell type in the dataset.
#'
#' @param cell_df A data frame containing information about cells, including "CellType", "gRNA", and "Gene" columns.
#' @param nGenes Number of unique genes to modify.
#' @param fraction_grnas Fraction of gRNAs for each selected gene to modify.
#' @param y The fold change in the number of cells of the specified cell type.
#' @param celltype The target cell type to modify (default: "Dopamine neurons").
#' @param verbose Logical; if \code{TRUE}, detailed information will be printed.
#' @return A modified data frame where the specified cell type's frequency is adjusted for selected genes.
#' @examples
#' modified_df <- modify_gene_effect(cell_df, nGenes = 5, fraction_grnas = 0.8, y = 1.5)
modify_gene_effect <- function(cell_df, nGenes = 5, fraction_grnas = .8, y = 1.5,
                               celltype = "Dopamine neurons", verbose = TRUE) {
  # Check if the necessary columns are present
  if (!all(c("CellType", "gRNA", "Gene") %in% colnames(cell_df))) {
    stop("The data frame must contain 'CellType', 'gRNA', and 'Gene' columns.")
  }

  # Check if the provided cell type is a valid cell type in the data
  if (!celltype %in% cell_df$CellType) {
    stop(paste(celltype, "must be a valid cell type in the data frame."))
  }

  # Get unique genes
  unique_genes <- unique(cell_df$Gene)

  # Check if there are enough genes to modify
  if (length(unique_genes) < nGenes) {
    stop("Not enough unique genes to modify.")
  }

  # Randomly select nGenes genes
  selected_genes <- sample(unique_genes, nGenes)

  if (verbose) {
    cat("Modified genes:\n")
    cat(paste(selected_genes, collapse = " ", sep = " "), "\n")
  }

  # Loop through each selected gene to modify the gRNA effect
  for (gene in selected_genes) {
    # Get rows corresponding to the current gene
    gene_rows <- which(cell_df$Gene == gene)
    current_df <- cell_df[gene_rows, ]

    # Get unique gRNAs for this gene
    unique_grnas <- unique(current_df$gRNA)

    # Calculate how many gRNAs to modify based on the fraction_grnas
    num_grnas_to_modify <- round(length(unique_grnas) * fraction_grnas)

    # Randomly select a fraction of gRNAs to modify
    selected_grnas <- sample(unique_grnas, num_grnas_to_modify)

    if (verbose) {
      cat(paste("Modified gRNAs for gene", gene, ":\n"))
      cat(paste(selected_grnas, collapse = " ", sep = " "), "\n")
    }

    # Loop through each selected gRNA to modify the cell frequencies
    for (grna in selected_grnas) {
      # Get rows corresponding to the current gRNA
      grna_rows <- which(cell_df$gRNA == grna)
      current_df_grna <- cell_df[grna_rows, ]

      # Calculate the number of target cell types to add
      current_celltype_count <- sum(current_df_grna$CellType == celltype)
      additional_cells <- round(current_celltype_count * (y))

      # Increase the number of specified cell type
      non_target_rows <- which(current_df_grna$CellType != celltype)

      if (length(non_target_rows) > 0) {
        # Change the cell type of non-target rows to the specified celltype
        rows_to_change <- sample(non_target_rows, additional_cells, replace = TRUE)

        if (verbose) {
          cat("Rows to change for gRNA:\n")
          cat(grna, "\n")
          cat(grna_rows[rows_to_change], "\n")
        }

        cell_df$CellType[grna_rows[rows_to_change]] <- celltype
      }
    }
  }

  return(cell_df)
}

#' Test Gene Influence Using Mantel-Haenszel Test
#'
#' Performs the Mantel-Haenszel test to evaluate the influence of specific genes on cell type distributions
#' across samples. This test assesses whether the distribution of cell types is associated with the presence
#' of specific genes.
#'
#' @param cell_df A data frame containing columns "Sample", "CellType", and "Gene".
#' @return A data frame with results of the Mantel-Haenszel test for each gene, including p-values and adjusted p-values.
#' @examples
#' test_results <- test_gene_influence(cell_df)
test_gene_influence <- function(cell_df) {
  # Check if the necessary columns are present
  if (!all(c("Sample", "CellType", "Gene") %in% colnames(cell_df))) {
    stop("The data frame must contain 'Sample', 'CellType', and 'Gene' columns.")
  }

  # Remove cells with missing gRNA information
  cell_df <- cell_df[!is.na(cell_df$Gene), ]

  # Get the unique genes and cell types
  unique_genes <- unique(cell_df$Gene)
  unique_cell_types <- unique(cell_df$CellType)

  # Initialize a data frame to store the results with dynamic columns for cell types
  results <- data.frame(
    Gene = character(),
    p_value = numeric(),
    adjusted_p_value = numeric(),
    stringsAsFactors = FALSE
  )

  # Add columns for each cell type to the results data frame
  for (cell_type in unique_cell_types) {
    results[[cell_type]] <- numeric()
  }

  # Loop through each gene and perform the Mantel-Haenszel test
  for (gene in unique_genes) {
    # Create a binary variable for presence of the current gene
    cell_df$Gene_present <- ifelse(cell_df$Gene == gene, 1, 0)

    # Create a 3D contingency table: CellType x Gene_present x Sample
    contingency_table <- table(cell_df$CellType, cell_df$Gene_present, cell_df$Sample)

    # Calculate cell type frequencies for this gene
    cell_type_freq <- prop.table(table(cell_df$CellType[cell_df$Gene == gene])) * 100

    # Prepare a row for this gene's results
    result_row <- data.frame(
      Gene = gene,
      p_value = NA,  # Placeholder for now
      adjusted_p_value = NA,  # Placeholder for now
      stringsAsFactors = FALSE
    )

    # Add frequencies for each cell type to the result row
    for (cell_type in unique_cell_types) {
      result_row[[cell_type]] <- ifelse(cell_type %in% names(cell_type_freq),
                                        round(cell_type_freq[cell_type], 2),
                                        0)
    }

    # Check if the table has sufficient data for testing
    if (all(dim(contingency_table) > 1)) {  # At least 2 rows, 2 columns, and 2 strata needed
      # Perform Mantel-Haenszel test
      test_result <- mantelhaen.test(contingency_table)
      result_row$p_value <- test_result$p.value
    }

    # Append the result row to the results data frame
    results <- rbind(results, result_row)
  }

  # Adjust p-values for multiple testing (Benjamini-Hochberg method)
  results$adjusted_p_value <- p.adjust(results$p_value, method = "BH")

  return(results)
}


#' Run Power Simulation
#'
#' This function runs power simulations for cell-based experiments with specified parameters.
#'
#' @param n_samples Number of samples to simulate.
#' @param n_cells Number of cells per sample.
#' @param nGenes_vals Vector of gene counts to modify.
#' @param fraction_grnas_vals Vector of fractions of gRNAs to modify in each gene.
#' @param fraction_cell_proportion_change_vals Vector of effect sizes (cell type proportion changes).
#' @param n_sims Number of simulations to run.
#' @param significance_level The significance level for hypothesis testing.
#' @param seurat_obj_metadata Metadata from a Seurat object
#' @param gRNA_count_data Path to gRNA data
#'
#' @return A data frame of sensitivity and specificity across simulations.
#' @export
power_simulation <- function(n_samples, n_cells, nGenes_vals, fraction_grnas_vals, fraction_cell_proportion_change_vals, n_sims = 10,
                             significance_level = 0.05, seurat_obj_metadata = seurat_obj.metadata,
                             annotation = "seurat_annotation",
                             celltype = "B",
                             gRNA_count_data) {

  # Load libraries
  requireNamespace("dplyr")
  requireNamespace("foreach")
  requireNamespace("doParallel")

  # Placeholder to store results
  power_results <- data.frame()
  all_simulation_results <- list()  # List to store results of each simulation

  # Counter for storing in the list
  sim_index <- 1

  # Loop over all combinations of nGenes, fraction_grnas, and fraction_cell_proportion_change
  for (nGenes in nGenes_vals) {
    for (fraction_grnas in fraction_grnas_vals) {
      for (fraction_cell_proportion_change in fraction_cell_proportion_change_vals) {

        # Run simulations in parallel
        simulation_counts <- foreach(sim = 1:n_sims, .combine = 'rbind', .packages = c('dplyr'),
                                     .export = c('n_cells', 'nGenes_vals', 'fraction_grnas_vals',
                                                 'fraction_cell_proportion_change_vals',
                                                 'gRNA_counts', 'seurat_obj.metadata', 'n_sims',
                                                 'generate_celltype_df_with_variability', 'assign_gRNA',
                                                 'modify_gene_effect','test_gene_influence')) %dopar% {

          # Generate initial cell type data
          cell_df <- generate_celltype_df_with_variability(n_cells  = 5000,
                                                           seurat_obj_metadata = seurat_obj.metadata,
                                                           n_samples = n_samples)

          # Assign gRNAs and genes to cells
          cell_df <- assign_gRNA(cell_df, gRNA_count_data = gRNA_count_data, missing_prob = 0.2)

          # Apply gene modification effect
          modified_df <- modify_gene_effect(cell_df, nGenes = nGenes, fraction_grnas = fraction_grnas,
                                            y = fraction_cell_proportion_change,
                                            celltype = celltype, verbose = FALSE)

          modified_genes <- anti_join(cell_df, modified_df) %>% distinct(Gene) %>% pull()

          # Run the Mantel-Haenszel test to check for gene influence
          test_results <- test_gene_influence(modified_df)

          true_pos <- test_results %>% filter(Gene %in% modified_genes) %>% filter(adjusted_p_value < significance_level) %>% nrow()
          false_pos <- test_results %>% filter(!(Gene %in% modified_genes)) %>% filter(adjusted_p_value < significance_level) %>% nrow()
          true_neg <- test_results %>% filter(!(Gene %in% modified_genes)) %>% filter(adjusted_p_value >= significance_level) %>% nrow()
          false_neg <- test_results %>% filter(Gene %in% modified_genes) %>% filter(adjusted_p_value >= significance_level) %>% nrow()

          sensitivity <- ifelse((true_pos + false_neg) > 0, true_pos / (true_pos + false_neg), 0)
          specificity <- ifelse((true_neg + false_pos) > 0, true_neg / (true_neg + false_pos), 0)

          # Return sensitivity and specificity for each simulation run
          return(c(sensitivity, specificity))
        }

        # Calculate average sensitivity and specificity across simulations
        avg_sensitivity <- mean(simulation_counts[, 1], na.rm = TRUE)
        avg_specificity <- mean(simulation_counts[, 2], na.rm = TRUE)

        # Calculate sd sensitivity and specificity across simulations
        sd_sensitivity <- sd(simulation_counts[, 1], na.rm = TRUE)
        sd_specificity <- sd(simulation_counts[, 2], na.rm = TRUE)

        # Store the results for each combination of nGenes, fraction_grnas, and fraction_cell_proportion_change
        power_results <- rbind(power_results, data.frame(
          nGenes = nGenes,
          fraction_grnas = fraction_grnas,
          fraction_cell_proportion_change = fraction_cell_proportion_change,
          sensitivity = avg_sensitivity,
          specificity = avg_specificity,
          sd_sensitivity = sd_sensitivity,
          sd_specificity = sd_specificity,
          n= nrow(simulation_counts)
        ))

        sim_index <- sim_index + 1
        if(sim_index %% 1000 == 0){
          print("Simulation: ")
          print(sim_index)
        }
      }
    }
  }

  # Return the overall power results and the detailed simulation results
  return(power_results)
}

#' Generate Cell Type Data with Variability
#'
#' Generates synthetic cell type data for a specified number of samples, introducing variability in cell type frequencies.
#' Each sample is created based on the frequencies of cell types in the provided Seurat object metadata.
#'
#' @param n_cells Total number of cells to generate.
#' @param seurat_obj_metadata Metadata from a Seurat object.
#' @param n_samples Number of samples to generate.
#' @param sample_col Column name in the Seurat metadata that contains sample identifiers (default: "orig.ident").
#' @param celltype_col Column name in the Seurat metadata that contains cell type labels (default: "NamedClusters").
#' @param enforce_celltype Cell type to enforce with a minimum number of cells (default: "Dopamine neurons").
#' @param min_cells Minimum number of cells of the enforced cell type to include in each sample (default: 50).
#' @return A data frame with generated cell types for the specified number of samples.
#' @examples
#' cell_df <- generate_celltype_df_with_variability(n_cells = 5000, seurat_obj_metadata = seurat_obj.metadata, n_samples = 5)
generate_celltype_df_with_variability <- function(n_cells, seurat_obj_metadata = seurat_obj.metadata,
                                                  n_samples, sample_col = "orig.ident",
                                                  celltype_col = "seurat_annotations",
                                                  enforce_celltype = "Dopamine neurons",
                                                  min_cells = 50) {
  # Get the metadata from the Seurat object
  metadata <- seurat_obj_metadata

  # Check if the sample column and cell type column exist
  if (!(sample_col %in% colnames(metadata))) {
    stop(paste("Column", sample_col, "not found in the Seurat object's metadata."))
  }

  if (!(celltype_col %in% colnames(metadata))) {
    stop(paste("Column", celltype_col, "not found in the Seurat object's metadata."))
  }

  # Create an empty data frame to store the results
  result_df <- data.frame()

  # Get unique cell types and sample identifiers
  unique_celltypes <- unique(metadata[[celltype_col]])
  unique_samples <- unique(metadata[[sample_col]])

  # Calculate mean frequencies and standard deviations for each cell type across samples
  freq_list <- list()
  for (celltype in unique_celltypes) {
    celltype_freqs <- sapply(unique_samples, function(sample) {
      # Get counts for the current cell type in the current sample
      count <- sum(metadata[[celltype_col]] == celltype & metadata[[sample_col]] == sample)
      total_cells <- sum(metadata[[sample_col]] == sample)
      if (total_cells == 0) return(0)
      return(count / total_cells)  # Frequency of the cell type in the sample
    })

    # Store the frequencies
    freq_list[[celltype]] <- celltype_freqs
  }

  # Convert the list to a data frame for easier processing
  freq_df <- do.call(rbind, freq_list)
  rownames(freq_df) <- unique_celltypes

  # Calculate mean and standard deviation for each cell type
  mean_freqs <- rowMeans(freq_df, na.rm = TRUE)
  sd_freqs <- apply(freq_df, 1, sd, na.rm = TRUE)

  # Loop to generate n_samples new samples
  for (i in 1:n_samples) {
    # Generate cell type frequencies for the new sample
    new_celltype_freqs <- rnorm(length(mean_freqs), mean = mean_freqs, sd = sd_freqs)
    new_celltype_freqs[new_celltype_freqs < 0] <- 0  # Ensure non-negative frequencies
    new_celltype_freqs <- new_celltype_freqs / sum(new_celltype_freqs)  # Normalize to sum to 1

    # Sample cell types according to the generated frequencies
    cell_types <- sample(
      names(mean_freqs),
      size = n_cells - min_cells,  # Reduce the total by the number of enforced cells
      replace = TRUE,
      prob = new_celltype_freqs
    )

    # Add the enforced cell type ('Dopamine neurons') with a minimum number of cells
    enforced_cell_types <- rep(enforce_celltype, min_cells)

    # Combine the sampled cell types with the enforced cell types
    cell_types <- c(cell_types, enforced_cell_types)

    # Shuffle the combined cell types to mix them
    cell_types <- sample(cell_types)

    # Create a data frame for this new sample
    sample_df <- data.frame(
      CellID = paste0("Sample", i, "_", sprintf("%05d", 1:n_cells - 1)),
      Sample = paste("NewSample", i),
      CellType = cell_types
    )

    # Append to the result data frame
    result_df <- rbind(result_df, sample_df)
  }

  return(result_df)
}


#' Assign gRNAs to Cells Based on Frequency
#'
#' Assigns gRNAs to cells from a provided list of gRNAs with associated genes,
#' allowing some cells to have missing gRNA information.
#'
#' @param cell_df A data frame representing the cells to which gRNAs will be assigned.
#' @param gRNA_count_data gRNA_count_data data frame, including "ID", "target_gene_symbol", and "count" columns.
#' @param missing_prob Probability of a cell having missing gRNA information (default: 0.2).
#' @return A modified data frame with gRNA and gene assignments, and a column indicating gRNA detection status.
#' @examples
#' assigned_df <- assign_gRNA(cell_df, gRNA_count_data = gRNA_count_data, missing_prob = 0.2)
assign_gRNA <- function(cell_df, gRNA_count_data, missing_prob = 0.2) {
  # Validate the missing probability
  if (!is.numeric(missing_prob) || missing_prob < 0 || missing_prob > 1) {
    stop("The 'missing_prob' must be a number between 0 and 1.")
  }

  # Calculate frequencies (probabilities) from counts
  gRNA_count_data$Frequency <- gRNA_count_data$count / sum(gRNA_count_data$count)

  # Check for required columns
  if (!all(c("ID", "target_gene_symbol", "count") %in% colnames(gRNA_count_data))) {
    stop("The gRNA data file must contain 'ID', 'target_gene_symbol', and 'count' columns.")
  }

  # Determine which cells will have missing gRNA
  missing_status <- runif(nrow(cell_df)) < missing_prob

  # Assign gRNAs to cells based on calculated frequencies
  assigned_gRNAs <- sample(gRNA_count_data$ID, size = nrow(cell_df), replace = TRUE, prob = gRNA_count_data$Frequency)

  # Map the assigned gRNAs to corresponding genes
  assigned_genes <- gRNA_count_data$target_gene_symbol[match(assigned_gRNAs, gRNA_count_data$ID)]

  # Add gRNA and Gene information to the cell data frame
  cell_df$gRNA <- assigned_gRNAs
  cell_df$Gene <- assigned_genes

  # Set gRNA and Gene to NA for missing cells
  cell_df$gRNA[missing_status] <- NA
  cell_df$Gene[missing_status] <- NA

  # Add a column indicating whether a gRNA was detected
  cell_df$gRNA_present <- ifelse(missing_status, "missing", "detected")

  return(cell_df)
}

#' Custom Color Palette
#'
#' A custom color palette designed for data visualization, providing a wide variety of vibrant, deep,
#' and pastel colors for different elements in a plot. This palette is useful when visualizing
#' multiple categories and ensuring that the colors are distinguishable and visually appealing.
#'
#' The color descriptions below are for reference to help understand the tone and usage of each color.
#'
#' @return A vector of hexadecimal color codes representing the custom color palette.
#' @examples
#' palette <- get_custom_palette()
#' barplot(1:15, col = palette, main = "Barplot with Custom Palette")
get_custom_palette <- function() {
  palette <- c(
    "#E63946",  # Vibrant Red
    "#F77F00",  # Bright Orange
    "#F1C40F",  # Bright Yellow
    "#457B9D",  # Deep Blue
    "#1D3557",  # Darkest Blue (Navy)
    "#5D9CEC",  # Medium Blue (Slightly softer but not too light)
    "#00BFFF",  # Deep Sky Blue
    "#D50032",  # Crimson
    "#FF6F61",  # Coral
    "#8D3DAF",  # Purple
    "#28B463",  # Vibrant Green
    "#2ECC71",  # Emerald Green
    "#3498DB",  # Light Blue
    "#9B59B6",  # Amethyst (Purple tone)
    "#E74C3C",  # Strong Red
    "#F39C12",  # Vivid Orange-Yellow
    "#27AE60",  # Darker Green
    "#16A085",  # Dark Teal
    "#2980B9",  # Rich Blue
    "#D35400",  # Pumpkin Orange
    "#BDC3C7",  # Silver
    "#7F8C8D",  # Medium Gray
    "#E0E0E0",  # Light Gray
    "#F8C471",  # Soft Yellow
    "#A93226",  # Dark Red
    "#7D3C98"   # Dark Purple
  )
  return(palette)
}


#' Plot Sensitivity vs. Fraction Cell Proportion Change from Simulation Results
#'
#' This function generates a plot visualizing the relationship between sensitivity and fraction cell proportion change for a simulation study.
#' The plot shows how sensitivity changes across different sample sizes, while allowing the user to customize the color palette.
#'
#' @param results A data frame containing simulation results. The data frame must have the following columns:
#' \itemize{
#'   \item \code{fraction_cell_proportion_change}: Numeric values representing the change in the fraction of cell proportions.
#'   \item \code{sensitivity}: Numeric values representing the sensitivity of the simulation.
#'   \item \code{n_samples}: Integer values representing the number of samples used in the simulation.
#'   \item \code{sd_sensitivity}: Numeric values representing the standard deviation of sensitivity.
#'   \item \code{nGenes}: Numeric values or factors representing the number of genes in the simulation.
#'   \item \code{fraction_grnas}: Numeric values or factors representing the fraction of gRNAs in the simulation.
#' }
#' @param palette A vector of colors to be used in the plot. Default is \code{get_custom_palette()}.
#'
#' @return A ggplot object showing the sensitivity vs. fraction cell proportion change plot, faceted by the number of genes and fraction of gRNAs.
#'
#' @details
#' This function creates a line and scatter plot using \code{ggplot2} to visualize the sensitivity of the simulation as a function of the fraction of cell proportion change.
#' The plot is grouped by the number of samples and faceted by the number of genes and fraction of gRNAs. The point size corresponds to the standard deviation of the sensitivity.
#'
#' @examples
#' \dontrun{
#'   # Assuming `simulation_results` is a data frame containing the required columns:
#'   plot_power_simulation(simulation_results)
#' }
#'
#' @import ggplot2
#' @import cowplot
#' @import scales
#' @export
plot_power_simulation <- function(results,palette = get_custom_palette()){

  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("cowplot", quietly = TRUE) ||
      !requireNamespace("scales", quietly = TRUE)) {
    stop("Please install ggplot2, cowplot, and scales packages to use this function.")
  }

  ggplot2::ggplot(results, aes(x = fraction_cell_proportion_change, y = sensitivity)) +
    geom_line(aes(group = n_samples, color = as.factor(n_samples))) +
    geom_point(aes(group = n_samples, ,size = sd_sensitivity,color = as.factor(n_samples))) +
    facet_grid(nGenes ~ fraction_grnas) +
    labs(
      title = "Sensitivity vs. Fraction Cell Proportion Change",
      x = "Fraction Cell Proportion Change",
      y = "Sensitivity",
      color = "Number of Samples"
    ) +
    cowplot::theme_cowplot()+scale_color_manual(values = palette) +
    ylim(c(0,1))+
    scale_x_continuous(labels = scales::percent_format(scale = 100)) +
    scale_y_continuous(labels = scales::percent_format(scale = 100)) +
    theme(
      panel.grid.major = element_line(color = "grey80", size = 0.5, linetype = "dashed"),  # Major grid lines
      panel.grid.minor = element_line(color = "grey90", size = 0.25, linetype = "dotted")  # Minor grid lines
    )
}
