requireNamespace(testthat)
requireNamespace(PerturbPower)

test_that("Power simulation runs without errors", {
  load(file = "data/gRNA_counts.rda")
  load(file = "data/seurat_obj_metadata.rda")

  cl <- setup_parallel()
  clusterExport(cl, varlist = c("generate_celltype_df_with_variability", "assign_gRNA",
                                "modify_gene_effect", "test_gene_influence", "power_simulation"))
  clusterExport(cl, varlist = c("seurat_obj.metadata", "gRNA_counts"))


  result <- power_simulation(n_samples = 5, n_cells = 5000, nGenes_vals = c(8),
                             fraction_grnas_vals = c(0.8),
                             fraction_cell_proportion_change_vals = seq(1, 2, by = 0.2),
                             seurat_obj_metadata = seurat_obj.metadata,
                          gRNA_count_data = gRNA_counts, n_sims = 10)
  stop_parallel(cl)

  plot_power_simulation(result)

  expect_s3_class(result, "data.frame")
})
