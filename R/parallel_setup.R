#' Setup Parallelization
#'
#' Sets up parallel processing with the specified number of cores.
#'
#' @param n_cores Number of cores to use. If NULL, uses detectCores() - 1.
#' @return The cluster object created for parallel processing.
#' @export
#'
setup_parallel <- function(n_cores = NULL) {
  requireNamespace("foreach")
  requireNamespace("doParallel")
  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores() - 1
  }
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  print(paste("Parallel processing set up with", n_cores, "cores"))
  return(cl)
}

#' Stop Parallelization
#'
#' Stops the parallel backend by shutting down the cluster.
#'
#' @param cl The cluster object to stop.
#' @export
stop_parallel <- function(cl) {
  parallel::stopCluster(cl)
}
