# Multi-run data require different pre-processing

#' Run the basehit model one protein at a time
#'
#' @inheritParams run_model
#' @param num_threads
#' @details The count file should have the first row specifying proteins, the second specifying
#'   barcodes, and all others after that specifying the output counts for each strain counts for
#'   each barcode (i.e. wide format, strain x barcode)\
#'
#'  Parallelized with furrr, so run a call to plan() that's appropriate for your system
#'
#' @export
model_proteins_separately = function(count_path,
                          out_dir = 'outputs/bh_out/',
                          out_name = 'results.xlsx',
                          model_path = 'src/stan_files/single_protein.stan',
                          algorithm = 'vb') {

}
