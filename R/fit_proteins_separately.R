# Huge, multi-run data require different pre-processing and fitting the model
# one protein at a time.

melt_join = function(count_df, prot_to_bc) {
  melted = data.table::melt(count_df,
                            id.vars = c('sample_id', 'run_id'),
                            variable.name = 'barcode',
                            value.name = 'count')

  prot_to_bc[melted, on = 'barcode'][order(sample_id, protein)]
}

filter_melt_join = function(pattern, inv, res, prot_to_bc) {
  if (inv) {
    melt_join(res[!grepl(pattern, sample_id)], prot_to_bc)
  } else {
    melt_join(res[grepl(pattern, sample_id)], prot_to_bc)
  }
}

read_pre = function(count_file, prot_to_bc) {
  res = fread(count_file, skip = 1) %>%
    dplyr::rename(sample_id = Barcode)

  res$run_id = stringr::str_extract(tail(strsplit(count_file, '/|\\\\')[[1]], 1), '[0-9]+')

  filter_melt_join(pattern = "[Pp]re", inv = FALSE,
                   res = res, prot_to_bc = prot_to_bc)
}

read_outputs_one = function(count_file, prot_to_bc, bh_wr_input) {
  # Read in the whole file
  res = fread(count_file, skip = 1) %>%
    dplyr::rename(sample_id = Barcode)

  # Get the run id from the file name
  run_i = stringr::str_extract(tail(strsplit(count_file, '/|\\\\')[[1]], 1), '[0-9]+')
  res$run_id = run_i

  # Pick out the barcodes that were well-represented IN THAT RUN
  run_wr = bh_wr_input[run_id == run_i]
  select_cols = c('sample_id', 'run_id', run_wr$barcode)

  patterns = c('[Bb]eads', '[Pp]re|[Bb]eads')
  inv = c(FALSE, TRUE)

  # Produce a list of two dataframes for beads and outputs.
  # Only include barcodes that were well-represented in that run.
  count_lists = purrr::map2(patterns, inv,
                            filter_melt_join,
                            res = res[,..select_cols],
                            prot_to_bc = prot_to_bc)

  names(count_lists) = c('beads', 'output')

  return(count_lists)
}

get_prot_bc_map = function(count_file) {
  data.table::transpose(fread(as.character(count_file),
                              nrows = 2,
                              header = FALSE))[-1,] %>%
    rlang::set_names(c('protein', 'barcode'))
}


# Read in the csvs from count path. Split into three for Pre, beads, and output
read_multirun = function(count_path,
                         count_threshold = 4,
                         verbose = TRUE){


  count_files = list.files(count_path,
                           pattern = '*.csv', full.names = TRUE)
  if (verbose) {
    message(paste0("* Reading data from: ", count_path))
    message(paste0("* Detected ", length(count_files), " count files"))
  }

  prot_to_bc = get_prot_bc_map(count_files[1])

  if (verbose) message(paste0("* Detected ", nrow(prot_to_bc), " barcodes allocated across ",
                              dplyr::n_distinct(prot_to_bc$protein),
                              " unique proteins in the design"))

  bh_pre = data.table::rbindlist(furrr::future_map(count_files, read_pre, prot_to_bc = prot_to_bc))

  if (verbose) {
    frac_wr = bh_pre[, .(frac_wr = round(mean(count > count_threshold), digits = 3)), by = run_id]
    message("Fraction of barcodes well-represented in Pre sample by run: ")
    message(paste0(capture.output(frac_wr), collapse = '\n'))
  }

  bh_wr_input = bh_pre[count > count_threshold,
                       .(barcode, run_id, pre_count = count)]

  count_dfs = purrr::map(count_files, read_outputs_one,
                                prot_to_bc = prot_to_bc,
                                bh_wr_input = bh_wr_input) %>%
    purrr::transpose() %>%
    purrr::map(data.table::rbindlist)

  count_list = list(pre = bh_pre,
                    wr_pre = bh_wr_input,
                    beads = count_dfs[[1]],
                    wr_output = count_dfs[[2]])

  return(count_list)
}

filter_multirun = function(count_list,
                           out_dir = NULL) {

  prot_samp_df = unique(count_list$wr_output[,.(sample_id, run_id, protein)])

  sample_id_df = unique(prot_samp_df[,.(sample_id)]) %>%
    separate(sample_id,
             into = c('strain', 'repl', 'plate'),
             sep = '-',
             remove = FALSE) %>%
    mutate(across(.fns = as.factor))

  id_df = sample_id_df[prot_samp_df, on = 'sample_id'] %>%
    mutate(across(.fns = as.factor))

  id_df$ps = factor(paste(id_df$protein, id_df$strain, sep = ':'))
  id_df = id_df[, `:=`(p_i = as.numeric(protein),
                       s_i = as.numeric(strain),
                       ps_i = as.numeric(ps))][]

  samp_strain_df = sample_id_df[,.(sample_id, strain)] %>%
    unique

  # n_nz = "number non-zero"
  bh_n_nz = samp_strain_df[count_list$wr_output, on = 'sample_id'] %>%
    .[, .(n_nz = sum(count != 0),
          p_nz = mean(count != 0)),
      by = .(strain, protein)]

  ps_df = id_df[, .(sample_id, strain, protein, ps)] %>%
    unique

  bh_wr_ixns = ps_df[bh_n_nz[n_nz >= 3 | p_nz > .5], on = .(strain, protein)] # well represented
  bh_wr_output = ps_df[count_list$wr_output, on = .(sample_id, protein)][ps %in% bh_wr_ixns$ps]

  bc_input_by_run = bh_pre[, bc_run := paste(barcode, run_id, sep = '_')][bc_run %in% bh_wr_input$wr_bc_runs] %>%
    .[,.(barcode, run_id, count)] %>%
    dplyr::rename(pre_count = count)
}

fit_models = function (splits_dir,
                       out_dir) {

}

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
                                     algorithm = 'vb',
                                     save_split = TRUE,
                                     save_fits = FALSE) {

  if (!is.null(out_dir) & !dir.exists(out_dir)) {
    if (dir.exists(out_dir)){
      stop('output directory already exists')
    }
    dir.create(out_dir)
  }

  count_list = read_multirun(count_path) # a list of four count dataframes: pre, WR pre, beads, and WR output

  filtered_data = filter_multirun(count_list)

  write_splits(filtered_data)

  identify_bead_binders(filtered_outputs) # some proteins bind the beads strongly

  fit_models(algorithm = algorithm,
             save_fits = save_fits,
             num_threads = num_threads)

  if (!save_splits) {
    unlink(split_data_dir)
  }


}
