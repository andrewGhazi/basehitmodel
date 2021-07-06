# Huge, multi-run data require different pre-processing and fitting the model
# one protein at a time.

#' Stats function to calculate the probability of a type S error from a set of
#' posterior samples
p_type_s = function(x){
  est_sign = sign(mean(x))

  1 - mean(sign(x) == est_sign)
}

#' Convert a table of basehit counts from wide to long and join on the header
#' line giving the mapping from barcodes to proteins
melt_join = function(count_df, prot_to_bc) {
  melted = data.table::melt(count_df,
                            id.vars = c('sample_id', 'run_id'),
                            variable.name = 'barcode',
                            value.name = 'count')

  prot_to_bc[melted, on = 'barcode'][order(sample_id, protein)]
}

#' melt_join either including or excluding a certain pattern
filter_melt_join = function(pattern, inv, res, prot_to_bc) {
  if (inv) {
    melt_join(res[!grepl(pattern, sample_id)], prot_to_bc)
  } else {
    melt_join(res[grepl(pattern, sample_id)], prot_to_bc)
  }
}

#' Read the pre counts
read_pre = function(count_file, prot_to_bc) {
  res = fread(count_file, skip = 1) %>%
    dplyr::rename(sample_id = Barcode)

  res$run_id = stringr::str_extract(tail(strsplit(count_file, '/|\\\\')[[1]], 1), '[0-9]+')

  filter_melt_join(pattern = "[Pp]re", inv = FALSE,
                   res = res, prot_to_bc = prot_to_bc)
}

#' Read the count file for one run
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

#' get the protein to barcode mapping
get_prot_bc_map = function(count_file) {
  data.table::transpose(fread(as.character(count_file),
                              nrows = 2,
                              header = FALSE))[-1,] %>%
    rlang::set_names(c('protein', 'barcode'))
}


#' Read in the csvs from count path. Split into three for Pre, beads, and
#' output.
#' @details Drop barcodes with less than or equal to count_threshold counts in
#'   the pre sample
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

  bh_wr_input = bh_pre[count > count_threshold, # this drops poorly represented barcodes
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

#' Filter multirun data
#' @description Filter multirun data to ixns with at least min_n_nz nonzero
#'   output counts or at least a proportion of nonzero counts > min_frac_nz
filter_multirun = function(count_list,
                           cache_dir = NULL,
                           save_outputs = TRUE,
                           check_cache = TRUE,
                           min_n_nz = 3,
                           min_frac_nz = .5,
                           drop_multirun_strains = TRUE) {


  if (check_cache) {
    if (file.exists(paste0(cache_dir, 'model_inputs.RData'))) {
      message(paste0("Loading cached model inputs from ", cache_dir,
                     " . Delete those and re-run if that's not what you want to use."))
      load(paste0(cache_dir, 'model_inputs.RData'))
      return(list(bh_input = bh_input, bh_eta = bh_eta))
    }
  }

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

  bh_wr_ixns = ps_df[bh_n_nz[n_nz >= min_n_nz | p_nz > min_frac_nz], on = .(strain, protein)] # well represented
  bh_wr_output = ps_df[count_list$wr_output, on = .(sample_id, protein)][ps %in% bh_wr_ixns$ps]

  bc_input_by_run = bh_pre[, bc_run := paste(barcode, run_id, sep = '_')][bc_run %in% bh_wr_input$wr_bc_runs] %>%
    .[,.(barcode, run_id, count)] %>%
    dplyr::rename(pre_count = count)

  io_by_run = bc_input_by_run[bh_wr_output, on = .(barcode, run_id)]

  if (drop_multirun_strains) {

    samples_to_drop = unique(bh_input[,.(sample_id, strain)]) %>%
      .[strain %in% c('AB9', 'AIEC', 'RG151') & !grepl('plt[1234]$', sample_id)]
    # ^ TODO generalize this for other strains/formats if necessary

    io_by_run = io_by_run[!(sample_id %in% samples_to_drop$sample_id)]
  }

  bh_input = io_by_run %>%
    .[, sd := sum(count) / 1e5, by = strain]

  bh_input$barcode = droplevels(barcode)
  bh_input$protein = droplevels(protein)
  bh_input$sample_id = droplevels(sample_id)
  bh_input$strain = droplevels(strain)
  bh_input$ps = droplevels(ps)

  bh_input = bh_input[,`:=`(s_i = as.integer(strain),
                            p_i = as.integer(protein),
                            ps_i = as.integer(ps))]

  bh_input$eta_i = as.integer(factor(paste(bh_input$pre_count, bh_input$strain, bh_input$protein,
                                           sep = ':')))

  # These are the unique levels of eta in the negative binomial regression. We pass them to the stan
  # model so that each only has to be computed once as a way to reduce redundancy in the computation.
  bh_eta = unique(bh_input[,.(pre_count, strain, sd, protein, p_i, s_i, ps_i, eta_i)])

  if (save_outputs) {
    save(bh_input, bh_eta,
         file = paste0(cache_dir, 'model_inputs.RData'))
  }

  return(list(bh_input = bh_input,
              bh_eta = bh_eta))
}

#' Write split data
#' @description Write out the data for each protein separately. It's written out
#'   to files named by index because sometimes proteins have weird characters or
#'   slashes in their names.
write_splits = function(filtered_data,
                        split_data_dir,
                        verbose = TRUE) {

  if (verbose) message("Writing out split data by protein...")

  bh_input = filtered_data$bh_input
  bh_eta = filtered_data$bh_eta

  prots = unique(bh_input$protein)
  n_write = length(prots)

  for (i in 1:n_write){
    fwrite(bh_input[protein == prots[i]],
           file = paste0(split_data_dir, i, '.csv.gz'))
    fwrite(bh_eta[protein == prots[i]],
           file = paste0(split_data_dir, i, '_eta.csv.gz'))
  }

  invisible()
}

get_entropy = function(values, offset = 0){
  p = (values + offset) / sum(values + offset)
  lp = log(p)
  nz = p != 0
  -sum(p[nz]*lp[nz])
}

identify_bead_binders = function(bh_input,
                                 drop_multirun_strains = TRUE,
                                 binding_threshold = 1) {

  s_id = bh_input[,.(sample_id)] %>%
    unique %>%
    tidyr::separate(sample_id,
             into = c('strain_ex', 'repl', 'plt'),
             sep = '-',
             remove = FALSE) %>%
    dplyr::select(sample_id, repl)

  count_input = bh_input[,.(protein, strain, sample_id, sd, pre_count, count, run_id)]

  if (drop_multirun_strains) {
    count_input = count_input[!(strain %in% c("AB9", 'AIEC', 'RG151') & (run_id != 1))]
  }

  adj_counts = count_input[, adj_count := count / pre_count][]

  concordance_scores = s_id[adj_counts, on = 'sample_id'] %>%
    .[,.(sum_adj = sum(adj_count)), by = .(protein, strain, repl)] %>%
    .[,.(concordance = get_entropy(sum_adj)), by = .(protein, strain)]

  concordance_scores$high_bead_binder = concordance_scores$concordance > binding_threshold

  return(concordance_scores)
}

read_split = function(i, split_data_dir) {
  p1 = fread(paste0(split_data_dir, i, '.csv.gz'))
  p1_eta = fread(paste0(split_data_dir, i, '_eta.csv.gz'))
  return(list(p1 = p1, p1_eta = p1_eta))
}

#' @param protein the protein
#' @param p1_eta the subset of bh_eta for the protein in question,
#' @param stan_model
fit_one_protein = function(protein,
                           split_data_dir,
                           save_fits,
                           prots,
                           ixn_prior_width,
                           algorithm,
                           out_dir) {

  i = which(prots == protein)

  split_data = read_split(i, split_data_dir)

  p1 = split_data$p1
  p1_eta = split_data$p1_eta

  if (prots[i] != protein | prots[i] != bh_input$protein[1]) {
    stop("Protein indexing error")
  }

  data_list = list(n_strain = nlevels(p1$strain),
                   n_bc = n_distinct(p1$barcode),
                   n_ps = nlevels(p1$strain), # because we're only looking at one protein
                   N = nrow(p1),
                   strain_i = p1$s_i,
                   n_eta = nrow(p1_eta),
                   eta_ps_i = p1_eta$ps_i,
                   eta_presum = log(p1_eta$sd) + log(p1_eta$pre_count),
                   eta_i = p1$eta_i,
                   ixn_prior_width = ixn_prior_width,
                   count_obs = p1$count)

  if (algorithm == 'variational') {
    # TODO expose more of the stan parameters to the user
    protein_fit = protein_model$variational(data = data_list,
                                            seed = 1234,
                                            output_samples = 5000)
  } else {
    protein_fit = protein_model$sample(data = data_list,
                                       seed = 1234,
                                       iter_sampling = 5000)
  }

  if (save_fits) {
    protein_fit$save_object(paste0(out_dir, '/protein_fit.RDS'))
  }

  protein_summary = protein_fit$summary('mean' = mean,
                                        'se' = posterior::mcse_mean,
                                        'med' = median,
                                        'qs' = ~quantile(.x, probs = c(.0025, .005, .025, .05, .1, .25,
                                                                       .75, .9, .95, .975, .995, .9975)),
                                        'conv' = posterior::default_convergence_measures(),
                                        'p_type_s' = p_type_s)
  protein_summary$protein = prot_name

  protein_summary = data.table::as.data.table(protein_summary)

  save(protein_summary, p1_eta,
       file = paste0(out_dir, paste0(prot_name, '.RData')))

  rm(protein_fit)
  gc()

  return(protein_summary)
}

fit_models = function (algorithm = algorithm,
                       split_data_dir = split_data_dir,
                       save_fits = save_fits,
                       ixn_prior_width,
                       out_dir,
                       seed) {

  protein_model = cmdstan_model('src/stan_files/single_protein.stan')

  prots = unique(bh_input$protein)

  summaries = furrr::future_map(.x = prots,
                                .f = fit_one_protein,
                                split_data_dir = split_data_dir,
                                save_fits = save_fits,
                                ixn_prior_width = ixn_prior_width,
                                prots = prots,
                                algorithm = algorithm,
                                out_dir = out_dir,
                                .options = furrr_options(seed = seed),
                                .progress = TRUE)

  res = data.table(protein = prots,
                   summary = summaries)

  return(res)
}

#' Run the basehit model one protein at a time
#'
#' @inheritParams run_model
#' @param num_threads
#' @details The count file should have the first row specifying proteins, the
#'   second specifying barcodes, and all others after that specifying the output
#'   counts for each strain counts for each barcode (i.e. wide format, strain x
#'   barcode)
#'
#'   Implemented with furrr, so run a call to plan() that's appropriate for your
#'   system in order to parallelize
#'
#' @export
model_proteins_separately = function(count_path,
                                     out_dir = 'outputs/bh_out/',
                                     cache_dir = 'outputs/bh_cache/',
                                     out_name = 'results.xlsx',
                                     model_path = 'src/stan_files/single_protein.stan',
                                     ixn_prior_width = .15,
                                     algorithm = 'variational',
                                     save_split = TRUE,
                                     save_fits = FALSE,
                                     seed = 1234) {

  if (!is.null(out_dir) & !dir.exists(out_dir)) {
    if (dir.exists(out_dir)){
      stop('output directory already exists')
    }

    #TODO add a check for cache_dir
    dir.create(out_dir)
    split_data_dir = file.path(gsub('\\/$', '', out_dir ), 'splits')
    dir.create(split_data_dir)
  }

  count_list = read_multirun(count_path) # a list of four count dataframes: pre, WR pre, beads, and WR output

  filtered_data = filter_multirun(count_list,
                                  cache_dir = cache_dir,
                                  save_outputs = TRUE)

  write_splits(filtered_data,
               split_data_dir)

  bb = identify_bead_binders(filtered_outputs) # some proteins bind the beads strongly

  model_fits = fit_models(algorithm = algorithm,
                          split_data_dir = split_data_dir,
                          save_fits = save_fits,
                          num_threads = num_threads,
                          ixn_prior_width = ixn_prior_width,
                          out_dir = out_dir,
                          seed = seed)

  write_summary(bb, model_fits)

  if (!save_splits) {
    unlink(split_data_dir)
  }

}
