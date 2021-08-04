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
                         cache_dir = NULL,
                         verbose = TRUE) {

  if (!missing(cache_dir) && file.exists(paste0(cache_dir, 'count_list.RData'))) {
    if (verbose) message("* loading cached count data")
    load(paste0(cache_dir, 'count_list.RData'))
    return(count_list)
  }

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

  if (verbose) message('* Reading counts from Pre samples...')
  bh_pre = data.table::rbindlist(furrr::future_map(count_files, read_pre, prot_to_bc = prot_to_bc))

  if (verbose) {
    frac_wr = bh_pre[, .(frac_wr = round(mean(count > count_threshold), digits = 3)), by = run_id]
    message("* Fraction of barcodes well-represented in Pre sample by run: ")
    message(paste0(capture.output(frac_wr), collapse = '\n'))
  }

  bh_wr_input = bh_pre[count > count_threshold, # this drops poorly represented barcodes
                       .(barcode, run_id, pre_count = count)][, wr_bc_runs := paste(barcode, run_id, sep = '_')]

  if (verbose) message("* Reading output counts... ")
  count_dfs = purrr::map(count_files, read_outputs_one,
                                prot_to_bc = prot_to_bc,
                                bh_wr_input = bh_wr_input) %>%
    purrr::transpose() %>%
    purrr::map(data.table::rbindlist)

  count_list = list(pre = bh_pre,
                    wr_pre = bh_wr_input,
                    beads = count_dfs[[1]],
                    wr_output = count_dfs[[2]],
                    prot_to_bc = prot_to_bc)

  if (!missing(cache_dir)) {
    save(count_list,
         file = paste0(cache_dir, 'count_list.RData'))
  }

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
                           drop_multirun_strains = TRUE,
                           verbose = TRUE) {


  cache_dir = end_with_slash(cache_dir)

  if (check_cache) {
    if (!dir.exists(cache_dir)) {
      if (verbose) message("* Cache directory doesn't exist. Creating it...")
      dir.create(cache_dir)
    }

    if (file.exists(paste0(cache_dir, 'model_inputs.RData'))) {
      message(paste0("* Loading cached model inputs from ", cache_dir,
                     " . Delete those and re-run if that's not what you want to use."))
      load(paste0(cache_dir, 'model_inputs.RData'))
      return(list(bh_input = bh_input, bh_eta = bh_eta))
    }
  }

  bh_pre = count_list$pre
  bh_wr_input = count_list$wr_pre

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
    dplyr::rename("pre_count" = "count")

  io_by_run = bc_input_by_run[bh_wr_output, on = .(barcode, run_id)]

  if (drop_multirun_strains) {

    samples_to_drop = unique(io_by_run[,.(sample_id, strain)]) %>%
      .[strain %in% c('AB9', 'AIEC', 'RG151') & !grepl('plt[1234]$', sample_id)]
    # ^ TODO generalize this for other strains/formats if necessary

    io_by_run = io_by_run[!(sample_id %in% samples_to_drop$sample_id)]
  }

  bh_input = io_by_run %>%
    .[, sd := sum(count) / 1e5, by = strain]

  # I know some of these are already factors but not which...
  bh_input$barcode = droplevels(as.factor(bh_input$barcode))
  bh_input$protein = droplevels(as.factor(bh_input$protein))
  bh_input$sample_id = droplevels(as.factor(bh_input$sample_id))
  bh_input$strain = droplevels(as.factor(bh_input$strain))
  bh_input$ps = droplevels(as.factor(bh_input$ps))

  bh_input = bh_input[,`:=`(s_i = as.integer(strain),
                            p_i = as.integer(protein),
                            ps_i = as.integer(ps))]

  bh_input$eta_i = as.integer(factor(paste(bh_input$pre_count, bh_input$strain, bh_input$protein,
                                           sep = ':')))

  # These are the unique levels of eta in the negative binomial regression. We pass them to the stan
  # model so that each only has to be computed once as a way to reduce redundancy in the computation.
  bh_eta = unique(bh_input[,.(pre_count, strain, sd, protein, p_i, s_i, ps_i, eta_i)])

  if (save_outputs) {
    if (verbose) message(paste0("* Saving filtered inputs to: ", paste0(cache_dir, 'model_inputs.RData')))
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

  bh_input = filtered_data$bh_input
  bh_eta = filtered_data$bh_eta

  proteins = unique(bh_input$protein)
  n_write = length(proteins)

  existing_splits = list.files(split_data_dir,
                               pattern = ".csv.gz")
  files_to_write = c(paste0(1:n_write, ".csv.gz"),
                     paste0(1:n_write, "_eta.csv.gz"))
  if (all(files_to_write %in% existing_splits)) {
    message("* Split data already in place.")
    return(invisible())
  }
  for (i in 1:n_write){
    fwrite(bh_input[protein == proteins[i]],
           file = paste0(split_data_dir, i, '.csv.gz'))
    fwrite(bh_eta[protein == proteins[i]],
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

identify_bead_binders = function(wr_pre, prot_to_bc,
                                 bh_beads, bh_output,
                                 binding_threshold,
                                 save_bead_binders,
                                 out_dir,
                                 verbose = TRUE) {
  beads_avg = wr_pre[,.(barcode, pre_count, run_id)][bh_beads, on = .(barcode, run_id)][pre_count > 4][,.(mean_bead_count = mean(count / pre_count)), by = barcode]

  adj_output = wr_pre[,.(barcode, pre_count, run_id)][bh_output, on = .(barcode, run_id)][pre_count > 4][, pre_adj_count := count / pre_count][, .(mean_prot_output = mean(pre_adj_count)), by = 'protein']

  bead_output_by_prot = prot_to_bc[beads_avg, on = 'barcode'][, .(prot_mean_bead_out = mean(mean_bead_count)), by = 'protein']

  high_bead_binding = adj_output[bead_output_by_prot, on = 'protein'][prot_mean_bead_out > binding_threshold]

  if (verbose) message(paste0("* Identified ", nrow(high_bead_binding), " proteins with high output counts in bead samples..."))

  if (save_bead_binders){
    if (verbose) message(paste0("* Saving high bead binding protein data frame to: ", paste0(out_dir, 'high_bead_binding.RData')))
    save(high_bead_binding,
         file = paste0(out_dir, 'high_bead_binding.RData'))
  }

  return(high_bead_binding)
}

get_concordance = function(bh_input,
                           drop_multirun_strains = TRUE) {

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

  return(concordance_scores)
}

read_split = function(i, split_data_dir) {
  p1 = fread(paste0(split_data_dir, i, '.csv.gz'))
  p1_eta = fread(paste0(split_data_dir, i, '_eta.csv.gz'))
  return(list(p1 = p1, p1_eta = p1_eta))
}

#' Fit the model to one protein
#' @param protein the protein
#' @param p1_eta the subset of bh_eta for the protein in question,
#' @param stan_model protein model to use
fit_one_protein = function(protein,
                           # bh_input,
                           # protein_model,
                           split_data_dir,
                           save_fits,
                           proteins,
                           ixn_prior_width,
                           algorithm,
                           out_dir) {

  i = which(proteins == protein)
  model_path = system.file("stan", "single_protein_by_ixn_width.stan",
                           package = 'bhm',
                           mustWork = TRUE)
  protein_model = cmdstan_model(stan_file = model_path, quiet = TRUE)

  split_data = read_split(i, split_data_dir)

  p1 = split_data$p1
  p1_eta = split_data$p1_eta

  eta_map = p1_eta[,.(eta_i,
                      new_i = 1:nrow(p1_eta))]

  p1 = eta_map[p1, on = 'eta_i'] %>%
    select(-eta_i) %>%
    dplyr::rename(eta_i = new_i)

  p1_eta = eta_map[p1_eta, on = 'eta_i'] %>%
    select(-eta_i) %>%
    dplyr::rename(eta_i = new_i)

  p1[, ps := factor(paste(protein, strain, sep = ':'))]

  p1_eta[, ps := factor(paste(protein, strain, sep = ':'),
                        levels = levels(p1$ps))]
  p1_eta$eta_i = 1:nrow(p1_eta)

  p1$strain = factor(p1$strain)
  p1$s_i = as.numeric(p1$strain)

  p1_eta$strain = factor(p1_eta$strain,
                         levels = levels(p1$strain))
  p1_eta$s_i = as.numeric(p1_eta$strain)

  p1_eta$ps = droplevels(p1_eta$ps)
  p1_eta$ps_i = as.numeric(p1_eta$ps)

  # if (proteins[i] != protein | proteins[i] != bh_input$protein[1]) {
  #   stop("Protein indexing error")
  # }

  data_list = list(n_strain = dplyr::n_distinct(p1$strain),
                   n_bc = dplyr::n_distinct(p1$barcode),
                   n_ps = dplyr::n_distinct(p1$strain), # because we're only looking at one protein
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
                                            output_samples = 5000)
  } else {
    protein_fit = protein_model$sample(data = data_list,
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

  protein_summary$protein = protein
  protein_summary = data.table::as.data.table(protein_summary)
  protein_summary[, s_i := as.numeric(stringr::str_extract(variable, '[0-9]+'))]

  strain_i_df =  unique(p1[,.(strain, s_i)])

  protein_summary = strain_i_df[protein_summary, on = 's_i']

  # save(protein_summary, p1_eta,
  #      file = paste0(out_dir, paste0(prot_name, '.RData')))

  rm(protein_fit)
  gc()

  return(protein_summary)
}

fit_safely = purrr::safely(fit_one_protein)

fit_models = function (algorithm = algorithm,
                       split_data_dir = split_data_dir,
                       save_fits = save_fits,
                       bh_input,
                       proteins,
                       ixn_prior_width,
                       out_dir,
                       seed) {

  model_path = system.file("stan", "single_protein_by_ixn_width.stan",
                           package = 'bhm',
                           mustWork = TRUE)
  protein_model = cmdstan_model(stan_file = model_path, quiet = TRUE)
  # compile it once at first so the first iterations of the future_map don't try to all compile it together

  summaries = furrr::future_map(.x = proteins,
                                .f = fit_safely,
                                split_data_dir = split_data_dir,
                                save_fits = save_fits,
                                ixn_prior_width = ixn_prior_width,
                                # bh_input = bh_input,
                                proteins = proteins,
                                algorithm = algorithm,
                                out_dir = out_dir,
                                .options = furrr_options(seed = seed),
                                .progress = TRUE)

  res = data.table(protein = proteins,
                   summary = summaries)

  return(res)
}

get_summary = function(fit_summaries, bead_binders, concordance_scores,
                       weak_score_threshold = .5,
                       strong_score_threshold = 1,
                       weak_concordance = .75,
                       strong_concordance = .95,
                       out_dir,
                       verbose = TRUE) {

  worked = fit_summaries %>%
    dplyr::filter(purrr::map_lgl(summary, ~is.null(.x$error))) %>%
    dplyr::mutate(i = 1:n()) %>%
    as.data.table

  errored = fit_summaries %>%
    dplyr::filter(purrr::map_lgl(summary, ~!is.null(.x$error)))

  if (verbose) message(paste0("* ", nrow(worked), ' out of ', nrow(fit_summaries),
                              " (", round(100*nrow(worked)/nrow(fit_summaries), digits = 2), '%) protein fits ran without error.'))

  parameters = worked %>%
    dplyr::mutate(res = purrr::map(summary, 1)) %>%
    dplyr::select(-summary, -protein) %>%
    tidyr::unnest_legacy(res) %>%
    dplyr::select(-i, -s_i)

  ixn_scores =   parameters[grepl('prot_str', variable)] %>%
    dplyr::select(protein, strain, ixn_score = mean, se, p_type_s, `0.25%`:`99.75%`) %>%
    dplyr::arrange(desc(ixn_score))

  other_params = parameters[!grepl('prot_str', variable)] %>%
    dplyr::select(protein, strain, variable, posterior_mean = mean, se, p_type_s, `0.25%`:`99.75%`)

  with_concord = concordance_scores[ixn_scores, on = .(protein, strain)]

  with_concord$note = c('', "This protein showed normalized bead output above the specified bead binding enrichment threshold (default 1)")[(with_concord$protein %in% bead_binders$protein) + 1]

  with_concord = with_concord %>%
    dplyr::mutate(strong_hit = (ixn_score > strong_score_threshold) & !(`0.5%` < 0 & `99.5%` > 0) & (concordance > strong_concordance),
                  weak_hit = (ixn_score > weak_score_threshold) & !(`2.5%` < 0 & `97.5%` > 0) & (concordance > weak_concordance)) %>%
    dplyr::select(protein, strain, ixn_score, strong_hit, weak_hit, concordance, p_type_s, `0.25%`:`99.75%`, note)

  list(wc = with_concord,
       op = other_params)
}

#' Write out fit summary files with notes on bead binding
write_summary = function(fit_summary,
                         out_dir,
                         verbose = TRUE) {

  summs = fit_summary

  with_concord = summs$wc
  other_params = summs$op

  if (verbose & nrow(with_concord) > 1048576) {
    message("More scores than Excel can handle. Only writing out the first 10,000 to the Excel file. The full results will also be written to a tsv")
    openxlsx::write.xlsx(x = list('ixn_scores' = with_concord[1:1e4]),
                         file = paste0(out_dir, 'ixn_scores.xlsx'))

    data.table::fwrite(with_concord,
                       file = paste0(out_dir, "ixn_scores.tsv"),
                       sep = '\t')
    data.table::fwrite(other_params,
                       file = paste0(out_dir, "other_parameters.tsv"),
                       sep = '\t')

  } else {
    openxlsx::write.xlsx(x = list('ixn_estimates' = with_concord,
                                  'other_parameters' = other_params),
                         file = paste0(out_dir, 'ixn_scores.xlsx'))
    data.table::fwrite(with_concord,
                       file = paste0(out_dir, "ixn_scores.tsv"),
                       sep = '\t')
    data.table::fwrite(other_params,
                       file = paste0(out_dir, "other_parameters.tsv"),
                       sep = '\t')
  }

  NULL
}

end_with_slash = function(dir_path){
  if (!grepl('\\/$', dir_path)){
    dir_path = paste0(dir_path, '/')
  }
  return(dir_path)
}

check_out_dir = function(out_dir) {

  out_dir = end_with_slash(out_dir)

  if (!is.null(out_dir) & !dir.exists(out_dir)) {
    if (dir.exists(out_dir)){
      stop('output directory already exists')
    } else {
      #TODO add a check for cache_dir
      dir.create(out_dir)
      # split_data_dir = file.path(gsub('\\/$', '', out_dir ), 'splits/')
      # dir.create(split_data_dir)
    }
  }

  return(out_dir)
}

#' Run the basehit model one protein at a time
#'
#' @inheritParams run_model
#' @param bead_binding_threshold proteins with enrichment in the beads above this threshold get noted in the output
#' @details The count file should have the first row specifying proteins, the
#'   second specifying barcodes, and all others after that specifying the output
#'   counts for each strain counts for each barcode (i.e. wide format, strain x
#'   barcode)
#'
#'   Implemented with furrr, so run a call to plan() that's appropriate for your
#'   system in order to parallelize
#' @export
model_proteins_separately = function(count_path,
                                     out_dir = 'outputs/bh_out/',
                                     cache_dir = 'outputs/bh_cache/',
                                     split_data_dir = NULL,
                                     ixn_prior_width = .15,
                                     algorithm = 'variational',
                                     save_split = TRUE,
                                     save_fits = FALSE,
                                     bead_binding_threshold = 1,
                                     save_bead_binders = TRUE,
                                     min_n_nz = 3,
                                     min_frac_nz = .5,
                                     verbose = TRUE,
                                     seed = 1234) {

  if (verbose) message("Step 1/8 Checking output directory...")
  out_dir = check_out_dir(out_dir)

  if (verbose) message("Step 2/8 Reading data...")
  count_list = read_multirun(count_path,
                             cache_dir = cache_dir) # a list of four 5 dataframes: pre, WR pre, beads, and WR output. Also the prot_to_bc header lines

  if (verbose) message("Step 3/8 Filtering data...")
  filtered_data = filter_multirun(count_list,
                                  cache_dir = cache_dir,
                                  min_n_nz = min_n_nz,
                                  min_frac_nz = min_frac_nz,
                                  save_outputs = TRUE,
                                  verbose = verbose)

  if (verbose) message("Step 4/8 Writing out data split by protein...")
  if (missing(split_data_dir)){
    message("* Using split data directory under cache directory...")
    split_data_dir = paste0(cache_dir, 'splits/')
  }
  if (!dir.exists(split_data_dir)) {
    message("* Split data directory doesn't exist. Creating it...")
    dir.create(split_data_dir)
  }
  write_splits(filtered_data,
               split_data_dir)

  if (verbose) message("Step 5/8 Identifying bead binding proteins...")
  bead_binding = identify_bead_binders(wr_pre = count_list$wr_pre,
                                       prot_to_bc = count_list$prot_to_bc,
                                       bh_beads = count_list$beads,
                                       bh_output = count_list$wr_output,
                                       binding_threshold = bead_binding_threshold,
                                       save_bead_binders = save_bead_binders,
                                       out_dir = out_dir,
                                       verbose = verbose)

  if (verbose) message("Step 6/8 Computing concordance scores...")
  concordance = get_concordance(filtered_data$bh_input,
                                drop_multirun_strains = TRUE)

  if (verbose) message("Step 7/8 Fitting model by protein...")
  model_fits = fit_models(algorithm = algorithm,
                          split_data_dir = split_data_dir,
                          bh_input = filtered_data$bh_input,
                          proteins = unique(filtered_data$bh_input$protein),
                          save_fits = save_fits,
                          ixn_prior_width = ixn_prior_width,
                          out_dir = out_dir,
                          seed = seed)

  if (verbose) message("Step 8/8 Writing result tables...")
  save(model_fits, bead_binding, concordance,
       file = paste0(out_dir, 'all_outputs.RData'))

  fit_summary = get_summary(model_fits, bead_binding, concordance)

  write_summary(fit_summary, out_dir, verbose)

  if (!save_split) {
    if (verbose) message("Removing split data directory...")
    unlink(split_data_dir)
  }

  fit_summary$wc
}
