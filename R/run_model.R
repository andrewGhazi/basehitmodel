#' Run the basehit model
#'
#' @param count_path The path to a Mapped_bcs.csv file
#' @param out_dir path to output, must end with a slash
#' @param out_name what you'd like the summary xlsx file to be called
#' @param model_path path to the stan model file to use
#' @param vb_threads The number of threads to use for the vb fit
#' @details The count file should have the first row specifying proteins, the second specifying
#'   barcodes, and all others after that specifying the output counts for each strain counts for
#'   each barcode (i.e. wide format, strain x barcode)
#'
#'   Changing the stan model will almost certainly break things if it has a different
#'   parameterization as the default.
run_model = function(count_path,
                     out_dir = 'outputs/bh_out/',
                     out_name = 'results.xlsx',
                     model_path = 'src/stan_files/basehit_reduce_redundancy_rs_nonsquare.stan',
                     vb_threads = 4) {

  if (dir.exists(out_dir)){
    stop("Please provide an output directory that doesn't already exist.")
  }
  bh_header = fread(count_path,
                    header = FALSE,
                    nrows = 2) %>%
    data.table::transpose() %>%
    .[-1,] %>%
    set_colnames(c('protein', 'barcode'))

  bh_data = fread(count_path,
                  skip = 1) %>%
    rename(sample = Barcode) %>%
    melt(id.vars = c('sample'),
         variable.name = 'barcode',
         value.name = 'count') %>%
    left_join(bh_header, by = 'barcode') %>%
    .[order(sample, protein)]


  bh_output = bh_data[!grepl('Beads|Pre', sample)][, ps := paste(sample, protein, sep = ':')][]
  bh_pre = bh_data[grepl('Pre', sample)][, ps := paste(sample, protein, sep = ':')][]
  bh_beads = bh_data[grepl('beads', tolower(sample))][, ps := paste(sample, protein, sep = ':')][]

  # mean(bh_output$count == 0) # 89.6%

  # WR = "well-represented"
  bh_wr_input = bh_pre[count>4] %>%
    rename(pre_count = count)

  # n_nz = "number non-zero"
  bh_n_nz = bh_output[barcode %in% bh_wr_input$barcode] %>%
    .[, .(n_nz = sum(count != 0),
          ps = paste(sample, protein, sep = ':')), by = .(sample, protein)]

  bh_wr = bh_n_nz[n_nz >= 3] # well represented
  bh_wr_output = bh_output[ps %in% bh_wr$ps & (barcode %in% bh_wr_input$barcode)]

  bh_input = bh_wr_input[,.(barcode, pre_count)] %>%
    .[bh_wr_output, on = 'barcode'] %>%
    .[,sd := sum(count) / 1e5,by = sample] %>%
    .[,`:=`(barcode = factor(barcode),
            protein = factor(protein),
            sample = factor(sample),
            ps = factor(ps))] %>%
    .[,`:=`(s_i = as.integer(sample),
            p_i = as.integer(protein),
            ps_i = as.integer(ps))] %>%
    .[,`:=`(eta_i = as.integer(factor(paste(pre_count, sample, protein, sep = ':'))))] %>%
    .[]

  # These are the unique levels of eta in the negative binomial regression. We pass them to the stan
  # model so that each only has tobe computed once as a way to reduce redundancy in the computation.
  bh_eta = bh_input[,.(pre_count, sample, sd, protein, p_i, s_i, ps_i, eta_i)] %>%
    unique %>%
    .[]

  #### Fit the model ----
  mod = cmdstan_model(model_path, cpp_options = list(stan_threads = TRUE))

  data_list = list(n_prot = nlevels(bh_input$protein),
                   n_strain = nlevels(bh_input$sample),
                   n_bc = nlevels(bh_input$barcode),
                   n_ps = bh_input$ps_i %>% max,
                   n_eta = nrow(bh_eta),
                   N = nrow(bh_input),
                   strain_i = bh_input$s_i,
                   eta_p_i = bh_eta$p_i,
                   eta_ps_i = bh_eta$ps_i,
                   eta_presum = log(bh_eta$sd) + log(bh_eta$pre_count),
                   eta_i = bh_input$eta_i,
                   count_obs = bh_input$count)

  refresh_freq = 1

  # out_dir = paste0(getwd(), '/outputs/bh_out')
  if (!dir.exists(out_dir)) dir.create(out_dir)

  bh_fit = mod$variational(data = data_list,
                           seed = 123,
                           threads = vb_threads,
                           algorithm = 'meanfield',
                           output_samples = 4000,
                           output_dir = out_dir,
                           refresh = refresh_freq)

  bh_fit$save_object(paste0(out_dir, '/bh_fit.RDS'))

  bh_summary = bh_fit$summary('mean' = mean, 'se' = posterior::mcse_mean, 'med' = median, 'qs' = ~quantile(.x, probs = c(.0025, .005, .025, .05, .1, .25,
                                                                                                                         .75, .9, .95, .975, .995, .9975)), 'conv' = posterior::default_convergence_measures())

  save(bh_summary,
       file = paste0(out_dir, '/bh_summary.RData'))

  #### compute intervals + concordance ----

  i_to_ps = bh_eta[,.(ps_i, sample, protein)] %>%
    unique %>%
    mutate(ps_i = as.character(ps_i))

  bh_groups = bh_input[,.(sample)] %>%
    unique %>%
    separate(sample,
             into = c('group', 'replicate'),
             sep = '-(?=[AB]$)',
             remove = FALSE)

  summary_df = bh_summary %>%
    filter(grepl('strain', variable)) %>%
    mutate(ps_i = str_extract(variable, '[0-9]+')) %>%
    left_join(i_to_ps, by = 'ps_i') %>%
    left_join(bh_groups %>% select(sample, group),
              by = 'sample') %>%
    select(group,  sample, protein, everything()) %>%
    filter(!is.na(group)) %>% #some samples not in the metadata :/
    arrange(group, protein, sample)

  check_hit = function(quantiles, levels, level = "95%"){

    switch(level,
           "50%" = !(quantiles[levels == "25%"] < 0 & quantiles[levels == "75%"] > 0),
           "80%" = !(quantiles[levels == "10%"] < 0 & quantiles[levels == "90%"] > 0),
           "90%" = !(quantiles[levels == "5%"] < 0 & quantiles[levels == "95%"] > 0),
           "95%" = !(quantiles[levels == "2.5%"] < 0 & quantiles[levels == "97.5%"] > 0),
           "99%" = !(quantiles[levels == "0.5%"] < 0 & quantiles[levels == "99.5%"] > 0),
           "99.5%" = !(quantiles[levels == "0.25%"] < 0 & quantiles[levels == "99.75%"] > 0))
  }

  hit_checks = summary_df %>%
    select(group, sample, protein, mean, `0.25%`:`99.75%`) %>%
    pivot_longer(`0.25%`:`99.75%`, names_to = 'quantile_level', values_to = 'quantile') %>%
    arrange(group, sample, protein) %>%
    as.data.table %>%
    .[,.(`50%_hit` = check_hit(quantile,quantile_level, "50%"),
         `80%_hit` = check_hit(quantile,quantile_level, "80%"),
         `90%_hit` = check_hit(quantile,quantile_level, "90%"),
         `95%_hit` = check_hit(quantile,quantile_level, "95%"),
         `99%_hit` = check_hit(quantile,quantile_level, "99%"),
         `99.5%_hit` = check_hit(quantile,quantile_level, "99.5%")), by = .(group, protein, sample)] %>%
    .[order(group, protein, sample)]

  hit_counts = hit_checks[, n := .N, by = .(group, protein)] %>%
    as_tibble %>%
    pivot_longer(`50%_hit`:`99.5%_hit`,
                 names_to = 'level',
                 values_to = 'is_hit') %>%
    group_by(level) %>%
    summarise(n_hits = sum(is_hit))

  concord = hit_checks[, n := .N, by = .(group, protein)] %>%
    as_tibble %>%
    pivot_longer(`50%_hit`:`99.5%_hit`,
                 names_to = 'level',
                 values_to = 'is_hit') %>%
    group_by(group, protein, level) %>%
    summarise(n = n[1],
              n_agree_hit = sum(is_hit),
              p_agree_hit = mean(is_hit)) %>%
    ungroup %>%
    arrange(-n_agree_hit)

  concord_wide = concord %>% select(-p_agree_hit) %>%
    pivot_wider(names_from = 'level',
                values_from = 'n_agree_hit') %>%
    arrange(-`99.5%_hit`)

  concord_p = concord %>% select(-n_agree_hit) %>%
    pivot_wider(names_from = 'level',
                values_from = 'p_agree_hit')

  ordered_concord_p = concord_wide %>%
    select(1:2) %>%
    left_join(concord_p, by = c('group', 'protein'))

  zf = bh_summary %>%
    filter(grepl('theta', variable)) %>%
    mutate(s_i = str_extract(variable, '[0-9]+')) %>%
    left_join(unique(bh_eta[,.(sample, s_i = as.character(s_i))]), by = 's_i') %>%
    select(sample, everything()) %>%
    select(-variable)

  ixns =  bh_summary %>%
    filter(grepl('prot_strain', variable)) %>%
    mutate(ps_i = str_extract(variable, '[0-9]+')) %>%
    left_join(unique(bh_eta[,.(protein, sample, ps_i = as.character(ps_i))]), by = 'ps_i') %>%
    select(protein, sample, everything()) %>%
    select(-variable, -ps_i) %>%
    arrange(protein, sample)

  prots = bh_summary %>%
    filter(grepl('protein', variable)) %>%
    mutate(p_i = str_extract(variable, '[0-9]+')) %>%
    left_join(unique(bh_eta[,.(protein, p_i = as.character(p_i))]), by = 'p_i') %>%
    select(protein, everything()) %>%
    select(-variable)

  op = bh_summary %>%
    filter(!grepl('protein|strain|theta', variable))

  openxlsx::write.xlsx(x = list('ixn_estimates' = ixns,
                                'protein_estimates' = prots,
                                'zero_fractions' = zf,
                                'other_parameters' = op,
                                'hit_checks' = hit_checks,
                                'hit_counts' = hit_counts,
                                'concordance_count' = concord_wide,
                                'concordance_fraction' = ordered_concord_p),
                       file = paste0(out_dir, out_name))

  return(bh_summary)
}
