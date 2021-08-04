#' @title Make a hit calling plot
#' @description Show both the posterior intervals and placement in the concordance distribution
#' @return a list of three plots
hit_calling_plot = function(protein_, strain_,
                            bead_binding,
                            concordances,
                            model_outputs,
                            bh_input,
                            bh_eta,
                            tags = c("A", 'B"'),
                            intervals = c(.95, .99),
                            algorithm = 'variational',
                            seed = 1234) {

  int_cols = paste0(100 * sort(c((1 - intervals) / 2,
                    (1 - intervals) / 2 + intervals)), "%")

  interval_df = model_fits[protein == protein_]$summary[[1]]$result[strain == strain_ & grepl('prot_st', variable)][, ..int_cols] %>%
    melt(variable.name = 'q',
         value.name = 'interaction_score')


  p1 = bh_input[protein == protein_]
  p1_eta = bh_eta[protein == protein_]
  prot_name = str_replace_all(pattern = '-|[:space:]|/|:', replacement = '_', string = p1$protein[1])
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
  p1_eta$ps = droplevels(p1_eta$ps)
  p1_eta$s_i = as.numeric(p1_eta$strain)
  p1_eta$ps_i = as.numeric(p1_eta$ps)

  message("Rerunning sampler...")

  data_list = list(n_strain = nlevels(p1$strain),
                   n_bc = n_distinct(p1$barcode),
                   n_ps = nlevels(p1$strain),
                   N = nrow(p1),
                   strain_i = p1$s_i,
                   n_eta = nrow(p1_eta),
                   eta_ps_i = p1_eta$ps_i,
                   eta_presum = log(p1_eta$sd) + log(p1_eta$pre_count),
                   eta_i = p1$eta_i,
                   count_obs = p1$count)


  protein_fit = protein_model$variational(data = data_list,
                                          seed = 1234,
                                          output_dir = 'outputs/IFNA17',
                                          output_samples = 5000)
  if (algorithm == 'variational') {
    # TODO expose more of the stan parameters to the user
    protein_fit = protein_model$variational(data = data_list,
                                            output_samples = 5000)
  } else {
    protein_fit = protein_model$sample(data = data_list,
                                       iter_sampling = 5000)
  }

  protein_summary = protein_fit$summary('mean' = mean,
                                        'se' = posterior::mcse_mean,
                                        'med' = median,
                                        'qs' = ~quantile(.x, probs = sort(c((1 - intervals) / 2,
                                                                            (1 - intervals) / 2 + intervals))),
                                        'conv' = posterior::default_convergence_measures(),
                                        'p_type_s' = p_type_s)

  int_df = as.data.table(protein_summary)[strain == strain_]
  interval_plot = ggplot(post_samples,
                         aes(x = ixn_score)) +
    geom_histogram(aes(y = ..density..)) +
    geom_vline(xintercept = c(.5, 1),
               color = 'grey60') +
    geom_vline(lty = 2,
               color = 'grey20',
               xintercept = 0) +
    geom_segment(data = int_df,
                 aes(x = x, xend = xend, y = y, yend = yend,
                     color = colors),
                 lwd = 1.5, show.legend = FALSE) +
    labs(y = 'posterior density',
         x = paste0(protein, ":", strain, " interaction score"),
         color = NULL) +
    scale_color_brewer(palette = 'Set1') +
    expand_limits(x = 0) +
    theme(text = element_text(size = 7)) +
    theme_light()

  both_plot = (interval_plot + labs(tag = tags[1])) + (entropy_plot + labs(tag = tags[2]))

  print(both_plot)
  plot_list = list(a = interval_plot,
                   b = entropy_plot,
                   ab = both_plot)
  return(plot_list)
}
