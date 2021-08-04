tick = function(x){
  paste0("`", x, "`")
}

#' @title Make a hit calling plot
#' @description Show both the posterior intervals and placement in the concordance distribution
#' @return a list of three plots
hit_calling_plot = function(protein_, strain_,
                            bead_binding,
                            concordances,
                            fit_summary, # TODO make it not need redundant arguments. this has the same information as the previous arg
                            weak_score_threshold = .5,
                            strong_score_threshold = 1,
                            weak_concordance = .75,
                            strong_concordance = .95,
                            tags = TRUE,
                            intervals = c(.95, .99),
                            algorithm = 'variational',
                            seed = 1234) {

  int_cols = c(paste0(100 * sort(c((1 - intervals) / 2,
                    (1 - intervals) / 2 + intervals)), "%"), "ixn_score")

  ixn_summary = fit_summary[protein == protein_ & strain == strain_]
  interval_df = ixn_summary[, ..int_cols]

  interval_plot = ggplot(interval_df,
         aes(x = interaction_score)) +
    geom_segment(lwd = 1,
                 color = 'black',
                 aes_string(x = tick(int_cols[1]),
                            xend = tick(int_cols[4]),
                            y = "0",
                            yend = "0")) +
    geom_segment(lwd = 3,
                 color = 'grey55',
                 aes_string(x = tick(int_cols[2]),
                            xend = tick(int_cols[3]),
                            y = "0",
                            yend = "0")) +
    geom_point(aes(y = 0,
                   x = ixn_score)) +
    geom_vline(xintercept = 0,
               lty = 2,
               color = "grey") +
    theme_light() +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank()) +
    labs(x = paste0(protein_, ":", strain_,
                    " interaction score"))

  ixn_summary = fit_summary[protein == protein_ & strain == strain_]

  # V These are interactions that would be weak hits by the posterior interval criterion alone
  entropy_plot = fit_summary[!(fit_summary[,which(names(fit_summary) == int_cols[2]), with = FALSE][[1]] < 0 &
                                 fit_summary[,which(names(fit_summary) == int_cols[3]), with = FALSE][[1]] > 0)] %>%
    ggplot(aes(concordance)) +
    geom_histogram(boundary = 0, bins = 100) +
    geom_vline(lty = 2,
               xintercept = c(weak_concordance, strong_concordance)) +
    geom_vline(color = 'red',
               xintercept = c(ixn_summary$concordance)) +
    theme_light() +
    theme(text = element_text(size = 7),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(y = NULL)

  layout_str = "
  12
  #2"

  if (tags) {
    tl = "A"
  } else {
    tl = NULL
  }

  both_plot = (interval_plot) + (entropy_plot ) +
    plot_layout(design = layout_str,
                tag_level = "keep") +
    plot_annotation(tag_levels = tl) &
    theme(plot.tag = element_text(size = 12))

  print(both_plot)
  plot_list = list(a = interval_plot,
                   b = entropy_plot,
                   ab = both_plot)
  return(plot_list)
}
