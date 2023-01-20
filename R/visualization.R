tick = function(x) {
  paste0("`", x, "`")
}

#' @title Make a hit calling plot
#' @description Show both the posterior intervals and placement in the
#'   concordance distribution
#' @details The concordance distribution shows the concordance values for all
#'   interactions that exclude 0 as a credible value by 95% interval
#' @return a list of three plots
#' @export
hit_calling_plot = function(protein_, strain_,
                            bead_binding,
                            concordances,
                            fit_summary,
                            weak_score_threshold = .5,
                            strong_score_threshold = 1,
                            weak_concordance = .75,
                            strong_concordance = .95,
                            tags = TRUE,
                            intervals = c(.95, .99)) {

  int_cols = c(paste0(100 * sort(c((1 - intervals) / 2,
                    (1 - intervals) / 2 + intervals)), "%"), "ixn_score")

  ixn_summary = fit_summary[protein == protein_ & strain == strain_]
  interval_df = ixn_summary[, ..int_cols]

  interval_plot = ggplot(interval_df,
         aes(x = interaction_score)) +
    geom_vline(xintercept = c(weak_score_threshold, strong_score_threshold),
               color = 'grey60') +
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
    theme(axis.text.y = element_blank(),
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

#' Show hit-calling metrics for multiple interactions
#' @description Plot both the posterior intervals of interaction score and
#'   position in the concordance distribution
#' @param ixns a character string giving multiple interactions as protein:strain
#' @param bead_binding bead-binding enrichment dataframe (saved as part of
#'   all_outputs.RData)
#' @param concordances concordance dataframe (saved as part of
#'   all_outputs.RData)
#' @param fit_summary the interaction score summary dataframe (result from
#'   basehitmodel::model_proteins_separately())
#' @param name_df a dataframe giving the strain identifiers and taxonomic names in two columns: "strain" and "strain_name"
#' @param weak_score_threshold score threshold for calling weak hits
#' @param strong_score_threshold score threshold for calling strong hits
#' @param weak_concordance concordance threshold for weak hits
#' @param strong_concordance concordance threshold for calling strong hits
#' @param global_entropy logical indicating whether to show the global entropy
#'   histogram of all profiled interactions (TRUE) or only interactions with
#'   scores that exclude 0 by the narrower interval (default FALSE)
#' @param tags logical: show A & B on the two panels?
#' @param intervals amount of probability mass to include in the two intervals
#' @return a list of three plots: the interval plot, the entropy plot, and the two stuck together
#' @export
multi_hit_calling_plot = function(ixns = c("CD55:AIEC", "CEACAM1:AIEC", "CD7:AIEC", "CEACAM1:NWP04"),
                                  bead_binding,
                                  concordances,
                                  fit_summary,
                                  name_df = NULL,
                                  weak_score_threshold = .5,
                                  strong_score_threshold = 1,
                                  weak_concordance = .75,
                                  strong_concordance = .95,
                                  global_entropy = FALSE,
                                  tags = TRUE,
                                  intervals = c(.95, .99)) {

  int_cols = c(paste0(100 * sort(c((1 - intervals) / 2,
                                   (1 - intervals) / 2 + intervals)), "%"), "ixn_score", 'protein', 'strain', "concordance")
  if (!missing(name_df)) {
    int_cols = c(int_cols, "strain_name")
  }
  ixn_df = data.table(ixn = ixns) %>%
    tidyr::separate(ixn, into = c('protein_', 'strain_'), remove = FALSE, sep = ":")

  if (!missing(name_df)) {
    name_df = as.data.table(name_df)
    fit_summary = fit_summary[name_df, on = 'strain'][, strain_name := paste(strain, '_', strain_name, sep = '')]
    interval_df = fit_summary[protein %in% ixn_df$protein_ | strain %in% ixn_df$strain_] %>%
      mutate(ixn = paste(protein, strain, sep = ":")) %>%
      .[ixn %in% ixn_df$ixn] %>%
      .[, ..int_cols] %>%
      mutate(interaction = forcats::fct_reorder(factor(paste(protein, gsub(' ', '_', strain_name), sep = ":")),
                                                .x = ixn_score))
  } else {
    interval_df = fit_summary[protein %in% ixn_df$protein_ | strain %in% ixn_df$strain_] %>%
      mutate(ixn = paste(protein, strain, sep = ":")) %>%
      .[ixn %in% ixn_df$ixn] %>%
      .[, ..int_cols] %>%
      mutate(interaction = forcats::fct_reorder(factor(paste(protein, gsub(' ', '_', strain), sep = ":")),
                                                .x = ixn_score))
  }



  interval_plot = ggplot(interval_df,
                         aes(x = ixn_score)) +
    geom_vline(xintercept = c(weak_score_threshold, strong_score_threshold),
               color = 'grey60') +
    geom_segment(lwd = 1,
                 color = 'black',
                 aes_string(x = tick(int_cols[1]),
                            xend = tick(int_cols[4]),
                            y = "interaction",
                            yend = "interaction")) +
    geom_segment(lwd = 3,
                 aes_string(x = tick(int_cols[2]),
                            xend = tick(int_cols[3]),
                            y = "interaction",
                            yend = "interaction",
                            color = "interaction")) +
    geom_point(aes(y = interaction,
                   x = ixn_score)) +
    geom_vline(xintercept = 0,
               lty = 2,
               color = "grey") +
    theme_light() +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank()) +
    labs(x = "interaction score") +
    guides(color = guide_legend(reverse = TRUE))+
    scale_color_brewer(palette = "Set1")

  ixn_summary = interval_df


  if (global_entropy){
    entropy_input = fit_summary
  } else {
    # V These are interactions that would be weak hits by the posterior interval criterion alone
    entropy_input = fit_summary[!(fit_summary[,which(names(fit_summary) == int_cols[2]), with = FALSE][[1]] < 0 &
                                    fit_summary[,which(names(fit_summary) == int_cols[3]), with = FALSE][[1]] > 0)]
  }

  entropy_plot = entropy_input %>%
    ggplot(aes(concordance)) +
    geom_histogram(boundary = 0, bins = 100) +
    geom_vline(lty = 2,
               xintercept = c(weak_concordance, strong_concordance)) +
    geom_vline(data = interval_df,
               aes(xintercept = concordance,
                   color = interaction),
               key_glyph = 'rect') +
    theme_light() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(y = NULL) +
    guides(color = guide_legend(reverse = TRUE)) +
    scale_color_brewer(palette = "Set1")

  layout_str = "
  12
  12"

  if (tags) {
    tl = "A"
  } else {
    tl = NULL
  }

  both_plot = (interval_plot+ theme(legend.position = 'none')) + (entropy_plot ) +
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

#' Plot barcode level model inputs
#'
#' @param bh_input a data frame of filtered model inputs
#' @param barcodes a character vector of barcodes to filter to
#' @param proteins a character vector of proteins to filter to
#' @param strains  a character vector of strains to filter to
#' @param force if TRUE, override the check preventing gigantic plots
#' @param log10_counts if TRUE, log10 the values in both panels
#' @details Grey cells in the top panel correspond to zeros in the original input (if log10_counts =
#'   TRUE). Empty cells (i.e. where you can see the underlying grid lines) are non-present in the
#'   input, likely due to being filtered out.
#' @examples
#' \dontrun{
#' bh_input = data.table::fread("~/Desktop/tmp/cache/bh_input.tsv.gz")
#' plot_model_inputs(bh_input,
#'                   strains = c("AB1", "AB10", "AB12"),
#'                   proteins = c("LSAMP", "THSD1_Epitope-1", "LRTM1"))
#' }
#' @export
plot_model_inputs = function(bh_input,
                             barcodes = NULL,
                             proteins = NULL,
                             strains  = NULL,
                             force = FALSE,
                             log10_counts = TRUE) {

  if (!data.table::is.data.table(bh_input)) bh_input = data.table::as.data.table(bh_input)
  if (!is.null(barcodes)) bh_input = bh_input[barcode %in% barcodes]
  if (!is.null(proteins)) bh_input = bh_input[protein %in% proteins]
  if (!is.null( strains)) bh_input = bh_input[strain  %in%  strains]

  if (nrow(bh_input) > 5000 && !force) stop("Too many count observations to plot. Set force = TRUE to override.")

  if (log10_counts) trans_fun = log10 else trans_fun = identity

  plot_input = bh_input |>
    dplyr::mutate(normalized_output = trans_fun(count / pre_count),
                  pre_count = trans_fun(pre_count))

  pre_plot = plot_input |>
    ggplot(aes(barcode, sample_id)) +
    geom_tile(aes(fill  = pre_count,
                  color = pre_count)) +
    facet_wrap(vars(protein), nrow = 1,
               scales = 'free_x') +
    scale_fill_viridis_c(option = "E") +
    scale_color_viridis_c(option = "E") +
    coord_cartesian(expand = FALSE) +
    theme_bw() +
    theme(axis.text.x  = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_blank(),
          axis.title.y = element_blank())

  out_plot = plot_input |>
    ggplot(aes(barcode, sample_id)) +
    geom_tile(aes(fill = normalized_output,
                  color = normalized_output)) +
    facet_wrap(vars(protein), nrow = 1,
               scales = "free_x") +
    scale_fill_viridis_c(option = "D") +
    scale_color_viridis_c(option = "D") +
    coord_cartesian(expand = FALSE) +
    theme_bw() +
    theme(axis.text.x = element_text(vjust  = .5,
                                     angle  = 90,
                                     family = "mono"),
          axis.title.x = element_blank(),
          panel.border = element_blank(),
          axis.title.y = element_blank())

  out_plot / pre_plot + patchwork::plot_layout(heights = c(6, 1), guides = 'collect')
}

#' Venn diagram of hit-calling thresholds
#'
#' @description This function draws a Venn diagram of the counts of interactions that pass the three
#'   hit-calling criteria for the provided thresholds.
#' @export
hit_calling_venn = function(score_df,
                            score_threshold = .5,
                            concordance_threshold = .75,
                            interval_width = .95) {
  interval_outer = 1 - interval_width

  interval_ends = (c(interval_outer/2, 1-interval_outer/2) * 100) |>
    as.character() |>
    paste('%', sep = '')

  if (!all(interval_ends %in% names(score_df))) {
    stop("score_df doesn't contain the interval endpoints for the given interval_width")
  }

  N = nrow(score_df)

  circle = function(cx,cy) {
    data.table(x = cx + cos(seq(0, 2*pi, length.out = 300)),
               y = cy + sin(seq(0, 2*pi, length.out = 300)))
  }

  a = 1
  h = sqrt(3)/2*a

  circle_df = data.table(cx = c(    0, -a/2,  a/2),
                         cy = c(2*h/3, -h/3, -h/3),
                         id = 1:3) |>
    dplyr::mutate(point_df = purrr::map2(cx,cy,
                                         circle)) |>
    tidyr::unnest_legacy()

  call_df = score_df |>
    dplyr::transmute(passes_interval = !(score_df[[interval_ends[1]]] < 0 & score_df[[interval_ends[2]]] > 0),
                     passes_size = ixn_score > score_threshold,
                     passes_concordance = concordance > concordance_threshold)

  call_dt = table(call_df) |>
    as.data.table() |>
    arrange(passes_interval, passes_size, passes_concordance) |>
    mutate(x = c(-1.4, 1, -1, 0,
                 0, .6, -.6, 0),
           y = c(1.2, -.4, -.4, -.7,
                 1, .4, .4, 0))

  lab_df = data.table(label = c("passes\ninterval",
                                "passes\nsize",
                                "passes\nconcordance"),
                      x = c(1.3, -1.9, 1.9),
                      y = c(1.3, -.9, -.9))

  ggplot(lab_df, aes(x,y)) +
    geom_path(data = circle_df,
              aes(group = id)) +
    geom_text(data = call_dt,
              aes(label = N)) +
    geom_text(data = lab_df,
              aes(label = label)) +
    coord_equal() +
    theme_void() +
    xlim(-2.5,2.5) +
    ylim(-1.5, 2)
}
