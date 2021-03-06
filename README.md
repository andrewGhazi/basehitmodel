
## Installation

This package requires
[cmdstanr](https://mc-stan.org/cmdstanr/articles/cmdstanr.html), and so
can’t go on CRAN (nor Bioc I think). Install `cmdstanr` first, then use
`devtools::install_github("andrewGhazi/basehitmodel")`.

## Example

Point `model_proteins_separately` to a directory of mapped basehit
counts and it will fit the ZINB model per-protein and save the results
in the specified output directory. An optional cache directory can be
used to save the pre-processing steps so that future model runs with
different statistical parameters (e.g. priors) faster. Help
documentation on the other arguments can be displayed with
`?model_proteins_separately`.

``` r
model_proteins_separately(count_path = "~/basehit/data/BASEHIT_rescreen/",
                          out_dir = "/path/to/out_dir/",
                          cache_dir = "/path/to/outputs/bh_cache/",
                          split_data_dir = NULL,
                          ixn_prior_width = 0.15,
                          algorithm = "variational",
                          iter_sampling = 5000,
                          iter_warmup = 1000,
                          save_split = TRUE,
                          save_fits = FALSE,
                          save_summaries = TRUE,
                          bead_binding_threshold = 1,
                          save_bead_binders = TRUE,
                          min_n_nz = 3,
                          min_frac_nz = 0.5,
                          verbose = TRUE,
                          seed = 1234)
```
