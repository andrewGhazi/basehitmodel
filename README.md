
## Installation

This package requires
[cmdstanr](https://mc-stan.org/cmdstanr/articles/cmdstanr.html), and so
can’t go on CRAN (nor Bioc I think). Install `cmdstanr` first, then use
it to install `CmdStan`:

``` r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/",
                                       getOption("repos")))

library(cmdstanr)
check_cmdstan_toolchain()
install_cmdstan(cores = 2)
```

Note that C++ compilers are also required (`check_cmdstan_toolchain()`
will check for this). Installing `CmdStan` can take around ~5 minutes,
depending on how fast your system can compile it. You can set
`Sys.setenv(MAKEFLAGS = "-j4")` before the install command to make it
compile in parallel if you wish.

Then use:

``` r
remotes::install_github("andrewGhazi/basehitmodel")
```

to install the remaining dependencies (listed in the DESCRIPTION file)
and `basehitmodel` itself. You can add `Ncpus=3` to that command to make
the dependency installation run in parallel, though you will probably
need to restart your R session if you set the parallel makeflag above.

Installing the R packages will take a ~20 seconds if your OS gets
precompiled binaries from CRAN, otherwise it will be ~2-5 minutes on a
normal laptop.

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
                          split_data_dir = "/path/to/outputs/splits/",
                          ixn_prior_width = 0.15,
                          algorithm = "variational",
                          iter_sampling = 2000,
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

## Software Versions

This package has been tested on R 4.1.0 and the following package
versions (though anything from ~2020 forwards *should* work):

| pkg        | ver    |
|------------|--------|
| dplyr      | 1.0.7  |
| data.table | 1.14.2 |
| magrittr   | 2.0.2  |
| openxlsx   | 4.2.4  |
| tibble     | 3.1.6  |
| tidyr      | 1.2.0  |
| cmdstanr   | 0.5.1  |
| readxl     | 1.3.1  |
| furrr      | 0.2.3  |
| stringr    | 1.4.0  |
| ggplot2    | 3.3.5  |
| patchwork  | 1.1.1  |
| forcats    | 0.5.1  |
| purrr      | 0.3.4  |
| posterior  | 1.2.0  |
| rlang      | 1.0.0  |
| stats      | 4.1.0  |
| utils      | 4.1.0  |
| readr      | 2.1.2  |
| progressr  | 0.10.0 |

Any OS that can install R and Stan should be sufficient. We have tested
the package most thoroughly on CentOS Linux release 7.6.1810 (Core), as
well as various recent versions of macOS and Windows 10.
