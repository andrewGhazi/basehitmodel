% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{multi_hit_calling_plot}
\alias{multi_hit_calling_plot}
\title{Show hit-calling metrics for multiple interactions}
\usage{
multi_hit_calling_plot(
  ixns = c("CD55:AIEC", "CEACAM1:AIEC", "CD7:AIEC", "CEACAM1:NWP04"),
  bead_binding,
  concordances,
  fit_summary,
  name_df = NULL,
  weak_score_threshold = 0.5,
  strong_score_threshold = 1,
  weak_concordance = 0.75,
  strong_concordance = 0.95,
  global_entropy = FALSE,
  tags = TRUE,
  intervals = c(0.95, 0.99)
)
}
\arguments{
\item{ixns}{a character string giving multiple interactions as protein:strain}

\item{bead_binding}{bead-binding enrichment dataframe (saved as part of
all_outputs.RData)}

\item{concordances}{concordance dataframe (saved as part of
all_outputs.RData)}

\item{fit_summary}{the interaction score summary dataframe (result from
basehitmodel::model_proteins_separately())}

\item{name_df}{a dataframe giving the strain identifiers and taxonomic names in two columns: "strain" and "strain_name"}

\item{weak_score_threshold}{score threshold for calling weak hits}

\item{strong_score_threshold}{score threshold for calling strong hits}

\item{weak_concordance}{concordance threshold for weak hits}

\item{strong_concordance}{concordance threshold for calling strong hits}

\item{global_entropy}{logical indicating whether to show the global entropy
histogram of all profiled interactions (TRUE) or only interactions with
scores that exclude 0 by the narrower interval (default FALSE)}

\item{tags}{logical: show A & B on the two panels?}

\item{intervals}{amount of probability mass to include in the two intervals}
}
\value{
a list of three plots: the interval plot, the entropy plot, and the two stuck together
}
\description{
Plot both the posterior intervals of interaction score and
  position in the concordance distribution
}
