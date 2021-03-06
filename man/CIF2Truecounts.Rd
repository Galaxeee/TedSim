% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GeneExpression_sim.R
\name{CIF2Truecounts}
\alias{CIF2Truecounts}
\title{Simulate true count matrix given simulated CIF values}
\usage{
CIF2Truecounts(
  ngenes,
  ncif,
  ge_prob,
  ncells,
  cif_res,
  prop_hge = 0.015,
  mean_hge = 5,
  scale_s = 1
)
}
\arguments{
\item{ngenes}{Number of genes}

\item{ncif}{number of cifs simulated}

\item{ge_prob}{ge probability}

\item{ncells}{Number of cells simulated}

\item{cif_res}{the CIFs simulated}

\item{mean_hge}{the mean of hge, default is 5}

\item{scale_s}{transcription rate scaler, or a vector to specify cell-type specific scale_s}

\item{prob_hge}{the probability of hge, default is 0.015}
}
\description{
Simulate true count matrix given simulated CIF values
}
