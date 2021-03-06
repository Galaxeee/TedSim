% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GeneExpression_sim.R
\name{SimulateCellStates}
\alias{SimulateCellStates}
\title{Simulate cell states based on the state tree}
\usage{
SimulateCellStates(
  par,
  cell_edges,
  state_edges,
  sif_mean,
  S,
  p_a,
  p_edge = NULL,
  max_walk = 2
)
}
\arguments{
\item{par}{current parent node}

\item{cell_edges}{edge list of the cell lineage tree}

\item{state_edges}{edge list of the cell state tree}

\item{sif_mean}{State Identity Factors}

\item{S}{current state of the parent node}

\item{p_a}{asymmetric division rate}

\item{p_edge}{branching possibilities}

\item{max_walk}{maximum walk distance on the state tree for one asymmetric division}
}
\description{
Simulate cell states based on the state tree
}
