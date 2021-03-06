% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GeneExpression_sim.R
\name{WalkTree}
\alias{WalkTree}
\title{Find the child state based on parent state and step size}
\usage{
WalkTree(S_parent, child, state_edges, sif_mean, p_edge, step)
}
\arguments{
\item{S_parent}{current state of the parent node}

\item{child}{child node label}

\item{state_edges}{edge list of the cell state tree}

\item{sif_mean}{edge list of the cell lineage}
}
\description{
Find the child state based on parent state and step size
}
