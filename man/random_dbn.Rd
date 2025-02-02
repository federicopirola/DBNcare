% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbn_random.R
\name{random_dbn}
\alias{random_dbn}
\title{random_dbn}
\usage{
random_dbn(nodes_names, is_same, prob_edge_intraslice, prob_edge_interslice)
}
\arguments{
\item{nodes_names}{a list of strings representing the names of nodes in the
network}

\item{is_same}{TRUE if you want to keep the same edges defined in G_0
in G_transition. FALSE if you DO NOT want G_transition to have all the relationships
defined in G_0.}

\item{prob_edge_intraslice}{a double parameter representing the probability of an edge to be present in the intraslice.}

\item{prob_edge_interslice}{a double parameter representing the probability of an edge to be present in the interslice.}
}
\value{
a DBN object.
}
\description{
Generates a random structure of a Dynamic Bayesian Network using the order
of the variable nodes_names
}
\examples{
random_dbn(c("A","B","C"), TRUE, 0.6, TRUE, 2)
}
