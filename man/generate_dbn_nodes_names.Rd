% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbn_random.R
\name{generate_dbn_nodes_names}
\alias{generate_dbn_nodes_names}
\title{generate_dbn_nodes_names}
\usage{
generate_dbn_nodes_names(nodes_names)
}
\arguments{
\item{nodes_names}{a vector containing the names of nodes.}
}
\value{
a list containing the nodes names for each time slice.
}
\description{
the function returns the nodes names for g_0 and g_transition.
}
\examples{
generate_dbn_nodes_names(c("A", "B", "C"))
}
