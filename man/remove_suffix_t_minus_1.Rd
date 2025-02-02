% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbn_visualization.R
\name{remove_suffix_t_minus_1}
\alias{remove_suffix_t_minus_1}
\title{Removes the suffix _t-1 from the node's name}
\usage{
remove_suffix_t_minus_1(node_name)
}
\arguments{
\item{node_name}{a string representing the name of the node}
}
\value{
the modified string
}
\description{
Removes the suffix _t-1 from the node's name
}
\examples{
remove_suffix_t_minus_1("A_t-1")
}
