% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbn_sampling.R
\name{filter_cpt}
\alias{filter_cpt}
\title{Filter the cpt of node n given the value of its parents}
\usage{
filter_cpt(df_cpt_n, parents_values)
}
\arguments{
\item{df_cpt_n}{a dataframe representing the cpt}

\item{parents_values}{a list with the parents of n and their values}
}
\value{
the filtered cpt
}
\description{
Filter the cpt of node n given the value of its parents
}
\examples{
filter_cpt(df_cp_n, parents_values)
}
