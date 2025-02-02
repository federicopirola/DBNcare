% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbn_initialization.R
\name{summary_DBN}
\alias{summary_DBN}
\title{Function returning a summary of input Dynamic Bayesian Network}
\usage{
summary_DBN(DBN, test = FALSE)
}
\arguments{
\item{DBN}{object of class 'DBN'}
}
\description{
Function returning a summary of input Dynamic Bayesian Network
}
\details{
The function returns from terminal:
\itemize{
\item the model string corresponding to the DAG at G_0
\item the model string corresponding to the DAG at G_transition
\item the Markovian Order of the process
\item the number of nodes, divided between Dynamic nodes and Static nodes
\item the number of arcs, divided between arcs of G_0, inner-slice arcs and intra-slice arcs of G_transition
\item the average Markov Blanket size of the nodes
\item the average neighboorhood size of the nodes
}
}
\examples{
summary(DBN_example)
}
