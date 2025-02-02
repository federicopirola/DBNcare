% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbn_initialization.R
\name{add_node_DBN}
\alias{add_node_DBN}
\title{Function for node addiction in Dynamic Bayesian Network}
\usage{
add_node_DBN(DBN, node, type = "Dynamic")
}
\arguments{
\item{DBN}{object of class 'DBN'}

\item{node}{object of class 'character': name of the node to add}

\item{type}{'Dynamic' (default) or 'Static' node}
}
\value{
object of class 'DBN'
}
\description{
Function for node addiction in Dynamic Bayesian Network
}
\examples{
DBN_example <- add_node_DBN(DBN=DBN_example, node="F", type='Dynamic')
DBN_example <- add_node_DBN(DBN=DBN_example, node="K", type='Static')
}
