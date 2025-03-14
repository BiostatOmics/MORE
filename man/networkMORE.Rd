% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output_analysis.R
\name{networkMORE}
\alias{networkMORE}
\title{networkMORE}
\usage{
networkMORE(
  outputRegpcond,
  cytoscape = TRUE,
  group1 = NULL,
  group2 = NULL,
  pc = 0,
  pathway = NULL,
  annotation = NULL,
  save = FALSE
)
}
\arguments{
\item{outputRegpcond}{Output object of RegulationPerCondition applied to MORE main function.}

\item{cytoscape}{TRUE for plotting the network in Cytoscape. FALSE to plot the network in R.}

\item{group1}{Name of the group to take as reference in the differential network creation. It also can be used for creating networks of a specific group. If it is not provided the networks of all conditions will be plotted. By default, NULL.}

\item{group2}{Name of the group to compare to the reference in the differential network creation. By default, NULL.}

\item{pc}{Percentile to consider to plot the most affecting regulators into the target omic. It must be a value comprissed between 0 and 1. By default, 0.}

\item{pathway}{If provided, the function will print the regulatory network involved in the specified pathway instead of the entire regulatory network. By default, NULL.}

\item{annotation}{Annotation matrix with target features in the first column, GO terms in the second and GO term description in the third. Only necessary when a specific pathway has to be plotted. By default, NULL.}

\item{save}{If TRUE a gml extension network is saved when cytoscape = FALSE. By default, FALSE.}
}
\value{
Plot of the network induced from more.
}
\description{
\code{networkMORE} Function to be applied to RegulationPerConidtion function output.
}
