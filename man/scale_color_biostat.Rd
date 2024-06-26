% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biostatcolors.R
\name{scale_color_biostat}
\alias{scale_color_biostat}
\title{Color scale constructor for BiostatOmics colors}
\usage{
scale_color_biostat(palette = "main", continuous = FALSE, reverse = FALSE, ...)
}
\arguments{
\item{palette}{A character string specifying the name of the desired palette from the `biostat_palettes`. Available options include: "main", "cblindfriendly", "sunshine", "hot", "warm" and "cold" (Default: "main").}

\item{continuous}{A logical value indicating whether the color aesthetic represents continuous data (Default: `FALSE`).}

\item{reverse}{A logical value indicating whether the colors in the selected palette should be reversed (Default: `FALSE`).}

\item{...}{Additional arguments passed either to `discrete_scale` or `scale_color_gradientn` from the `ggplot2` package, depending on the value of the `continuous` parameter.}
}
\value{
A `ggplot2` scale function suitable for adding to a `ggplot2` object.
}
\description{
The `scale_color_biostat` function provides a mechanism to integrate BioStatOmics group's curated color palettes into `ggplot2` visualizations.
}
\details{
The `scale_color_biostat` function acts as a bridge between these palettes and the `ggplot2` package, allowing users to apply the palettes to their plots Depending on the nature of the data (continuous or discrete), the function intelligently selects the appropriate scale from `ggplot2` to render the colors.

When the `continuous` parameter is set to `TRUE`, the function employs the `scale_color_gradientn` function from `ggplot2` to generate a continuous color scale. Conversely, for discrete data, the `discrete_scale` function is utilized. This ensures that the chosen palette is optimally represented in the plot, irrespective of the data type.
}
\examples{
library(ggplot2)
data("iris")
g <- ggplot(iris, aes(Sepal.Width, Sepal.Length, color = Species))
g <- g + geom_point(size = 4)
g <- g + scale_color_biostat(palette = "main")

}
\author{
Pedro Salguero Garcia, Maider Aguerralde Martin. Maintainer: magumar2@upv.edu.es
}
