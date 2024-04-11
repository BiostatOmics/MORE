#' @import ggplot2

biostat_colors <- c(`green` = "#3bdfa2",
                   `red` = "#d11d56",
                   `yellow` = "#d0ef3f",
                   `purple` = "#7332a7",
                   `violet` = "#d6b0e9",
                   `blue` = "#7480ca",
                   `b_dgreen` = "#126D6D",
                   `b_lgreen` = "#6CD0D0",
                   `b_yellow`="#FDDB73",
                   `b_orange` = "#FFB000")

#' Function to extract biostat colors as hex codes
#'
#' @param ... Character names of biostat_colors

biostat_cols <- function(...) {
  cols <- c(...)
  
  if (is.null(cols))
    return (biostat_colors)
  
  biostat_colors[cols]
}

#' Biostatomics color palette
biostat_palettes <- list(
  `main`  = biostat_cols("purple", "violet", "green","yellow","red","blue"),
  `cblindfriendly`= biostat_cols("b_dgreen", "b_lgreen","violet","b_yellow","b_orange"),
  `sunshine`= biostat_cols("b_orange", "blue"),
  
  `hot`   = biostat_cols("red", "b_orange"),
  `warm`  = biostat_cols("blue", "violet"),
  `grass`  = biostat_cols("green", "b_dgreen")
  )

#' Return Biostatomics Color Palettes
#'
#' @description
#' The `getbiostatPalettes` function retrieves a collection of color palettes, specifically designed for scientific visualizations. These palettes are part of the biostat collection.
#'
#' @details
#' By using the `getbiostatPalettes` function, users can access these palettes and incorporate them into their visualizations, ensuring that their plots and graphs are both informative and visually appealing.
#'
#' @return
#' A list containing the various color palettes from the BioStatOmics collection. Each palette in the list is represented as a vector of color values.
#'
#' @author
#' Pedro Salguero Garcia, Maider Aguerralde Martin. Maintainer: magumar2@upv.edu.es
#'
#' @export
#'
#' @examples
#' getbiostatPalettes()
#'
getbiostatPalettes <- function(){
  return(biostat_palettes)
}


#' Retrieve BioStatOmics Main Color Set
#'
#' @description
#' The `getbiostatColors` function provides access to a curated set of colors that are part of the BiostatOmicsColors package. These colors have been specifically chosen for their utility in scientific visualizations.
#'
#' @details
#' When using the `getbiostatColors` function, users can seamlessly integrate these colors into their R visualizations, benefiting from the expertise embedded in the BiostatOmics color selection.
#'
#' @return
#' A list containing the primary colors from the BioStatOmics collection. Each color in the list is represented as a hexadecimal color value.
#'
#' @author
#' Pedro Salguero Garcia, Maider Aguerralde Martin. Maintainer: magumar2@upv.edu.es
#'
#' @export
#'
#' @examples
#' getbiostatColors()
#'
getbiostatColors <- function(){
  return(biostat_colors)
}

#' Interpolate a BioStatOmics Color Palette
#'
#' @description
#' The `biostat_pal` function offers a flexible way to interpolate colors from the BiostatOmics color palettes. This function provides an interface to generate a range of colors based on the selected BioStatOmics palette, allowing for enhanced customization in scientific visualizations.
#'
#' @details
#' The BioStatOmics color palettes, available in the BioStatOmicsColors package, have been specifically curated for scientific visualizations. The `biostat_pal` function leverages the `colorRampPalette` function from the `grDevices` package to interpolate between the colors of the chosen BioStatOmics palette. This interpolation capability ensures that users can generate a continuous range of colors, suitable for representing a wide variety of data types and scales. Whether visualizing continuous data gradients or categorical distinctions, the interpolated BiostatOmics palettes can provide clarity and aesthetic appeal to the visual representation.
#'
#' @param palette
#' A character string specifying the name of the desired palette from the `biostat_palettes`. Available options include: "main", "cblindfriendly", "sunshine", "hot", "warm" and "cold".
#' @param reverse
#' A logical value indicating whether the colors in the selected palette should be reversed. Default is `FALSE`.
#' @param ...
#' Additional arguments to be passed to the `colorRampPalette` function from the `grDevices` package.
#'
#' @author
#' Pedro Salguero Garcia, Maider Aguerralde Martin. Maintainer: magumar2@upv.edu.es
biostat_pal <- function(palette = "main", reverse = FALSE, ...) {
  pal <- biostat_palettes[[palette]]
  
  if (reverse) pal <- rev(pal)
  
  grDevices::colorRampPalette(pal, ...)
}

#' Retrieve Colors from BioStatOmics group's Palettes
#'
#' @description
#' The `colorbiostat` function facilitates the extraction of a specified number of colors from the BioStatOmics group's curated color palettes. This function is designed to obtain a set of colors for their scientific visualizations.
#'
#' @details
#' BioStatOmics group's color palettes, available within the package, are tailored for scientific data visualization. The `colorbiostat` function is built upon these palettes, offering flexibility in color selection based on the user's requirements. It integrates with the `palette` argument to choose the color thematic.
#'
#' It's essential to note that if the requested number of colors (`n`) is less than or equal to the size of the chosen palette, the function will directly extract the colors without interpolation. However, if `n` surpasses the palette size, interpolation is employed to generate the required colors.
#'
#' @param n
#' An integer specifying the number of colors to be extracted from the chosen palette.
#' @param reverse
#' A logical value indicating whether the colors in the selected palette should be reversed (Default: `FALSE`).
#' @param palette
#' A character string specifying the name of the desired palette from the `biostat_palettes`. Available options include: "main", "cblindfriendly", "sunshine", "hot", "warm" and "cold" (Default: "main").
#'
#' @return
#' A character vector of colors corresponding to the specified number and palette.
#'
#' @author
#' Pedro Salguero Garcia, Maider Aguerralde Martin. Maintainer: magumar2@upv.edu.es
#'
#' @export
#'
#' @examples
#' library(ggplot2)
#' data("iris")
#' colorSpecies <- colorbiostat(3, palette = "cold")
#' plot(x = iris$Sepal.Length, y = iris$Sepal.Width, col = colorSpecies[iris$Species], pch = 16)
#'
colorbiostat <- function(n, reverse = FALSE, palette = "main"){
  if(!palette %in% names(biostat_palettes)){
    message("Palette musst be some of the following palettes: ", paste(names(biostat_palettes),sep= " ", collapse= ", "))
  }
  
  pal <- biostat_pal(palette = palette, reverse = reverse)
  
  if(n<=length(getbiostatPalettes()[[palette]])){
    colors <- getbiostatPalettes()[[palette]][1:n]
    names(colors) <- NULL
    return(colors)
  }
  
  return(pal(n))
}

#' Color scale constructor for BiostatOmics colors
#' @description
#' The `scale_color_biostat` function provides a mechanism to integrate BioStatOmics group's curated color palettes into `ggplot2` visualizations.
#'
#' @details
#' The `scale_color_biostat` function acts as a bridge between these palettes and the `ggplot2` package, allowing users to apply the palettes to their plots Depending on the nature of the data (continuous or discrete), the function intelligently selects the appropriate scale from `ggplot2` to render the colors.
#'
#' When the `continuous` parameter is set to `TRUE`, the function employs the `scale_color_gradientn` function from `ggplot2` to generate a continuous color scale. Conversely, for discrete data, the `discrete_scale` function is utilized. This ensures that the chosen palette is optimally represented in the plot, irrespective of the data type.
#'
#' @param palette
#' A character string specifying the name of the desired palette from the `biostat_palettes`. Available options include: "main", "cblindfriendly", "sunshine", "hot", "warm" and "cold" (Default: "main").
#' @param continuous
#' A logical value indicating whether the color aesthetic represents continuous data (Default: `FALSE`).
#' @param reverse
#' A logical value indicating whether the colors in the selected palette should be reversed (Default: `FALSE`).
#' @param ...
#' Additional arguments passed either to `discrete_scale` or `scale_color_gradientn` from the `ggplot2` package, depending on the value of the `continuous` parameter.
#' @return
#' A `ggplot2` scale function suitable for adding to a `ggplot2` object.
#'
#' @author
#' Pedro Salguero Garcia, Maider Aguerralde Martin. Maintainer: magumar2@upv.edu.es
#'
#' @examples
#' library(ggplot2)
#' data("iris")
#' g <- ggplot(iris, aes(Sepal.Width, Sepal.Length, color = Species))
#' g <- g + geom_point(size = 4)
#' g <- g + scale_color_biostat(palette = "main")
#'
#' @export
scale_color_biostat <- function(palette = "main", continuous = FALSE, reverse = FALSE, ...) {
  if(!palette %in% names(biostat_palettes)){
    message("Palette must be some of the following palettes: ", paste(names(biostat_palettes),sep=" ", collapse=", "))
    stop()
  }
  
  pal <- biostat_pal(palette = palette, reverse = reverse)
  
  if (!continuous) {
    ggplot2::discrete_scale("colour", paste0("biostat_", palette), palette = pal, ...)
  } else {
    ggplot2::scale_color_gradientn(colours = pal(256), ...)
  }
}

#' Fill scale constructor for biostatomics colors
#' @description
#' The `scale_fill_biostat` function provides a mechanism to integrate BioStatOmics group's curated color palettes into `ggplot2` visualizations.
#'
#' @details
#' The `scale_fill_biostat` function acts as a bridge between these palettes and the `ggplot2` package, allowing users to apply the palettes to their plots Depending on the nature of the data (continuous or discrete), the function intelligently selects the appropriate scale from `ggplot2` to render the colors.
#'
#' When the `continuous` parameter is set to `TRUE`, the function employs the `scale_fill_gradientn` function from `ggplot2` to generate a continuous color scale. Conversely, for discrete data, the `discrete_scale` function is utilized. This ensures that the chosen palette is optimally represented in the plot, irrespective of the data type.
#'
#' @param palette
#' A character string specifying the name of the desired palette from the `biostat_palettes`. Available options include: "main", "nature", "sunshine", "hot", "warm", "cold", and "complete" (Default: "main").
#' @param continuous
#' A logical value indicating whether the color aesthetic represents continuous data (Default: `FALSE`).
#' @param reverse
#' A logical value indicating whether the colors in the selected palette should be reversed (Default: `FALSE`).
#' @param ...
#' Additional arguments passed either to `discrete_scale` or `scale_fill_gradientn` from the `ggplot2` package, depending on the value of the `continuous` parameter.
#'
#' @return
#' A `ggplot2` scale function suitable for adding to a `ggplot2` object.
#'
#' @author
#' Pedro Salguero Garcia, Maider Aguerralde Martin. Maintainer: magumar2@upv.edu.es
#'
#' @examples
#' library(ggplot2)
#' data("iris")
#' g <- ggplot(iris, aes(x = Sepal.Width, fill = Species))
#' g <- g + geom_histogram(binwidth = 0.2, alpha = 0.8)
#' g <- g + labs(title = "Histogram of Sepal Width", x = "Sepal Width", y = "Frequency")
#' g <- g + scale_fill_biostat(palette = "main")
#'
#' @export
scale_fill_biostat <- function(palette = "main", continuous = FALSE, reverse = FALSE, ...) {
  if(!palette %in% names(biostat_palettes)){
    message("Palette must be some of the following palettes: ", paste(names(biostat_palettes),sep=" ", collapse=", "))
    stop()
  }
  
  pal <- biostat_pal(palette = palette, reverse = reverse)
  
  if(!continuous){
    ggplot2::discrete_scale("fill", paste0("biostat_", palette), palette = pal, ...)
  }else{
    ggplot2::scale_fill_gradientn(colours = pal(256), ...)
  }
  
}

