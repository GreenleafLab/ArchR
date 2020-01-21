##########################################################################################
# ggPlot Wrapper Methods For Easy Plotting
##########################################################################################

#' A ggplot-based dot plot wrapper function
#'
#' This function is a wrapper around ggplot geom_point to allow for a more intuitive plotting of ArchR data.
#'
#' @param x A numeric vector containing the x-axis values for each point.
#' @param y A numeric vector containing the y-axis values for each point.
#' @param color A numeric/categorical vector containing coloring information for each point.
#' @param discrete A boolean value indicating whether the supplied data is discrete (`TRUE`) or continuous (`FALSE`).
#' @param discreteSet The name of a custom palette from `ArchRPalettes` to use for categorical/discrete color.
#' @param continuousSet The name of a custom palette from `ArchRPalettes` to use for numeric color.
#' @param labelMeans A boolean value indicating whether the mean of each categorical/discrete color should be labeled.
#' @param pal QQQ UNCLEAR WHY THIS EXISTS OR WHAT IT IS FOR. IS THIS FOR A custom palette used to override discreteSet/continuousSet for coloring vector.
#' @param defaultColor The default color for points that do not have another color applied (i.e. `NA` values).
#' @param colorDensity A boolean value indicating whether the density of points on the plot should be indicated colorimetrically. If `TRUE`, continuousSet is used as the color palette.
#' @param size The numeric size of the points to be plotted.
#' @param xlim A numeric vector of two values indicating the lower and upper bounds of the x-axis on the plot.
#' @param ylim A numeric vector of two values indicating the lower and upper bounds of the y-axis on the plot.
#' @param extend A numeric value indicating the fraction to extend the x-axis and y-axis beyond the maximum and minimum values if `xlim` and `ylim` are not provided. For example, 0.05 will extend the x-axis and y-axis by 5 percent on each end.
#' @param xlabel The label to plot for the x-axis.
#' @param ylabel The label to plot for the y-axis.
#' @param title The title of the plot.
#' @param randomize A boolean value indicating whether to randomize the order of the points when plotting.
#' @param seed A numeric seed number for use in randomization.
#' @param colorTitle QQQ
#' @param colorOrder QQQ
#' @param alpha A number indicating the transparency to use for each point. See `ggplot2` for more details.
#' @param baseSize The base font size to use in the plot.
#' @param ratioYX The aspect ratio of the x and y axes on the plot.
#' @param labelType A string indicating how to label the points on the plot. Options include "ggrepel", "shadowtext".
#' @param fgColor The foreground color of the plot.
#' @param bgColor The background color of the plot.
#' @param labelSize The numeric font size of labels.
#' @param addFit A string indicating if a fit/regression line (see `geom_smooth` methods) should be included in the plot and what method to use for this fit.
#' @param rastr A boolean value that indicates whether the plot should be rasterized using `ggrastr`. This does not rasterize lines and labels, just the internal portions of the plot.
#' @param dpi The resolution in dots per inch to use for the plot.
#' @export
ggPoint <- function(
    x = NULL, 
    y = NULL, 
    color = NULL, 
    discrete = TRUE, 
    discreteSet = "stallion",
    continuousSet = "solar_extra", 
    labelMeans = FALSE,  
    pal = NULL, 
    defaultColor = "lightGrey",
    colorDensity = FALSE,
    size = 1, 
    xlim = NULL, 
    ylim = NULL, 
    extend = 0.05, 
    xlabel = "x", 
    ylabel = "y", 
    title = "", 
    randomize = FALSE, 
    seed = 1,
    colorTitle = NULL, 
    colorOrder = NULL, 
    alpha = 1, 
    baseSize = 6, 
    ratioYX = 1, 
    labelType = "ggrepel", 
    fgColor = NULL, 
    bgColor = "white", 
    labelSize = 1.5,
    addFit = NULL, 
    rastr = FALSE, 
    dpi = 300,
    ...
    ){

    .validInput(input = x, name = "x", valid = c("numeric"))
    .validInput(input = y, name = "y", valid = c("numeric"))
    .validInput(input = color, name = "color", valid = c("numeric", "character", "null"))
    .validInput(input = discrete, name = "discrete", valid = c("boolean"))
    .validInput(input = discreteSet, name = "discreteSet", valid = c("character"))
    .validInput(input = continuousSet, name = "continuousSet", valid = c("character"))
    .validInput(input = labelMeans, name = "labelMeans", valid = c("boolean"))
    .validInput(input = pal, name = "pal", valid = c("character"))
    .validInput(input = defaultColor, name = "defaultColor", valid = c("character"))
    .validInput(input = size, name = "size", valid = c("numeric"))
    .validInput(input = xlim, name = "xlim", valid = c("numeric", "null"))
    .validInput(input = ylim, name = "ylim", valid = c("numeric", "null"))
    .validInput(input = extend, name = "extend", valid = c("numeric"))
    .validInput(input = xlabel, name = "xlabel", valid = c("character"))
    .validInput(input = ylabel, name = "ylabel", valid = c("character"))
    .validInput(input = title, name = "title", valid = c("character"))
    .validInput(input = randomize, name = "randomize", valid = c("boolean"))
    .validInput(input = seed, name = "seed", valid = c("integer"))
    .validInput(input = colorTitle, name = "colorTitle", valid = c("character", "null"))
    .validInput(input = colorOrder, name = "colorOrder", valid = c("character", "null"))
    .validInput(input = alpha, name = "alpha", valid = c("numeric"))
    .validInput(input = baseSize, name = "baseSize", valid = c("numeric"))
    .validInput(input = ratioYX, name = "ratioYX", valid = c("numeric"))
    .validInput(input = labelType, name = "labelType", valid = c("character"))
    .validInput(input = fgColor, name = "fgColor", valid = c("character", "null"))
    .validInput(input = bgColor, name = "bgColor", valid = c("character", "null"))
    .validInput(input = labelSize, name = "labelSize", valid = c("numeric"))
    .validInput(input = addFit, name = "addFit", valid = c("character", "null"))
    .validInput(input = rastr, name = "rastr", valid = c("boolean"))
    .validInput(input = dpi, name = "dpi", valid = c("numeric"))

    stopifnot(length(y) == length(x))

    if(randomize){
      set.seed(seed)
      idx <- sample(seq_along(x), length(x))
    }else{
      idx <- seq_along(x)
    }

    df <- data.frame(x = x, y = y)
    include <- which(is.finite(x) & is.finite(y))
    
    if(length(include) != length(x)){
      message("Some values are not finite! Excluding these points!")
      df <- df[include,]
      x <- x[include]
      y <- y[include]
      if(!is.null(color)){
        color <- color[include]
      }
    }

    if(is.null(xlim)){
        xlim <- range(df$x) %>% extendrange(f = extend)
    }

    if(is.null(ylim)){
        ylim <- range(df$y) %>% extendrange(f = extend)
    }

    ratioXY <- ratioYX * diff(xlim)/diff(ylim)

    #Plot
    .requirePackage("ggplot2")

    if (is.null(color) & !colorDensity) {

      p <- ggplot(df[idx,], aes(x = x, y = y)) + 
          coord_equal(ratio = ratioXY, xlim = xlim, ylim = ylim, expand = F) + 
          xlab(xlabel) + ylab(ylabel) + 
          ggtitle(title) + 
          theme_ArchR(baseSize = baseSize)

      if(rastr){
        if(!requireNamespace("ggrastr", quietly = TRUE)){
          message("ggrastr is not available for rastr of points, continuing without rastr!")
          p <- p + geom_point(size = size, alpha = alpha, color = defaultColor)
        }else{
          .requirePackage("ggrastr")
          p <- p + geom_point_rast(
              size = size, raster.dpi = dpi, alpha = alpha, color = defaultColor)
        }
      }else{
        p <- p + geom_point(size = size, alpha = alpha, color = defaultColor)
      }
        
    }else {

        if(colorDensity){
          
          discrete <- FALSE
          df <- .getDensity(x, y, n = 100, sample = NULL) #change
          df <- df[order(df$density), ,drop=FALSE]
          df$color <- df$density
          
          if(is.null(colorTitle)){
            colorTitle <- "density"
          }

        }else if(discrete){

          if(!is.null(colorOrder)){
            if(!all(color %in% colorOrder)){
              stop("Not all colors are in colorOrder!")
            }
          }else{
            colorOrder <- gtools::mixedsort(unique(color))
          }

          if(is.null(colorTitle)){
            colorTitle <- "color"
          }

          stopifnot(length(color) == nrow(df))
          df$color <- factor(color, levels = colorOrder)
          
        }else{
          stopifnot(length(color) == nrow(df))
          df$color <- color
        }

        p <- ggplot(df[idx,], aes(x = x, y = y, color = color)) +  
              coord_equal(ratio = ratioXY, xlim = xlim, ylim = ylim, expand = F) + 
              xlab(xlabel) + ylab(ylabel) + 
              ggtitle(title) + theme_ArchR(baseSize = baseSize) +
              theme(legend.direction = "horizontal", legend.box.background = element_rect(color = NA)) +
              labs(color = colorTitle)

      if(rastr){
        
        if(!requireNamespace("ggrastr", quietly = TRUE)){
          message("ggrastr is not available for rastr of points, continuing without rastr!")
          p <- p + geom_point(size = size, alpha = alpha)
        }else{
          .requirePackage("ggrastr")
          p <- p + geom_point_rast(size = size, raster.dpi = dpi, alpha = alpha, 
            raster.width=par('fin')[1], raster.height = (ratioYX * par('fin')[2]))
        }
      
      }else{

          p <- p + geom_point(size = size, alpha = alpha)

      }

      if (discrete) {
          
          if (!is.null(pal)) {
              p <- p + scale_color_manual(values = pal)
          }else {
              p <- p + scale_color_manual(values = paletteDiscrete(set = discreteSet, values = colorOrder))
          }

          if (labelMeans) {
              
              dfMean <- split(df, df$color) %>% lapply(., function(x) {
                data.frame(x = mean(x[, 1]), y = mean(x[, 2]), color = x[1, 3])
              }) %>% Reduce("rbind", .)

              #Check Packages!
              if(tolower(labelType) == "repel" | tolower(labelType) == "ggrepel"){
                .requirePackage("ggrepel")
              }else if(tolower(labelType) == "shadow" | tolower(labelType) == "shadowtext"){
                .requirePackage("shadowtext")
              }

              if(tolower(labelType) == "repel" | tolower(labelType) == "ggrepel"){

                if(!is.null(fgColor)){
                  p <- p + ggrepel::geom_label_repel(data = dfMean, aes(x, y, label = color), color = fgColor, size = labelSize)
                }else{
                  p <- p + ggrepel::geom_label_repel(data = dfMean, aes(x, y, label = color), size = labelSize)
                }

              }else if(tolower(labelType) == "shadow" | tolower(labelType) == "shadowtext"){

                if(!is.null(fgColor)){
                  p <- p + shadowtext::geom_shadowtext(data = dfMean, aes(x, y, label = color), color = fgColor, bg.colour = bgColor, size = labelSize)
                }else{
                  p <- p + shadowtext::geom_shadowtext(data = dfMean, aes(x, y, label = color), bg.colour = bgColor, size = labelSize)
                }
        
              }else{
                stop("Error unrecognized label type!")
              }
          }

      }else{

          if (!is.null(pal)) {
              p <- p + scale_colour_gradientn(colors = pal)
          }else {
              p <- p + scale_colour_gradientn(colors = paletteContinuous(set = continuousSet))
          }
      }

    }

    if (!is.null(addFit)) {
        p <- p + geom_smooth(data = df, aes(color = NULL), method = addFit, color = "black") + 
          ggtitle(paste0(title, "\nPearson = ", round(cor(df$x, df$y), 3), "\nSpearman = ", round(cor(df$x, df$y, method = "spearman"), 3)))
    }

    p <- p + theme(legend.position = "bottom", legend.key = element_rect(size = 2))#, legend.spacing.x = unit(0.1, 'cm'), legend.spacing.y = unit(0.1, 'cm'))

    if(!is.null(ratioYX)){
      attr(p, "ratioYX") <- ratioYX
    }

    return(p)

}

#' A ggplot-based one-to-one dot plot wrapper function
#'
#' This function is a wrapper around ggplot geom_point to allow for plotting one-to-one sample comparisons in ArchR.
#'
#' @param x A numeric vector containing the x-axis values for each point.
#' @param y A numeric vector containing the y-axis values for each point.
#' @param size The numeric size of the points to plot.
#' @param alpha A number indicating the transparency to use for each point. See `ggplot2` for more details.
#' @param xlabel The label to plot for the x-axis.
#' @param ylabel The label to plot for the y-axis.
#' @param title The title of the plot.
#' @param min QQQ UNCLEAR. ALSO SHOULDNT THIS APPLY TO BOTH AXES? x-limits min quantile [0,1]
#' @param max QQQ UNCLEAR. ALSO SHOULDNT THIS APPLY TO BOTH AXES? x-limits max quantile [0,1]
#' @param nPlot The number of points to plot. When this value is less than the total points, the `sample` function is used to extract random data points to be plotted.
#' @param nKernel QQQ NO IDEA WHAT KDE2D IS OR WHAT THIS IS USED FOR. The value of `n` to use for the QQQ `kde2d` from the `MASS` package.
#' @param densityMax The quantile that should be represented by the maximum color on the continuous scale designated by `pal`. Values above `densityMax` will be thresholded to the maximum color on the color scale.
#' @param extend QQQ DOUBLE CHECK. ALSO, HOW DOES THIS VALUE RELATE TO min AND max?? A numeric value indicating the fraction to extend the x-axis and y-axis beyond the maximum value on either axis. For example, 0.05 will extend the x-axis and y-axis by 5 percent on each end beyond `max(c(x,y))` and `min(c(x,y))`.
#' @param baseSize The base font size to use in the plot.
#' @param rastr A boolean value that indicates whether the plot should be rasterized. This does not rasterize lines and labels, just the internal portions of the plot.
#' @param pal A custom palette from `ArchRPalettes` used to display the density of points on the plot.
#' @param ... Additional params to be supplied to ggPoint
#' @export
ggOneToOne <- function (
  x = NULL,
  y = NULL,
  size = 2, 
  alpha = 1,
  xlabel = "x", 
  ylabel = "y", 
  title = "Correlation",
  min = 0.05, 
  max = 0.9999, 
  nPlot = 100 * 10^3, 
  nKernel = 100, 
  densityMax = 0.95, 
  extend = 0.05, 
  baseSize = 6, 
  rastr = TRUE,
  pal = paletteContinuous(set = "blue_yellow"),
  ...
  ){

  .validInput(input = x, name = "x", valid = c("numeric"))
  .validInput(input = y, name = "y", valid = c("numeric"))
  .validInput(input = size, name = "size", valid = c("numeric"))
  .validInput(input = alpha, name = "alpha", valid = c("numeric"))
  .validInput(input = xlabel, name = "xlabel", valid = c("character"))
  .validInput(input = ylabel, name = "ylabel", valid = c("character"))
  .validInput(input = title, name = "title", valid = c("character"))
  .validInput(input = min, name = "min", valid = c("numeric"))
  .validInput(input = max, name = "max", valid = c("numeric"))
  .validInput(input = nPlot, name = "nPlot", valid = c("integer"))
  .validInput(input = nKernel, name = "nKernel", valid = c("numeric"))
  .validInput(input = densityMax, name = "densityMax", valid = c("numeric"))
  .validInput(input = extend, name = "extend", valid = c("numeric"))
  .validInput(input = baseSize, name = "baseSize", valid = c("numeric"))
  .validInput(input = rastr, name = "rastr", valid = c("boolean"))
  .validInput(input = pal, name = "pal", valid = c("character"))
  
  #Check for NA
  idx <- which(!is.na(x) & !is.na(y) & !is.infinite(x) & !is.infinite(y))
  x <- x[idx]
  y <- y[idx]
  
  #Ratio X/Y
  lim <- quantile(c(x, y), c(min, max)) %>% extendrange(f = extend)
  ratioXY <- diff(lim)/diff(lim)
  
  #Calculate Correlations
  pearson <- round(cor(x, y, method = "pearson", use = "complete"), 3)
  spearman <- round(cor(x, y, method = "spearman", use = "complete"), 3)
  title <- sprintf("%s \nPearson = %s , Spearman = %s", title, pearson, spearman)
  
  #Get Density
  message("adding denisty...")
  df <- .getDensity(x, y, n = nKernel, sample = nPlot) #change
  df <- df[order(df[, "density"]), ]
  
  #GGPlot
  message("plotting...")
  gg <- ggPoint(
      x = df$x, 
      y = df$y, 
      color = df$density, 
      pal = pal,
      xlabel = xlabel,
      ylabel = ylabel,
      discrete = FALSE, 
      colorTitle = "density",
      xlim = lim, 
      ylim = lim, 
      size = size, 
      alpha = alpha, 
      title = title, 
      baseSize = baseSize,
      rastr = rastr
    ) + geom_abline(slope = 1, intercept = 0, lty = "dashed")

  return(gg)

}

.getDensity <- function(x, y, n = 100, sample = NULL, densityMax = 0.95){
  #modified from http://slowkow.com/notes/ggplot2-color-by-density/
  df <- data.frame(x=x,y=y)
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  df$density <- dens$z[ii]
  df$density[df$density > quantile(unique(df$density),densityMax)] <- quantile(unique(df$density),densityMax) #make sure the higher end doesnt bias colors
  if(!is.null(sample)){
    df <- df[sample(nrow(df), min(sample,nrow(df))),]
  }
  return(df)
}

#' A ggplot-based violin plot wrapper function
#'
#' This function is a wrapper around ggplot geom_violin to allow for plotting violin plots in ArchR.
#' 
#' @param x A character vector containing the categorical x-axis values for each y-axis value.
#' @param y A numeric vector containing the y-axis values for each point.
#' @param xlabel The label to plot for the x-axis.
#' @param ylabel The label to plot for the y-axis.
#' @param xOrder A character vector indicating a custom order for plotting x-axis categorical values. Should contain all possible values of `x` in the desired order.
#' @param addPoints A boolean value indicating whether individual points should be added to the plot using `geom_quasirandom`.
#' @param size QQQ ARE BOXPLOTS REQUIRED ON THE VIOLIN? SEEMS ODD The line width for boxplot lines.
#' @param baseSize The base font size to use in the plot.
#' @param ratioYX The aspect ratio of the x and y axes on the plot.
#' @param sampleRatio QQQ The fraction of the total number of points to be displayed over violins. A value of 0.1 would plot 10 percent of the total data points over the violin.
#' @param title The title of the plot.
#' @param pal QQQ IS IT POSSIBLE TO PASS A USER-DEFINED PALETTE HERE? IF SO WE SHOULD UPDATE ALL PALETTE RELATED PARAMS. A custom palette from `ArchRPalettes` for discrete coloring.
#' @export
ggViolin <- function(
  x = NULL, 
  y = NULL, 
  xlabel = NULL, 
  ylabel = NULL, 
  xOrder = NULL,
  addPoints = FALSE,
  size = 1,  
  baseSize = 6, 
  ratioYX = NULL,
  sampleRatio = 0.1, 
  title = "", 
  pal = paletteDiscrete(values=x, set = "stallion"),
  ...
  ){

  .validInput(input = x, name = "x", valid = c("character"))
  .validInput(input = y, name = "y", valid = c("numeric"))
  .validInput(input = xlabel, name = "xlabel", valid = c("character", "null"))
  .validInput(input = ylabel, name = "ylabel", valid = c("character", "null"))
  .validInput(input = xOrder, name = "xOrder", valid = c("character", "null"))
  .validInput(input = addPoints, name = "addPoints", valid = c("boolean"))
  .validInput(input = size, name = "size", valid = c("numeric"))
  .validInput(input = baseSize, name = "baseSize", valid = c("numeric"))
  .validInput(input = ratioYX, name = "ratioYX", valid = c("numeric", "null"))
  .validInput(input = sampleRatio, name = "sampleRatio", valid = c("numeric"))
  .validInput(input = pal, name = "pal", valid = c("character"))
  
  names(y) <- x
  me = round(mean(stats::aggregate(y ~ names(y), FUN = mean)[, 2]), 2)
  sd = round(sd(stats::aggregate(y ~ names(y), FUN = mean)[, 2]), 2)
  min = round(min(y), 2)
  max = round(max(y), 2)
  df <- data.frame(x, y)

  if(!is.null(xOrder)){
    if(!all(x %in% xOrder)){
      stop("Not all x values are present in xOrder!")
    }
  }else{
    xOrder <- gtools::mixedsort(unique(x))
  }

  df$x <- factor(df$x, xOrder)
  
  p <- ggplot(df, aes_string(x = "x", y = "y", color = "x")) + 
    geom_violin(aes_string(fill="x"), alpha = 0.35) +
    geom_boxplot(size = size, outlier.size = 0, outlier.stroke = 0, fill = NA) + 
    scale_color_manual(values = pal, guide = FALSE) + 
    scale_fill_manual(values = pal, guide = FALSE) + 
    theme_ArchR(xText90 = TRUE, baseSize = baseSize) +
    ggtitle(title)

  if(!is.null(ratioYX)){
    p <- p + coord_fixed(ratioYX, expand = TRUE)
  }

  if(addPoints){

    if(requireNamespace("ggrastr", quietly = TRUE)){
      .requirePackage("ggrastr")
      p <- p + ggrastr::geom_quasirandom_rast(data = df[sample(seq_len(nrow(df)), floor(nrow(df) * sampleRatio)),], alpha = 1, 
            aes(x = x, y = y, color = x, fill = x), 
            size = 0.5, dodge.width=1)
    }else{
      message("ggrastr is not available for rastr of points, continuing without points!")
    }

  }

  if (!is.null(xlabel)) {
    p <- p + xlab(xlabel)
  }
  
  if (!is.null(ylabel)) {
    p <- p + ylab(ylabel)
  }
  
  p <- p + theme(legend.position = "bottom")

  if(!is.null(ratioYX)){
    attr(p, "ratioYX") <- ratioYX
  }
  
  return(p)

}


#' A ggplot-based Hexplot wrapper function summary of points in a standardized manner
#'
#' This function will plot x,y coordinate values summarized in hexagons in a standardized manner
#'
#' @param x A numeric vector containing the x-axis values for each point.
#' @param y A numeric vector containing the y-axis values for each point.
#' @param color A numeric/categorical vector containing coloring information for each point.
#' @param pal A custom continuous palette from `ArchRPalettes` for coloration of hexes.
#' @param bins The number of bins to be used for plotting the hexplot. `bins` indicates the total number of hexagons that will fit within the surface area of the plot. 
#' @param xlim A numeric vector of two values indicating the lower and upper bounds of the x-axis on the plot.
#' @param ylim A numeric vector of two values indicating the lower and upper bounds of the y-axis on the plot.
#' @param extend A numeric value indicating the fraction to extend the x-axis and y-axis beyond the maximum and minimum values if `xlim` and `ylim` are not provided. For example, 0.05 will extend the x-axis and y-axis by 5 percent on each end.
#' @param xlabel The label to plot for the x-axis.
#' @param ylabel The label to plot for the y-axis.
#' @param title The title of the plot.
#' @param colorTitle The label to use for the legend corresponding to `color`.
#' @param baseSize The base font size to use in the plot.
#' @param ratioYX The aspect ratio of the x and y axes on the plot.
#' @param FUN The function to use for summarizing data into hexagons. Typically "mean" or something similar.
#' @param quantCut If this is not null, a quantile cut is performed to threshold the top and bottom of the distribution. This prevents skewed color scales caused by strong outliers. The format of this should be c(a,b) where a is the upper threshold and b is the lower threshold. For example, quantileCut = c(0.025,0.975) will take the top and bottom 2.5 percent of values and set them to the value of the 97.5th and 2.5th percentile values respectively.
#' @param addPoints QQQ WHY WOULD THIS EVER BE USED? ARE POINTS PLOTTED INSTEAD OF HEXES OR ON TOP OF HEXES? A boolean value indicating whether individual points should be shown on the hexplot.
#' @param ... additional params to pass
#' @export
ggHex <- function(
  x = NULL, 
  y = NULL, 
  color = NULL, 
  pal = paletteContinuous(set = "solar_extra"), 
  bins = 200,
  xlim = NULL, 
  ylim = NULL, 
  extend = 0.05, 
  xlabel = "x", 
  ylabel = "y",
  title = "", 
  colorTitle = "values", 
  baseSize = 6,
  ratioYX = 1, 
  FUN = "mean", 
  quantCut = c(0.01, 0.99),
  addPoints = FALSE,
  ...
  ){

    .validInput(input = x, name = "x", valid = c("numeric"))
    .validInput(input = y, name = "y", valid = c("numeric"))
    .validInput(input = color, name = "color", valid = c("numeric"))
    .validInput(input = pal, name = "pal", valid = c("character"))
    .validInput(input = bins, name = "bins", valid = c("integer"))
    .validInput(input = xlim, name = "xlim", valid = c("numeric", "null"))
    .validInput(input = ylim, name = "ylim", valid = c("numeric", "null"))
    .validInput(input = xlabel, name = "xlabel", valid = c("character"))
    .validInput(input = ylabel, name = "ylabel", valid = c("character"))
    .validInput(input = title, name = "title", valid = c("character"))
    .validInput(input = colorTitle, name = "colorTitle", valid = c("character", "null"))
    .validInput(input = baseSize, name = "baseSize", valid = c("numeric"))
    .validInput(input = ratioYX, name = "ratioYX", valid = c("numeric"))
    .validInput(input = FUN, name = "FUN", valid = c("character"))
    .validInput(input = quantCut, name = "quantCut", valid = c("numeric"))
    .validInput(input = addPoints, name = "addPoints", valid = c("boolean"))

    df <- data.frame(x = x, y = y)
    include <- which(is.finite(x) & is.finite(y))

    if(length(include) != length(x)){
      message("Some values are not finite! Excluding these points!")
      df <- df[include,]
      if(!is.null(color)){
        color <- color[include]
      }
    }
    df$color <- color

    if (is.null(xlim)) {
        xlim <- range(df$x) %>% extendrange(f = extend)
    }
    if (is.null(ylim)) {
        ylim <- range(df$y) %>% extendrange(f = extend)
    }
    ratioXY <- ratioYX * diff(xlim)/diff(ylim)

    p <- ggplot()

    if(addPoints){
      if(requireNamespace("ggrastr", quietly = TRUE)){
        .requirePackage("ggrastr")
        p <- p + geom_point_rast(data = df, aes(x=x,y=y), color = "lightgrey")
      }else{
        message("ggrastr is not available for rastr of points, continuing without points!")
      }
    }

    values <- ggplot_build(p + stat_summary_hex(data = df, aes(x=x,y=y,z=color), fun = FUN, bins = bins, color = NA))$data[[1]]$value

    limits <- quantile(values, c(min(quantCut), max(quantCut)), na.rm=TRUE)

    p <- p + stat_summary_hex(data = df, aes(x=x,y=y,z=color), fun = FUN, bins = bins, color = NA) +
        scale_fill_gradientn(
          colors = pal,
          limits = limits, 
          oob = scales::squish
        ) +
        xlab(xlabel) + 
        ylab(ylabel) + 
        ggtitle(title) +
        theme_ArchR(baseSize = baseSize) +
        coord_equal(ratio = ratioXY, xlim = xlim, ylim = ylim, expand = FALSE) +
        theme(legend.direction="horizontal", legend.box.background = element_rect(color = NA)) +
        labs(color = colorTitle)

    p <- p + theme(legend.position = "bottom")
    
    if(!is.null(ratioYX)){
      attr(p, "ratioYX") <- ratioYX
    }

    p

}

#' Align ggplot plots vertically or horizontally
#'
#' This function aligns ggplots vertically or horizontally
#'
#' @param ... All additional arguments will be interpreted as `ggplot2` plot objects and used if and only if `plotList` is `NULL`
#' @param plotList A list of `ggplot2` plot objects to be aligned.
#' @param sizes QQQ HOW DOES IT KNOW WHICH VALUE TO MATCH WITH EACH PLOT? A numeric vector or list of values indicating the relative size for each of the objects in `plotList` or supplied in `...`.
#' @param type A string indicating wheter vertical ("v") or horizontal ("h") alignment should be used for the multi-plot layout.
#' @param draw A boolean value indicating whether to draw the plot(s) (`TRUE`) or return a graphical object (`FALSE`).
#' @export
ggAlignPlots <- function(
  ..., 
  plotList = NULL, 
  sizes = NULL, 
  type = "v",  
  draw = TRUE
  ){
  
  .validInput(input = plotList, name = "plotList", valid = c("list", "null"))
  .validInput(input = sizes, name = "sizes", valid = c("numeric"))
  .validInput(input = type, name = "type", valid = c("character"))
  .validInput(input = draw, name = "draw", valid = c("boolean"))
  if(type %ni% c("v", "h")){
    stop("type must be v (vertical) or h (horizontal)!")
  }

  #http://stackoverflow.com/a/21503904

  .requirePackage("gtable")

  if(is.null(plotList)){
    plotList <- list(...)
  }

  ## test that only passing plots
  stopifnot(do.call(all, lapply(plotList, inherits, "gg")))

  gl <- lapply(plotList, ggplotGrob)

  #if ncols do not match fill with empty gtables_add_cols
  if(type == "v" | type == "vertical"){
    maxCol <- max(unlist(lapply(gl, ncol)))
    gl <- lapply(gl, function(x){
      while(ncol(x) < max(maxCol)){
        x <- gtable::gtable_add_cols(x, unit(1, "null"))
      }
      return(x)
    })
  }

  combined <- Reduce(function(x, y)
    if(type == "v" | type == "vertical"){
      gtable:::rbind_gtable(x,y,"first")
    }else{
      gtable:::cbind_gtable(x,y,"first")
    }, gl[-1], gl[[1]])

  if(type == "v" | type == "vertical"){
    combined$widths <- do.call(grid::unit.pmax, lapply(gl, "[[", "widths"))
    #remove vertical spaces from background layout
    combined$heights[combined$layout$t[grepl("background", combined$layout$name)][-1]] <- grid::unit(rep(0,length(combined$heights[combined$layout$t[grepl("background", combined$layout$name)][-1]])), "cm")
    if(!missing(sizes)){
      sList <- lapply(seq_along(gl), function(x){
        orig <- gl[[x]]$heights[gl[[x]]$layout$t[grepl("panel", gl[[x]]$layout$name)]]
        new <- rep(sizes[[x]]/length(orig),length(orig))
        return(new)
      })
      s <- grid::unit(unlist(sList), "null")
      combined$heights[combined$layout$t[grepl("panel", combined$layout$name)]] <- s
    }
  }else if(type == "h" | type == "horizontal"){
    combined$heights <- do.call(grid::unit.pmax, lapply(gl, "[[", "heights"))
    if(!missing(sizes)){
      sList <- lapply(seq_along(gl), function(x){
        orig <- gl[[x]]$widths[gl[[x]]$layout$l[grepl("panel", gl[[x]]$layout$name)]]
        new <- rep(sizes[[x]]/length(orig),length(orig))
        return(new)
      })
      s <- grid::unit(unlist(sList), "null")
      combined$widths[combined$layout$l[grepl("panel", combined$layout$name)]] <- s
    }
  }else{
    stop("Unrecognized type ", type)
  }

  if(draw){
    grid::grid.newpage()
    grid::grid.draw(combined)
  }else{
    combined    
  }
  
}

#' ggplot2 default theme for ArchR
#'
#' This function returns a ggplot2 theme that is black borded with black font.
#' 
#' @param color The color to be used for text, lines, ticks, etc for the plot.
#' @param baseSize The base font size to use in the plot.
#' @param baseLineSize The base line width (in points) to be used throughout the plot.
#' @param baseRectSize The base line width (in points) to use for rectangular boxes throughout the plot.
#' @param plotMarginCm The width in centimeters of the whitespace margin around the plot.
#' @param legendPosition The location to put the legend. Valid options are "bottom", "top", "left", and "right.
#' @param legendTextSize QQQ The base text size (in points) for the legend text.
#' @param axisTickCm The length in centimeters to be used for the axis ticks.
#' @param xText90 A boolean value indicating whether the x-axis text should be rotated 90 degrees counterclockwise.
#' @param yText90 A boolean value indicating whether the y-axis text should be rotated 90 degrees counterclockwise.
#' @export
theme_ArchR <- function(
  color = "black",
  baseSize = 6, 
  baseLineSize = 0.5,
  baseRectSize = 0.5,
  plotMarginCm = 1,
  legendPosition = "bottom",
  legendTextSize = 5,
  axisTickCm = 0.1,
  xText90 = FALSE,
  yText90 = FALSE,
  ...
  ){

  .validInput(input = color, name = "color", valid = c("character"))
  .validInput(input = baseSize, name = "baseSize", valid = c("numeric"))
  .validInput(input = baseLineSize, name = "baseLineSize", valid = c("numeric"))
  .validInput(input = baseRectSize, name = "baseRectSize", valid = c("numeric"))
  .validInput(input = plotMarginCm, name = "plotMarginCm", valid = c("numeric"))
  .validInput(input = legendPosition, name = "legendPosition", valid = c("character"))
  .validInput(input = legendTextSize, name = "legendTextSize", valid = c("numeric"))
  .validInput(input = axisTickCm, name = "axisTickCm", valid = c("numeric"))
  .validInput(input = xText90, name = "xText90", valid = c("boolean"))
  .validInput(input = yText90, name = "yText90", valid = c("boolean"))

  theme <- theme_bw() + theme(
      axis.text = element_text(color = color, size = baseSize), 
      axis.title = element_text(color = color, size = baseSize),
      title = element_text(color = color, size = baseSize),
      plot.margin = unit(c(plotMarginCm, plotMarginCm, plotMarginCm, plotMarginCm), "cm"),
      panel.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = color, size = (4/3) * baseRectSize * as.numeric(grid::convertX(grid::unit(1, "points"), "mm"))),
      axis.ticks.length = unit(axisTickCm, "cm"), 
      axis.ticks = element_line(color = color, size = baseLineSize * (4/3) * as.numeric(grid::convertX(grid::unit(1, "points"), "mm"))),
      legend.key = element_rect(fill = "transparent", colour = NA),
      legend.text = element_text(color = color, size = legendTextSize),
      legend.background = element_rect(fill = "transparent"),
      legend.box.background = element_rect(fill = "transparent"),
      legend.position = legendPosition,
      strip.text = element_text(size = baseSize, color="black"),
      plot.background = element_rect(fill = "transparent", color = NA)
    )
  
  if(xText90){
    theme <- theme %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }

  if(yText90){
    theme <- theme %+replace% theme(axis.text.y = element_text(angle = 90, vjust = 1))
  }

  return(theme)

}


