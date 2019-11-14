#' Ggplot Points in a standardized manner
#'
#' This function will plot x,y coordinates in a standardized manner
#'
#' @param x x vector of data to be plot
#' @param y y vector of data to be plot
#' @param color color vector of data to be plot (must be same length as x,y)
#' @param discrete discrete is color discrete?
#' @param discreteSet default color ArchR_palette for discrete colors
#' @param continuousSet default color ArchR_palette for continuous colors
#' @param pal custom palette option
#' @param defaultColor default color 
#' @param colorDensity color x,y coordinates by relative density
#' @param size size of points
#' @param xlim xlimits for plot
#' @param ylim ylimits for plot
#' @param extend extend limits by this proportion if not set by xlim or ylim
#' @param xlabel label for x-axis
#' @param ylabel label for y-axis
#' @param title title of plot
#' @param randomize randomize x,y plotting order
#' @param seed seed for randomizing points
#' @param colorTitle title for legend corresponding to color
#' @param colorOrder order of colors for plotting/palettes
#' @param alpha alpha transperancy for points
#' @param baseSize baseSize for fonts
#' @param labelMeans label group means by color
#' @param labelType label type to add (ggrepel or shadowtext)
#' @param fgColor foreground color of labels
#' @param bgColor background color of labels
#' @param labelSize label size to be created
#' @param addFit add an x,y fit line using geom_smooth
#' @param ratioYX ratio of y to x in plot
#' @param rastr rastr points
#' @param dpi rastr resolution
#' @param ... additional params to pass
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
    ...){
    
    stopifnot(is.numeric(x))
    stopifnot(is.numeric(y))
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

    p <- p + theme(legend.position = "bottom")

    return(p)

}

#' GG Plot One to One Heatscatter
#'
#' @param x x vector of data to be plot
#' @param y y vector of data to be plot
#' @param size plot point size
#' @param alpha alpha transperancy of points
#' @param xlabel xlabel
#' @param ylabel ylabel
#' @param title ggtitle
#' @param min x/y min quantile [0,1]
#' @param max x/y max quantile [0,1]
#' @param nPlot number of points to plot (correlations computed prior)
#' @param nKernel n for MASS::kde2d default = 100
#' @param baseSize base size of fonts in plot
#' @param pal color palette to use for density
#' @export
#'
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
  pal = paletteContinuous(set = "viridis"),
  ...){
  
  #Check is Numeric
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(y))
  
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

#' GG Violin Plot
#' 
#' @param x categorical values to each y value
#' @param y numeric values
#' @param xlabel xlabel
#' @param ylabel ylabel
#' @param xOrder custom order of x for plotting
#' @param points add points using ggrastr geom_quasirandom?
#' @param size size of barplot lines
#' @param baseSize base size of fonts in plot
#' @param pal color palette see paletteDiscrete for examples
#' @export
#'
ggViolin <- function(
  x = NULL, 
  y = NULL, 
  xlabel = NULL, 
  ylabel = NULL, 
  xOrder = NULL,
  points = FALSE,
  size = 1,  
  baseSize = 6, 
  ratioYX = NULL,
  sampleRatio = 0.1, 
  title = "", 
  pal = paletteDiscrete(values=x, set = "stallion"),
  ...){

  stopifnot(!is.numeric(x))
  stopifnot(is.numeric(y))
  
  names(y) <- x
  me = round(mean(stats::aggregate(y ~ names(y), FUN = mean)[, 2]), 2)
  sd = round(sd(stats::aggregate(y ~ names(y), FUN = mean)[, 2]), 2)
  min = round(min(y), 2)
  max = round(max(y), 2)
  df <- data.frame(x, y)

  if(!is.null(xOrder)){
    if(!all(x %in% xOrder)){
      stop("Not x colors are in xOrder!")
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

  if(points){

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

  return(p)

}


#' Ggplot Hexplot summary of points in a standardized manner
#'
#' This function will plot x,y coordinates values summarized in hexagons in a standardized manner
#'
#' @param x x vector of data to be plot
#' @param y y vector of data to be plot
#' @param color color vector of values to be plot (must be same length as x,y)
#' @param pal custom palette option
#' @param bins number of bins for hexplot
#' @param xlim xlimits for plot
#' @param ylim ylimits for plot
#' @param extend extend limits by this proportion if not set by xlim or ylim
#' @param xlabel label for x-axis
#' @param ylabel label for y-axis
#' @param title title of plot
#' @param colorTitle title for legend corresponding to color
#' @param baseSize baseSize for fonts
#' @param ratioYX ratio of y to x in plot
#' @param FUN function for summarizing hexagons
#' @param ... additional params to pass
#' @export
ggHex <- function(
  x = NULL, 
  y = NULL, 
  color = NULL, 
  pal = paletteContinuous(set = "solar_extra"), 
  bins = 150,
  xlim = NULL, 
  ylim = NULL, 
  extend = 0.05, 
  xlabel = "x", 
  ylabel = "y",
  title = "", 
  colorTitle = "values", 
  baseSize = 6,
  ratioYX = 1, 
  FUN = "median", 
  addPoints = FALSE,
  ...){

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

    p <- p + stat_summary_hex(data = df, aes(x=x,y=y,z=color), fun = FUN, bins = bins, color = NA) +
        scale_fill_gradientn(colors = pal) +
        xlab(xlabel) + 
        ylab(ylabel) + 
        ggtitle(title) +
        theme_ArchR(baseSize = baseSize) +
        coord_equal(ratio = ratioXY, xlim = xlim, ylim = ylim, expand = FALSE) +
        theme(legend.direction="horizontal", legend.box.background = element_rect(color = NA)) +
        labs(color = colorTitle)

    p <- p + theme(legend.position = "bottom")
    
    p

}

#' Align Ggplots vertically or horizontally
#'
#' This function aligns ggplots vertically or horizontally
#'
#' @param ... ggplots
#' @param sizes sizes are a vector or list of values for each ggplot ie c(1,1) for two plots
#' @param type v,vertical or h,horizontal
#' @param plotList add a list of plots to be aligned
#' @param grobList add a list of grobs to be aligned
#' @export
#'
ggAlignPlots <- function(..., sizes, type = "v", plotList = NULL, grobList = NULL, draw = TRUE){
  
  #http://stackoverflow.com/a/21503904

  .requirePackage("gtable")

  if(is.null(grobList)){

    if(is.null(plotList)){
      plotList <- list(...)
    }

    ## test that only passing plots
    stopifnot(do.call(all, lapply(plotList, inherits, "gg")))

    gl <- lapply(plotList, ggplotGrob)

  }else{
    
    gl <- grobList
    rm(grobList)
    gc()

  }

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
#' @param color color of theme
#' @param baseSize is the size of the font for the axis text and title
#' @param baseFamily is family for font
#' @param baseLineSize is the size of line
#' @param baseRectSize is the size of rectangle boxes
#' @param plotMarginCm plot margin in cm
#' @param legendPosition where is the legend default bottom
#' @param legendTextSize 0.75*base_size
#' @param axisTickCm axis tick length in cm
#' @param xText90 rotate x axis text 90 degrees
#' @param yText90 rotate y axis text 90 degrees
#' @export
theme_ArchR <- function(
  color = "black",
  baseSize = 6, 
  baseFamily = "",
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





# #' GG Plot One to One Heatscatter
# #'
# #' @param x x
# #' @param y y
# #' @param size geom_point size
# #' @param alpha geom_point alpha
# #' @param xlabel xlabel
# #' @param ylabel ylabel
# #' @param title ggtitle
# #' @param min xmin quantile [0,1]
# #' @param max xmax quantile [0,1]
# #' @param plot_n number of points to plot
# #' @param kernel_n n for MASS::kde2d default = 100
# #' @param plot_n number of points to plot
# #' @param baseSize base_font size
# #' @param pal continuous color palette to use
# #' @export
# ggLine <- function(x, y, color = NULL, discrete = TRUE, discreteSet = "stallion", 
#     continuousSet = "solar_extra", pal = NULL, size = 1, xlim = NULL, ylim = NULL, 
#     extend = 0.05, xlabel = "x", ylabel = "y", title = "", 
#     alpha = 1, baseSize = 6, ratioYX = 1, 
#     nullColor = "lightGrey"){
    
#     stopifnot(is.numeric(x))
#     stopifnot(is.numeric(y))
#     stopifnot(length(y)==length(x))

#     df <- data.frame(x = x, y = y)
#     include <- which(is.finite(x) & is.finite(y))
#     if(length(include) != length(x)){
#       message("Some values are not finite! Excluding these points!")
#       df <- df[include,]
#       x <- x[include]
#       y <- y[include]
#       if(!is.null(color)){
#         color <- color[include]
#       }
#     }
#     if (is.null(xlim)) {
#         xlim <- range(df$x) %>% extendrange(f = extend)
#     }
#     if (is.null(ylim)) {
#         ylim <- range(df$y) %>% extendrange(f = extend)
#     }
#     ratioXY <- ratioYX * diff(xlim)/diff(ylim)

#     #Plot
#     library(ggplot2)

#     if (is.null(color)) {

#       p <- ggplot(df, aes(x = x, y = y)) + coord_equal(ratio = ratioXY, xlim = xlim, 
#               ylim = ylim, expand = F) + xlab(xlabel) + ylab(ylabel) + 
#               ggtitle(title) + theme_ArchR(baseSize = baseSize)

#       p <- p + geom_line(size = size, alpha = alpha, color = nullColor)
        
#     }else {

#         if(discrete){
#           stopifnot(length(color) == nrow(df))
#             df$color <- factor(color, levels = sort(unique(color)))
#         }else {
#           stopifnot(length(color) == nrow(df))
#             df$color <- color
#         }
#         p <- ggplot(df, aes(x = x, y = y, color = color)) +  
#             coord_equal(ratio = ratioXY, xlim = xlim, 
#             ylim = ylim, expand = F) + xlab(xlabel) + ylab(ylabel) + 
#             ggtitle(title) + theme_ArchR(baseSize = baseSize) +
#             theme(legend.direction="horizontal" , legend.box.background = element_rect(color = NA))

#         p <- p + geom_line(size = size, alpha = alpha)

#         if (discrete) {
#             if (!is.null(pal)) {
#                 p <- p + scale_color_manual(values = pal)
#             }else {
#                 p <- p + scale_color_manual(values = paletteDiscrete(set = discreteSet, 
#                   values = sort(unique(color))))
#             }

#         }else {
#             if (!is.null(pal)) {
#                 p <- p + scale_colour_gradientn(colors = pal)
#             }else {
#                 p <- p + scale_colour_gradientn(colors = paletteContinuous(set = continuousSet))
#             }
#         }
#     }

#     return(p)
# }



