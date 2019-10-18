#' GG Plot One to One Heatscatter
#'
#' @param x x
#' @param y y
#' @param size geom_point size
#' @param alpha geom_point alpha
#' @param xlabel xlabel
#' @param ylabel ylabel
#' @param title ggtitle
#' @param min xmin quantile [0,1]
#' @param max xmax quantile [0,1]
#' @param plot_n number of points to plot
#' @param kernel_n n for MASS::kde2d default = 100
#' @param plot_n number of points to plot
#' @param baseSize base_font size
#' @param pal continuous color palette to use
#' @export
ggPoint <- function(x, y, color = NULL, discrete = TRUE, discreteSet = "stallion", 
    labelMeans = FALSE, continuousSet = "solar_extra", pal = NULL, colorDensity = FALSE,
    size = 1, xlim = NULL, ylim = NULL, extend = 0.05, xlabel = "x", randomize = FALSE, seed = 1,
    ylabel = "y", title = "", alpha = 1, baseSize = 6, ratioYX = 1, 
    labelType = "ggrepel", bgColor = "white", fgColor = NULL, labelSize = 1.5,
    addFit = NULL, nullColor = "lightGrey", rastr = FALSE, dpi = 300){
    
    stopifnot(is.numeric(x))
    stopifnot(is.numeric(y))
    stopifnot(length(y)==length(x))

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
    if (is.null(xlim)) {
        xlim <- range(df$x) %>% extendrange(f = extend)
    }
    if (is.null(ylim)) {
        ylim <- range(df$y) %>% extendrange(f = extend)
    }
    ratioXY <- ratioYX * diff(xlim)/diff(ylim)

    #Plot
    library(ggplot2)

    if (is.null(color) & !colorDensity) {

      p <- ggplot(df[idx,], aes(x = x, y = y)) + coord_equal(ratio = ratioXY, xlim = xlim, 
              ylim = ylim, expand = F) + xlab(xlabel) + ylab(ylabel) + 
              ggtitle(title) + theme_ArchR(baseSize = baseSize)

      if(rastr){
        if(!requireNamespace("ggrastr", quietly = TRUE)){
          message("ggrastr is not available for rastr of points, continuing without rastr!")
          p <- p + geom_point(size = size, alpha = alpha, color = nullColor)
        }else{
          .requirePackage("ggrastr")
          p <- p + geom_point_rast(size = size, raster.dpi = dpi, alpha = alpha, color = nullColor)
        }
      }else{
        p <- p + geom_point(size = size, alpha = alpha, color = nullColor)
      }
        
    }else {

        if(colorDensity){
          discrete <- FALSE
          df <- getDensity(x, y, n = 100, sample = NULL) #change
          df <- df[order(df$density), ,drop=FALSE]
          df$color <- df$density
        }else if(discrete){
          stopifnot(length(color) == nrow(df))
            df$color <- factor(color, levels = sort(unique(color)))
        }else {
          stopifnot(length(color) == nrow(df))
            df$color <- color
        }
        p <- ggplot(df[idx,], aes(x = x, y = y, color = color)) +  coord_equal(ratio = ratioXY, xlim = xlim, 
            ylim = ylim, expand = F) + xlab(xlabel) + ylab(ylabel) + 
            ggtitle(title) + theme_ArchR(baseSize = baseSize) +
            theme(legend.direction="horizontal" , legend.box.background = element_rect(color = NA))

      if(rastr){
        
        if(!requireNamespace("ggrastr", quietly = TRUE)){
          message("ggrastr is not available for rastr of points, continuing without rastr!")
          p <- p + geom_point(size = size, alpha = alpha)
        }else{
          .requirePackage("ggrastr")
          p <- p + geom_point_rast(size = size, raster.dpi = dpi, alpha = alpha)
        }          
      
      }else{

          p <- p + geom_point(size = size, alpha = alpha)

      }

      if (discrete) {
          
          if (!is.null(pal)) {
              p <- p + scale_color_manual(values = pal)
          }else {
              p <- p + scale_color_manual(values = paletteDiscrete(set = discreteSet, 
                values = sort(unique(color))))
          }

          if (labelMeans) {
              dfMean <- split(df, df$color) %>% lapply(., function(x) {
                data.frame(x = mean(x[, 1]), y = mean(x[, 2]), 
                  color = x[1, 3])
              }) %>% Reduce("rbind", .)

              if(tolower(labelType) == "repel" | tolower(labelType) == "ggrepel"){

                if(!is.null(fgColor)){
                  p <- p + ggrepel::geom_label_repel(data = dfMean, aes(x, y, label = color), color = fgColor, size = labelSize)
                }else{
                  p <- p + ggrepel::geom_label_repel(data = dfMean, aes(x, y, label = color), size = labelSize)
                }

              }else if(tolower(labelType) == "shadow" | tolower(labelType) == "shadowtext"){

                if(!is.null(fgColor)){
                  p <- p + geom_shadowtext(data = dfMean, aes(x, y, label = color), color = fgColor, bg.colour = bgColor, size = labelSize)
                }else{
                  p <- p + geom_shadowtext(data = dfMean, aes(x, y, label = color), bg.colour = bgColor, size = labelSize)
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
        p <- p + geom_smooth(data = df, aes(color = NULL), method = addFit, 
            color = "black") + ggtitle(paste0(title, "\nPearson = ", 
            round(cor(df$x, df$y), 3), "\nSpearman = ", round(cor(df$x, 
                df$y, method = "spearman"), 3)))
    }
    return(p)
}

#' GG Plot One to One Heatscatter
#'
#' @param x x
#' @param y y
#' @param size geom_point size
#' @param alpha geom_point alpha
#' @param xlabel xlabel
#' @param ylabel ylabel
#' @param title ggtitle
#' @param min xmin quantile [0,1]
#' @param max xmax quantile [0,1]
#' @param nPlot number of points to plot
#' @param nKernel n for MASS::kde2d default = 100
#' @param baseSize base_font size default is 12
#' @param pal continuous color palette to use
#' @export
#'
ggOneToOne <- function (x, y, nPlot = 100 * 10^3, 
          nKernel = 100, size = 2, 
          xlabel = "x", ylabel = "y", title = "Sample Correlation", 
          min = 0.1, max = 0.9999, 
          densityMax = 0.95, extend = 0.05, 
          alpha = 1, baseSize = 12, 
          pal = paletteContinuous(set = "viridis")){
  
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
  df <- getDensity(x, y, n = nKernel, sample = nPlot) #change
  df <- df[order(df[, "density"]), ]
  
  #GGPlot
  message("plotting...")
  gg <- ggPlotPoint(
      x = df$x, 
      y = df$y, 
      color = df$density, 
      pal = pal,
      xlabel = xlabel,
      ylabel = ylabel,
      discrete = FALSE, 
      xlim = lim, 
      ylim = lim, 
      size = size, 
      alpha = alpha, 
      title = title, 
      baseSize = baseSize
    ) + geom_abline(slope = 1, intercept = 0, lty = "dashed")
  return(gg)
}

#modified from http://slowkow.com/notes/ggplot2-color-by-density/
getDensity <- function(x, y, n = 100, sample = NULL, densityMax = 0.95){
  df <- data.frame(x=x,y=y)
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  df$density <- dens$z[ii]
  df$density[df$density > quantile(unique(df$density),densityMax)] <- quantile(unique(df$density),densityMax) #make sure the higher end doesnt bias colors
  if(!is.null(sample)){
    df <- nSample(df,sample,type="r")
  }
  return(df)
}

#' GG Violin Plot
#' 
#' @param x categorical values to each y value
#' @param y numeric values
#' @param xlabel xlabel
#' @param ylabel ylabel
#' @param base_size base_size of theme
#' @param size size of barplot lines
#' @param pal color palette see paletteDiscrete for examples
#' @export
#'
ggViolin <- function (x, y, base_size = 12, xlabel = NULL, ylabel = NULL, points = FALSE, baseSize = 6, ratioYX = 1,
  sampleRatio = 0.1, size = 1, title = "", pal = paletteDiscrete(values=x, set = "stallion")) {
  stopifnot(!is.numeric(x))
  stopifnot(is.numeric(y))
  names(y) <- x
  me = round(mean(stats::aggregate(y ~ names(y), FUN = mean)[, 2]), 2)
  sd = round(sd(stats::aggregate(y ~ names(y), FUN = mean)[, 2]), 2)
  min = round(min(y), 2)
  max = round(max(y), 2)
  df <- data.frame(x, y)
  df$x <- factor(df$x, gtools::mixedsort(unique(paste0(df$x))))
  p <- ggplot(df, aes_string(x = "x", y = "y", color = "x")) + coord_fixed(ratioYX, expand = TRUE) +
    geom_violin(aes_string(fill="x"), alpha = 0.35) +
    geom_boxplot(size = size, outlier.size = 0, outlier.stroke = 0, fill = NA) + 
    scale_color_manual(values = pal, guide = FALSE) + 
    scale_fill_manual(values = pal, guide = FALSE) + 
    theme_ArchR(xText90 = TRUE, baseSize = baseSize) +
    ggtitle(title)

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
  return(p)
}


#' GG Violin Plot
#' 
#' @param x categorical values to each y value
#' @param y numeric values
#' @param xlabel xlabel
#' @param ylabel ylabel
#' @param base_size base_size of theme
#' @param size size of barplot lines
#' @param pal color palette see paletteDiscrete for examples
#' @export
#'
ggHex <- function(x, y, color, extend = 0.05, ratioYX = 1, xlim = NULL, ylim = NULL,
  bins = 150, pal = paletteContinuous(set = "solar_extra"), title = "", baseSize = 12,
    xlabel = "x" , ylabel ="y", fun = "median", ...){

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

    p <- ggplot() +
        stat_summary_hex(data = df, aes(x=x,y=y,z=color), fun = fun, bins = bins, color = NA) +
        scale_fill_gradientn(colors = pal) +
        xlab(xlabel) + 
        ylab(ylabel) + 
        ggtitle(title) +
        theme_ArchR(baseSize = baseSize) +
        coord_equal(ratio = ratioXY, xlim = xlim, ylim = ylim, expand = FALSE) +
        theme(legend.direction="horizontal", legend.box.background = element_rect(color = NA))
    
    p

}


#' GG Plot One to One Heatscatter
#'
#' @param x x
#' @param y y
#' @param size geom_point size
#' @param alpha geom_point alpha
#' @param xlabel xlabel
#' @param ylabel ylabel
#' @param title ggtitle
#' @param min xmin quantile [0,1]
#' @param max xmax quantile [0,1]
#' @param plot_n number of points to plot
#' @param kernel_n n for MASS::kde2d default = 100
#' @param plot_n number of points to plot
#' @param baseSize base_font size
#' @param pal continuous color palette to use
#' @export
ggLine <- function(x, y, color = NULL, discrete = TRUE, discreteSet = "stallion", 
    continuousSet = "solar_extra", pal = NULL, size = 1, xlim = NULL, ylim = NULL, 
    extend = 0.05, xlabel = "x", ylabel = "y", title = "", 
    alpha = 1, baseSize = 6, ratioYX = 1, 
    nullColor = "lightGrey"){
    
    stopifnot(is.numeric(x))
    stopifnot(is.numeric(y))
    stopifnot(length(y)==length(x))

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
    if (is.null(xlim)) {
        xlim <- range(df$x) %>% extendrange(f = extend)
    }
    if (is.null(ylim)) {
        ylim <- range(df$y) %>% extendrange(f = extend)
    }
    ratioXY <- ratioYX * diff(xlim)/diff(ylim)

    #Plot
    library(ggplot2)

    if (is.null(color)) {

      p <- ggplot(df, aes(x = x, y = y)) + coord_equal(ratio = ratioXY, xlim = xlim, 
              ylim = ylim, expand = F) + xlab(xlabel) + ylab(ylabel) + 
              ggtitle(title) + theme_ArchR(baseSize = baseSize)

      p <- p + geom_line(size = size, alpha = alpha, color = nullColor)
        
    }else {

        if(discrete){
          stopifnot(length(color) == nrow(df))
            df$color <- factor(color, levels = sort(unique(color)))
        }else {
          stopifnot(length(color) == nrow(df))
            df$color <- color
        }
        p <- ggplot(df, aes(x = x, y = y, color = color)) +  
            coord_equal(ratio = ratioXY, xlim = xlim, 
            ylim = ylim, expand = F) + xlab(xlabel) + ylab(ylabel) + 
            ggtitle(title) + theme_ArchR(baseSize = baseSize) +
            theme(legend.direction="horizontal" , legend.box.background = element_rect(color = NA))

        p <- p + geom_line(size = size, alpha = alpha)

        if (discrete) {
            if (!is.null(pal)) {
                p <- p + scale_color_manual(values = pal)
            }else {
                p <- p + scale_color_manual(values = paletteDiscrete(set = discreteSet, 
                  values = sort(unique(color))))
            }

        }else {
            if (!is.null(pal)) {
                p <- p + scale_colour_gradientn(colors = pal)
            }else {
                p <- p + scale_colour_gradientn(colors = paletteContinuous(set = continuousSet))
            }
        }
    }

    return(p)
}

#' Align GG Plots
#' @param ... ggplots
#' @param sizes sizes are a vector or list of values for each ggplot ie c(1,1) for two plots
#' @param type v,vertical or h,horizontal
#' @export
#'
ggAlignPlots <- function(..., sizes, type = "v", plotList = NULL, grobList = NULL){
  
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
  grid::grid.newpage()
  grid::grid.draw(combined)

}

#' ggplot2 default theme for ArchR
#'
#' This function returns a ggplot2 theme that is black borded with black font.
#' 
#' @param color color of theme
#' @param base_size is the size of the font for the axis text and title
#' @param base_family is family for font
#' @param base_line_size is the size of line
#' @param base_rect_size is the size of rectangle boxes
#' @param plot_margin_cm plot margin in cm
#' @param legend_position where is the legend default bottom
#' @param legend_text_size 0.75*base_size
#' @param axis_tick_length_cm axis tick length in cm
#' @param rotate_x_axis_text_90 rotate x axis text 90 degrees
#' @param rotate_y_axis_text_90 rotate y axis text 90 degrees
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







