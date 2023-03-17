
shinyServer <- function(input,output, session){
  
  
  # EMBEDS ------------------------------------------------------------------------------------
  
  plot1 <- shiny::reactive({
    
    # availableMatrices <- getAvailableMatrices(ArchRProj)
    
    if(input$matrix_EMBED1_forComparison %in% availableMatrices){
      mat <- availableMatrices[ availableMatrices  %in% input$matrix_EMBED1_forComparison]
      
      p_empty <-  ggplot2::ggplot()  +
        ggplot2::xlab("Dimension 1") + ggplot2::ylab("Dimension 2") + ggplot2::theme_bw(base_size=10)+
        ggplot2::ggtitle(base::paste0("EMBED of IterativeLSI colored by \n", mat,": ",input$EMBED1_forComparison)) +
        ggplot2::theme(
          panel.background = ggplot2::element_rect(fill='transparent'), #transparent panel bg
          plot.background = ggplot2::element_rect(fill='transparent', color=NA), #transparent plot bg
          panel.grid.major = ggplot2::element_blank(), #remove major gridlines
          panel.grid.minor = ggplot2::element_blank(), #remove minor gridlines
          legend.background = ggplot2::element_rect(fill='transparent'), #transparent legend bg
          legend.box.background = ggplot2::element_rect(fill='transparent'),  axis.title=ggplot2::element_text(size=14),           
          plot.title = ggplot2::element_text(size=16)
        )
      
      plotfunctions::emptyPlot(0,0, axes=FALSE)
      legend_plot <- plotfunctions::gradientLegend(valRange=base::c(base::scale()[[mat]][,input$EMBED1_forComparison][1],base::scale()[[mat]][,input$EMBED1_forComparison][2]), 
                                    color = raster::color()[[mat]], pos=.5, side=1)
      
      
      p <- rhdf5::h5read(paste0(subOutDir,"/",mat,"_plotBlank72.h5"), base::paste0(mat, "/", input$EMBED1_forComparison))
      temp_jpg <- t(base::matrix(farver::decode_native(p), nrow = 216))
      
      last_plot <- cowplot::ggdraw() + cowplot::draw_image(temp_jpg, x = 0, y = 0, scale = 0.85) + 
        cowplot::draw_plot(p_empty, scale = 0.8) +
        cowplot::draw_text("Log2(+1)", x = 0.35, y = 0.05, size = 12)
      
      base::print(last_plot, vp=grid::viewport(0.5, 0.6, 1, 1))        
    }else{
      p_empty <- ggplot2::ggplot()  +
        ggplot2::xlab("Dimension 1") + ggplot2::ylab("Dimension 2") + ggplot2::theme_bw(base_size=10)+
        ggplot2::ggtitle(input$matrix_EMBED1_forComparison) +
        ggplot2::theme(
          panel.background = ggplot2::element_rect(fill='transparent'), #transparent panel bg
          plot.background = ggplot2::element_rect(fill='transparent', color=NA), #transparent plot bg
          panel.grid.major = ggplot2::element_blank(), #remove major gridlines
          panel.grid.minor = ggplot2::element_blank(), #remove minor gridlines
          legend.background = ggplot2::element_rect(fill='transparent'), #transparent legend bg
          legend.box.background = ggplot2::element_rect(fill='transparent'),  axis.title=ggplot2::element_text(size=14),
          plot.title = ggplot2::element_text(size=16)
        )
      
      plotfunctions::emptyPlot(0,0, axes=FALSE)
      
      graphics::legend('bottom', legend=embed_legend[[1]],
             pch=15, col = color_embeddings[[1]],
             horiz = FALSE, x.intersp = 1, text.width=0.35,
             cex = 0.7, bty="n", ncol = 4)
      
      p <- rhdf5::h5read(base::paste0(subOutDir,"/mainEmbeds.h5"), input$matrix_EMBED1_forComparison)# input$EMBED1_forComparison))
      temp_jpg <- t(base::matrix(farver::decode_native(p), nrow = 216))
      
      last_plot <- cowplot::ggdraw() + cowplot::draw_image(temp_jpg, x = 0, y = 0, scale = 0.7) +
        cowplot::draw_plot(p_empty, scale = 0.7) +
        cowplot::draw_text("color", x = 0.1, y = 0.135, size = 12)
      
      base::print(last_plot, vp=grid::viewport(0.5, 0.6, 1, 1))
    }
    
  })
  
  plot2 <- shiny::reactive({
    
    # availableMatrices <- getAvailableMatrices(ArchRProj)
    
    if(input$matrix_EMBED2_forComparison %in% availableMatrices){
      mat <- availableMatrices[ availableMatrices  %in% input$matrix_EMBED2_forComparison]
      
      p_empty <-  ggplot2::ggplot()  +
        ggplot2::xlab("Dimension 1") + ggplot2::ylab("Dimension 2") + ggplot2::theme_bw(base_size=10)+
        ggplot2::ggtitle(base::paste0("EMBED of IterativeLSI colored by \n", mat,": ",input$EMBED2_forComparison)) +
        ggplot2::theme(
          panel.background = ggplot2::element_rect(fill='transparent'), #transparent panel bg
          plot.background = ggplot2::element_rect(fill='transparent', color=NA), #transparent plot bg
          panel.grid.major = ggplot2::element_blank(), #remove major gridlines
          panel.grid.minor = ggplot2::element_blank(), #remove minor gridlines
          legend.background = ggplot2::element_rect(fill='transparent'), #transparent legend bg
          legend.box.background = ggplot2::element_rect(fill='transparent'),  axis.title=ggplot2::element_text(size=14),           
          plot.title = ggplot2::element_text(size=16)
        )
      
      plotfunctions::emptyPlot(0,0, axes=FALSE)
      legend_plot <- plotfunctions::gradientLegend(valRange=base::c(base::scale()[[mat]][,input$EMBED2_forComparison][1],base::scale()[[mat]][,input$EMBED2_forComparison][2]), 
                                    color = raster::color()[[mat]], pos=.5, side=1)
      
      p <- h5read(base::paste0(subOutDir,"/",mat,"_plotBlank72.h5"), base::paste0(mat, "/", input$EMBED2_forComparison))
      temp_jpg <- t(base::matrix(farver::decode_native(p), nrow = 216))
      last_plot <- cowplot::ggdraw() + cowplot::draw_image(temp_jpg, x = 0, y = 0, scale = 0.85) + 
        cowplot::draw_plot(p_empty, scale = 0.8) +
        cowplot::draw_text("Log2(+1)", x = 0.35, y = 0.05, size = 12)
      
      base::print(last_plot, vp=grid::viewport(0.5, 0.6, 1, 1))        
    }else{
      
      
      p_empty <- ggplot2::ggplot()  +
        ggplot2::xlab("Dimension 1") + ggplot2::ylab("Dimension 2") + ggplot2::theme_bw(base_size=10)+
        ggplot2::ggtitle(input$matrix_EMBED2_forComparison) +
        ggplot2::theme(
          panel.background = ggplot2::element_rect(fill='transparent'), #transparent panel bg
          plot.background = ggplot2::element_rect(fill='transparent', color=NA), #transparent plot bg
          panel.grid.major = ggplot2::element_blank(), #remove major gridlines
          panel.grid.minor = ggplot2::element_blank(), #remove minor gridlines
          legend.background = ggplot2::element_rect(fill='transparent'), #transparent legend bg
          legend.box.background = ggplot2::element_rect(fill='transparent'),  axis.title=ggplot2::element_text(size=14),
          plot.title = ggplot2::element_text(size=16)
        )
      
      plotfunctions::emptyPlot(0,0, axes=FALSE)
      
      legend('bottom', legend=embed_legend[[1]],
             pch=15, col = color_embeddings[[1]],
             horiz = FALSE, x.intersp = 1, text.width=0.35,
             cex = 0.7, bty="n", ncol = 4)
      
      p <- h5read(base::paste0(subOutDir,"/mainEmbeds.h5"), input$matrix_EMBED2_forComparison)# input$EMBED2_forComparison))
      temp_jpg <- t(base::matrix(farver::decode_native(p), nrow = 216))
      
      last_plot <- cowplot::ggdraw() + cowplot::draw_image(temp_jpg, x = 0, y = 0, scale = 0.7) +
        cowplot::draw_plot(p_empty, scale = 0.7) +
        cowplot::draw_text("color", x = 0.1, y = 0.135, size = 12)
      
      base::print(last_plot, vp=grid::viewport(0.5, 0.6, 1, 1))
      
    }
    
  })
  
  
  #Output Handler: Downloads EMBEDS
  output$download_EMBED1 <- shiny::downloadHandler(
    filename <- function(){
      base::paste0("EMBED-",base::paste(input$matrix_EMBED1_forComparison,input$EMBED1_forComparison,sep="-"),input$plot_choice_download_EMBED1)
    },
    content = function(file){
      if(input$plot_choice_download_EMBED1==".pdf")
      {grDevices::pdf(file = file,onefile=FALSE, width = input$EMBED1_plot_width, height = input$EMBED1_plot_height)}

      else if(input$plot_choice_download_EMBED1==".png")
      {grDevices::png(file = file, width = input$EMBED1_plot_width, height = input$EMBED1_plot_height,units="in",res=1000)}

      else
      {grDevices::tiff(file = file, width = input$EMBED1_plot_width, height = input$EMBED1_plot_height,units="in",res=1000)}
      
      plot1 = plot1()


      gridExtra::grid.arrange(plot1)
      grDevices::dev.off()
    }
  )
  
  output$download_EMBED2 <- shiny::downloadHandler(
    filename <- function(){
      base::paste0("EMBED-",base::paste(input$matrix_EMBED2_forComparison,input$EMBED2_forComparison,sep="-"),input$plot_choice_download_EMBED2)
    },
    content = function(file){

      if(input$plot_choice_download_EMBED2==".pdf")
      {grDevices::pdf(file = file,onefile=FALSE, width = input$EMBED2_plot_width, height = input$EMBED2_plot_height)}

      else if(input$plot_choice_download_EMBED2==".png")
      {grDevices::png(file = file, width = input$EMBED2_plot_width, height = input$EMBED2_plot_height,units="in",res=1000)}

      else
      {grDevices::tiff(file = file, width = input$EMBED2_plot_width, height = input$EMBED2_plot_height,units="in",res=1000)}


        plot2 <- plot2()

      gridExtra::grid.arrange(plot2)
      grDevices::dev.off()
    }
  )
  
  output$EMBED_plot_1 <- DT::renderDT(NULL)
  output$EMBED_plot_2 <- DT::renderDT(NULL)
  
  color <- shiny::reactive({base::readRDS(base::paste0(subOutDir,"/pal.rds"))})
  scale <- shiny::reactive({base::readRDS(base::paste0(subOutDir,"/scale.rds"))})
  
  #plot EMBED1
  output$EMBED_plot_1<- shiny::renderPlot({
    
    plot1()
    
  }, height = 450,width=450)
  
  # #plot EMBED2
  output$EMBED_plot_2<- shiny::renderPlot({
    
    plot2()
    
  } ,height = 450,width=450)

  #update EMBED dropdown based on selected Matrix--------------------------------
  
  #Update dropdown for EMBED1
  featureNames1 <- shiny::reactive({
    
    if(!(input$matrix_EMBED1_forComparison %in% groupBy)){
      # availableMatrices <- getAvailableMatrices(ArchRProj)
      matName <- availableMatrices[ availableMatrices  %in% input$matrix_EMBED1_forComparison]
      featureNames <- rhdf5::h5read(file = base::paste0(subOutDir, "/", matName, "_plotBlank72.h5"),
                             name = matName)
      Feature_dropdown1 = base::names(featureNames)
      return(Feature_dropdown1)
    }
    
  })

  shiny::observeEvent(input$matrix_EMBED1_forComparison,{
    
    if(!(input$matrix_EMBED1_forComparison %in% groupBy)){
  shiny::updateSelectizeInput(session, 'EMBED1_forComparison', label = 'Feature Name',
                           choices = base::sort(featureNames1()),
                           server = TRUE, selected = base::sort(featureNames1())[1])
      }
    })

  # })
  
  
  featureNames2 <- shiny::reactive({
    
    if(!(input$matrix_EMBED2_forComparison %in% groupBy)){
          
        # availableMatrices <- getAvailableMatrices(ArchRProj)
        matName <- availableMatrices[ availableMatrices  %in% input$matrix_EMBED2_forComparison]
        featureNames <- rdhf5::h5read(file = base::paste0(subOutDir, "/", matName, "_plotBlank72.h5"),
               name = matName)
                             
        Feature_dropdown2 = base::names(featureNames)
        return(Feature_dropdown2)
    
    }
    
  })
  
  shiny::observeEvent(input$matrix_EMBED2_forComparison,{
    if(!(input$matrix_EMBED2_forComparison %in% groupBy)){
    shiny::updateSelectizeInput(session, 'EMBED2_forComparison', label = 'Feature Name',
                         choices = base::sort(featureNames2()),
                         server = TRUE, selected = base::sort(featureNames2())[1])
    }
    })

  #Update dropdown for EMBED2
  # shiny::observeEvent(input$matrix_EMBED2_forComparison,{
  #   if(shiny::isolate(input$matrix_EMBED2_forComparison)=="Motif Matrix")
  #   {
  
  #     shiny::updateSelectizeInput(session, 'EMBED2_forComparison', label = 'Feature Name',
  #                          choices = base::sort(MM_dropdown),
  #                          server = TRUE, selected = base::sort(MM_dropdown)[2])
  #   }
  
  #   else if(shiny::isolate(input$matrix_EMBED2_forComparison)=="Gene Score Matrix")
  #   {
  #     shiny::updateSelectizeInput(session, 'EMBED2_forComparison', label = 'Feature Name',
  #                          choices = base::sort(GSM_dropdown),
  #                          server = TRUE, selected = base::sort(GSM_dropdown)[2])
  #   }
  #   else if(shiny::isolate(input$matrix_EMBED2_forComparison)=="Gene Integration Matrix")
  #   {
  
  #     shiny::updateSelectizeInput(session, 'EMBED2_forComparison', label = 'Feature Name',
  #                          choices = base::sort(GIM_dropdown),
  #                          server = TRUE,selected = base::sort(GIM_dropdown)[2])
  #   }
  
  # })

  # Plot Browser ----------------------------------------------------------------

  # Observe the inputs for ATAC-Seq Explorer
  shiny::observeEvent(input$range_min, {
    shiny::updateSliderInput(session, "range",
                      value = base::c(input$range_min, base::max(input$range)))
  })
  
  shiny::observeEvent(input$range_max, {
    shiny::updateSliderInput(session, "range",
                      value = base::c(input$range_min,input$range_max))
  })
  
  shiny::observeEvent(input$range , {
    
    shiny::updateNumericInput(session, "range_min", value = base::min(input$range))
    shiny::updateNumericInput(session, "range_max", value = base::max(input$range))
    
  }, priority = 200)
  
  # Output Handler:downloads file
  output$down <- shiny::downloadHandler(
    filename <- function(){
      base::paste0("ArchRBrowser-",input$gene_name,input$plot_choice_download_peakBrowser)
    },
    content = function(file){
      
      if(input$plot_choice_download_peakBrowser==".pdf")
      {grDevices::pdf(file = file,onefile=FALSE, width = input$plot_width, height = input$plot_height)}
      
      else if(input$plot_choice_download_peakBrowser==".png")
      {grDevices::png(file = file, width = input$plot_width, height = input$plot_height,units="in",res=1000)}
      
      else
      {grDevices::tiff(file = file, width = input$plot_width, height = input$plot_height,units="in",res=1000)}
      
      
        p_browser_atacClusters<- ArchR::plotBrowserTrack(
          ArchRProj = ArchRProj,
          ShinyArchR = TRUE,
          plotSummary = base::c("bulkTrack", input$selectPlotSummary),
          baseSize = 11,
          facetbaseSize = 11,
          groupBy = input$browserContent,
          geneSymbol = shiny::isolate(input$gene_name),
          upstream = -base::min(shiny::isolate(input$range))*1000,
          downstream = base::max(shiny::isolate(input$range))*1000,
          tileSize = shiny::isolate(input$tile_size),
          ylim =  base::c(0, shiny::isolate(input$ymax)),
          loops = ArchR::getCoAccessibility(ArchRProj)
          
        )[[input$gene_name]]
     
      
      gridExtra::grid.arrange(p_browser_atacClusters)  
      
      grDevices::dev.off()
    }
  )
  output$browser_atacClusters <- DT::renderDT(NULL)
  
  #handles error
  restartFN <- shiny::observeEvent(input$restartButton, {
    if (shiny::isolate(input$gene_name) == ""){
      
      output$browser_atacClusters <- shiny::renderPlot({
        p <- ggplot2::ggplot() +
          ggplot2::xlim(base::c(-5,5)) + ggplot2::ylim(base::c(-5,5)) +
          ggplot2::geom_text(size=20, ggplot2::aes(x = 0, y = 0, label = "Please supply\na valid gene name!")) + ggplot2::theme_void()
        base::print(p)
      })
    }else{
      
      # Plots scATACSeq clusters
      output$browser_atacClusters<- shiny::renderPlot({
        grid::grid.newpage()

          p_browser_atacClusters<- ArchR::plotBrowserTrack(
            ArchRProj = ArchRProj,
            ShinyArchR = TRUE,
            plotSummary = base::c("bulkTrack", input$selectPlotSummary),
            baseSize = 11,
            facetbaseSize = 11,
            groupBy =  input$browserContent,
            geneSymbol = shiny::isolate(input$gene_name),
            upstream = -base::min(shiny::isolate(input$range))*1000,
            downstream = base::max(shiny::isolate(input$range))*1000,
            tileSize = shiny::isolate(input$tile_size),
            ylim =  base::c(0, shiny::isolate(input$ymax)),
            loops = ArchR::getCoAccessibility(ArchRProj)
            
          )[[input$gene_name]]
        
        
        grid::grid.draw(p_browser_atacClusters)
        
      },height = 900)
      
    }
  })
}
z