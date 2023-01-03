
shinyServer <- function(input,output, session){
  
  
  # EMBEDS ------------------------------------------------------------------------------------
  
  plot1 <- reactive({
    
    availableMatrices <- getAvailableMatrices(ArchRProj)
    
    if(input$matrix_EMBED1_forComparison %in% availableMatrices){
      mat <- availableMatrices[ availableMatrices  %in% input$matrix_EMBED1_forComparison]
      
      p_empty <-  ggplot()  +
        xlab("Dimension 1") + ylab("Dimension 2") + theme_bw(base_size=10)+
        ggtitle(paste0("EMBED of IterativeLSI colored by \n", mat,": ",input$EMBED1_forComparison)) +
        theme(
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          panel.grid.major = element_blank(), #remove major gridlines
          panel.grid.minor = element_blank(), #remove minor gridlines
          legend.background = element_rect(fill='transparent'), #transparent legend bg
          legend.box.background = element_rect(fill='transparent'),  axis.title=element_text(size=14),           
          plot.title = element_text(size=16)
        )
      
      emptyPlot(0,0, axes=FALSE)
      legend_plot <- gradientLegend(valRange=c(scale()[[mat]][,input$EMBED1_forComparison][1],scale()[[mat]][,input$EMBED1_forComparison][2]), 
                                    color = color()[[mat]], pos=.5, side=1)
      
      
      p <- h5read(paste0("./",subOutputDir,"/",mat,"_plotBlank72.h5"),  paste0(mat, "/", input$EMBED1_forComparison))
      temp_jpg <- t(matrix(decode_native(p), nrow = 216))
      
      last_plot <- ggdraw() + draw_image(temp_jpg, x = 0, y = 0, scale = 0.85) + 
        draw_plot(p_empty, scale = 0.8) +
        draw_text("Log2(+1)", x = 0.35, y = 0.05, size = 12)
      
      print(last_plot, vp=viewport(0.5, 0.6, 1, 1))        
    }else{
      p_empty <-  ggplot()  +
        xlab("Dimension 1") + ylab("Dimension 2") + theme_bw(base_size=10)+
        ggtitle(input$matrix_EMBED1_forComparison) +
        theme(
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          panel.grid.major = element_blank(), #remove major gridlines
          panel.grid.minor = element_blank(), #remove minor gridlines
          legend.background = element_rect(fill='transparent'), #transparent legend bg
          legend.box.background = element_rect(fill='transparent'),  axis.title=element_text(size=14),
          plot.title = element_text(size=16)
        )
      
      emptyPlot(0,0, axes=FALSE)
      
      legend('bottom', legend=embed_legend[[1]],
             pch=15, col = color_embeddings[[1]],
             horiz = FALSE, x.intersp = 1, text.width=0.35,
             cex = 0.7, bty="n", ncol = 4)
      
      p <- h5read(paste0("./",subOutputDir,"/mainEmbeds.h5"), input$matrix_EMBED1_forComparison)# input$EMBED1_forComparison))
      temp_jpg <- t(matrix(decode_native(p), nrow = 216))
      
      last_plot <- ggdraw() + draw_image(temp_jpg, x = 0, y = 0, scale = 0.7) +
        draw_plot(p_empty, scale = 0.7) +
        draw_text("color", x = 0.1, y = 0.135, size = 12)
      
      print(last_plot, vp=viewport(0.5, 0.6, 1, 1))
    }
    
  })
  
  plot2 <- reactive({
    
    availableMatrices <- getAvailableMatrices(ArchRProj)
    
    if(input$matrix_EMBED2_forComparison %in% availableMatrices){
      mat <- availableMatrices[ availableMatrices  %in% input$matrix_EMBED2_forComparison]
      
      p_empty <-  ggplot()  +
        xlab("Dimension 1") + ylab("Dimension 2") + theme_bw(base_size=10)+
        ggtitle(paste0("EMBED of IterativeLSI colored by \n", mat,": ",input$EMBED2_forComparison)) +
        theme(
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          panel.grid.major = element_blank(), #remove major gridlines
          panel.grid.minor = element_blank(), #remove minor gridlines
          legend.background = element_rect(fill='transparent'), #transparent legend bg
          legend.box.background = element_rect(fill='transparent'),  axis.title=element_text(size=14),           
          plot.title = element_text(size=16)
        )
      
      emptyPlot(0,0, axes=FALSE)
      legend_plot <- gradientLegend(valRange=c(scale()[[mat]][,input$EMBED2_forComparison][1],scale()[[mat]][,input$EMBED2_forComparison][2]), 
                                    color = color()[[mat]], pos=.5, side=1)
      
      p <- h5read(paste0("./",subOutputDir,"/",mat,"_plotBlank72.h5"),  paste0(mat, "/", input$EMBED2_forComparison))
      temp_jpg <- t(matrix(decode_native(p), nrow = 216))
      last_plot <- ggdraw() + draw_image(temp_jpg, x = 0, y = 0, scale = 0.85) + 
        draw_plot(p_empty, scale = 0.8) +
        draw_text("Log2(+1)", x = 0.35, y = 0.05, size = 12)
      
      print(last_plot, vp=viewport(0.5, 0.6, 1, 1))        
    }else{
      
      
      p_empty <-  ggplot()  +
        xlab("Dimension 1") + ylab("Dimension 2") + theme_bw(base_size=10)+
        ggtitle(input$matrix_EMBED2_forComparison) +
        theme(
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          panel.grid.major = element_blank(), #remove major gridlines
          panel.grid.minor = element_blank(), #remove minor gridlines
          legend.background = element_rect(fill='transparent'), #transparent legend bg
          legend.box.background = element_rect(fill='transparent'),  axis.title=element_text(size=14),
          plot.title = element_text(size=16)
        )
      
      emptyPlot(0,0, axes=FALSE)
      
      legend('bottom', legend=embed_legend[[1]],
             pch=15, col = color_embeddings[[1]],
             horiz = FALSE, x.intersp = 1, text.width=0.35,
             cex = 0.7, bty="n", ncol = 4)
      
      p <- h5read(paste0("./",subOutputDir,"/mainEmbeds.h5"), input$matrix_EMBED2_forComparison)# input$EMBED2_forComparison))
      temp_jpg <- t(matrix(decode_native(p), nrow = 216))
      
      last_plot <- ggdraw() + draw_image(temp_jpg, x = 0, y = 0, scale = 0.7) +
        draw_plot(p_empty, scale = 0.7) +
        draw_text("color", x = 0.1, y = 0.135, size = 12)
      
      print(last_plot, vp=viewport(0.5, 0.6, 1, 1))
      
    }
    
  })
  
  
  #Output Handler: Downloads EMBEDS
  output$download_EMBED1<-downloadHandler(
    filename <- function(){
      paste0("EMBED-",paste(input$matrix_EMBED1_forComparison,input$EMBED1_forComparison,sep="-"),input$plot_choice_download_EMBED1)
    },
    content = function(file){
      if(input$plot_choice_download_EMBED1==".pdf")
      {pdf(file = file,onefile=FALSE, width = input$EMBED1_plot_width, height = input$EMBED1_plot_height)}

      else if(input$plot_choice_download_EMBED1==".png")
      {png(file = file, width = input$EMBED1_plot_width, height = input$EMBED1_plot_height,units="in",res=1000)}

      else
      {tiff(file = file, width = input$EMBED1_plot_width, height = input$EMBED1_plot_height,units="in",res=1000)}
      
      plot1 = plot1()


      grid.arrange(plot1)
      dev.off()
    }
  )
  
  output$download_EMBED2<-downloadHandler(
    filename <- function(){
      paste0("EMBED-",paste(input$matrix_EMBED2_forComparison,input$EMBED2_forComparison,sep="-"),input$plot_choice_download_EMBED2)
    },
    content = function(file){

      if(input$plot_choice_download_EMBED2==".pdf")
      {pdf(file = file,onefile=FALSE, width = input$EMBED2_plot_width, height = input$EMBED2_plot_height)}

      else if(input$plot_choice_download_EMBED2==".png")
      {png(file = file, width = input$EMBED2_plot_width, height = input$EMBED2_plot_height,units="in",res=1000)}

      else
      {tiff(file = file, width = input$EMBED2_plot_width, height = input$EMBED2_plot_height,units="in",res=1000)}


        plot2 <- plot2()

      grid.arrange(plot2)
      dev.off()
    }
  )
  
  output$EMBED_plot_1 <- DT::renderDT(NULL)
  output$EMBED_plot_2 <- DT::renderDT(NULL)
  
  color <- reactive({readRDS(paste0("./",subOutputDir,"/pal.rds"))})
  scale <- reactive({readRDS(paste0("./",subOutputDir,"/scale.rds"))})
  
  #plot EMBED1
  output$EMBED_plot_1<- renderPlot({
    
    plot1()
    
  }, height = 450,width=450)
  
  # #plot EMBED2
  output$EMBED_plot_2<- renderPlot({
    
    plot2()
    
  } ,height = 450,width=450)

  #update EMBED dropdown based on selected Matrix--------------------------------
  
  #Update dropdown for EMBED1
  featureNames1 <- reactive({
    
    if(!(input$matrix_EMBED1_forComparison %in% groupBy)){
      availableMatrices <- getAvailableMatrices(ArchRProj)
      matName <- availableMatrices[ availableMatrices  %in% input$matrix_EMBED1_forComparison]
      featureNames <- h5read(file = paste0(getOutputDirectory(ArchRProj), "/",outputDir, "/", subOutputDir, "/", matName, "_plotBlank72.h5"),
                             name = matName)
      Feature_dropdown1 = names(featureNames)
      return(Feature_dropdown1)
    }
    
  })

  observeEvent(input$matrix_EMBED1_forComparison,{
    
    if(!(input$matrix_EMBED1_forComparison %in% groupBy)){
  updateSelectizeInput(session, 'EMBED1_forComparison', label = 'Feature Name',
                           choices = sort(featureNames1()),
                           server = TRUE,selected =sort(featureNames1())[1])
      }
    })

  # })
  
  
  featureNames2 <- reactive({
    
    if(!(input$matrix_EMBED2_forComparison %in% groupBy)){
          
        availableMatrices <- getAvailableMatrices(ArchRProj)
        matName <- availableMatrices[ availableMatrices  %in% input$matrix_EMBED2_forComparison]
        featureNames <- h5read(file = paste0(getOutputDirectory(ArchRProj), "/",outputDir, "/", subOutputDir, "/", matName, "_plotBlank72.h5"),
               name = matName)
                             
        Feature_dropdown2 = names(featureNames)
        return(Feature_dropdown2)
    
    }
    
  })
  
  observeEvent(input$matrix_EMBED2_forComparison,{
    if(!(input$matrix_EMBED2_forComparison %in% groupBy)){
    updateSelectizeInput(session, 'EMBED2_forComparison', label = 'Feature Name',
                         choices = sort(featureNames2()),
                         server = TRUE,selected =sort(featureNames2())[1])
    }
    })

  #Update dropdown for EMBED2
  # observeEvent(input$matrix_EMBED2_forComparison,{
  #   if(isolate(input$matrix_EMBED2_forComparison)=="Motif Matrix")
  #   {
  # 
  #     updateSelectizeInput(session, 'EMBED2_forComparison', label = 'Feature Name',
  #                          choices = sort(MM_dropdown),
  #                          server = TRUE,selected =sort(MM_dropdown)[2])
  #   }
  # 
  #   else if(isolate(input$matrix_EMBED2_forComparison)=="Gene Score Matrix")
  #   {
  #     updateSelectizeInput(session, 'EMBED2_forComparison',  label = 'Feature Name',
  #                          choices = sort(GSM_dropdown),
  #                          server = TRUE,selected =sort(GSM_dropdown)[2])
  #   }
  #   else if(isolate(input$matrix_EMBED2_forComparison)=="Gene Integration Matrix")
  #   {
  # 
  #     updateSelectizeInput(session, 'EMBED2_forComparison',  label = 'Feature Name',
  #                          choices = sort(GIM_dropdown),
  #                          server = TRUE,selected =sort(GIM_dropdown)[2])
  #   }
  # 
  # })

  # Plot Browser ----------------------------------------------------------------

  # Observe the inputs for ATAC-Seq Explorer
  observeEvent(input$range_min, {
    updateSliderInput(session, "range",
                      value = c(input$range_min,max(input$range)))
  })
  
  observeEvent(input$range_max, {
    updateSliderInput(session, "range",
                      value = c(input$range_min,input$range_max))
  })
  
  observeEvent(input$range , {
    
    updateNumericInput(session, "range_min", value = min(input$range))
    updateNumericInput(session, "range_max", value = max(input$range))
    
  }, priority = 200)
  
  # Output Handler:downloads file
  output$down<-downloadHandler(
    filename <- function(){
      paste0("ArchRBrowser-",input$gene_name,input$plot_choice_download_peakBrowser)
    },
    content = function(file){
      
      if(input$plot_choice_download_peakBrowser==".pdf")
      {pdf(file = file,onefile=FALSE, width = input$plot_width, height = input$plot_height)}
      
      else if(input$plot_choice_download_peakBrowser==".png")
      {png(file = file, width = input$plot_width, height = input$plot_height,units="in",res=1000)}
      
      else
      {tiff(file = file, width = input$plot_width, height = input$plot_height,units="in",res=1000)}
      

        p_browser_atacClusters<- plotBrowserTrack(
          ArchRProj = ArchRProj,
          # ShinyArchR = ShinyArchR,
          plotSummary = c("bulkTrack", input$selectPlotSummary),
          baseSize = 11,
          facetbaseSize = 11,
          groupBy = input$browserContent,
          geneSymbol = isolate(input$gene_name),
          upstream = -min(isolate(input$range))*1000,
          downstream = max(isolate(input$range))*1000,
          tileSize = isolate(input$tile_size),
          ylim =  c(0, isolate(input$ymax)),
          loops = getCoAccessibility(ArchRProj)
          
        )[[input$gene_name]]
     
      
      grid.arrange(p_browser_atacClusters)  
      
      dev.off()
    }
  )
  output$browser_atacClusters <- DT::renderDT(NULL)
  
  #handles error
  restartFN <- observeEvent(input$restartButton, {
    if (isolate(input$gene_name) == ""){
      
      output$browser_atacClusters <- renderPlot({
        p <- ggplot() +
          xlim(c(-5,5)) + ylim(c(-5,5)) +
          geom_text(size=20, aes(x = 0, y = 0, label = "Please supply\na valid gene name!")) + theme_void()
        print(p)
      })
    }else{
      
      # Plots scATACSeq clusters
      output$browser_atacClusters<- renderPlot({
        grid::grid.newpage()

          p_browser_atacClusters<- plotBrowserTrack(
            ArchRProj = ArchRProj,
            # ShinyArchR = ShinyArchR,
            plotSummary = c("bulkTrack", input$selectPlotSummary),
            baseSize = 11,
            facetbaseSize = 11,
            groupBy =  input$browserContent,
            geneSymbol = isolate(input$gene_name),
            upstream = -min(isolate(input$range))*1000,
            downstream = max(isolate(input$range))*1000,
            tileSize = isolate(input$tile_size),
            ylim =  c(0, isolate(input$ymax)),
            loops = getCoAccessibility(ArchRProj)
            
          )[[input$gene_name]]
        
        
        grid::grid.draw(p_browser_atacClusters)
        
      },height = 900)
      
    }
  })
}
