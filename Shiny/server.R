
shinyServer <- function(input,output, session){
  
  
  # UMAPS ------------------------------------------------------------------------------------
  
  #Output Handler: Downloads UMAPS
  output$download_UMAP1<-downloadHandler(
    filename <- function(){
      paste0("UMAP-",paste(input$matrix_UMAP1_forComparison,input$UMAP1_forComparison,sep="-"),input$plot_choice_download_UMAP1)
    },
    content = function(file){
      if(input$plot_choice_download_UMAP1==".pdf")
      {pdf(file = file,onefile=FALSE, width = input$UMAP1_plot_width, height = input$UMAP1_plot_height)}
      
      else if(input$plot_choice_download_UMAP1==".png")
      {png(file = file, width = input$UMAP1_plot_width, height = input$UMAP1_plot_height,units="in",res=1000)}
      
      else
      {tiff(file = file, width = input$UMAP1_plot_width, height = input$UMAP1_plot_height,units="in",res=1000)}
      
      if((input$matrix_UMAP1_forComparison)=="Gene Score Matrix")
      {
        p_empty <-  ggplot()  +
          xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") + theme_bw(base_size=10)+
          ggtitle(paste0("UMAP of IterativeLSI colored by \n GeneScoreMatrix: ",input$UMAP1_forComparison)) +
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
        legend_plot <- gradientLegend(valRange=c(scale()$gsm[,input$UMAP1_forComparison][1],scale()$gsm[,input$UMAP1_forComparison][2]), 
                                      color = color()$gsm, pos=.5, side=1)
        
        
        p <- h5read("./inputData/plotBlank72.h5", paste0("GSM/", input$UMAP1_forComparison))
        temp_jpg <- t(matrix(decode_native(p), nrow = 216))
        
        
        last_plot <- ggdraw() + draw_image(temp_jpg, x = 0, y = 0, scale = 0.85) + 
          draw_plot(p_empty, scale = 0.8) +
          draw_text("Log2(+1)", x = 0.35, y = 0.05, size = 12)
        
        plot1 = print(last_plot, vp=viewport(0.5, 0.6, 1, 1))        
        
      }
      
      else if((input$matrix_UMAP1_forComparison)=="Gene Integration Matrix")
      {
        p_empty <-  ggplot()  +
          xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") + theme_bw(base_size=10)+
          ggtitle(paste0("UMAP of IterativeLSI colored by \n GeneIntegrationMatrix: ",input$UMAP1_forComparison)) +
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
        legend_plot <- gradientLegend(valRange=c(scale()$gim[,input$UMAP1_forComparison][1],scale()$gim[,input$UMAP1_forComparison][2]), 
                                      color = color()$gim, pos=.5, side=1)
        
        p <- h5read("./inputData/plotBlank72.h5", paste0("GIM/", input$UMAP1_forComparison))
        temp_jpg <- t(matrix(decode_native(p), nrow = 216))
        
        
        last_plot <- ggdraw() + draw_image(temp_jpg, x = 0, y = 0, scale = 0.85) + 
          draw_plot(p_empty, scale = 0.8) +
          draw_text("Log2(+1)", x = 0.35, y = 0.05, size = 12)
        
        plot1 = print(last_plot, vp=viewport(0.5, 0.6, 1, 1))   
        
        
      }
      
      else if((input$matrix_UMAP1_forComparison)=="Motif Matrix")
      {
        p_empty <-  ggplot()  +
          xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") + theme_bw(base_size=10)+
          ggtitle(paste0("UMAP of IterativeLSI colored by \n MotifMatrix: ",input$UMAP1_forComparison)) +
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
        legend_plot <- gradientLegend(valRange=c(scale()$mm[,input$UMAP1_forComparison][1],scale()$mm[,input$UMAP1_forComparison][2]), 
                                      color = color()$mm, pos=.5, side=1)
        
        p <- h5read("./inputData/plotBlank72.h5", paste0("MM/", input$UMAP1_forComparison))
        
        temp_jpg <- t(matrix(decode_native(p), nrow = 216))
        
        last_plot <- ggdraw() + draw_image(temp_jpg, x = 0, y = 0, scale = 0.85) + 
          draw_plot(p_empty, scale = 0.8) +
          draw_text("Log2(+1)", x = 0.35, y = 0.05, size = 12)
        
        plot1 = print(last_plot, vp=viewport(0.5, 0.6, 1, 1))   
      }
      
      else
      {
        
        
        if((input$matrix_UMAP1_forComparison)=="Clusters"){
          
          title = "Colored by scATAC-seq clusters" 
          
        }
        
        if((input$matrix_UMAP1_forComparison)=="Constrained"){
          
          title = "UMAP: constrained integration" 
          
        }
        
        if((input$matrix_UMAP1_forComparison)=="Constrained remap"){
          
          title = "UMAP: Constrained remmaped clusters" 
          
        }
        
        if((input$matrix_UMAP1_forComparison)=="Sample"){
          
          title = "Colored by original identity" 
          
        }
        
        if((input$matrix_UMAP1_forComparison)=="Unconstrained"){
          
          title = "UMAP: unconstrained integration" 
          
        }
        
        p_empty <-  ggplot()  +
          xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") + theme_bw(base_size=10)+
          ggtitle(title) +
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
        
        if((input$matrix_UMAP1_forComparison)=="Clusters"){
          
          legend('bottom', legend=umap_legend_names$Clusters, 
                 pch=15, col = color_umaps$Clusters,
                 horiz = FALSE, x.intersp = 1, text.width=0.35,
                 cex = 0.7, bty="n", ncol = 4)
          
          
          
        }
        
        if((input$matrix_UMAP1_forComparison)=="Sample"){
          legend('bottom', legend=umap_legend_names$Sample, pch=15, 
                 col = color,
                 horiz = TRUE, x.intersp = 1, text.width=0.6,
                 cex = 0.7, bty="n")
        }
        
        if((input$matrix_UMAP1_forComparison)=="Constrained"){
          
          legend('bottom', legend=umap_legend_names$Constrained, 
                 pch=15, col = color_umaps$Constrained,
                 horiz = FALSE, x.intersp = 1, text.width=0.35,
                 cex = 0.5, bty="n", ncol = 5)
          
          
        }
        
        
        if((input$matrix_UMAP1_forComparison)=="Unconstrained"){
          legend('bottom', legend=umap_legend_names$Unconstrained, 
                 pch=15, col = color_umaps$unconstrained,
                 horiz = FALSE, x.intersp = 1, text.width=0.35,
                 cex = 0.5, bty="n", ncol = 5)
        }
        
        
        if((input$matrix_UMAP1_forComparison)=="Constrained remap"){
          legend('bottom', legend=umap_legend_names$`Constrained remap`, 
                 pch=15, col = color_umaps$`Constrained remap`,
                 horiz = FALSE, x.intersp = 1, text.width=0.35,
                 cex = 0.7, bty="n", ncol = 4)
        }
        
        p <- h5read("./inputData/mainUMAPs.h5", input$matrix_UMAP1_forComparison)
        temp_jpg <- t(matrix(decode_native(p), nrow = 216))
        
        last_plot <- ggdraw() + draw_image(temp_jpg, x = 0, y = 0, scale = 0.7) + 
          draw_plot(p_empty, scale = 0.7) +
          draw_text("color", x = 0.1, y = 0.135, size = 12)
        
        plot1 = print(last_plot, vp=viewport(0.5, 0.6, 1, 1)) 
      }
      
      
      grid.arrange(plot1)
      dev.off()
    }
  )
  
  output$download_UMAP2<-downloadHandler(
    filename <- function(){
      paste0("UMAP-",paste(input$matrix_UMAP2_forComparison,input$UMAP2_forComparison,sep="-"),input$plot_choice_download_UMAP2)
    },
    content = function(file){
      
      if(input$plot_choice_download_UMAP2==".pdf")
      {pdf(file = file,onefile=FALSE, width = input$UMAP2_plot_width, height = input$UMAP2_plot_height)}
      
      else if(input$plot_choice_download_UMAP2==".png")
      {png(file = file, width = input$UMAP2_plot_width, height = input$UMAP2_plot_height,units="in",res=1000)}
      
      else
      {tiff(file = file, width = input$UMAP2_plot_width, height = input$UMAP2_plot_height,units="in",res=1000)}
      
      
      
      
      if((input$matrix_UMAP2_forComparison)=="Gene Score Matrix")
      {
        p_empty <-  ggplot()  +
          xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") + theme_bw(base_size=10)+
          ggtitle(paste0("UMAP of IterativeLSI colored by \n GeneScoreMatrix: ",input$UMAP1_forComparison)) +
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
        legend_plot <- gradientLegend(valRange=c(scale()$gsm[,input$UMAP1_forComparison][1],scale()$gsm[,input$UMAP1_forComparison][2]), 
                                      color = color()$gsm, pos=.5, side=1)
        
        
        p <- h5read("./inputData/plotBlank72.h5", paste0("GSM/", input$UMAP1_forComparison))
        temp_jpg <- t(matrix(decode_native(p), nrow = 216))
        
        
        last_plot <- ggdraw() + draw_image(temp_jpg, x = 0, y = 0, scale = 0.85) + 
          draw_plot(p_empty, scale = 0.8) +
          draw_text("Log2(+1)", x = 0.35, y = 0.05, size = 12)
        
        plot2 = print(last_plot, vp=viewport(0.5, 0.6, 1, 1))        
        
      }
      
      else if((input$matrix_UMAP2_forComparison)=="Gene Integration Matrix")
      {
        p_empty <-  ggplot()  +
          xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") + theme_bw(base_size=10)+
          ggtitle(paste0("UMAP of IterativeLSI colored by \n GeneIntegrationMatrix: ",input$UMAP1_forComparison)) +
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
        legend_plot <- gradientLegend(valRange=c(scale()$gim[,input$UMAP1_forComparison][1],scale()$gim[,input$UMAP1_forComparison][2]), 
                                      color = color()$gim, pos=.5, side=1)
        
        p <- h5read("./inputData/plotBlank72.h5", paste0("GIM/", input$UMAP1_forComparison))
        temp_jpg <- t(matrix(decode_native(p), nrow = 216))
        
        
        last_plot <- ggdraw() + draw_image(temp_jpg, x = 0, y = 0, scale = 0.85) + 
          draw_plot(p_empty, scale = 0.8) +
          draw_text("Log2(+1)", x = 0.35, y = 0.05, size = 12)
        
        plot2 = print(last_plot, vp=viewport(0.5, 0.6, 1, 1))   
        
        
      }
      
      else if((input$matrix_UMAP2_forComparison)=="Motif Matrix")
      {
        p_empty <-  ggplot()  +
          xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") + theme_bw(base_size=10)+
          ggtitle(paste0("UMAP of IterativeLSI colored by \n MotifMatrix: ",input$UMAP1_forComparison)) +
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
        legend_plot <- gradientLegend(valRange=c(scale()$mm[,input$UMAP1_forComparison][1],scale()$mm[,input$UMAP1_forComparison][2]), 
                                      color = color()$mm, pos=.5, side=1)
        
        p <- h5read("./inputData/plotBlank72.h5", paste0("MM/", input$UMAP1_forComparison))
        
        temp_jpg <- t(matrix(decode_native(p), nrow = 216))
        
        last_plot <- ggdraw() + draw_image(temp_jpg, x = 0, y = 0, scale = 0.85) + 
          draw_plot(p_empty, scale = 0.8) +
          draw_text("Log2(+1)", x = 0.35, y = 0.05, size = 12)
        
        plot2 = print(last_plot, vp=viewport(0.5, 0.6, 1, 1))   
        
        
      }
      
      else
      {
        
        if((input$matrix_UMAP2_forComparison)=="Clusters"){
          
          title = "Colored by scATAC-seq clusters" 
          
        }
        
        if((input$matrix_UMAP2_forComparison)=="Constrained"){
          
          title = "UMAP: constrained integration" 
          
        }
        
        if((input$matrix_UMAP2_forComparison)=="Constrained remap"){
          
          title = "UMAP: Constrained remmaped clusters" 
          
        }
        
        if((input$matrix_UMAP2_forComparison)=="Sample"){
          
          title = "Colored by original identity" 
          
        }
        
        if((input$matrix_UMAP2_forComparison)=="Unconstrained"){
          
          title = "UMAP: unconstrained integration" 
          
        }
        
        p_empty <-  ggplot()  +
          xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") + theme_bw(base_size=10)+
          ggtitle(title) +
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
        
        if((input$matrix_UMAP2_forComparison)=="Clusters"){
          
          legend('bottom', legend=umap_legend_names$Clusters, 
                 pch=15, col = color_umaps$Clusters,
                 horiz = FALSE, x.intersp = 1, text.width=0.35,
                 cex = 0.7, bty="n", ncol = 4)
          
          
          
        }
        
        if((input$matrix_UMAP2_forComparison)=="Sample"){
          legend('bottom', legend=umap_legend_names$Sample, pch=15, 
                 col = color_umaps$Sample,
                 horiz = TRUE, x.intersp = 1, text.width=0.6,
                 cex = 0.7, bty="n")
        }
        
        if((input$matrix_UMAP2_forComparison)=="Constrained"){
          
          legend('bottom', legend=umap_legend_names$Constrained, 
                 pch=15, col = color_umaps$Constrained,
                 horiz = FALSE, x.intersp = 1, text.width=0.35,
                 cex = 0.5, bty="n", ncol = 5)
          
          
        }
        
        
        if((input$matrix_UMAP2_forComparison)=="Unconstrained"){
          legend('bottom', legend=umap_legend_names$Unconstrained, 
                 pch=15, col = color_umaps$unconstrained,
                 horiz = FALSE, x.intersp = 1, text.width=0.35,
                 cex = 0.5, bty="n", ncol = 5)
        }
        
        
        if((input$matrix_UMAP2_forComparison)=="Constrained remap"){
          legend('bottom', legend=umap_legend_names$`Constrained remap`, 
                 pch=15, col = color_umaps$`Constrained remap`,
                 horiz = FALSE, x.intersp = 1, text.width=0.35,
                 cex = 0.7, bty="n", ncol = 4)
        }
        
        p <- h5read("./inputData/mainUMAPs.h5", input$matrix_UMAP2_forComparison)
        temp_jpg <- t(matrix(decode_native(p), nrow = 216))
        
        last_plot <- ggdraw() + draw_image(temp_jpg, x = 0, y = 0, scale = 0.7) + 
          draw_plot(p_empty, scale = 0.7) +
          draw_text("color", x = 0.1, y = 0.135, size = 12)
        
        plot2 = print(last_plot, vp=viewport(0.5, 0.6, 1, 1)) 
      }
      
      
      grid.arrange(plot2)
      dev.off()
    }
  )
  
  output$UMAP_plot_1 <- DT::renderDT(NULL)
  output$UMAP_plot_2 <- DT::renderDT(NULL)
  
  color <- reactive({readRDS("./inputData/pal.rds")})
  scale <- reactive({readRDS("./inputData/scale.rds")})
  
  #plot UMAP1
  output$UMAP_plot_1<- renderPlot({
    
    
    if((input$matrix_UMAP1_forComparison)=="Gene Score Matrix")
    {
     
      p_empty <-  ggplot()  +
        xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") + theme_bw(base_size=10)+
        ggtitle(paste0("UMAP of IterativeLSI colored by \n GeneScoreMatrix: ",input$UMAP1_forComparison)) +
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
      legend_plot <- gradientLegend(valRange=c(scale()$gsm[,input$UMAP1_forComparison][1],scale()$gsm[,input$UMAP1_forComparison][2]), 
                                    color = color()$gsm, pos=.5, side=1)
      
      
      p <- h5read("./inputData/plotBlank72.h5", paste0("GSM/", input$UMAP1_forComparison))
      temp_jpg <- t(matrix(decode_native(p), nrow = 216))
      
      
      last_plot <- ggdraw() + draw_image(temp_jpg, x = 0, y = 0, scale = 0.85) + 
        draw_plot(p_empty, scale = 0.8) +
        draw_text("Log2(+1)", x = 0.35, y = 0.05, size = 12)
      
      print(last_plot, vp=viewport(0.5, 0.6, 1, 1))        
      
    }
    
    else if((input$matrix_UMAP1_forComparison)=="Gene Integration Matrix")
    {
      # getUmap(input$UMAP1_forComparison,GSM_Umaps_data_fileIndexer,"GSM_Umaps_data","plot_scaffold_GSM",isolate(input$matrix_UMAP1_forComparison))
      
      
      p_empty <-  ggplot()  +
        xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") + theme_bw(base_size=10)+
        ggtitle(paste0("UMAP of IterativeLSI colored by \n GeneIntegrationMatrix: ",input$UMAP1_forComparison)) +
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
      legend_plot <- gradientLegend(valRange=c(scale()$gim[,input$UMAP1_forComparison][1],scale()$gim[,input$UMAP1_forComparison][2]), 
                                    color = color()$gim, pos=.5, side=1)
      
      p <- h5read("./inputData/plotBlank72.h5", paste0("GIM/", input$UMAP1_forComparison))
      temp_jpg <- t(matrix(decode_native(p), nrow = 216))
      
      
      last_plot <- ggdraw() + draw_image(temp_jpg, x = 0, y = 0, scale = 0.85) + 
        draw_plot(p_empty, scale = 0.8) +
        draw_text("Log2(+1)", x = 0.35, y = 0.05, size = 12)
      
      print(last_plot, vp=viewport(0.5, 0.6, 1, 1))   
      
      
    }
    
    else if((input$matrix_UMAP1_forComparison)=="Motif Matrix")
    {
      p_empty <-  ggplot()  +
        xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") + theme_bw(base_size=10)+
        ggtitle(paste0("UMAP of IterativeLSI colored by \n MotifMatrix: ",input$UMAP1_forComparison)) +
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
      legend_plot <- gradientLegend(valRange=c(scale()$mm[,input$UMAP1_forComparison][1],scale()$mm[,input$UMAP1_forComparison][2]), 
                                    color = color()$mm, pos=.5, side=1)
      
      p <- h5read("./inputData/plotBlank72.h5", paste0("MM/", input$UMAP1_forComparison))
      
      temp_jpg <- t(matrix(decode_native(p), nrow = 216))
      
      last_plot <- ggdraw() + draw_image(temp_jpg, x = 0, y = 0, scale = 0.85) + 
        draw_plot(p_empty, scale = 0.8) +
        draw_text("Log2(+1)", x = 0.35, y = 0.05, size = 12)
      
      print(last_plot, vp=viewport(0.5, 0.6, 1, 1))   
      
      
    }
    
    else
    {
      
      
      if((input$matrix_UMAP1_forComparison)=="Clusters"){
        
        title = "Colored by scATAC-seq clusters" 
        
      }
      
      if((input$matrix_UMAP1_forComparison)=="Constrained"){
        
        title = "UMAP: constrained integration" 
        
      }
      
      if((input$matrix_UMAP1_forComparison)=="Constrained remap"){
        
        title = "UMAP: Constrained remmaped clusters" 
        
      }
      
      if((input$matrix_UMAP1_forComparison)=="Sample"){
        
        title = "Colored by original identity" 
        
      }
      
      if((input$matrix_UMAP1_forComparison)=="Unconstrained"){
        
        title = "UMAP: unconstrained integration" 
        
      }
      
      p_empty <-  ggplot()  +
        xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") + theme_bw(base_size=10)+
        ggtitle(title) +
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
      
      if((input$matrix_UMAP1_forComparison)=="Clusters"){
        
        legend('bottom', legend=umap_legend_names$Clusters, 
               pch=15, col = color_umaps$Clusters,
               horiz = FALSE, x.intersp = 1, text.width=0.35,
               cex = 0.7, bty="n", ncol = 4)
        
        
        
      }
      
      if((input$matrix_UMAP1_forComparison)=="Sample"){
        legend('bottom', legend=umap_legend_names$Sample, pch=15, 
               col = color_umaps$Sample,
               horiz = TRUE, x.intersp = 1, text.width=0.6,
               cex = 0.7, bty="n")
      }
      
      if((input$matrix_UMAP1_forComparison)=="Constrained"){
        
        legend('bottom', legend=umap_legend_names$Constrained, 
               pch=15, col = color_umaps$Constrained,
               horiz = FALSE, x.intersp = 1, text.width=0.35,
               cex = 0.5, bty="n", ncol = 5)
        
        
      }
      
      
      if((input$matrix_UMAP1_forComparison)=="Unconstrained"){
        legend('bottom', legend=umap_legend_names$Unconstrained, 
               pch=15, col = color_umaps$unconstrained,
               horiz = FALSE, x.intersp = 1, text.width=0.35,
               cex = 0.5, bty="n", ncol = 5)
      }
      
      
      if((input$matrix_UMAP1_forComparison)=="Constrained remap"){
        legend('bottom', legend=umap_legend_names$`Constrained remap`, 
               pch=15, col = color_umaps$`Constrained remap`,
               horiz = FALSE, x.intersp = 1, text.width=0.35,
               cex = 0.7, bty="n", ncol = 4)
      }
      
      p <- h5read("./inputData/mainUMAPs.h5", input$matrix_UMAP1_forComparison)# input$UMAP1_forComparison))
      temp_jpg <- t(matrix(decode_native(p), nrow = 216))
      
      last_plot <- ggdraw() + draw_image(temp_jpg, x = 0, y = 0, scale = 0.7) + 
        draw_plot(p_empty, scale = 0.7) +
        draw_text("color", x = 0.1, y = 0.135, size = 12)
      
      print(last_plot, vp=viewport(0.5, 0.6, 1, 1)) 
      
      # umaps[input$matrix_UMAP1_forComparison]
      
      
    }
    
  }, height = 450,width=450)
  
  # #plot UMAP2
  output$UMAP_plot_2<- renderPlot({
    if((input$matrix_UMAP2_forComparison)=="Gene Score Matrix")
    {
      
      p_empty <-  ggplot()  +
        xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") + theme_bw(base_size=10)+
        ggtitle(paste0("UMAP of IterativeLSI colored by \n GeneScoreMatrix: ",input$UMAP2_forComparison)) +
        theme(
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          panel.grid.major = element_blank(), #remove major gridlines
          panel.grid.minor = element_blank(), #remove minor gridlines
          legend.background = element_rect(fill='transparent'), #transparent legend bg
          legend.box.background = element_rect(fill='transparent'),  axis.title=element_text(size=14),           plot.title = element_text(size=16)
        )
      
      emptyPlot(0,0, axes=FALSE)
      legend_plot <- gradientLegend(valRange=c(scale()$gsm[,input$UMAP2_forComparison][1],scale()$gsm[,input$UMAP2_forComparison][2]), 
                                    color = color()$gsm, pos=.5, side=1)
      
      p <- h5read("./inputData/plotBlank72.h5", paste0("GSM/", input$UMAP2_forComparison))# input$UMAP1_forComparison))
      temp_jpg <- t(matrix(decode_native(p), nrow = 216))
      
      last_plot <- ggdraw() + draw_image(temp_jpg, x = 0, y = 0, scale = 0.85) + 
        draw_plot(p_empty, scale = 0.8) +
        draw_text("Log2(+1)", x = 0.35, y = 0.05, size = 12)
      
      print(last_plot, vp=viewport(0.5, 0.6, 1, 1))        
      
    }
    
    else if((input$matrix_UMAP2_forComparison)=="Gene Integration Matrix")
    {
      # getUmap(input$UMAP2_forComparison,GSM_Umaps_data_fileIndexer,"GSM_Umaps_data","plot_scaffold_GSM",isolate(input$matrix_UMAP2_forComparison))
      
      p_empty <-  ggplot()  +
        xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") + theme_bw(base_size=10)+
        ggtitle(paste0("UMAP of IterativeLSI colored by \n GeneIntegrationMatrix: ",input$UMAP2_forComparison)) +
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
      legend_plot <- gradientLegend(valRange=c(scale()$gim[,input$UMAP2_forComparison][1],scale()$gim[,input$UMAP2_forComparison][2]), 
                                    color = color()$gim, pos=.5, side=1)
      
      p <- h5read("./inputData/plotBlank72.h5", paste0("GIM/", input$UMAP2_forComparison))# input$UMAP1_forComparison))
      temp_jpg <- t(matrix(decode_native(p), nrow = 216))
      
      last_plot <- ggdraw() + draw_image(temp_jpg, x = 0, y = 0, scale = 0.85) +       
        draw_plot(p_empty, scale = 0.8) +
        draw_text("Log2(+1)", x = 0.35, y = 0.05, size = 12)
      
      print(last_plot, vp=viewport(0.5, 0.6, 1, 1)) 
      
      
    }
    
    else if((input$matrix_UMAP2_forComparison)=="Motif Matrix")
    {
      # getUmap(input$UMAP2_forComparison,MM_Umaps_data_fileIndexer,"MM_Umaps_data","plot_scaffold_MM",isolate(input$matrix_UMAP2_forComparison))
      
      p_empty <-  ggplot()  +
        xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") + theme_bw(base_size=10)+
        ggtitle(paste0("UMAP of IterativeLSI colored by \n MotifMatrix: ",input$UMAP2_forComparison)) +
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
      legend_plot <- gradientLegend(valRange=c(scale()$mm[,input$UMAP2_forComparison][1],scale()$mm[,input$UMAP2_forComparison][2]), 
                                    color = color()$mm, pos=.5, side=1)
      
      p <- h5read("./inputData/plotBlank72.h5", paste0("MM/", input$UMAP2_forComparison))# input$UMAP1_forComparison))
      temp_jpg <- t(matrix(decode_native(p), nrow = 216))
      
      
      last_plot <- ggdraw() + draw_image(temp_jpg, x = 0, y = 0, scale = 0.85) + 
        draw_plot(p_empty, scale = 0.8) +
        draw_text("Log2(+1)", x = 0.35, y = 0.05, size = 12)
      
      print(last_plot, vp=viewport(0.5, 0.6, 1, 1)) 
      
      
    }
    
    else
    {
      # umaps[input$matrix_UMAP2_forComparison]
      
      
      if((input$matrix_UMAP2_forComparison)=="Clusters"){
        
        title = "Colored by scATAC-seq clusters" 
        
      }
      
      if((input$matrix_UMAP2_forComparison)=="Constrained"){
        
        title = "UMAP: constrained integration" 
        
      }
      
      if((input$matrix_UMAP2_forComparison)=="Constrained remap"){
        
        title = "UMAP: Constrained remmaped clusters" 
        
      }
      
      if((input$matrix_UMAP2_forComparison)=="Sample"){
        
        title = "Colored by original identity" 
        
      }
      
      if((input$matrix_UMAP2_forComparison)=="Unconstrained"){
        
        title = "UMAP: unconstrained integration" 
        
      }
      
      p_empty <-  ggplot()  +
        xlab("UMAP Dimension 1") + ylab("UMAP Dimension 2") + theme_bw(base_size=10)+
        ggtitle(title) +
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
      
      if((input$matrix_UMAP2_forComparison)=="Clusters"){
        
        legend('bottom', legend=umap_legend_names$Clusters, 
               pch=15, col = color_umaps$Clusters,
               horiz = FALSE, x.intersp = 1, text.width=0.35,
               cex = 0.7, bty="n", ncol = 4)
        
        
        
      }
      
      if((input$matrix_UMAP2_forComparison)=="Sample"){
        legend('bottom', legend=umap_legend_names$Sample, pch=15, 
               col = color_umaps$Sample,
               horiz = TRUE, x.intersp = 1, text.width=0.6,
               cex = 0.7, bty="n")
      }
      
      if((input$matrix_UMAP2_forComparison)=="Constrained"){
        
        legend('bottom', legend=umap_legend_names$Constrained, 
               pch=15, col = color_umaps$Constrained,
               horiz = FALSE, x.intersp = 1, text.width=0.35,
               cex = 0.5, bty="n", ncol = 5)
        
        
      }
      
      
      if((input$matrix_UMAP2_forComparison)=="Unconstrained"){
        legend('bottom', legend=umap_legend_names$Unconstrained, 
               pch=15, col = color_umaps$unconstrained,
               horiz = FALSE, x.intersp = 1, text.width=0.35,
               cex = 0.5, bty="n", ncol = 5)
      }
      
      
      if((input$matrix_UMAP2_forComparison)=="Constrained remap"){
        legend('bottom', legend=umap_legend_names$`Constrained remap`, 
               pch=15, col = color_umaps$`Constrained remap`,
               horiz = FALSE, x.intersp = 1, text.width=0.35,
               cex = 0.7, bty="n", ncol = 4)
      }
      
      p <- h5read("./inputData/mainUMAPs.h5", input$matrix_UMAP2_forComparison)
      temp_jpg <- t(matrix(decode_native(p), nrow = 216))
      
      last_plot <- ggdraw() + draw_image(temp_jpg, x = 0, y = 0, scale = 0.7) + 
        draw_plot(p_empty, scale = 0.7) +
        draw_text("color", x = 0.1, y = 0.135, size = 12)
      
      print(last_plot, vp=viewport(0.5, 0.6, 1, 1)) 
      
    }
    
  } ,height = 450,width=450)

  #update Umap dropdown based on selected Matrix--------------------------------
  
  #Update dropdown for UMAP1
  observeEvent(input$matrix_UMAP1_forComparison,{
    if(isolate(input$matrix_UMAP1_forComparison)=="Motif Matrix")
    {
      updateSelectizeInput(session, 'UMAP1_forComparison', label = 'Feature Name',
                           choices = sort(MM_dropdown),
                           server = TRUE,selected =sort(MM_dropdown)[1])
    }
    
    else if(isolate(input$matrix_UMAP1_forComparison)=="Gene Score Matrix")
    {
      updateSelectizeInput(session, 'UMAP1_forComparison', label = 'Feature Name',
                           choices = sort(GSM_dropdown),
                           server = TRUE,selected =sort(GSM_dropdown)[1])
    }
    else if(isolate(input$matrix_UMAP1_forComparison)=="Gene Integration Matrix")
    {
      updateSelectizeInput(session, 'UMAP1_forComparison', label = 'Feature Name',
                           choices = sort(GIM_dropdown),
                           server = TRUE,selected =sort(GIM_dropdown)[1])
    }
    
  })
  
  #Update dropdown for UMAP2
  observeEvent(input$matrix_UMAP2_forComparison,{
    if(isolate(input$matrix_UMAP2_forComparison)=="Motif Matrix")
    {
      
      updateSelectizeInput(session, 'UMAP2_forComparison', label = 'Feature Name',
                           choices = sort(MM_dropdown),
                           server = TRUE,selected =sort(MM_dropdown)[2])
    }
    
    else if(isolate(input$matrix_UMAP2_forComparison)=="Gene Score Matrix")
    {
      updateSelectizeInput(session, 'UMAP2_forComparison',  label = 'Feature Name',
                           choices = sort(GSM_dropdown),
                           server = TRUE,selected =sort(GSM_dropdown)[2])
    }
    else if(isolate(input$matrix_UMAP2_forComparison)=="Gene Integration Matrix")
    {
      
      updateSelectizeInput(session, 'UMAP2_forComparison',  label = 'Feature Name',
                           choices = sort(GIM_dropdown),
                           server = TRUE,selected =sort(GIM_dropdown)[2])
    }
   
  })
 
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
      
      
      if(isolate(input$browserContent)=="Unconstrained")
      {
        p_browser_atacClusters<- plotBrowserTrack(
          ArchRProj = ArchRProj,
          ShinyArchR = TRUE,
          plotSummary = c("bulkTrack", input$selectPlotSummary),
          baseSize = 11,
          facetbaseSize = 11,
          groupBy = "Clusters",
          geneSymbol = isolate(input$gene_name),
          upstream = -min(isolate(input$range))*1000,
          downstream = max(isolate(input$range))*1000,
          tileSize = isolate(input$tile_size),
          ylim =  c(0, isolate(input$ymax)),
          loops = getCoAccessibility(ArchRProj)
          
        )[[input$gene_name]]
      }
      else
      {
        
        p_browser_atacClusters <- plotBrowserTrack(
          ArchRProj = ArchRProj,
          ShinyArchR = TRUE,
          plotSummary = c("bulkTrack", input$selectPlotSummary),
          groupBy = "Clusters",
          baseSize = 11,
          facetbaseSize = 11,
          geneSymbol = isolate(input$gene_name),
          upstream =-min(isolate(input$range))*1000 ,
          downstream = max(isolate(input$range))*1000,
          tileSize = isolate(input$tile_size),
          ylim =  c(0, isolate(input$ymax)),
          loops = getPeak2GeneLinks(ArchRProj)
        )[[input$gene_name]]
      }
      
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
        
        if(isolate(input$browserContent)=="Unconstrained")
        {
          p_browser_atacClusters<- plotBrowserTrack(
            ArchRProj = ArchRProj,
            ShinyArchR = TRUE,
            plotSummary = c("bulkTrack", input$selectPlotSummary),
            baseSize = 11,
            facetbaseSize = 11,
            groupBy = "Clusters",
            geneSymbol = isolate(input$gene_name),
            upstream = -min(isolate(input$range))*1000,
            downstream = max(isolate(input$range))*1000,
            tileSize = isolate(input$tile_size),
            ylim =  c(0, isolate(input$ymax)),
            loops = getCoAccessibility(ArchRProj)
            
          )[[input$gene_name]]
        }
        else
        {
          p_browser_atacClusters <- plotBrowserTrack(
            ArchRProj = ArchRProj,
            ShinyArchR = TRUE,
            plotSummary = c("bulkTrack", input$selectPlotSummary),
            groupBy = "Clusters",
            baseSize = 11,
            facetbaseSize = 11,
            geneSymbol = isolate(input$gene_name),
            upstream =-min(isolate(input$range))*1000 ,
            downstream = max(isolate(input$range))*1000,
            tileSize = isolate(input$tile_size),
            ylim =  c(0, isolate(input$ymax)),
            loops = getPeak2GeneLinks(ArchRProj)
          )[[input$gene_name]]
        }
        
        
        
        grid::grid.draw(p_browser_atacClusters)
        
      },height = 900)
      
    }
  })
}
