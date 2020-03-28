library(shiny)
library(shinyWidgets)
require(tidyverse)
require(cowplot)
require(shinyjs)
library("ggplot2")
library(shinythemes)
library(annotables)
library(rlist)
library(BuenColors)
library(GenomicAlignments)
options(repos = BiocManager::repositories())
require(dplyr)
require(data.table)
require(readr)
require(stringr)
require(tidyr)

genes <-annotables::grch38 %>%
  dplyr::select(symbol)
erythroid_line=readRDS('erythroid_line.rds')
hsc_differentiation_pathway=readRDS('hsc_differentiation_pathway.rds')

server <- function(input, output, session) {
  
  erythroid_line_fcn <- function() {
    
    df <- erythroid_line %>% dplyr::filter(symbol==input$genes)
    min_val <- df %>%
      gather(cell, pct, -(ensg:transcript_id)) %>%
      group_by(cell) %>%
      mutate(max = max(pct)) 
    
    unique_t <- min_val %>% dplyr::filter(pct==max) %>% pull(transcript_id) %>% unique
    
    max<-df %>%
      gather(cell, pct, -(ensg:transcript_id)) %>%
      group_by(cell) %>% dplyr::filter(transcript_id %in% unique_t) %>%
      ungroup() %>%
      group_by(cell) %>%
      mutate(other = 100 - sum(pct)) %>%
      spread(transcript_id, pct) %>%
      gather(transcript_id, pct, -(ensg:cell))
    max = max %>% dplyr::filter(cell !='hs_BFU') %>% dplyr::filter(cell != 'hs_CFU') %>% dplyr::filter(cell != 'hs_CD34')
    max$type = case_when(max$cell %in% c('CB_BFU', 'CB_CD34', 'CB_CFU', 'hs_early_basophilic',
                                         'hs_late_basophilic', 'hs_orthochromatic', 'hs_polychromatic', 
                                         'hs_proerythroblast') ~ 'Cord Blood',
                         
                         max$cell %in% c('hs_PB_BFU', 'hs_PB_CD34', 'hs_PB_CFU', 
                                         'hs_PB_early_basophilic', 'hs_PB_late_basophilic',
                                         'hs_PB_orthochromatic', 'hs_PB_polychromatic', 
                                         'hs_PB_proerythroblast') ~ 'Peripheral Blood')
    
    max$cell = case_when(max$cell %in% c('hs_PB_BFU', 'CB_BFU') ~ 'BFUE', 
                         max$cell %in% c('hs_PB_CD34', 'CB_CD34') ~ 'CD34',
                         max$cell %in% c('hs_PB_CFU', 'CB_CFU') ~ 'CFUE',
                         max$cell %in% c('hs_PB_early_basophilic', 'hs_early_basophilic') ~ 'eBaso',
                         max$cell %in% c('hs_PB_late_basophilic', 'hs_late_basophilic') ~ 'lBaso',
                         max$cell %in% c('hs_PB_orthochromatic', 'hs_orthochromatic') ~ 'Ortho',
                         max$cell %in% c('hs_PB_polychromatic', 'hs_polychromatic') ~ 'Poly',
                         max$cell %in% c('hs_PB_proerythroblast', 'hs_proerythroblast') ~ 'Pro',)
    
    max$num = case_when(max$cell == "CD34" ~ 8,
                        max$cell == "BFUE" ~ 7,
                        max$cell == "CFUE" ~ 6,
                        max$cell == "Pro" ~ 5,
                        max$cell == "eBaso" ~ 4,
                        max$cell == "lBaso" ~ 3,
                        max$cell == "Poly" ~ 2,
                        max$cell == "Ortho" ~ 1)
    max <- max[order(max$num),]
    
    max$cell = with(max, reorder(cell,num))
    max_pb_cb=max
    # hsc differentiation pathways
    df <- hsc_differentiation_pathway %>% dplyr::filter(symbol==input$genes)
    
    min_val <- df %>%
      gather(cell, pct, -(ensg:transcript_id)) %>%
      group_by(cell) %>%
      mutate(max = max(pct)) 
    
    unique_t <- min_val %>% dplyr::filter(pct==max) %>% pull(transcript_id) %>% unique
    
    max<-df %>%
      gather(cell, pct, -(ensg:transcript_id)) %>%
      group_by(cell) %>% dplyr::filter(transcript_id %in% unique_t) %>%
      ungroup() %>%
      group_by(cell) %>%
      mutate(other = 100 - sum(pct)) %>%
      spread(transcript_id, pct) %>%
      gather(transcript_id, pct, -(ensg:cell))
    
    max$num = case_when(max$cell == "HSC" ~ 23,
                        max$cell == "LSC" ~ 22,
                        max$cell == "pHSC" ~ 21,
                        max$cell == 'Blast' ~ 20,
                        max$cell == "MPP" ~ 19,
                        max$cell == "LMPP" ~ 18,
                        max$cell == "CMP" ~ 17,
                        max$cell == "CLP" ~ 16,
                        max$cell == "GMP-A" ~ 15,
                        max$cell == "GMP-B" ~ 14,
                        max$cell == "GMP-C" ~ 13,
                        max$cell == "MEP" ~ 12,
                        max$cell == "NKcell" ~ 11,
                        max$cell == "CD4Tcell" ~ 10,
                        max$cell == "CD8Tcell" ~ 9,
                        max$cell == "Bcell" ~ 8,
                        max$cell == "PDC" ~ 7,
                        max$cell == "Mono" ~ 6,
                        max$cell == "Ery" ~ 5,
                        max$cell == "MDC" ~ 4,
                        max$cell == "Neutrophils" ~ 3,
                        max$cell == 'Basophils' ~2,
                        max$cell == "platelet" ~ 1)
    max$cell = with(max, reorder(cell,num))
    
    max$type='x'
    max$plot='plot2'
    max_pb_cb$plot='plot1'
    max_bind <- rbind(max_pb_cb,max)
    
    myColors <- colors()[1:length(unique(max_bind$transcript_id))]
    names(myColors) <- unique(max_bind$transcript_id)
    colScale <- scale_colour_manual(name = "transcript_id",values = myColors)

    max_bind$transcript_id <- as.factor(max_bind$transcript_id)
  
    p1=( ggplot(subset(max_bind, plot == "plot1"), aes(y=pct, x=cell,fill=transcript_id))
      + geom_bar(stat="identity",lwd=0.5, aes(color=transcript_id))
      + scale_colour_hue(drop=F, aesthetics = c("colour", "fill")) + coord_flip() + pretty_plot() + facet_wrap(~type)
    )

    p2=( ggplot(subset(max_bind, plot == "plot2"), aes( y=pct, x=cell,fill=transcript_id))
      + geom_bar(stat="identity",lwd=0.5, aes(color=transcript_id))
      + scale_colour_hue(drop=F, aesthetics = c("colour", "fill")) + coord_flip() + pretty_plot()
    )
    p1

  }
  hsc_differentiation_pathway_fcn <- function() {
    
    df <- erythroid_line %>% dplyr::filter(symbol==input$genes)
    min_val <- df %>%
      gather(cell, pct, -(ensg:transcript_id)) %>%
      group_by(cell) %>%
      mutate(max = max(pct)) 
    
    unique_t <- min_val %>% dplyr::filter(pct==max) %>% pull(transcript_id) %>% unique
    
    max<-df %>%
      gather(cell, pct, -(ensg:transcript_id)) %>%
      group_by(cell) %>% dplyr::filter(transcript_id %in% unique_t) %>%
      ungroup() %>%
      group_by(cell) %>%
      mutate(other = 100 - sum(pct)) %>%
      spread(transcript_id, pct) %>%
      gather(transcript_id, pct, -(ensg:cell))
    max = max %>% dplyr::filter(cell !='hs_BFU') %>% dplyr::filter(cell != 'hs_CFU') %>% dplyr::filter(cell != 'hs_CD34')
    max$type = case_when(max$cell %in% c('CB_BFU', 'CB_CD34', 'CB_CFU', 'hs_early_basophilic',
                                         'hs_late_basophilic', 'hs_orthochromatic', 'hs_polychromatic', 
                                         'hs_proerythroblast') ~ 'Cord Blood',
                         
                         max$cell %in% c('hs_PB_BFU', 'hs_PB_CD34', 'hs_PB_CFU', 
                                         'hs_PB_early_basophilic', 'hs_PB_late_basophilic',
                                         'hs_PB_orthochromatic', 'hs_PB_polychromatic', 
                                         'hs_PB_proerythroblast') ~ 'Peripheral Blood')
    
    max$cell = case_when(max$cell %in% c('hs_PB_BFU', 'CB_BFU') ~ 'BFUE', 
                         max$cell %in% c('hs_PB_CD34', 'CB_CD34') ~ 'CD34',
                         max$cell %in% c('hs_PB_CFU', 'CB_CFU') ~ 'CFUE',
                         max$cell %in% c('hs_PB_early_basophilic', 'hs_early_basophilic') ~ 'eBaso',
                         max$cell %in% c('hs_PB_late_basophilic', 'hs_late_basophilic') ~ 'lBaso',
                         max$cell %in% c('hs_PB_orthochromatic', 'hs_orthochromatic') ~ 'Ortho',
                         max$cell %in% c('hs_PB_polychromatic', 'hs_polychromatic') ~ 'Poly',
                         max$cell %in% c('hs_PB_proerythroblast', 'hs_proerythroblast') ~ 'Pro',)
    
    max$num = case_when(max$cell == "CD34" ~ 8,
                        max$cell == "BFUE" ~ 7,
                        max$cell == "CFUE" ~ 6,
                        max$cell == "Pro" ~ 5,
                        max$cell == "eBaso" ~ 4,
                        max$cell == "lBaso" ~ 3,
                        max$cell == "Poly" ~ 2,
                        max$cell == "Ortho" ~ 1)
    max <- max[order(max$num),]
    
    max$cell = with(max, reorder(cell,num))
    max_pb_cb=max
    # hsc differentiation pathways
    df <- hsc_differentiation_pathway %>% dplyr::filter(symbol==input$genes)
    
    min_val <- df %>%
      gather(cell, pct, -(ensg:transcript_id)) %>%
      group_by(cell) %>%
      mutate(max = max(pct)) 
    
    unique_t <- min_val %>% dplyr::filter(pct==max) %>% pull(transcript_id) %>% unique
    
    max<-df %>%
      gather(cell, pct, -(ensg:transcript_id)) %>%
      group_by(cell) %>% dplyr::filter(transcript_id %in% unique_t) %>%
      ungroup() %>%
      group_by(cell) %>%
      mutate(other = 100 - sum(pct)) %>%
      spread(transcript_id, pct) %>%
      gather(transcript_id, pct, -(ensg:cell))
    
    max$num = case_when(max$cell == "HSC" ~ 23,
                        max$cell == "LSC" ~ 22,
                        max$cell == "pHSC" ~ 21,
                        max$cell == 'Blast' ~ 20,
                        max$cell == "MPP" ~ 19,
                        max$cell == "LMPP" ~ 18,
                        max$cell == "CMP" ~ 17,
                        max$cell == "CLP" ~ 16,
                        max$cell == "GMP-A" ~ 15,
                        max$cell == "GMP-B" ~ 14,
                        max$cell == "GMP-C" ~ 13,
                        max$cell == "MEP" ~ 12,
                        max$cell == "NKcell" ~ 11,
                        max$cell == "CD4Tcell" ~ 10,
                        max$cell == "CD8Tcell" ~ 9,
                        max$cell == "Bcell" ~ 8,
                        max$cell == "PDC" ~ 7,
                        max$cell == "Mono" ~ 6,
                        max$cell == "Ery" ~ 5,
                        max$cell == "MDC" ~ 4,
                        max$cell == "Neutrophils" ~ 3,
                        max$cell == 'Basophils' ~2,
                        max$cell == "platelet" ~ 1)
    max$cell = with(max, reorder(cell,num))
    
    max$type='x'
    max$plot='plot2'
    max_pb_cb$plot='plot1'
    max_bind <- rbind(max_pb_cb,max)
    
    myColors <- colors()[1:length(unique(max_bind$transcript_id))]
    names(myColors) <- unique(max_bind$transcript_id)
    colScale <- scale_colour_manual(name = "transcript_id",values = myColors)
    
    max_bind$transcript_id <- as.factor(max_bind$transcript_id)
    
    p1=( ggplot(subset(max_bind, plot == "plot1"), aes(y=pct, x=cell,fill=transcript_id))
         + geom_bar(stat="identity",lwd=0.5, aes(color=transcript_id))
         + scale_colour_hue(drop=F, aesthetics = c("colour", "fill")) + coord_flip() + pretty_plot() + facet_wrap(~type)
    )
    
    p2=( ggplot(subset(max_bind, plot == "plot2"), aes( y=pct, x=cell,fill=transcript_id))
         + geom_bar(stat="identity",lwd=0.5, aes(color=transcript_id))
         + scale_colour_hue(drop=F, aesthetics = c("colour", "fill")) + coord_flip() + pretty_plot()
    )
    p2
    
  }
  
  output$download <- downloadHandler(
    filename = "output.pdf", 
    content = function(file) {
      ggsave(filename= file, width = 7, height = 5)
      print(erythroid_line_fcn())
      print(erythroid_line_fcn())
      dev.off()
    })
  output$download2 <- downloadHandler(
    filename = "output.pdf", 
    content = function(file) {
      ggsave(filename= file, width = 7, height = 5)
      print(erythroid_line_fcn())
      print(erythroid_line_fcn())
      dev.off()
    })
  output$plot <- renderPlot({
    erythroid_line_fcn()
  })
  output$plot2 <- renderPlot({
    hsc_differentiation_pathway_fcn()
  })
}
ui <- fluidPage(
  useShinyjs(),
  title = "gene_isoform",
  hr(),
  div(id = "form",
      selectInput(inputId = 'genes',
                  label = 'Select Gene:',
                  choices = genes, selected = 'BCL11A'),
      downloadButton(outputId = 'download', label = 'Save'),
      plotOutput('plot2'),
      downloadButton(outputId = 'download2', label = 'Save'),   
      plotOutput('plot')
      

      ))

shinyApp(ui = ui, server = server)
