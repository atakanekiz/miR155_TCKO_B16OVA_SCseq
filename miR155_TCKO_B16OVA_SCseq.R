library(Seurat)
library(dplyr)
library(shiny)
  

combined <- readRDS("clustered_combined_w_tsne.rds")
all_genes <- sort(rownames(combined@data))
all_samples <- levels(combined@meta.data$cond)

d9wt_cells <- WhichCells(combined, subset.name = "cond", accept.value = "WT(d9)")
d9ko_cells <- WhichCells(combined, subset.name = "cond", accept.value = "miR-155 TCKO(d9)")
d12wt_cells <- WhichCells(combined, subset.name = "cond", accept.value = "WT(d12)")
d12ko_cells <- WhichCells(combined, subset.name = "cond", accept.value = "miR-155 TCKO(d12)")


########################################################################################################################
##############################################        UI         #######################################################
########################################################################################################################

ui <- fluidPage(
  
  titlePanel("SCseq of tumor infiltrating immune cells in murine melanoma"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      selectInput(inputId = "gene",
                  label = "Select a gene name",
                  choices = all_genes, 
                  selected = "Mir155hg",
                  multiple = F),
      
      checkboxGroupInput(inputId = "sample",
                         label = "Select samples to include in plot",
                         choices = all_samples,
                         selected = all_samples),
      
      sliderInput("expr_range", "Change color scaling",
                 min = 0, max = 0,
                 value = c(0, 0))
      
    ),
    
    mainPanel(
      
      tabsetPanel(
        
        tabPanel("tSNE Plots",
                 
                 fluidRow(
                 column(8, plotOutput(outputId = "featureplot",width = 600, height = 500)),
                 column(8, imageOutput(outputId = "clusters", width = 500, height = 500))
                 )
               
                 
        ),
        
        tabPanel("Violin Plots",
                 
                 plotOutput(outputId = "vlnplot")
                 
        )
      )
      
    )
    
  )
  
)


########################################################################################################################
##############################################      Server       #######################################################
########################################################################################################################

server <- function(input, output, session){
  
  
  # min_expr <- reactive({
  #   
  #   req(input$gene)
  #   
  #   round(min(combined@data[input$gene,]),2)
  #   
  # })
  # 
  # 
  # max_expr <- reactive({
  #   
  #   req(input$gene)
  #   
  #   round(max(combined@data[input$gene,]),2)
  #   
  # })
  
  observe({
    
    min_expr <- round(min(combined@data[input$gene,]),2)
    max_expr <- round(max(combined@data[input$gene,]),2)
    
    updateSliderInput(session, inputId = "expr_range", 
                      value = c(min_expr, max_expr), min = min_expr, max = max_expr)
    
  })
 
  

  output$featureplot <- renderPlot({
    
    validate(need(input$sample, message = "Select one or more samples to plot"))
    
    
    cells <- WhichCells(combined, subset.name = "cond", accept.value = input$sample)
    

    FeaturePlot(combined, features.plot = input$gene, cells.use = cells, pt.size = 2, 
                min.cutoff = input$expr_range[1], max.cutoff = input$expr_range[2], no.legend = F,
                cols.use = c("gray90","darkblue"))
    
  })
  
  output$clusters <- renderImage(
    
    {filename <- "Cluster_definitions.png"
    list(src=filename, width="auto", height = "auto")
    
    }, deleteFile = F
  
  )

  output$vlnplot <- renderPlot({
    
    cells <- WhichCells(combined, subset.name = "cond", accept.value = input$sample)
    
    VlnPlot(combined, features.plot = input$gene, cells.use = cells, x.lab.rot = T, group.by = "ident")
    
    
    
  })
  
  
  
}

shinyApp(ui=ui, server=server)






