
library(shiny)

# Define UI for slider demo application
shinyUI(pageWithSidebar(
  
  #  Application title
  headerPanel("Rapa Nui Demography"),
  
  # Sidebar with sliders that demonstrate various available options
  sidebarPanel(
    
    tags$h3("Model Parameters"),
    
    sliderInput("N0", "Initial Population (Ratio of Carrying Capacity):", 
                min=0.001, max=1, value=0.05, step=0.01),
    
    sliderInput("r", "Instrinsic Growth Rate:", 
                min=0.001, max=0.1, value=0.01, step=0.001),
    
    sliderInput("b1", "Pollen Cover Parameter", 
                min=-0.03, max=0.03, value=0, step=0.001),
    
    sliderInput("b2", "SOI Cover Parameter:", 
                min=-1, max=1, value=0, step=0.01),

  ),
  
  # Show a table summarizing the values entered
  mainPanel(
    plotOutput("graph1"),
    plotOutput("graph2"),
    plotOutput("graph3"),
  )
))