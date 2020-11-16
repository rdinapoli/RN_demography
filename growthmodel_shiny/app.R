library(shiny)
library(here)
load(here('R_imagefiles','variables.RData')) 

# Define UI for application that draws a histogram
ui <- fluidPage(
  pageWithSidebar(
    #  Application title
    headerPanel("Rapa Nui Demography Model"),
    
    # Sidebar with sliders that demonstrate various available options
    sidebarPanel(
      
      tags$h3("Model Parameters"),
      
      sliderInput("N0", "Initial Population (Ratio of Carrying Capacity):", 
                  min=0.001, max=0.5, value=0.05, step=0.01),
      
      sliderInput("r", "Instrinsic Growth Rate:", 
                  min=0.001, max=0.1, value=0.01, step=0.001),
      
      sliderInput("a", "Carrying Capacity Intercept:", 
                  min=0.5, max=1.5, value=1, step=0.001),
      
      sliderInput("b1", "Pollen Cover Parameter", 
                  min=-0.03, max=0.03, value=0, step=0.001),
      
      sliderInput("b2", "SOI Cover Parameter:", 
                  min=-0.5, max=0.5, value=0, step=0.01),
      
      checkboxInput("showK",label = "Show K",value=FALSE),
      
      checkboxInput("showSOI","Show SOI",value=FALSE),
      
      checkboxInput("showPalm","Show Palm Pollen",value=FALSE)
    ),
    
    # Show a table summarizing the values entered
    mainPanel(
      plotOutput("graph")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output){
  
  mydata <- reactive({
    # Model Parameters:
    N0 <- input$N0  
    r <- input$r      
    a  <- input$a     
    b1 <- input$b1  
    b2 <- input$b2  
    showK <- input$showK
    showSOI <- input$showSOI
    showPalm <- input$showPalm
    timeRange <- c(800,150)
    x1 <- palm$PollenPerc
    x2 <- soi$SOIpr
    CalBP = timeRange[1]:timeRange[2]
    tt = 0:length(CalBP)
    pop = numeric(length=length(tt))
    K = numeric(length=length(tt))
    pop[1] = N0
    K[1] = NA
    for (i in 2:length(tt))
    {
      pop[i] = pop[i-1] * exp(r*(1-pop[i-1]/(a+b1*x1[i-1]+b2*x2[i-1])))
      K[i] = (a+b1*x1[i-1]+b2*x2[i-1])
    }
    pop = pop[-1]
    K = K[-1]
    result = data.frame(CalBP=CalBP,Pop=pop,SOI=soi$SOIpr,Palm=palm$PollenPerc,K=K)
    list(result=result,showK=showK,showSOI=showSOI,showPalm=showPalm)
  })
  
  output$graph <- renderPlot({
    if (any(mydata()[[1]]$K<0))
    {
      plot(0,0,xlim=c(800,150),xlab='Cal BP',ylab = 'Population Size',type='l',lwd=2.5,ylim=c(0,1),axes=FALSE)
      axis(1)
      axis(2)
      text(x=500,y=0.5,labels = 'Negative K!',cex=4)
    } else {
    par(mar=c(5,4,2,7))
    plot(mydata()[[1]]$CalBP, mydata()[[1]]$Pop,xlim=c(800,150),xlab='Cal BP',ylab = 'Population Size',type='l',lwd=2.5,ylim=c(0,max(c(1,mydata()[[1]]$Pop,mydata()[[1]]$K))),axes=FALSE)
    axis(1)
    axis(2)
    if (mydata()[[2]]==TRUE){lines(mydata()[[1]]$CalBP, mydata()[[1]]$K,lty=2)}
    if (mydata()[[3]]==TRUE&mydata()[[4]]==FALSE)
    {par(new=TRUE)
      plot(mydata()[[1]]$CalBP, mydata()[[1]]$SOI,xlim=c(800,150),xlab='',ylab = '',axes=FALSE,type='l',col='Darkorange')
      axis(4,col.ticks = 'Darkorange',col = 'Darkorange',col.axis='Darkorange')
      mtext('SOI',side = 4,line=2,col = 'Darkorange')
    }
    if (mydata()[[3]]==FALSE&mydata()[[4]]==TRUE)
    {par(new=TRUE)
      plot(mydata()[[1]]$CalBP, mydata()[[1]]$Palm,xlim=c(800,150),xlab='',ylab = '',axes=FALSE,type='l',col='Darkgreen')
      axis(4,col.ticks = 'Darkgreen',col = 'Darkgreen',col.axis='Darkgreen')
      mtext('Palm Pollen %',side = 4,line=2,col = 'Darkgreen')
    }
    
    
    if (mydata()[[3]]==TRUE&mydata()[[4]]==TRUE)
    {
      par(new=TRUE)
      plot(mydata()[[1]]$CalBP, mydata()[[1]]$SOI,xlim=c(800,150),xlab='',ylab = '',axes=FALSE,type='l',col='Darkorange')
      axis(4,col.ticks = 'Darkorange',col = 'Darkorange',col.axis='Darkorange')
      mtext('SOI',side = 4,line=2,col = 'Darkorange')
      par(new=TRUE)
      plot(mydata()[[1]]$CalBP, mydata()[[1]]$Palm,xlim=c(800,150),xlab='',ylab = '',axes=FALSE,type='l',col='Darkgreen')
      axis(4,col.ticks = 'Darkgreen',col = 'Darkgreen',col.axis='Darkgreen',line=3.5)
      mtext('Palm Pollen %',side = 4,line=5.5,col = 'Darkgreen')
    }
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
