library(shiny)
library(here)
load(here('R_imagefiles','variables.RData')) 

# Simulation and Shiny Application of Flue Season Dynamics
shinyServer(function(input, output) {
  
  mydata <- reactive({
    # Model Parameters:
    N0 <- input$N0  
    r <- input$r      
    a  <- 1   
    b1 <- input$b1  
    b2 <- input$b2  
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
    list(result=result)
    
  })
  
  output$datatable <- renderTable(mydata()[["wide"]])
  
  output$graph1 <- renderPlot({
    plot(mydata()[[1]]$CalBP, mydata()[[1]]$Pop,xlim=c(800,150),xlab='Cal BP',ylab = 'Population Size',type='l',lwd=2,ylim=c(0,max(c(1,mydata()[[1]]$Pop,mydata()[[1]]$K))))
    lines(mydata()[[1]]$CalBP, mydata()[[1]]$K,lty=2)
  })

  output$graph2 <- renderPlot({
    plot(mydata()[[1]]$CalBP, mydata()[[1]]$Palm,xlim=c(800,150),xlab='Cal BP',ylab = 'Palm Tree Coverage',type='l',col='Darkgreen')
  })
    
  output$graph3 <- renderPlot({
    plot(mydata()[[1]]$CalBP, mydata()[[1]]$SOI,xlim=c(800,150),xlab='Cal BP',ylab = 'SOI',type='l',lwd=2,col='Lightblue')
  })

})