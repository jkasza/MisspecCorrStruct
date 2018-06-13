################################################################### 
# This file contains the server for the Shiny app 
###################################################################

library(shiny)
library(ggplot2)
library(reshape2)
library(scales)

source("MisspecCorrStruct.R", local=TRUE)



shinyServer(function(input, output, session) {
  
output$ratiovarplot <- renderPlot({
 vartreat_versus3_plot(input$design, input$m, input$rho0, input$model)
})

mymatrix <- reactiveValues(data=NULL) 

observe({
  req(input$file1)
  mymatrix$data <- read.csv(input$file1$datapath)
})

observeEvent(input$reset, {
  mymatrix$data <- NULL
  reset('file1')
})


output$mymat_ratiovarplot <- renderPlotly({
  if(!is.null(mymatrix$data))   mymat <- mymatrix$data
 # else mymat <- SWdesmat(7)
  myplot <- vartreat_versus3_plotDESMAT(mymat, input$m, input$rho0, input$model)
  myplot
  })
  
output$Plotexplan <- renderUI({
        HTML(paste("The plot displays V&#770",
        tags$sub(1), "/V",tags$sub(3),  " or  V&#770",tags$sub(2), 
        "/V",tags$sub(3), ". V&#770", tags$sub(1), " is the variance of the 
        treatment effect estimator obtained using Model 1 if the expected values of 
        the variance components obtained using the Model 1 ANOVA formulas
        (with expectation under Model 3, the autoregressive model) 
        are used to estimate variance  components. V",tags$sub(3), " is the variance of 
        the treatment effect estimator under Model 3 using the correct Model 3 
        within-cluster correlation structure and the true value of the 
        decay parameters. Results are displayed for a range of
        study lengths and  a range of decay values.", sep=""))
})

output$text1 <- renderText({ 
  "An example of a stepped wedge design matrix:"
})
output$text2 <- renderText({ 
  "An example of a parallel design matrix:"
})
output$text3 <- renderText({ 
  "An example of a parallel w/ baseline design matrix:"
})
output$text4 <- renderText({ 
  "An example of a CRXO design matrix:"
})


output$SWxmat <- renderTable({
  head(SWdesmat(7), n=7)
},digits=0)
output$pllelxmat <- renderTable({
  head(plleldesmat(7), n=7)
},digits=0)
output$pllelbxmat <- renderTable({
  head(pllelBLdesmat(7), n=7)
},digits=0)
output$crxoxmat <- renderTable({
  head(CRXOdesmat(7), n=7)
},digits=0)



output$Contactdetails <- renderUI({
  HTML(paste("This Shiny app accompanies the paper &quot; Inference for the treatment 
             effect in multiple-period cluster randomised trials when random 
             effect correlation structure is misspecified &quot; by Jessica Kasza and 
            Andrew Forbes. For questions or comments, please contact Jessica Kasza: 
            jessica.kasza &quot; at &quot; monash.edu", sep=""))
})
                                
  
})