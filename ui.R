################################################################### 
# This file contains the UI for the Shiny app 
###################################################################

library(shiny)
library(plotly)

shinyUI(pageWithSidebar(
  h3("Estimating variance components when random effect correlation structure is misspecified."),
  sidebarPanel(
    selectInput("design", label = ("Which cluster randomised design would you like to consider?"), 
                choices = list("Stepped wedge" = 1, "CRXO" = 4, "Parallel" = 2, "Parallel with baseline" =5), selected = 1),
    selectInput("model", label = ("Compare Model 3 (correlation decay model) to:"), 
                choices = list("Model 1 (one-way model)" = 1, "Model 2 (two-way model)" = 2), selected = 1),
    
    fileInput('file1', 'Upload a design matrix instead:',
              accept=c('text/plain', '.txt', '.csv')),
    helpText("The file must be a comma separated .csv or .txt file consisting of 0s and 1s, with a column for each time period. Do not include row or column names."),
    actionButton('reset', 'Clear file'),
    
    numericInput("m",
                 "Number of subjects/cluster-period:",
                 min = 1,
                 max=1000,
                 step = 1,
                 value = 100),
    sliderInput("rho0",
                label = HTML(paste("Intra-cluster correlation, &rho;",tags$sub(0),":", sep="")),
                #"rho 0:"
                min = 0.001,
                max = 0.2,
                step = 0.001,
                value = 0.05),
    
    helpText(HTML(paste("Note: calculations assume &sigma;",tags$sup(2),tags$sub(HTML(paste("1&epsilon;"))), "+", "&sigma;",
                        tags$sup(2),tags$sub(HTML(paste("1&alpha;"))), " = 1",  " so  &rho; = &sigma;",tags$sup(2),tags$sub(HTML(paste("1&alpha;"))), "/(&sigma;",tags$sup(2),tags$sub(HTML(paste("1&alpha;"))), "+", "&sigma;",
                        tags$sup(2),tags$sub(HTML(paste("1&epsilon;"))), ")  = &sigma;",tags$sup(2),tags$sub(HTML(paste("1&alpha;"))), sep="")))
    
    
  ),
  mainPanel(
    tabsetPanel(
      tabPanel("Ratio of variances of treatment effect estimator", value=1, plotOutput("ratiovarplot")),
      htmlOutput("Plotexplan"),
      tabPanel("Design Schematics", value=2, textOutput("text1"),
               tableOutput("SWxmat"),
               textOutput("text4"),
               tableOutput("crxoxmat"),
               textOutput("text2"),
               tableOutput("pllelxmat"),
               textOutput("text3"),
               tableOutput("pllelbxmat")),
      tabPanel("Ratios for user-input design", value=3,  plotlyOutput("mymat_ratiovarplot")),
      tabPanel("Ref. and contact details", value=4, htmlOutput("Contactdetails")) )
  )
))

