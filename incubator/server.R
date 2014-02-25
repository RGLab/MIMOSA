require(rCharts)
library(shiny)
require(data.table)
require(xtable)

# Define server logic required to summarize and view the selected dataset
shinyServer(function(input, output,session) {
  data <- reactive({
    inFile <- input$file1
    if(is.null(inFile)){
      return(NULL)
    }
    fread(inFile$datapath)
  })
  rx <- reactive({
    inFile <- input$file2
    if(is.null(inFile)){
      return(NULL)
    }
    fread(inFile$datapath)
  })
  
  output$rxheaders<-renderUI({
    tbl<-rx()
    if(is.null(tbl)){
      tbl<-matrix(nrow=0,ncol=1);
      colnames(tbl)<-""
    }
    selectInput(inputId='rxheaders',label='Select Treatment Headers',choices=colnames(tbl),multiple=TRUE)
  })
  
  output$selectheader<-renderUI({
    tbl<-data()
    if(is.null(tbl)){
      tbl<-matrix(nrow=0,ncol=1);
      colnames(tbl)<-""
    }
    selectInput(inputId='selectheaders',label='Select Data Headers',choices=colnames(tbl),multiple=TRUE)   
  })
})