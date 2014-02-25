library(shiny)
require(data.table)
require(knitr)
require(rCharts)

# Define UI for dataset viewer application
shinyUI(pageWithSidebar(
  
  # 
  headerPanel("MIMOSA"),
  
  # Sidebar with controls to select a dataset and specify the number
  # of observations to view
  sidebarPanel(
    fileInput('file1','Choose Data File', accept=c("text/csv","text/comma-separated-values,text/plain",".csv")),
    tags$hr(),
    ### File selector for the metadata
    
    fileInput('file2','Choose Metadata File', accept=c("text/csv","text/comma-separated-values,text/plain",".csv")),
    tags$hr()    
  ),
  
  # Show a summary of the dataset and an HTML table with the requested
  # number of observations
  mainPanel(
    h3("Data Headers \n"),
    htmlOutput('selectheader'),
    h3("Treatment Headers \n"),
    htmlOutput('rxheaders')
      )
))