#########################################
# Title: gffcompareShiny
# Author: Harry Smith
# Description: This is a program for an 
# RShiny app that takes as input a gtf
# from gffcompare and summuarizes/
# visualizes the class codes using ggplot2
#########################################

library(shiny)

## Sect 1: UI
ui <- fluidPage()

## Sect 2: SERVER
server <- function(input, output) {}

## Sect 3: KNIT
shinyApp(ui = ui, server = server)


