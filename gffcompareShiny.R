#########################################
# Title: gffcompareShiny
# Author: Harry Smith
# Description: This is a program for an 
# RShiny app that takes as input a gtf
# from gffcompare and summuarizes/
# visualizes the class codes using ggplot2
#########################################

library(shiny)
library(tidyverse)
library(reshape2)
library(scales)
library(RColorBrewer)
options(stringsAsFactors = F)
options(dplyr.width = Inf)

## Sect 1: UI
ui <- fluidPage(
  ## Input functions
  selectInput(inputId = "ds", label = "Select a dataset (gtf)")
)

## Sect 2: SERVER
server <- function(input, output) {}

## Sect 3: KNIT
shinyApp(ui = ui, server = server)


