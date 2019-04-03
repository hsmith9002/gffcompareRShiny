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
library(rtracklayer) 
library(GenomicRanges)
options(stringsAsFactors = F)
options(dplyr.width = Inf)
'%!in%' <- function(x,y)!('%in%'(x,y))
options(shiny.maxRequestSize=500*1024^2)

## Sect 1: UI
ui <- fluidPage(
  ## Other
  headerPanel("Class Code Shiny Proof of Concept"),
  helpText("Note: The class codes summarized in the 
           table do not include 'intergenic'unknown' comparisons
           as defined by gffcompare. All other codes 
           are used exactly as defined on gffcompare manual."),
  ## Input functions
  fileInput(inputId = "file1", 
            label = "Upload gtf",
            accept = ".gtf"),
  actionButton("result1","Generate Result"),
  ## Output functions
  plotOutput(outputId = "bar"),
  tableOutput("table"),
  downloadButton("downloadplot", label = "Download plot")
)

## Sect 2: SERVER
server <- function(input, output, session) {
  ## generate output bar chart 
  
  ## prepare plot objects
  data <- reactiveValues()
  observeEvent(input$result1, {
    inFile <- input$file1
  if (is.null(inFile))
    return(NULL)
  path <- inFile$datapath
  
  ############################
  #import gtf
  ############################
  
  progress <- Progress$new(session, min=0, max=10)
  on.exit(progress$close())
  
  progress$set(message = 'Importing and parsing gtf file using RTracklayer',
               detail = 'This may take a few minutes...')   
  
  df <- import(path)
  df <- as.data.frame(df)
  
  ############################
  #Summarize class codes
  ############################
  
  ## Filter exon rows
  irfilt <- df %>% 
    filter(type == "transcript")
  ## Generate data sets for intergenic and nonintergenic transcripts
  interg <- irfilt %>% filter(class_code == "u")
  noninterg <- irfilt %>% filter(class_code != "u")
  if(sum(dim(interg)[1], dim(noninterg)[1]) != dim(irfilt)[1]) stop("Dimension of filtered data sets do not equal original")
  ## Summarize class codes for non-intergenic transcripts and return
  noninterg$class_code <- factor(noninterg$class_code, levels = c("=", "c", "k", "m", "n", "j", "e", "o", "s", "x", "i", "y", "p", "r", "u", "z"), 
                                 labels = c("Exact match", "Contained in reference", "Containment of reference", 
                                            "Full intron chain match", "partial intron chain match", "At least 1 exon junction match", 
                                            "single exon partially covering intron", "other same strand overlap", "intron match on opposite strand",
                                            "exon overlap on opposite strand", "fully contained within reference intron", "contains a reference within introns",
                                            "no overlap", "repeat", "intergenic", "z"))
  out <- as.data.frame(table(noninterg$class_code))
  levs <- factor(c("=", "c", "k", "m", "n", "j", "e", "o", "s", "x", "i", "y", "p", "r", "u", "z"),
                 levels = c("=", "c", "k", "m", "n", "j", "e", "o", "s", "x", "i", "y", "p", "r", "u", "z"), 
                 labels = c("Exact match", "Contained in reference", "Containment of reference", 
                            "Full intron chain match", "partial intron chain match", "At least 1 exon junction match", 
                            "single exon partially covering intron", "other same strand overlap", "intron match on opposite strand",
                            "exon overlap on opposite strand", "fully contained within reference intron", "contains a reference within introns",
                            "no overlap", "repeat", "intergenic", "z"))
  labs <- c("Exact match", "Contained in reference", "Containment of reference", 
            "Full intron chain match", "partial intron chain match", "At least 1 exon junction match", 
            "single exon partially covering intron", "other same strand overlap", "intron match on opposite strand",
            "exon overlap on opposite strand", "fully contained within reference intron", "contains a reference within introns",
            "no overlap", "repeat", "intergenic", "z")
  out.levs <- out$Var1[which(out$Var1 %in% levs)]
  out.labs <- out$Var1[which(out$Var1 %in% labs)]
  out$Var1 <- factor(out$Var1, levels = out.levs, labels = out.labs)
  out$prop <- (out$Freq/dim(irfilt)[1])
  out$Perc <- percent((out$Freq/dim(irfilt)[1]))
  ## Check to make sure the number of class codes in summary sums to the total number of non-intergenic transcripts
  if(sum(out$Freq) != dim(noninterg)[1]) stop("Total class codes do not sum to total non-intergenic txts")
  out <- as.data.frame(out)
  output$table <- renderTable(out[-length(out$Var1), ])
  ############################
  #Prepare for plotting
  ############################
  
  bnccc_man <- out %>% 
    filter(Var1 %in% c("Exact match", "Contained in reference", 
                       "Containment of reference", "At least 1 exon junction match")) %>%
    mutate(Var1 = factor(as.character(Var1), 
                         levels = levs, 
                         labels = labs)) 
  ccvec <- setdiff(c("Exact match", "Contained in reference", 
                     "Containment of reference", "At least 1 exon junction match"), bnccc_man$Var1)
  bnccc_tmp <- data.frame(Var1 = ccvec, 
                          Freq = as.numeric(rep(0, length(ccvec))), 
                          prop = as.numeric(rep(0, length(ccvec))),
                          Perc = percent(as.numeric(rep(0, length(ccvec)))))
  
  bnccc_other <- c("z", sum(out$Freq[which(out$Var1 %!in% c("Exact match", "Contained in reference", 
                                                            "Containment of reference", "At least 1 exon junction match"))]),
                   sum(out$prop[which(out$Var1 %!in% c("Exact match", "Contained in reference", 
                                                       "Containment of reference", "At least 1 exon junction match"))]),
                   percent(sum(out$prop[which(out$Var1 %!in% c("Exact match", "Contained in reference", 
                                                               "Containment of reference", "At least 1 exon junction match"))])))
  bnccc_man <- rbind(bnccc_man, bnccc_tmp, bnccc_other)
  bnccc_man <- bnccc_man %>%
    mutate(prop = as.numeric(prop))
  bnccc_test <- bnccc_man[order(levs[which(levs %in% bnccc_man$Var1)]), ]
  
  ############################
  #Generate plot
  ############################
  
  data$plot <- ggplot(bnccc_man, aes(x=Var1, y=prop, fill=Var1)) +
    geom_bar(stat="identity") +
    scale_y_continuous(labels=percent, 
                       limits=c(0,0.70)) + 
    scale_x_discrete(labels=c("Exact match", 
                              "Contained in \n reference", 
                              "Containment of \n reference", 
                              "At least 1 \n exon junction \n match",
                              "Ambiguous overlap")) +
    scale_fill_manual('leged', values = brewer.pal(n = 5, name = "Set2")) +
    theme(legend.position='none',
          plot.title = element_text(hjust=0, size = 20)) +
    geom_text(data=bnccc_man, 
              aes(label=Perc[order(bnccc_man$Var1)], 
                  y=prop[order(bnccc_man$Var1)]+0.015,
                  x = (seq(1,5,1))+0.085), 
              size=5.5) +
    xlab("Class Code") + 
    ylab(paste0("Number of Transcripts non-intergenic: ",
                prettyNum(sum(as.numeric(bnccc_man$Freq)), 
                          big.mark = ",", 
                          big.interval = 3))) 
  })
  
  output$bar <- renderPlot({
    data$plot
  }
  )
  ## Download plot
  output$downloadplot <- downloadHandler(
    filename = function(){
      paste("classCode", ".png", sep = "")
    },
    content = function(file) {
     ggsave(file, data$plot)
    }
  )
}

## Sect 3: KNIT
shinyApp(ui = ui, server = server)


