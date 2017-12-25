db <- read.csv('db.csv')
library(shinycssloaders) #loading spinners
library(ape) #GC, read.GenBank
library(plyr) #rbind.fill.matrix
library(seqinr) #count, comp, c2s
library(DT) #nmers
library(shinythemes) #for theme
library(shinydashboard) #base count bars
library(plotly) #pie chart
library(Biostrings) #matchPattern

findORFs <- function(sequence){
  # 1. find all start and stop codons
  start <- 'atg'
  stops <- c('taa','tag','tga')
  startCodons <-  attr(attr(matchPattern(start,sequence),'ranges'),'start') #locations of start codons
  stopCodons <- numeric()
  for (i in 1:3)
    stopCodons <- append(stopCodons, sort(attr(attr(matchPattern(stops[i],sequence),'ranges'),'start')))
  
  # 2. find start codons in 3 reading frames
  frameOneStarts = numeric();  frameTwoStarts = numeric();  frameThreeStarts = numeric()
  for (i in 1:length(startCodons)){
    if (startCodons[i] %% 3 == 1)
      frameOneStarts <- append(frameOneStarts,startCodons[i])
    else if (startCodons[i] %% 3 == 2)
      frameTwoStarts <- append(frameTwoStarts,startCodons[i])
    else if (startCodons[i] %% 3 == 0)
      frameThreeStarts <- append(frameThreeStarts,startCodons[i])
  }
  
  # 3. find stop codons in 3 reading frames
  frameOneStops = numeric();  frameTwoStops = numeric();  frameThreeStops = numeric()
  for (i in 1:length(stopCodons)){
    if (stopCodons[i] %% 3 == 1)
      frameOneStops <- append(frameOneStops,stopCodons[i])
    else if (stopCodons[i] %% 3 == 2)
      frameTwoStops <- append(frameTwoStops,stopCodons[i])
    else if (stopCodons[i] %% 3 == 0)
      frameThreeStops <- append(frameThreeStops,stopCodons[i])
  }
  
  # 4. determine potential open reading frames
  currentFrameStarts = numeric(); currentFrameStops = numeric(); 
  currentReadingFrames = numeric() #this vector will contain ORFs as pairs of start & stop locations
  openReadingFrames <- matrix(); #final matrix to be returned minus initial temporary row
  for (frame in 1:3){
    if (frame == 1) {currentFrameStarts<-frameOneStarts; currentFrameStops<-frameOneStops}
    else if (frame == 2) {currentFrameStarts<-frameTwoStarts; currentFrameStops<-frameTwoStops}
    else {currentFrameStarts<-frameThreeStarts; currentFrameStops<-frameThreeStops}
    
    #loop over start codons in current frame
    for (i in 1:length(currentFrameStarts)){
      #first check that the start codon is past the previous stop codon to avoid frame overlapping
      if(  (length(currentReadingFrames) == 0) || (currentFrameStarts[i] > tail(currentReadingFrames,1))  ){
        #now loop over the stop codons within the same frame
        for(j in 1:length(currentFrameStops)){
          #check that the stop codon is past the current start codon
          if(currentFrameStops[j] > currentFrameStarts[i]){
            currentReadingFrames <- append(currentReadingFrames,c(currentFrameStarts[i],currentFrameStops[j]))
            break
          }
        } # end of loop for stops
      }
    } # end of loop for starts
    openReadingFrames <- rbind.fill.matrix(openReadingFrames, t(currentReadingFrames))
    currentReadingFrames <- numeric()
  } # end of loop for frames
  return(openReadingFrames[-1,])
}

progressBar <- function(value = 0, label = FALSE, color = "aqua", size = NULL,
                        striped = FALSE, active = FALSE, vertical = FALSE) {
  stopifnot(is.numeric(value))
  if (value < 0 || value > 100)
    stop("'value' should be in the range from 0 to 100.", call. = FALSE)
  if (!(color %in% shinydashboard:::validColors || color %in% shinydashboard:::validStatuses))
    stop("'color' should be a valid status or color.", call. = FALSE)
  if (!is.null(size))
    size <- match.arg(size, c("sm", "xs", "xxs"))
  text_value <- paste0(value, "%")
  if (vertical)
    style <- htmltools::css(height = text_value, `min-height` = "2em")
  else
    style <- htmltools::css(width = text_value, `min-width` = "2em")
  tags$div(
    class = "progress",
    class = if (!is.null(size)) paste0("progress-", size),
    class = if (vertical) "vertical",
    class = if (active) "active",
    tags$div(
      class = "progress-bar",
      class = paste0("progress-bar-", color),
      class = if (striped) "progress-bar-striped",
      style = style,
      role = "progressbar",
      `aria-valuenow` = value,
      `aria-valuemin` = 0,
      `aria-valuemax` = 100,
      tags$span(class = if (!label) "sr-only", text_value)
    )
  )
}

progressGroup <- function(text, value, valueLabel, min = 0, max = value, color = "aqua", striped = T, active = T, label =T) {
  stopifnot(is.character(text))
  stopifnot(is.numeric(value))
  if (value < min || value > max)
    stop(sprintf("'value' should be in the range from %d to %d.", min, max), call. = FALSE)
  tags$div(
    class = "progress-group",
    tags$span(class = "progress-text", text),
    tags$span(class = "progress-number", sprintf("%s", valueLabel)),
    progressBar(round(value / max * 100), color = color, striped = striped, active = active, label = label)
  )
}
