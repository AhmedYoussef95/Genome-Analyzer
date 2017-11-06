server <- function(input, output){
  
  sequence <- reactive({
    acc <- db[db$organism == input$organism,2]
    read.GenBank(acc, as.character = TRUE)[[1]]
  })
  
  output$descrip <- renderUI({
    helpText(db[db$organism == input$organism,3])
  })
  
  output$`nmer table` <- DT::renderDataTable({
    t <- sort(count(sequence(),input$nmers), decreasing = T)
    names(t) <- toupper(names(t))
    for (i in 1:length(t))
      t[i] <- prettyNum(t[i],",")
    t <- data.frame(t)
    DT::datatable(t) %>% formatStyle(columns = c(1,2),color = 'black')
    }, options = list(columnDefs = list(list(targets=c(2),orderable=F))) )
  
  output$length <- renderUI({
    l <- paste("Genome Length = ",prettyNum(length(sequence()), ",")," bp")
    progressGroup(l,length(sequence()),"",1,length(sequence()))
    #p(strong(paste("Genome Length = ",prettyNum(length(sequence()), ",")," bp")),class = "text-center")
  })
  
  output$Acount <- renderUI({
    count <- count(sequence(),1)[1]
    progressGroup("A:", count, prettyNum(count,','), 1, length(sequence()))
  })
  
  output$Ccount <- renderUI({
    count <- count(sequence(),1)[2]
    progressGroup("C:", count, prettyNum(count,','), 1, length(sequence()))
    })
  
  output$Gcount <- renderUI({
    count <- count(sequence(),1)[3]
    progressGroup("G:", count, prettyNum(count,','), 1, length(sequence()))
  })
  
  output$Tcount <- renderUI({
    count <- count(sequence(),1)[4]
    progressGroup("T:", count, prettyNum(count,','), 1, length(sequence()))
  })
  
  output$gc <- renderPlot({
    window <- round(length(sequence())/(input$window))
    starts <- seq(from = 1, to = length(sequence())-window, by = window)
    n <- length(starts)
    chunkGCs <- numeric(n) # Make a vector of the same length as vector "starts", but just containing zeroes
    for (i in 1:n) {
      chunk <- sequence()[starts[i]:(starts[i]+window-1)]
      chunkGC <- GC(chunk)*100
      chunkGCs[i] <- chunkGC
    }
    plot(starts,chunkGCs,type="b",xlab="Nucleotide start position",ylab="GC content (%)",main = "GC Content")
    abline(h = mean(chunkGCs), col = "#EE7600", lty = 2)
    text(x = starts[2],y = mean(chunkGCs),pos=3, col = '#EE7600', labels = paste('mean',round(mean(chunkGCs),1),'%'))
  })
}