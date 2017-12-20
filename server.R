server <- function(input, output){
  
  sequence <- reactive({
    if (input$acc == "")
      acc <- db[db$organism == input$organism,2]
    else
      acc <- input$acc
    try(read.GenBank(acc, as.character = TRUE, species.names = F)[[1]])
  })
  
  output$descrip <- renderUI({
    if (input$acc == "")
      helpText(db[db$organism == input$organism,3])
    else{
      tryCatch({
        helpText(attr(read.GenBank(input$acc),"description"))
      }, error = function(e){helpText("Invalid accession number")})
    }
  })
  
  output$`stat label` <-renderUI({
    p(strong(input$tabs))
  })
  
  output$stat <- renderUI({
    if (input$tabs == 'Base Counts')
      helpText('All DNA sequences are composed of the four nucleotides: A, C, G, and T.
               The relative abundance of each nucleotide is of importance to many scientific experiments.')
    else if (input$tabs == 'GC Content')
      helpText('The GC content measures the relative abundance of the G-C base pairing along the DNA sequence.
               Protein-coding genes tend to be located at CG-rich regions, making it useful to identify such regions.
               The slider value determines the number of segments to divide the sequence into.')
    else if (input$tabs == 'N-mers')
      helpText('An N-mer is a group of N consecutive nucleotides in a sequence. 
               There are 4^N unique N-mers in any DNA sequence.
               Several N-mers of biological importance across species have been identified, such as the TATA promoter region.')
    
  })
  
  
  
  output$length <- renderUI({
    l <- paste("Genome Length = ",prettyNum(length(sequence()), ",")," bp")
    progressGroup(l,length(sequence()),"",1,length(sequence()))
  })
  
  output$Acount <- renderUI({
    tryCatch({
      count <- count(sequence(),1)[1]
      progressGroup("A:", count, prettyNum(count,','), 1, length(sequence()))
    }, error = function(e){helpText("Invalid accession number")})
  })
  
  output$Ccount <- renderUI({
    tryCatch({
      count <- count(sequence(),1)[2]
      progressGroup("C:", count, prettyNum(count,','), 1, length(sequence()))
    }, error = function(e){helpText("")})
  })
  
  output$Gcount <- renderUI({
    tryCatch({
      count <- count(sequence(),1)[3]
      progressGroup("G:", count, prettyNum(count,','), 1, length(sequence()))
    }, error = function(e){helpText("")})
  })
  
  output$Tcount <- renderUI({
    tryCatch({
      count <- count(sequence(),1)[4]
      progressGroup("T:", count, prettyNum(count,','), 1, length(sequence()))
    }, error = function(e){helpText("")})
  })
  
  output$pie <- renderPlotly({
    ntcounts <- c(count(sequence(),1)[1], count(sequence(),1)[2], count(sequence(),1)[3], count(sequence(),1)[4])
    labels <- c("A", "C", "G", "T")
    plot_ly(labels = labels, values = ntcounts,  textposition = 'inside',
            textinfo = 'label+percent', textfont = list(color = 'white', size = 14)) %>%
      add_pie(hole = 0.5) %>%
      layout( paper_bgcolor='transparent')  
    })
  
  output$gc <- renderPlot({
    tryCatch({
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
    }, error = function(e){helpText("Invalid accession number")})
  })
  
  output$`nmer table` <- DT::renderDataTable({
    t <- sort(count(sequence(),input$nmers), decreasing = T)
    names(t) <- toupper(names(t))
    for (i in 1:length(t))
      t[i] <- prettyNum(t[i],",")
    t <- data.frame(t)
    DT::datatable(t) %>% formatStyle(columns = c(1,2),color = 'black')
  }, options = list(columnDefs = list(list(targets=c(2),orderable=F))) )
}
