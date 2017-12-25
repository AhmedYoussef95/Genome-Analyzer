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
    else if (input$tabs == 'Open Reading Frames')
      helpText('An open reading frame is a potential protein-coding region within a DNA sequence. 
               Open reading frames begin with the start codon ATG and end with one of three different stop codons.
               Each DNA sequence has six different reading frames; three for the forward strand and likewise for the reverse strand.
               NOTE: Plotting open reading frames for larger genomes may take a while.')
    
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
        chunkGC <- GC(chunk)*100 #convert GC ratio to %
        chunkGCs[i] <- chunkGC
      }
      plot(starts,chunkGCs,type="b",xlab="Nucleotide position",ylab="GC content (%)",main = "GC Content")
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
  
  output$orf <- renderPlot({
    #findORFs is defined in global.R
    forwardStrand <- findORFs(c2s(sequence()))
    reverseStrand <- findORFs(c2s(rev(comp(sequence()))))
    openReadingFrames <- rbind.fill.matrix(forwardStrand,reverseStrand) #matrix that contains the 6 ORFs
    
    #set up plot with 6 segments
    plot(c(1,length(sequence())), c(0,0), ylim=c(0,6), type="l", axes=F, 
           xlab="Nucleotide Position", ylab = '', main="Predicted open reading frames")
     segments(1,1,length(sequence()),1)
     segments(1,2,length(sequence()),2)
     segments(1,3,length(sequence()),3)
     segments(1,4,length(sequence()),4)
     segments(1,5,length(sequence()),5)
     text(0,0.6,"Frame 1"); text(0,1.6,"Frame 2"); text(0,2.6,"Frame 3");
     text(0,3.6,"Frame 4"); text(0,4.6,"Frame 5"); text(0,5.6,"Frame 6");
     axis(1, pos=0)
    
    # determine lengths of open reading frames
     orflengths = numeric()
    for (frame in 1:6){
      v <- openReadingFrames[frame,]
      v <- v[!is.na(v)]
      for(i in seq(1,length(v)-1,2))
        orflengths <- append(orflengths,(v[i+1]-v[i]))
        #rect(v[i],frame-1,v[i+1],frame-0.2,col="orange",border="red")
    }
     
    # plot frames whose lengths fall in the 95th percentile (more probable candidates)
     for (frame in 1:6){
       v <- openReadingFrames[frame,]
       v <- v[!is.na(v)]
       for(i in seq(1,length(v)-1,2)){
         if(  (v[i+1]-v[i]) > quantile(orflengths, probs=c(0.95))  )
           rect(v[i],frame-1,v[i+1],frame-0.6,col="#EE7600",border="black")
       }
     }
       
      })
}
