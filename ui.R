shinyUI( 
  fluidPage(
  theme = shinytheme('superhero'),
  titlePanel(title=div(img(src="dna_icon.png", height = 40, width = 35),'Genome Analyzer'), windowTitle = 'Genome Analyzer'),
  sidebarLayout(
    sidebarPanel(
    #  helpText('Compute and visualize statistics for genomes of different organisms'), hr(),
      selectInput('organism','Organism', selectize = F,
                  list("Virus" = db$organism[1:2],"Human" = db$organism[3:4], "Bacteria" = db$organism[5:8], "Other" = db$organism[9:11])),
      textInput('acc','Or enter GenBank accession number'),hr(),
      p(strong('Description')),
      uiOutput('descrip'), hr(),
      uiOutput('stat label'),
      uiOutput('stat'), hr(),
      tags$div(class="header", checked=NA,
               tags$p("Source code:"),
               tags$a(href="https://github.com/AhmedYoussef95/Genome-Analyzer", "github.com/AhmedYoussef95")
      ), hr(),
      helpText(p(em('©️ Ahmed Youssef')))
      ),
    
    mainPanel(
      tabsetPanel(type = "tabs", id = "tabs",
                  tabPanel("Base Counts", box(uiOutput("length"),uiOutput('Acount'), uiOutput('Ccount'), uiOutput('Gcount'), uiOutput('Tcount')), box(plotlyOutput('pie'))),
                  tabPanel("GC Content", tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: #EE7600}")), sliderInput('window','Number of Windows', min = 5, max=300, value = 20,width = "100%"),withSpinner(plotOutput('gc', height = 500),type = getOption("spinner.type", default = 8))),
                  tabPanel("N-mers", numericInput('nmers','Number of consecutive nucleotides (N)',value = 2,min = 2,max = 6),DT::dataTableOutput('nmer table'),style = "background-color: #7e8b9a;"),
                  tabPanel("Open Reading Frames", withSpinner(plotOutput('orf', height = 600),type = getOption("spinner.type", default = 8)))
                  )
    )
  )
  
))
