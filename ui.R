shinyUI( 
  fluidPage(
  theme = shinytheme('superhero'),
  titlePanel(title=div(img(src="dna_icon.png", height = 40, width = 35),'Genome Analyzer'), windowTitle = 'Genome Analyzer'),
  sidebarLayout(
    sidebarPanel(
     # helpText('Display statistics for genomes of different organisms'), hr(),
      selectInput('organism','Organism',choices = db$organism), hr(),
      p(strong('Description')),
      uiOutput('descrip'), hr(),
      tags$div(class="header", checked=NA,
               tags$p("Source code:"),
               tags$a(href="https://github.com/AhmedYoussef95/Genome-Analyzer", "github.com/AhmedYoussef95")
      ), hr(),
      helpText(p(em('©️ Ahmed Youssef')))
      #textInput('nationality','Nationality:', value = 'All'),
      #checkboxGroupInput('foot','Preferred foot:',choices = c('Left','Right'),selected = c('Left','Right'))
    ),
    
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Base Counts", box(uiOutput("length"),uiOutput('Acount'), uiOutput('Ccount'), uiOutput('Gcount'), uiOutput('Tcount'))),
                  tabPanel("GC Content", tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: #EE7600}")), sliderInput('window','Window Size', min = 5, max=500, value = 20,width = "100%"),withSpinner(plotOutput('gc', height = 500),type = getOption("spinner.type", default = 8))),
                  tabPanel("N-mers", numericInput('nmers','N-mer',value = 2,min = 2,max = 6),DT::dataTableOutput('nmer table'),style = "background-color: #7e8b9a;")
                  )
    )
  )
  
))