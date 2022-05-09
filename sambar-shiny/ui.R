library(shiny)
library(shinydashboard)
library(dashboardthemes)
library(shinyjs)
library(shinyWidgets)
library(shinycssloaders)
library(dplyr)
library(shinyBS)
library(flexdashboard)
library(bsplus)
source("theme.R", local = T)$value

# css <- "
# #reverseSlider .irs-bar, .irs-bar-edge{
#     border-top: 1px solid #ddd;
#     border-bottom: 1px solid #ddd;
#     background: linear-gradient(to bottom, #DDD -50%, #FFF 150%);
# }
# #reverseSlider .irs-line, .irs-bar{
#     background: #db2c35;
#     border: 1px solid #db2c35;
# }
# "

header <- dashboardHeader(title = "SAMBA: Sampling Biomarker Analysis")

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Data", tabName = "data", icon = icon("database"))
    ,
    menuItem("Calculate & Plot", tabName = "calc", icon = icon("chart-area"))
    ,
    menuItem("About", tabName = "about", icon = icon("info"))
    )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "data",
            h2("Data set-up"),
            fluidPage(theme = "styles.css",
                      fluidRow(
                        column(6, 
                               fluidRow(
                                 column(10,fileInput(inputId = "sdatain", 
                                                     label = "Upload sampling data here"
                                                     , multiple = T) %>% shinyInput_label_embed(
                                                       icon("info-circle") %>%
                                                         bs_embed_tooltip(title = 
                                                                            "Sampling data consists of two csv files with columns as metabolites and rows as
                                                            samples. The WT file must end with _WT and the mutant/disease file must end with something that
                                                            isn't WT (such as _KO or _MUT). Can be .csv or .csv.gz", placement = "bottom")
                                                     ),
                                        )
                                 ),
                               materialSwitch("exchk", label = "Use example data", value = F, status="primary", 
                                              right = T),
                               materialSwitch("datachk", label = "Show table previews",value = F, status = "primary", 
                                              right = T),
                               useShinyjs(),
                               uiOutput("tables")),
                        column(6, fileInput(inputId = "metabdict", label = "Upload a metabolite dictionary", multiple = F)
                               %>% shinyInput_label_embed(
                                 icon("info-circle") %>%
                                   bs_embed_tooltip(title = 
                                   "A metabolite dictionary is a .tsv file generated using sambapy from the corresponding model. 
                                   It contains an ID column (corresponds to the IDs in the sampling results) and a Name column (
                                   the names you want plotted). It is purely for plotting purposes.", placement = "bottom")
                               ),
                               materialSwitch("reconchk", label = "Use Recon 2.2 dictionary", value = F, status="primary", 
                                              right = T),
                               materialSwitch("metabchk", label = "Show metabolite dictionary preview", value = F, 
                                              right = T, status = "primary"),
                               dataTableOutput("metabdicttable")
                               ))
                      ),
            )
    ,
    tabItem(tabName = "calc",
            h2("Calculate and plot differences between each WT and mut condition"),
            fluidPage(
              actionButton(inputId = "calcbutton", label = "Calculate and plot!"),
              fluidRow(column(5,
                              tags$style(
                                HTML(
                                  ".js-irs-0 .irs-single
                                 {
                                    background: #db2c35;
                                 }
                                .js-irs-0 .irs-bar, .js-irs-0 .irs-bar-edge{
                                    border-top: 1px solid #ddd;
                                    border-bottom: 1px solid #ddd;
                                    background: linear-gradient(to bottom, #DDD -50%, #FFF 150%);
                                }
                                .js-irs-0 .irs-line{
                                    background: #db2c35;
                                    border: 1px solid #db2c35;
                                }
                                ")
                              ),
                              sliderInput(inputId = "thrslider",label = "Threshold",
                                   min = 0, max = 10, value = 1, step = 0.1)%>% 
                                shinyInput_label_embed(
                                     icon("info-circle") %>%
                                       bs_embed_tooltip(title = "Threshold for z-score cut-off. A higher threshold will 
                                                        select only the strongest metabolite changes.", placement = "bottom")
                                   ))
                       # ,column(2,numericInput(inputId = "thr_numeric", label = "Threshold",value = 1, min =0,
                       #              max = 10))
                       )
              
              ,
              tabsetPanel(
                tabPanel(title = "Distributions", 
                         numericInput(inputId = "maxn", label = "Maximum number of metabolites to show", value = 8, min = 1, 
                                      step = 1),
                         plotOutput("distplots")%>% withSpinner(color="#0dc5c1")),
                tabPanel(title = "Scatterplot", 
                         fluidRow(
                           column(4,materialSwitch("autozoom", label = "Automatic zoom", value = F, right = T, status = "primary")),
                           column(4,numericInput(inputId = "outlierdist", label = "Outlier distance", value = 3, 
                                                        min = 0)%>% 
                                    shinyInput_label_embed(
                                      icon("info-circle") %>%
                                        bs_embed_tooltip(title = "The distance from the upper quantile over which a 
                                        metabolite is considered an outlier (for the autozoom).", placement = "bottom")
                                    )),
                           column(4,numericInput(inputId = "outlierct", label = "Outlier count", value = 0, 
                                                        min = 0)%>% 
                                    shinyInput_label_embed(
                                      icon("info-circle") %>%
                                        bs_embed_tooltip(title = "The minimum number of outliers there needs to be for 
                                                         autozoom to be added.", placement = "bottom")
                                    ))),
                         plotOutput("scatterplot", dblclick = "scatter_dblc", brush = brushOpts(
                           id = "scatter_brush",
                           resetOnNew = TRUE
                         ))%>% withSpinner(color="#0dc5c1"),
                         ),
                tabPanel(title = "Z-scores",
                         fixedRow(column(5, dataTableOutput("zscoretable")),
                                  column(5, sliderInput("binwidthslider", label = "Bin width",
                                                        min = 0.1, max = 5, value = 0.1, step = 0.1),
                                         materialSwitch("labelschk", label = "Show labels",value = F, 
                                                           status = "primary", right = T),
                                         plotOutput("zscoredistrib")))
                         )
                )
              )
            ),
    
    tabItem(tabName = "about",
            h2("About SAMBA"),
            tags$p("Code for SAMBApy and SAMBAR can be found here:"),
            tags$a(href = "https://forgemia.inra.fr/metexplore/cbm/samba", "https://forgemia.inra.fr/metexplore/cbm/samba"))
    )
  
  )

dashboardPage(skin = "purple", header, sidebar, body)