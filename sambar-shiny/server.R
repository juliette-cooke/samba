library(shiny)
library(tools)
library(stringr)
library(data.table)
library(sambar)
library(ggplot2)
library(dplyr)
library(DT)

options(shiny.maxRequestSize = 50*1024^2)

# Single table in
# function(input, output){
#   sdata = reactive({
#     req(input$sdatain)
#     dataFile = data.table::fread(input$sdatain$datapath)
#     return(dataFile)
#   })
#   output$tables = renderDataTable(sdata())
# }

# Multiple tables in
function(input, output, session){
  # Take in the data
  sdata = reactive({
    if (input$exchk)  {
      withProgress(message = "Loading...", value = 0,{
      sdata = load_sampling_results("/home/juliette/these/code/git/samba/data/", "xaa")
      print("done")
      incProgress(1/2)
      for (i in 1:length(sdata)){
        output[[paste0('T', i)]] = renderDataTable({sdata[[i]]}, options = list(
          pageLength = 5))
      }})
      sdata
    } else{
      req(input$sdatain)
      ntables = length(input$sdatain[, 1])
      sdata = list()
      names_l = c()
      for(intable in 1:ntables){
        name.split = unlist(str_split(file_path_sans_ext(input$sdatain$name[intable], compression = T), "_"))
        name = if (name.split[length(name.split)] != "WT") "MUT" else "WT"
        names_l = c(names_l, name)
        sdata[[intable]] = fread(input$sdatain[[intable, 'datapath']])
        output[[paste0('T', intable)]] = 
          renderDataTable({sdata[[intable]]}, options = list(
            pageLength = 5))
      }
      names(sdata) = names_l
      sdata
      print("Done")
    }
  })
  
  # Disable file upload if using example data
  observeEvent(input$exchk, {
    if (input$exchk) {
      shinyjs::disable("sdatain")
      sdata()
      } 
    else {
      shinyjs::enable("sdatain")}
  })
  
  # Necessary to load sdatain
  observe({
    req(input$sdatain)
    sdata()})
  
  # Datatables to show once you've uploaded files
  observeEvent(list(input$datachk, input$exchk), {
    if (input$datachk) show("tables") else hide("tables")
  })
  output$tables = renderUI({
    req(isTruthy(sdata()) || isTruthy(input$exchk))
    tagList(lapply(1:length(sdata()), function(i) {
      dataTableOutput(paste0('T', i))
    }))
  })
  
  # Metab dict
  m_dict = reactive({
    if (input$reconchk)  {
      read.csv("/home/juliette/these/code/git/samba/test_data/Recon-2_from_matlab_metab_id_name.tsv", sep = "\t")
    } else{
      read.csv(input$metabdict$datapath, sep = "\t")
    }
  })
  observeEvent(input$reconchk, {
    if (input$reconchk) shinyjs::disable("metabdict") else shinyjs::enable("metabdict")
  })
  observeEvent(input$metabchk, {
    if (input$metabchk) show("metabdicttable") else hide("metabdicttable")
  })
  output$metabdicttable = renderDataTable({
    req(isTruthy(input$metabdict) || isTruthy(input$reconchk))
    m_dict()
  })
  
  
  # Calculate diff and zscore
  diff_z = eventReactive(input$calcbutton, {
    diff = calc_diff(sdata())
    zscore = calc_zscore(diff)
    gc()
    list("diff"=diff, "zscore"=zscore)
    })

  # Zscore table and histogram
  output$zscoretable = renderDataTable(
    {
      metab_map = as.vector(m_dict()$Name)
      names(metab_map) = m_dict()$ID
      datatable(diff_z()$zscore %>% mutate(abszscore = round(abs(zscore), 2))%>%
                  mutate(zscore = round(zscore,2)) %>%
                  mutate(Metab = str_replace_all(string=diff_z()$zscore$Metab, pattern= fixed(metab_map))),
                options = list(
                  autoWidth = T,
                  scrollX=T
                ), rownames = F) %>%
      formatStyle(columns="zscore",
      backgroundColor = styleInterval(c(-input$thrslider, input$thrslider), c("red", "white", "red"))
    )
  })
  
  output$zscoredistrib = renderPlot({
    plot_zscore_distrib(diff_z()$zscore, m_dict(), input$thrslider, binwidth = input$binwidthslider,
                        labels = input$labelschk)
  })

  # Plot distribs
  distrib = eventReactive(list(input$calcbutton, input$thrslider, input$maxn), {
    plot_distrib(sdata(), diff_z()$zscore, thr=input$thrslider, max_n = input$maxn)
  })
  output$distplots = renderPlot({
    distrib()
    })

  # Plot scatterplot
  ranges <- reactiveValues(x = NULL, y = NULL)
  # scatter = eventReactive(list(input$calcbutton, input$thrslider, outlierd(), input$outlierct), {
  #   means = calc_sdata_means(sdata())
  #   plot_scatter(sdata(), diff_z()$zscore, metab_dict = m_dict(), means, thr=input$thrslider, 
  #                outlier.dist = outlierd(),outlier.count = input$outlierct)+
  #     coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand=F)
  # })
  # 
  # output$scatterplot = renderPlot({
  #   scatter()
  # })
  output$scatterplot = renderPlot({
    means = calc_sdata_means(sdata())
    plot_scatter(sdata(), diff_z()$zscore, metab_dict = m_dict(), means, thr=input$thrslider,
                   outlier.dist = outlierd(),outlier.count = input$outlierct)+
        coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand=T)
    })

  
  outlierd = reactiveVal(3)
  observeEvent(input$scatter_dblc, {
    brush = input$scatter_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })

  observeEvent(list(input$autozoom, input$outlierdist), {
    if (input$autozoom){ 
      show("outlierdist")
      show("outlierct")
      outlierd(input$outlierdist)
      # If autozoom is toggled on, reset the graph coord:
      ranges$x <- NULL
      ranges$y <- NULL
      }
    
    else{
      
      hide("outlierdist")
      hide("outlierct") 
      outlierd(10000)
    }
  })
  
  # Commented for now as it can create an infinite update loop
  # observeEvent(input$thr_numeric,{
  #   updateSliderInput(session,"thrslider", value = input$thr_numeric)
  # }, ignoreInit = T)
  # 
  # observeEvent(input$thrslider,{
  #   updateNumericInput(session,"thr_numeric", value = input$thrslider )
  # }, ignoreInit = T)

}
