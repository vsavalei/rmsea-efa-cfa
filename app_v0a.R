library(shiny)
library(plotly) 
library(DT)
library(dplyr)

source("functions_app.R")

#--------------------------------------------------------------------------------------------------------------------#
# UI code

# for adjusting the overlapping anchors on the slider 
js1 <- paste("function MRdoesOverlap() {",
             "   var $lastLabel = $('#sliderange2 .irs-grid-text:last');", # MR 
             "   var $prevLastLabel = $lastLabel.prevAll('.irs-grid-text').first();",
             "   return $lastLabel.offset().left < $prevLastLabel.offset().left + $prevLastLabel.width();",
             "}\n",
             "Shiny.addCustomMessageHandler('regrid', function(force) {",
             "   if (MRdoesOverlap() | force) {",
             "      console.log('Overlap detected - adjusting tick number');",
             "      var $sld = $('#range').data('ionRangeSlider');", #range is the input name 
             "      var ticks_n = $sld.options.grid_num;",
             "      $sld.update({grid_num: Math.round(ticks_n)});",
             "   }",
             "});",sep = "\n")


ui <- fluidPage(
  tags$head(tags$script(HTML(js1), type = "text/javascript")), #corrsponds js (for adjusting the overlapping anchors )


  
  h1("SEM Fit Indices' Sensitivity to cross-loadings"),
  h4("How sensitive are RMSEA, CFI, and SRMR to omitted cross-loadings? This app will generate population covariance matrices 
	from CFA models with 2 or 3 factors and with a varying number of crossloadings. The app will then fit the CFA model 
	with	no cross-loadings to these population matrices and compute the population RMSEA, CFI, and SRMR."),
  h6("This app was developed by", HTML(paste0(a("Victoria Savalei", href="http://ubcsemlab.com/"))), "and Lihan (Bill) Chen"),
  
  sidebarLayout(
    sidebarPanel(
      "Factor loadings are drawn from uniform distributions on [a,b] with mean L=(a+b)/2 and range R=(b-a). 
				 The user can set ranges (MR and CR) to zero to specify constant main loadings and constant cross-loadings. 
				 Variables' variances are assumed to be 1. For this reason, 
				 the maximum possible value of a crossloading is determined by the corresponding value of the main loading and the factor correlation. 
				 The upper limit on the choice of the average CL value below is computed from the average ML value specified above. This is done in 
				 an attempt to minimize impossible configurations (with negative error variances). However, 
				 impossible configurations can still be generated unless both MR and CR are set to zero. Check the printed residual 
				 variances for negative values. Plots will be omitted for all such configurations.",  
      br(),br(), 
      
      radioButtons("custom", "I would like to examine a ", 
                   c("Two factor model"="twoFactor","Three factor model"="threeFactor")),
      
      
      conditionalPanel(
        condition = "input.custom == 'twoFactor'",
        sliderInput("p2", "Total number of variables:", min=4, max=50, step=2, value=8), 
        
        sliderInput("fcor2", "Factor correlation in the true model:", min=0, max=1, value=.2), 
        sliderInput("aveloading2", "Average Main Loading (ML)", min=0, max=1,value=.7),
        uiOutput("sliderange2"), #MR
        uiOutput("slidemax_cross2"), #CL
        uiOutput("sliderange_cross2"), #CR
        
      ),
      
      conditionalPanel(
        condition = "input.custom == 'threeFactor'",
        sliderInput("p3", "Total number of variables:", min=6, max=51, step=3, value=9), 
        
        sliderInput("fcor3", "Factor correlation in the true model:", min=0, max=1, value=.2), 
        sliderInput("aveloading3", "Average Main Loading (ML)", min=0, max=1,value=.7),
        uiOutput("sliderange3"), #MR
        uiOutput("slidemax_cross3"), #CL
        uiOutput("sliderange_cross3"), #CR
        
      ),
      
      # conditionalPanel(
      #   condition = "input.custom == 'spec'",
      #   strong("Specify the main factor loading for each variable:"),br(),
      #   uiOutput("slideload"),actionButton("updateButtonSpec", "Plot!!")
      # ),
      
      actionButton("updateButton", "Compute fit index values!")
    ),
    
    mainPanel(  
        p(textOutput("plottext")),
        plotlyOutput("plots"),
        br(),
        uiOutput("tabletext"),
        br(),
        dataTableOutput("table")
      # conditionalPanel(
      #   condition = "input.custom == 'twoFactor'", 
      # 
      #   # h4(textOutput("printload")),h4(textOutput("printload_cross")),
      #   # h4(textOutput("printres")),h4(textOutput("printres2")),
      #   
      #   p(textOutput("plottext")), 
      #   plotlyOutput("plots"),
      # 
      #   uiOutput("tabletext"),
      #   dataTableOutput("table"),
      # ),
      
      # conditionalPanel(
      #   condition = "input.custom == 'threeFactor'",
      # 
      #   # h4(textOutput("printload")),h4(textOutput("printload_cross")),
      #   # h4(textOutput("printres")),h4(textOutput("printres2")),
      # 
      #   p(textOutput("plottext")),
      #   plotlyOutput("plots"),
      # 
      #   uiOutput("tabletext"),
      #   dataTableOutput("table"),
      # )
      
      # conditionalPanel(
      #   condition = "input.custom == 'spec'",
      #   h4(textOutput("printloadspec")),
      #   br(),
      #   h4(textOutput("plottextspec")),
      #   plotlyOutput("plotspec"),
      #   h4(textOutput("plottextspecadd"))
      # )
      
    )
  )
)

#--------------------------------------------------------------------------------------------------------------------#
# helper function for repeated identical text formatting
renderLoadingRecap <- function(loadtxt, loadings, p, digits=3){
  return(renderText({
    tnoend <- paste(sapply(loadings[1:(p-1)], function(x){
      sprintf(paste0('%.', digits, 'f'), x)}), collapse=", ")
    tend <- paste(sapply(loadings[p], function(x){
      sprintf(paste0('%.', digits, 'f'), x)}), collapse=", ")
    paste(loadtxt, tnoend," and ", tend, ".", sep="")
  }))
}

#--------------------------------------------------------------------------------------------------------------------#
# server code
server <- function(input, output, session) {
  
  session$onFlushed(function() {
    session$sendCustomMessage("regrid", FALSE);
  }, FALSE);
  
  # nsdl = 3 #mumber of sign digits
  # pm = 10 #cutoff for p to stop displaying some output to prevent clutter  
  
  #define input sliders for TWO factor model
  
  output$slidemax_cross2 <- renderUI({
    sliderInput("aveloading_cross2", "Average Cross Loading (CL)", min = 0, 
                max = round((sqrt(1-(input$aveloading2)^2+(input$fcor2*input$aveloading2)^2)-input$fcor2*input$aveloading2),2), 
                round = -2, step = 0.01, value = .2) 
  })
  
  #for randomly generated factor loadings:
  output$sliderange2 <- renderUI({
    sliderInput("range2", "Main Loadings Range (MR)", min = 0, max = round((min(2*input$aveloading2, 2*(1-input$aveloading2))),2), 
                value = min(.1,input$aveloading2, (1-input$aveloading2)), round = -3, step = 0.01) 
  })
  
  #for randomly generated factor cross-loadings:
  output$sliderange_cross2 <- renderUI({
    sliderInput("range_cross2", "Cross-Loadings Range (CR)", 
                min = 0, 
                max = round(min(1,2*input$aveloading_cross2,  
                                2*(sqrt(1-(input$aveloading2)^2+(input$fcor2*input$aveloading2)^2)-input$fcor2*input$aveloading2-input$aveloading_cross2)),2), 
                value = min(0,input$input$aveloading_cross2, (1-input$input$aveloading_cross2)), round = -2, step = 0.01) 
  })
  
  #define input sliders for THREE factor model
  output$slidemax_cross3 <- renderUI({
    sliderInput("aveloading_cross3", "Average Cross Loading (CL)", min = 0, 
                max = round((sqrt(1-(input$aveloading3)^2+(input$fcor3*input$aveloading3)^2)-input$fcor3*input$aveloading3),2), 
                round = -2, step = 0.01, value = .2) 
  })
  
  #for randomly generated factor loadings:
  output$sliderange3 <- renderUI({
    sliderInput("range3", "Main Loadings Range (MR)", min = 0, max = round((min(2*input$aveloading3, 2*(1-input$aveloading3))),2), 
                value = min(.1,input$aveloading3, (1-input$aveloading3)), round = -3, step = 0.01) 
  })
  
  #for randomly generated factor cross-loadings:
  output$sliderange_cross3 <- renderUI({
    sliderInput("range_cross3", "Cross-Loadings Range (CR)", 
                min = 0, 
                max = round(min(1,2*input$aveloading_cross3,  
                                2*(sqrt(1-(input$aveloading3)^2+(input$fcor3*input$aveloading3)^2)-input$fcor3*input$aveloading3-input$aveloading_cross3)),2), 
                value = min(0,input$input$aveloading_cross3, (1-input$input$aveloading_cross3)), round = -2, step = 0.01) 
  })
  
  # #for user defined factor loadings: 
  # numP <- reactive({input$p2})  # get number of loadings from ui slider 
  # 
  # output$slideload <- renderUI({   # create custom loading sliders
  #   lapply(1:(numP()), function(i) {
  #     sliderInput(inputId = paste0("load", i), # this creates input$load1, input$load2, etc.
  #                 label = paste0("Loading ", i),
  #                 min = 0, max = 1, value = .5)
  #   })
  # })
  
  # button pressed to update randomly generated loadings
  # code needed
  
  observeEvent(input$updateButton,{      
    
    pSwitch <- switch(input$custom,
                      twoFactor = input$p2,
                      threeFactor = input$p3)
    
    numCrossLoadingSwitch <- switch(input$custom,
                                    twoFactor = input$p2,
                                    threeFactor = input$p3 * 2)
    
    aveloadingSwitch <- switch(input$custom,
                               twoFactor = input$aveloading2,
                               threeFactor = input$aveloading3)
    
    rangeSwitch <- switch(input$custom,
                          twoFactor = input$range2,
                          threeFactor = input$range3)
    
    range_crossSwitch <- switch(input$custom,
                          twoFactor = input$range_cross2,
                          threeFactor = input$range_cross3)
    
    aveloading_crossSwitch <- switch(input$custom,
                          twoFactor = input$aveloading_cross2,  
                          threeFactor = input$aveloading_cross3 )         
    
    fcorSwitch <- switch(input$custom,
                         twoFactor = input$fcor2,
                         threeFactor = input$fcor3)
    

    # loadtxt <- "The true model is a 2-factor model with the randomly generated main loadings of " 
    genLoadingss <- runif(pSwitch, min=aveloadingSwitch-.5*rangeSwitch, max=aveloadingSwitch+.5*rangeSwitch)
    
    
    # loadtxt_cross <- "The randomly generated values of crossloadings to be added to the true model, one by one, are " 
    numCrossLoading <- runif(numCrossLoadingSwitch, min=aveloading_crossSwitch -.5*range_crossSwitch, max=aveloading_crossSwitch+.5*range_crossSwitch)  #### CHANGED
    
    
    # Allow switching the main function between two-factor model and three-factor model
    mainFunc <- switch(input$custom,
                       twoFactor = main.2f,
                       threeFactor = main.3f)
    #replace rmseas with results
    results <- as.data.frame(mainFunc(isolate(pSwitch),isolate(fcorSwitch),genLoadingss,numCrossLoading))
    #print(tail(results, n=1))
    
    print(select(tail(results, n=1),matches("x[0-9]{1,2}same")))
    
    print(select(tail(results, n=1),matches("x[0-9]{1,2}dif")))
    
    genResiduals_2f_seq <- as.vector(select(tail(results, n=1),matches("x[0-9]{1,2}same")))%>%
      unname()%>%
      unlist()
   
    genResiduals_2f_alt <- as.vector(select(tail(results, n=1),matches("x[0-9]{1,2}dif")))%>%
      unname()%>%
      unlist()
 
  
    
    # genResiduals_2f_seq <- 1-genLoadingss^2-numCrossLoading^2-2*genLoadingss*numCrossLoading*fcorSwitch # formula for 2f model 
    # 
    # genLoadingss_cross_2f_reordered <- c(numCrossLoading[c(TRUE, FALSE)], numCrossLoading[c(FALSE, TRUE)]) ###Is this a typo? The order of c(TRUE, FALSE)
    # 
    # genResiduals_2f_alt <- 1-genLoadingss^2-genLoadingss_cross_2f_reordered^2-2*genLoadingss*genLoadingss_cross_2f_reordered*fcorSwitch 
    
    
    # the table output
    table2f <- as.data.frame(rbind(
      #paste(LETTERS[1:8]),
      round(genLoadingss, 3),
      round(numCrossLoading, 3),
      round(genResiduals_2f_seq, 3),
      round(genResiduals_2f_alt, 3)
    ))
    
    row.names(table2f) <- c("Main loadings","Cross loadings","Residual vairances: to 1st factor, then 2nd","Residual vairances: to alternating factor")

    table3f <- as.data.frame(rbind(
      round(genLoadingss, 3),
      round(numCrossLoading, 3),
      round(genResiduals_2f_seq, 3),
      round(genResiduals_2f_alt, 3)
    ))
    row.names(table3f) <- c("Main loadings","Cross loadings","Same1 and Same2","Diff1 and Diff2")
    
    # To customize the table header for 3-factor model 
    sketch = htmltools::withTags(table(
      class = 'display',
      thead(
        tr(
          th(rowspan = 2, 'Order added'),
          th(class = 'dt-center', colspan = pSwitch, 'Same'),
          th(class = 'dt-center', colspan = pSwitch, 'Different')
        ),
        tr(
          lapply(paste0(rep("P", numCrossLoadingSwitch), c(1:numCrossLoadingSwitch)), th)
        )
      )
    ))
    
    if (input$custom =="threeFactor"){
      output$table <- renderDataTable(datatable(table3f,
                                                container = sketch,
                                                colnames = F,
                                                caption = "Parameter values when all cross-loadings are added",  
                                                options = list( scrollX = T, # to add a horizontal scroller in case of having a wide table 
                                                                dom = 't', # to hide the table's filer and search function 
                                                                ordering = F # to hide the ordering function 
                                                                
                                                                # # to hide the header of the table as it's irrevalent
                                                                # headerCallback = JS("function(thead, data, start, end, display){",
                                                                #                     "  $(thead).remove();",
                                                                #                     "}")
                                                )))
    }

    
    if (input$custom =="twoFactor"){
      output$table <- renderDataTable(
        datatable(
          table2f,
          colnames = paste0(rep("P", numCrossLoadingSwitch), c(1:numCrossLoadingSwitch)),
          caption = "Parameter values when all cross-loadings are added",
          options = list(
            scrollX = T,
            # to add a horizontal scroller in case of having a wide table
            dom = 't',
            # to hide the table's filer and search function
            ordering = F # to hide the ordering function
            
            # # to hide the header of the table as it's irrevalent
            # headerCallback = JS("function(thead, data, start, end, display){",
            #                     "  $(thead).remove();",
            #                     "}")
          )
        )
      )
    }
    
    output$tabletext <- renderUI({
      HTML(paste(strong(em("Main loadings"))," for a ", numText, " are randomly generated." ,
                 strong(em("Crossloadings"))," are randomly generated and added to the true modelm, one by one.",
                 strong(em("Residual vairances: to 1st factor, then 2nd")),"means that when these loadings are added first to one factor." ,
                 strong(em("Residual vairances: to alternating factor")), " means that when these loadings are added to alternating factors.")
      )
    })
    
    # store numerical value 2 or 3 
    numText <- switch(input$custom,
                       twoFactor = "2-factor model",
                       threeFactor = "3-factor model")
  
    #define text
    output$plottext <- renderText({ 
      paste("In the plots below, the number of crossloadings in the true model is on the x-axis; this number 
	          varies from 0 to ", isolate(pSwitch), " (the number of variables). The cross-loadings in the true 
	          model are being added either a) to the first factor first and then to the second factor
	          or b) to the alternating factors. Their exact values are given by the list above. (If some 
	          values are missing from the plots, one of the residual variances is negative or the model 
	          failed to converge). The fitted model is a ", numText ," with no cross-loadings. 
	          You can hover over the curve to get specific fit index values.", sep="") 
    })
    

    if (input$custom =="threeFactor"){
      output$plots <- renderPlotly({
        
        upperbound_rmsea = max(c(results$rmsea_same1_f,results$rmsea_dif1_f,results$rmsea_same2_f,results$rmsea_dif2_f,0.08)) + 0.005
        upperbound_srmr = max(c(results$srmr_same1_f,results$srmr_dif1_f,results$srmr_same2_f,results$srmr_dif2_f,0.08)) + 0.005
        lowerbound_cfi = min(c(results$cfi_same1_f,results$cfi_same2_f,results$cfi_dif1_f,results$cfi_dif2_f,0.9))-0.005

        plot1 <- ggplot(data=results,aes(x=number_crossloadings)) +
          geom_line(mapping = aes( y=rmsea_same1_f, color="Same1"))+
          suppressWarnings(geom_point(aes(y=rmsea_same1_f,
                         color="Same1",
                         text = paste0("# of cross Loadings: ", number_crossloadings, "<br>RMSEA: ", sprintf('%.3f', rmsea_same1_f)))))+
          
          geom_line(mapping = aes( y=rmsea_same2_f, color="Same2"))+
          suppressWarnings(geom_point(aes(y=rmsea_same2_f,
                         color="Same2",
                         text = paste0("# of cross Loadings: ", number_crossloadings, "<br>RMSEA: ", sprintf('%.3f', rmsea_same2_f)))))+
          
          geom_line(mapping = aes( y=rmsea_dif1_f, color="Diff1"))+
          suppressWarnings(geom_point(aes(y=rmsea_dif1_f,
                         color="Diff1",
                         text = paste0("# of cross Loadings: ", number_crossloadings, "<br>RMSEA: ", sprintf('%.3f', rmsea_dif1_f)))))+
          
          geom_line(mapping = aes( y=rmsea_dif2_f, color="Diff2"))+
          suppressWarnings(geom_point(aes(y=rmsea_dif2_f,
                         color="Diff2",
                         text = paste0("# of cross Loadings: ", number_crossloadings, "<br>RMSEA: ", sprintf('%.3f', rmsea_dif2_f)))))+
          
          #geom_abline(color="grey",slope=0, intercept=0.05) + 
          geom_abline(color="grey",slope=0, intercept=0.08) +  
          labs(color = "How crossloadings are added") + 
          xlab("Number of crossloadings in the true model")+
          ylab("RMSEA")+
          ylim(NA,upperbound_rmsea)
        p1 <- ggplotly(plot1,tooltip = c("text"))  %>% style(showlegend = FALSE)
        
        plot2 <- ggplot(data=results,aes(x=number_crossloadings)) +
          geom_line(mapping = aes( y=cfi_same1_f, color="Same1"))+
          suppressWarnings(geom_point(aes(y=cfi_same1_f,
                         color="Same1",
                         text = paste0("# of cross Loadings: ", number_crossloadings, "<br>CFI: ", sprintf('%.3f', cfi_same1_f)))))+
          
          geom_line(mapping = aes( y=cfi_same2_f, color="Same2"))+
          suppressWarnings(geom_point(aes(y=cfi_same2_f,
                         color="Same2",
                         text = paste0("# of cross Loadings: ", number_crossloadings, "<br>CFI: ", sprintf('%.3f', cfi_same2_f)))))+
          
          geom_line(mapping = aes( y=cfi_dif1_f, color="Diff1"))+
          suppressWarnings(geom_point(aes(y=cfi_dif1_f,
                         color="Diff1",
                         text = paste0("# of cross Loadings: ", number_crossloadings, "<br>CFI: ", sprintf('%.3f', cfi_dif1_f)))))+
          
          geom_line(mapping = aes( y=cfi_dif2_f, color="Diff2"))+
          suppressWarnings(geom_point(aes(y=cfi_dif2_f,
                         color="Diff2",
                         text = paste0("# of cross Loadings: ", number_crossloadings, "<br>CFI: ", sprintf('%.3f', cfi_dif2_f)))))+
          
          geom_abline(color="grey",slope=0, intercept=0.90) +
          labs(color = "How crossloadings are added") + 
          xlab("Number of crossloadings in the true model")+
          ylab("CFI")+
          ylim(lowerbound_cfi,NA)
        p2 <- ggplotly(plot2,tooltip = c("text"))  %>% style(showlegend = FALSE)
        
        plot3 <- ggplot(data=results,aes(x=number_crossloadings)) +
          geom_line(mapping = aes( y=srmr_same1_f, color="Same1"))+
          suppressWarnings(geom_point(aes(y=srmr_same1_f,
                         color="Same1",
                         text = paste0("# of cross Loadings: ", number_crossloadings, "<br>SRMR: ", sprintf('%.3f', srmr_same1_f)))))+
          
          geom_line(mapping = aes( y=srmr_same2_f, color="Same2"))+
          suppressWarnings(geom_point(aes(y=srmr_same2_f,
                         color="Same2",
                         text = paste0("# of cross Loadings: ", number_crossloadings, "<br>SRMR: ", sprintf('%.3f', srmr_same2_f)))))+
          
          geom_line(mapping = aes( y=srmr_dif1_f, color="Diff1"))+
          suppressWarnings(geom_point(aes(y=srmr_dif1_f,
                         color="Diff1",
                         text = paste0("# of cross Loadings: ", number_crossloadings, "<br>SRMR: ", sprintf('%.3f', srmr_dif1_f)))))+
          
          geom_line(mapping = aes( y=srmr_dif2_f, color="Diff2"))+
          suppressWarnings(geom_point(aes(y=srmr_dif2_f,
                         color="Diff2",
                         text = paste0("# of cross Loadings: ", number_crossloadings, "<br>SRMR: ", sprintf('%.3f', srmr_dif2_f)))))+
          
          geom_abline(color="grey",slope=0, intercept=0.90) +
          labs(color = "How crossloadings are added") + 
          xlab("Number of crossloadings in the true model")+
          ylab("SRMR")+
          ylim(NA,upperbound_srmr)
        p3 <- ggplotly(plot3,tooltip = c("text"))  #%>% style(showlegend = FALSE)
        
        subplot(
          p1,p2,p3,
          nrows = 1
        ) %>%
          layout(annotations = list( 
            list( 
              x = 0.151,
              y = 1.0,
              text = "RMSEA",  
              xref = "paper",  
              yref = "paper",  
              xanchor = "center",  
              yanchor = "bottom",  
              showarrow = FALSE 
            ),  
            list( 
              x = 0.5, 
              y = 1,  
              text = "CFI",  
              xref = "paper",  
              yref = "paper",  
              xanchor = "center",  
              yanchor = "bottom",  
              showarrow = FALSE 
            ),  
            list( 
              x = 0.85,  
              y = 1,  
              text = "SRMR",  
              xref = "paper",  
              yref = "paper",  
              xanchor = "center",  
              yanchor = "bottom",  
              showarrow = FALSE 
            ) ))
      })
    }
    
    if (input$custom =="twoFactor"){
      ##### THREE FACTOR MODEL SECTION 
      output$plots <- renderPlotly({
        # compute the range for the plots 
        upperbound_rmsea = max(c(results$rmsea_same_f,results$rmsea_dif_f,0.08)) + 0.005
        upperbound_srmr = max(c(results$srmr_same_f,results$srmr_dif_f,0.08)) + 0.005
        lowerbound_cfi = min(c(results$cfi_same_f,results$cfi_dif_f,0.9))-0.005
        
        print(results)
        
        plot1 <- ggplot(data=results,aes(x=number_crossloadings))+ 
          geom_line(aes(y=rmsea_same_f,
                        color="To 1st factor, then 2nd"))+ 
          suppressWarnings(geom_point(aes(y=rmsea_same_f,
                         color="To 1st factor, then 2nd",
                         text = paste0("# of cross Loadings: ", number_crossloadings, "<br>RMSEA: ", sprintf('%.3f', rmsea_same_f)))))+
          
          geom_line(aes(y=rmsea_dif_f,color ="To alternating factors"))+
          suppressWarnings(geom_point(aes(y=rmsea_dif_f,
                         color ="To alternating factors",
                         text = paste0("# of cross Loadings: ", number_crossloadings, "<br>RMSEA: ", sprintf('%.3f', rmsea_dif_f)))))+
          geom_abline(color="grey",slope=0, intercept=0.08) + labs(color = "How crossloadings are added") +
          xlab("Number of crossloadings in the true model")+
          ylim(NA,upperbound_rmsea) 
        
        p1 <- ggplotly(plot1,tooltip = c("text"))  %>% style(showlegend = FALSE)
        
        plot2 <- ggplot(data=results, aes(x=number_crossloadings))+ 
          geom_line(aes(y=cfi_same_f,
                        color="To 1st factor, then 2nd"))+ 
          suppressWarnings(geom_point(aes(y=cfi_same_f,
                         color="To 1st factor, then 2nd",
                         text = paste0("# of cross Loadings: ", number_crossloadings,
                                       "<br>CFI: ", sprintf('%.3f', cfi_same_f)))))+
          geom_line(aes(y=cfi_dif_f, 
                        color="To alternating factors"))+ 
          suppressWarnings(geom_point(aes(y=cfi_dif_f,
                         color="To alternating factors",
                         text = paste0("# of cross Loadings: ", number_crossloadings,
                                       "<br>CFI: ", sprintf('%.3f', cfi_dif_f)))))+
          geom_abline(color="grey",slope=0, intercept=0.90) + labs(color = "How crossloadings are added") +
          xlab("Number of crossloadings in the true model")+
          # ylab("CFI for the model with no crossloadings")+
          ylim(lowerbound_cfi,NA)
        p2 <- ggplotly(plot2,tooltip = c("text"))  %>% style(showlegend = FALSE)
        
        plot3 <- ggplot(data=results, aes(x=number_crossloadings))+ 
          geom_line(aes(y=srmr_same_f,
                        color="To 1st factor, then 2nd"))+ 
          geom_line(aes(y=srmr_dif_f,color="To alternating factors"))+ 
          suppressWarnings(geom_point(aes(y=srmr_same_f,
                         color="To 1st factor, then 2nd",
                         text = paste0("# of cross Loadings: ", number_crossloadings,
                                       "<br>SRMR: ", sprintf('%.3f', srmr_same_f)))))+
          suppressWarnings(geom_point(aes(y=srmr_dif_f,
                         color="To alternating factors",
                         text = paste0("# of cross Loadings: ", number_crossloadings,
                                       "<br>SRMR: ", sprintf('%.3f', srmr_dif_f)))))+
          geom_abline(color="grey",slope=0, intercept=0.08) + labs(color = "How crossloadings are added") +
          xlab("Number of crossloadings in the true model")+
          # ylab("SRMR for the model with no crossloadings")+
          ylim(NA,upperbound_srmr)
        p3 <- ggplotly(plot3,tooltip = c("text"))
        
        subplot(
          p1,p2,p3,
          nrows = 1
        ) %>%
          layout(annotations = list( 
            list( 
              x = 0.151,
              y = 1.0,
              text = "RMSEA",  
              xref = "paper",  
              yref = "paper",  
              xanchor = "center",  
              yanchor = "bottom",  
              showarrow = FALSE 
            ),  
            list( 
              x = 0.5, 
              y = 1,  
              text = "CFI",  
              xref = "paper",  
              yref = "paper",  
              xanchor = "center",  
              yanchor = "bottom",  
              showarrow = FALSE 
            ),  
            list( 
              x = 0.85,  
              y = 1,  
              text = "SRMR",  
              xref = "paper",  
              yref = "paper",  
              xanchor = "center",  
              yanchor = "bottom",  
              showarrow = FALSE 
            ) ))
        
      })
    }
    
    
    
  }) # observeEvent(input$updateButton
  
  
  
  # 	
  # 	# button pressed to update specified loadings
  # 	observeEvent(input$updateButtonSpec,{
  # 
  # 		loadtxtspec <- "The specified loadings are " 
  # 				
  #  		# specLoadings <- # defaults: we need this before numP() is reset
  # 		all <- rep(.5, numP())
  # 		# check for custom values
  # 		for(i in 1:(numP())){
  # 			# eval parse is ugly, but see comments on slideload
  # 			eval(parse(text=paste0("all[", i, "] <- input$load", i)))
  # 		}
  # 		
  # 		specLoadingss <- all # defaults: we need this before numP() is reset
  # 		
  # 		#relss <- alpha.omega.pop(specLoadingss)
  # 		
  # 		output$printloadspec <- renderText({
  # 				tnoend=paste(as.character(round(specLoadingss[1:(isolate(input$p2)-1)],3)), collapse=", ")
  # 				tend=paste(as.character(round(specLoadingss[isolate(input$p2)],3)), collapse=", ")
  # 				paste(loadtxtspec, tnoend," and ", tend, ".", sep="")
  # 			})
  # 		
  # 		#define text for plot
  # 		basictxt <- "The plot below varies the number of cross-loadings in the true model from 0 to p 
  # 		(the number of variables) on the x-axis. On the y-axis is the RMSEA when a 2-factor model 
  # 		with no cross-loadings is fit to the population covariance matrix 
  # 		generated by the true model. The blue curves show the RMSEA values when the cross-loadings 
  # 		are added to the first factor first (from 0 to p/2), and then to the 
  # 		second factor (p/2+1 to p). The red curves show the RMSEA values when the cross-loadings 
  # 		are adding to alternating factors."
  # 
  # 		output$plottextspec <- renderText({basictxt})
  # 
  # 		addtxt <- "Warning: these values of main loadings and cross-loadings are not very realistic!"
  # 		
  # 		if(input$p2<=pm){
  # 			addtxt = paste0("Single-click on a linetype within the legend to remove/re-add the corresponding plot. ",
  # 										"Double-click on a linetype to remove/re-add all other plots.")
  # 		}  
  # 
  # 		output$plottextspecadd <- renderText({ 
  # 			paste(duptxt, addtxt, sep="") 
  # 		})
  
  
  # 	#make plot
  # 	#relspecx <- as.data.frame(alpha.omega.pop.varyl(specLoadingss)) 
  # 	main.2f(p,fcor,l,cld)
  # 	
  # 	relspecx <- as.data.frame(alpha.omega.pop.varyl(specLoadingss)) 
  # 	
  # 	
  # 	#set up stuff for plotting
  # 	xsta <- relspecx$changing.loading
  # 	ysta <- relspecx$omega.minus.alpha
  # 	xend <- c(xsta[2:length(xsta)],NA)
  # 	yend <- c(ysta[2:length(ysta)],NA)
  # 
  # 	gs <- nrow(relspecx)/input$p2 #gridsize, as set in alpha.omega.pop.varyl
  # 	
  # 	xend[seq(gs,nrow(relspecx),by=gs)]=xend[seq(gs,nrow(relspecx),by=gs)-1] 
  # 	yend[seq(gs,nrow(relspecx),by=gs)]=yend[seq(gs,nrow(relspecx),by=gs)-1] 
  # 
  # 	dimrel <- dim(relspecx)[1]
  # 	
  # 	for(i in 1:(numP())){
  # 		addme=c(alpha.omega.pop(specLoadingss)$alpha,
  # 						alpha.omega.pop(specLoadingss)$omega,
  # 						t(specLoadingss),99,
  # 						alpha.omega.pop(specLoadingss)$omina,i,
  # 						t(specLoadingss[i]))
  # 		relspecx=rbind(relspecx,addme)    
  # 	}
  # 	dimrel2 <- dim(relspecx)[1]
  # 
  # 	#plot difference vs varing loading
  # 	output$plotspec <- renderPlotly({  
  # 
  # 		plt <- ggplot(data=relspecx[1:dimrel,],
  # 							aes(text = paste("omega = ",round(omega,digits=3),
  # 															"<br>alpha = ",round(alpha,digits=3),
  # 															paste0("<br>alpha underestimates omega by ", round(omega-alpha,digits=3),", or ",
  # 																ifelse(omega-alpha==0, 0, round((omega-alpha)*100/omega,digits=2)),"%"),
  # 															paste0("<br>changing loading value (loading ", which.l, ") = ",
  # 																round(changing.loading,digits=2))
  # 						))) +
  # 				geom_segment(data=relspecx[1:dimrel,],
  # 						aes(x=-1, xend=-1, y=0, yend=0, linetype=as.factor(which.l)), color='black') + xlim(0,1) + #hack for solving ggplotly legend color issue
  # 				geom_segment(data=relspecx[1:dimrel,],
  # 						aes(x=changing.loading, xend=xend, y=omega.minus.alpha, yend=yend,
  # 							linetype=as.factor(which.l), colour=omega)) +    
  # 				geom_point(data=relspecx[(dimrel+1):dimrel2,],
  # 						aes(x=changing.loading, y=omega.minus.alpha), colour= "blue") +
  # 				xlab("Changing loading value") + ylab("Difference between omega and alpha") +
  # 				scale_color_gradient2(name="omega",midpoint=.5,low="lightpink", mid="salmon",high="blue") +
  # 				theme(legend.title = element_blank())
  # 		
  # 		if (isolate(input$p2)<=pm){ # show linetype legend
  # 			plt <- ggplotly(plt,tooltip = c("text")) %>%
  # 					add_annotations(text="changing loading",  xref="paper", yref="paper",
  # 						x=0.15, xanchor="right", 
  # 						y=-.3, yanchor="bottom",    
  # 						legendtitle=TRUE, showarrow=FALSE) %>% 
  # 					layout(legend = list(orientation = "h", yanchor="bottom", xanchor="left",y = -.4, x =0.15))
  # 		} else { # hide linetype legend
  # 			plt <- ggplotly(plt,tooltip = c("text")) %>% 
  # 					layout(showlegend = FALSE) 
  # 		}
  # 	})
  # 
  
  #}) # closes observeEvent(input$updateButtonSpec  
  
  
  
} #end of server

# start the app
shinyApp(ui = ui, server = server)