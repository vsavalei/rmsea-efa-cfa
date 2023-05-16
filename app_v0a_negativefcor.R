library(shiny)
library(plotly) 
library(DT)
library(shinybusy)
library(Matrix)
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
  tags$head(tags$script(HTML(js1), type = "text/javascript"), #corresponds to js (for adjusting the overlapping anchors )
            
            tags$style(
              HTML(".shiny-notification {
             position:fixed;
             top: calc(50%);
             left: calc(50%);
             opacity: 1;
             }"))), 
  
  h1("SEM Fit Indices' Sensitivity to cross-loadings"),
  h4("How sensitive are RMSEA, CFI, and SRMR to omitted cross-loadings? This app will generate population covariance matrices 
	from CFA models with 2 or 3 factors and with a varying number of crossloadings. The app will then fit the CFA model 
	with	no cross-loadings to these population matrices and compute the population RMSEA, CFI, and SRMR."),
  h6("This app was developed by", HTML(paste0(a("Victoria Savalei", href="http://ubcsemlab.com/"))), "and Muhua (Karyn) Huang"),
  
  sidebarLayout(
    sidebarPanel(
      "Factor loadings are drawn from uniform distributions on [a,b] with mean L=(a+b)/2 and range R=(b-a). 
				 The user can set ranges (MR and CR) to zero to specify constant main loadings and constant cross-loadings. 
				 Variables' variances are assumed to be 1. For this reason, 
				 the maximum possible value of a crossloading is determined by the corresponding value of the main loading and the factor correlation. 
				 The upper limit on the choice of the average CL value below is computed from the average ML value specified above. This is done in 
				 an attempt to minimize impossible configurations (with negative error variances).", 
      strong("However, 
				 impossible configurations can still be generated unless both MR and CR are set to zero. Check the printed residual 
				 variances for negative values. Plots will be omitted for all such configurations."),  
      br(),br(), 
      radioButtons("custom", "I would like to examine a ", 
                   c("Two factor model"="twoFactor", "Three factor model"="threeFactor")),
      
      conditionalPanel(
        condition = "input.custom == 'twoFactor'",
        sliderInput("p2", "Total number of variables:", min=4, max=50, step=2, value=6), 
        uiOutput(outputId = "warningLargeP2"),
        #Vika added 2/5/2023: min=-1 i/0 0
        sliderInput("fcor2", "Factor correlation in the true model:", min=-1, max=1, value=.2, step=.05), 
        sliderInput("aveloading2", "Average Main Loading (ML):", min=0, max=1,value=.7),
        uiOutput("sliderange2"), #MR
        uiOutput("slidemax_cross2"), #CL
        uiOutput(outputId = "warningNegativeCL2"),
        uiOutput("sliderange_cross2"), #CR
        uiOutput(outputId = "warningCRgreaterCL2")
      ),
      
      conditionalPanel(
        condition = "input.custom == 'threeFactor'",
        sliderInput("p3", "Total number of variables:", min=6, max=51, step=3, value=9), 
        uiOutput(outputId = "warningLargeP3"),
        # Karyn May 6th: changed min from 0 to -1 
        sliderInput("fcor3", "Factor correlation in the true model:", min=-1, max=1, value=.2, step=.05), 
        sliderInput("aveloading3", "Average Main Loading (ML):", min=0, max=1,value=.7),
        uiOutput("sliderange3"), #MR
        uiOutput("slidemax_cross3"), #CL
        uiOutput(outputId = "warningNegativeCL3"),
        uiOutput("sliderange_cross3"), #CR
        uiOutput(outputId = "warningCRgreaterCL3")
      ),
      
      actionButton("updateButton", "Compute fit index values!")
    ),
    
    mainPanel(  
      p(textOutput("plottext")),
      plotlyOutput("plots"),
      br(),
      uiOutput("tabletext"),
      dataTableOutput("table"),
      dataTableOutput("table2"),
      conditionalPanel(
        condition = "input$custom = input$aveloading3",
        br(),
        uiOutput("matrix")
      )
    )
  )
)


#--------------------------------------------------------------------------------------------------------------------#
# server code
server <- function(input, output, session) {
  
  session$onFlushed(function() {
    session$sendCustomMessage("regrid", FALSE);
  }, FALSE);
  
  #define input sliders for TWO factor model
  #Vika change 2/40/2023: added min below (before was 0), changed max to not exceed 1
  output$slidemax_cross2 <- renderUI({
    sliderInput("aveloading_cross2", "Average Cross Loading (CL):", 
                min = round((-sqrt(1-(input$aveloading2)^2+(input$fcor2*input$aveloading2)^2)-input$fcor2*input$aveloading2),2), 
                max = round(
                  min(1, (sqrt(1-(input$aveloading2)^2+(input$fcor2*input$aveloading2)^2)-input$fcor2*input$aveloading2))
                  ,2), 
                round = -2, step = 0.01, value = .2) 
  })
  
  output$sliderange2 <- renderUI({
    sliderInput("range2", "Main Loadings Range (MR):", min = 0, max = round((min(2*input$aveloading2, 2*(1-input$aveloading2))),2), 
                value = min(.1,input$aveloading2, (1-input$aveloading2)), round = -3, step = 0.01) 
  })
  
  #vika 2/5/2023: In the min(1,2*input$aveloading_cross2,...), replace 1 with 2, 
  #and the second term with its generalized version 
  output$sliderange_cross2 <- renderUI({
    sliderInput("range_cross2", "Cross-Loadings Range (CR):", 
                min = 0, 
                max = round(min(2,
                                2*abs(-sqrt(1-(input$aveloading2)^2+(input$fcor2*input$aveloading2)^2)
                                      -input$fcor2*input$aveloading2-input$aveloading_cross2),  
                                2*(sqrt(1-(input$aveloading2)^2+(input$fcor2*input$aveloading2)^2)
                                   -input$fcor2*input$aveloading2-input$aveloading_cross2)),2), 
                value = min(0,input$input$aveloading_cross2, (1-input$input$aveloading_cross2)), round = -2, step = 0.01) 
  })
  
  
  #define input sliders for THREE factor model
  output$slidemax_cross3 <- renderUI({
    sliderInput("aveloading_cross3", "Average Cross Loading (CL):", 
                min = round((-sqrt(1-(input$aveloading3)^2+(input$fcor3*input$aveloading3)^2)-input$fcor3*input$aveloading3),2), 
                max = round((sqrt(1-(input$aveloading3)^2+(input$fcor3*input$aveloading3)^2)-input$fcor3*input$aveloading3),2), 
                round = -2, step = 0.01, value = .2)
  })
  
  #for randomly generated factor loadings:
  output$sliderange3 <- renderUI({
    sliderInput("range3", "Main Loadings Range (MR):", min = 0, max = round((min(2*input$aveloading3, 2*(1-input$aveloading3))),2), 
                value = min(.1,input$aveloading3, (1-input$aveloading3)), round = -3, step = 0.01) 
  })
  
  #for randomly generated factor cross-loadings:  # Karyn May 6th
  output$sliderange_cross3 <- renderUI({
    sliderInput("range_cross3", "Cross-Loadings Range (CR):", 
                min = 0, 
                max = round(min(2,
                                2*abs(-sqrt(1-(input$aveloading3)^2+(input$fcor3*input$aveloading3)^2)
                                      -input$fcor3*input$aveloading3-input$aveloading_cross3),  
                                2*(sqrt(1-(input$aveloading3)^2+(input$fcor3*input$aveloading3)^2)
                                   -input$fcor3*input$aveloading3-input$aveloading_cross3)),2), 
                value = min(0,input$input$aveloading_cross3, (1-input$input$aveloading_cross3)), round = -2, step = 0.01) 
  })
  
  # show warning message when p is large 
  output$warningLargeP2  <- renderUI({
    if(input$p2 > 20) {
      tagList(
        tags$p("Warning: For a model of this size, you may have to wait a while!", style = "color: red;")
      )
    }
  })
  output$warningLargeP3  <- renderUI({
    if(input$p3 > 18) {
      tagList(
        tags$p("Warning: For a model of this size, you may have to wait a while!", style = "color: red;")
      )
    }
  })
  
  # show warning message for potential negative variances
  # ATTENTION: suppress warning messages in the end 
  output$warningNegativeCL2  <- renderUI({
  if( input$aveloading_cross2 < 0) {
      tagList(
        tags$p("Warning: You’re allowing negative crossloadings! (CL < 0)", style = "color: red;")
      )
    }
  })
  output$warningCRgreaterCL2 <- renderUI({
    if((input$range_cross2 > 2* input$aveloading_cross2) && (input$aveloading_cross2 >= 0 )) {
      tagList(
        tags$p("Warning: You’re allowing negative crossloadings! (CR > 2*CL)", style = "color: red;")
      )
    }
  })
  output$warningNegativeCL3 <- renderUI({
    if(input$aveloading_cross3 < 0) {
      tagList(
        tags$p("Warning: You’re allowing negative crossloadings! (CL < 0)", style = "color: red;")
      )
    }
  })
  output$warningCRgreaterCL3 <- renderUI({
    if((input$range_cross3 > 2* input$aveloading_cross3) && (input$aveloading_cross3 >= 0 )) {
      tagList(
        tags$p("Warning: You’re allowing negative crossloadings! (CR > 2*CL)", style = "color: red;")
      )
    }
  })
  
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
                                     threeFactor = input$aveloading_cross3)         
    
    fcorSwitch <- switch(input$custom,
                         twoFactor = input$fcor2,
                         threeFactor = input$fcor3)
    
    # store numerical value 2 or 3 
    numText <- switch(input$custom,
                      twoFactor = "2-factor model",
                      threeFactor = "3-factor model")
    
    #set.seed(123) #enable to replicate figures in the paper
    
    genLoadingss <- runif(pSwitch, min=aveloadingSwitch-.5*rangeSwitch, max=aveloadingSwitch+.5*rangeSwitch) 
    
    numCrossLoading <- runif(numCrossLoadingSwitch, min=aveloading_crossSwitch -.5*range_crossSwitch, max=aveloading_crossSwitch+.5*range_crossSwitch)
    
    # Allow switching the main function between two-factor model and three-factor model
    mainFunc <- switch(input$custom,
                       twoFactor = main.2f,
                       threeFactor = main.3f)
    
    show_modal_spinner() # show the modal window
    
    resultsAndOrder <- mainFunc(isolate(pSwitch),isolate(fcorSwitch),genLoadingss,numCrossLoading)
    results <- as.data.frame(resultsAndOrder$results)
    orders <- resultsAndOrder$orders
    
    remove_modal_spinner() # remove it when done
    
    #send an alart for negative residual variances
    residuals <- select(tail(results, n=1),matches("x[0-9]{1,2}"))
    if (any(residuals < 0)){
      showNotification(paste("Nagative Variance Warning"), duration = 5, type = "error")
    }

    output$tabletext <- renderUI({
      HTML(paste(strong(em("Main loadings"))," for a ", numText, " are randomly generated." ,
                 strong(em("Crossloadings"))," are randomly generated and added to the true model, one by one."))
    })
    
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
    
    if (input$custom =="twoFactor"){
      
      ### 3f tables ====
      genResiduals_seq <- as.vector(select(residuals,matches("x[0-9]{1,2}same")))%>%
        unname()%>%
        unlist()
      
      genResiduals_alt <- as.vector(select(residuals,matches("x[0-9]{1,2}dif")))%>%
        unname()%>%
        unlist()
      
      # the table output
      table2f1 <- as.data.frame(rbind(
        round(genLoadingss, 3),
        round(numCrossLoading, 3),
        round(genResiduals_seq, 3)
      ))
      
      # the table output
      table2f2 <- as.data.frame(rbind(
        round(genLoadingss, 3),
        round(c(numCrossLoading[c(TRUE, FALSE)],numCrossLoading[c(FALSE, TRUE)]),3),
        round(genResiduals_alt, 3)
      ))
      
      row.names(table2f1) <- c("Main loadings","Cross loadings","Residual variances")
      row.names(table2f2) <- c("Main loadings","Cross loadings","Residual variances")
      
      # Render the first data table 
      output$table <- renderDataTable(
        datatable(
          table2f1,
          colnames = paste0(rep("item", numCrossLoadingSwitch), c(1:numCrossLoadingSwitch)),
          caption = "Parameter values when all cross-loadings are added in sequential order",
          options = list(
            scrollX = T,# to add a horizontal scroller in case of having a wide table
            dom = 't',# to hide the table's filer and search function
            ordering = F # to hide the ordering function
          )
        )
      )
      
      # Render the second data table 
      output$table2 <- renderDataTable(
        datatable(
          table2f2,
          colnames = paste0(rep("item", numCrossLoadingSwitch), c(1:numCrossLoadingSwitch)),
          caption = "Parameter values when all cross-loadings are added in alternating order",
          options = list(
            scrollX = T,# to add a horizontal scroller in case of having a wide table
            dom = 't',# to hide the table's filer and search function
            ordering = F # to hide the ordering function
          )
        )
      )
      
      ### 3f figures =====
      ColorblindnessFriendlyValues2 <- c("Same" = "#648FFF", "Alt" = "#FFB000")
      
      output$plots <- renderPlotly({
        # compute the range for the plots 
        upperbound_rmsea = max(c(results$rmsea_same_f,results$rmsea_dif_f,0.08)) + 0.005
        upperbound_srmr = max(c(results$srmr_same_f,results$srmr_dif_f,0.08)) + 0.005
        lowerbound_cfi = min(c(results$cfi_same_f,results$cfi_dif_f,0.9)) - 0.005
        
        # ensure the x-axis only has whole number 
        if ((max(results$number_crossloadings) < 6)){
          step_tick <- 1
        } else if (max(results$number_crossloadings) %% 4 == 0) {
          step_tick <- floor(max(results$number_crossloadings)/ 4)
        } else {
          step_tick <- floor(max(results$number_crossloadings)/ 3)
        } 
        
        plot1 <- ggplot(data=results,aes(x=number_crossloadings))+ 
          scale_x_continuous(breaks = seq(min(results$number_crossloadings), max(results$number_crossloadings), step_tick)) +
          scale_y_continuous(limits = c(NA,upperbound_rmsea), labels = function(x) sprintf("%.2f", x)) +
          geom_line(aes(y=rmsea_same_f,
                        color="Same"))+ 
          suppressWarnings(geom_point(aes(y=rmsea_same_f,
                                          color="Same",
                                          # shape = "Same",
                                          text = paste0("# of cross Loadings: ", number_crossloadings, "<br>RMSEA: ", sprintf('%.3f', rmsea_same_f)))))+
          geom_line(aes(y=rmsea_dif_f,color ="Alt"))+
          suppressWarnings(geom_point(aes(y=rmsea_dif_f,
                                          color ="Alt",
                                          # shape = "Alt",
                                          text = paste0("# of cross Loadings: ", number_crossloadings, "<br>RMSEA: ", sprintf('%.3f', rmsea_dif_f)))))+
          scale_color_manual(values = ColorblindnessFriendlyValues2, labels = c("Same", "Alt")) +
          # scale_shape_manual(values = c("Same" = 16, "Alt" = 17), labels = c("Same", "Alt")) +
          geom_abline(color="grey",slope=0, intercept=0.08) 
        
        p1 <- ggplotly(plot1,tooltip = c("text"))  %>% style(showlegend = FALSE)
        
        plot2 <- ggplot(data=results, aes(x=number_crossloadings))+ 
          scale_x_continuous(breaks = seq(min(results$number_crossloadings), max(results$number_crossloadings), step_tick)) +
          scale_y_continuous(limits = c(lowerbound_cfi, NA), labels = function(x) sprintf("%.2f", x)) +
          geom_line(aes(y=cfi_same_f,
                        color="Same"))+ 
          suppressWarnings(geom_point(aes(y=cfi_same_f,
                                          color="Same",
                                          # shape = "Same",
                                          text = paste0("# of cross Loadings: ", number_crossloadings,
                                                        "<br>CFI: ", sprintf('%.3f', cfi_same_f)))))+
          geom_line(aes(y=cfi_dif_f, 
                        color="Alt"))+ 
          suppressWarnings(geom_point(aes(y=cfi_dif_f,
                                          color="Alt",
                                          # shape = "Alt",
                                          text = paste0("# of cross Loadings: ", number_crossloadings,
                                                        "<br>CFI: ", sprintf('%.3f', cfi_dif_f)))))+
          scale_color_manual(values = ColorblindnessFriendlyValues2, labels = c("Same", "Alt")) +   
          # scale_shape_manual(values = c("Same" = 16, "Alt" = 17), labels = c("Same", "Alt")) +
          geom_abline(color="grey",slope=0, intercept=0.95)
        
        p2 <- ggplotly(plot2,tooltip = c("text"))  %>% style(showlegend = FALSE)
        
        plot3 <- ggplot(data=results, aes(x=number_crossloadings))+ 
          scale_x_continuous(breaks = seq(min(results$number_crossloadings), max(results$number_crossloadings), step_tick)) +
          scale_y_continuous(limits = c(NA,upperbound_srmr) , labels = function(x) sprintf("%.2f", x)) +
          geom_line(aes(y=srmr_same_f,
                        color="Same"))+ 
          geom_line(aes(y=srmr_dif_f,color="Alt"))+ 
          suppressWarnings(geom_point(aes(y=srmr_same_f,
                                          color="Same",
                                          # shape = "Same",
                                          text = paste0("# of cross Loadings: ", number_crossloadings,
                                                        "<br>SRMR: ", sprintf('%.3f', srmr_same_f)))))+
          suppressWarnings(geom_point(aes(y=srmr_dif_f,
                                          color="Alt",
                                          # shape = "Alt",
                                          text = paste0("# of cross Loadings: ", number_crossloadings,
                                                        "<br>SRMR: ", sprintf('%.3f', srmr_dif_f)))))+
          scale_color_manual(values = ColorblindnessFriendlyValues2, labels = c("Same", "Alt")) +   
          # scale_shape_manual(values = c("Same" = 16, "Alt" = 17), labels = c("Same", "Alt")) +
          geom_abline(color="grey",slope=0, intercept=0.08) + 
          labs(color = "Order") + theme(legend.position = c(0.8, 0.2))
        
        p3 <- ggplotly(plot3,tooltip = c("text"))
        
        subplot(
          p1,p2,p3,
          nrows = 1
        ) %>%
          layout(annotations = layout_3_figs)# %>%
          # layout(legend = list(orientation="h",
          #                      entrywidth=70,
          #                      yanchor="bottom",
          #                      y=4,
          #                      xanchor="right",
          #                      x=1))            
      })
    }
    
    if (input$custom =="threeFactor"){
      ### 3f tables ====
      genResiduals_same1 <- as.vector(select(residuals,matches("x[0-9]{1,2}same1")))%>%
        unname()%>%
        unlist()
      genResiduals_same2 <- as.vector(select(residuals,matches("x[0-9]{1,2}same2")))%>%
        unname()%>%
        unlist()
      genResiduals_dif1 <- as.vector(select(residuals,matches("x[0-9]{1,2}dif1")))%>%
        unname()%>%
        unlist()
      genResiduals_dif2 <- as.vector(select(residuals,matches("x[0-9]{1,2}dif2")))%>%
        unname()%>%
        unlist()
      
      table3f1 <- as.data.frame(rbind(
        round(genLoadingss, 3),
        round(genResiduals_same1, 3),
        round(genResiduals_same2, 3),
        round(genResiduals_dif1, 3),
        round(genResiduals_dif2, 3)
      ))
      row.names(table3f1) <- c("Main loadings","Residual Variance: Same1","Residual Variance: Same2","Residual Variance: Alt1","Residual Variance: Alt2")
      
      table3f2 <- as.data.frame(t(round(numCrossLoading, 3)))
      row.names(table3f2) <- c("Cross loading values")
      
      # number of items (added this to avoid )
      numItem<- numCrossLoadingSwitch/2
      
      output$table <- renderDataTable(datatable(table3f1,
                                                colnames = paste0(rep("item", numItem), c(1: numItem)),
                                                caption = "Parameter values when all cross-loadings are added",
                                                options = list( scrollX = T, # to add a horizontal scroller in case of having a wide table
                                                                dom = 't', # to hide the table's filer and search function
                                                                ordering = F # to hide the ordering function
                                                )))
      output$table2 <- renderDataTable(datatable(table3f2,
                                                 colnames = paste0(rep("CL", numCrossLoadingSwitch), c(1:numCrossLoadingSwitch)),
                                                 caption = "All the cross-loadings",
                                                 options = list( scrollX = T, # to add a horizontal scroller in case of having a wide table
                                                                 dom = 't', # to hide the table's filer and search function
                                                                 ordering = F # to hide the ordering function
                                                 )))
      
      # Print out different values for crossloadings 
      print("Same1 Cross loading values")
      print(sparseMatrix(j = orders[,"col.same1"], i = orders[,"row.same1"], x = orders[,"cld"]))
      print("Same2 Cross loading values")
      print(sparseMatrix(j = orders[,"col.same2"], i = orders[,"row.same2"], x = orders[,"cld"]))
      print("Alt1 Cross loading values")
      print(sparseMatrix(j = orders[,"col.dif1"], i = orders[,"row.dif1"], x = orders[,"cld"]))
      print("Alt2 Cross loading values")
      print(sparseMatrix(j = orders[,"col.dif2"], i = orders[,"row.dif2"], x = orders[,"cld"]))
      
      # render the matrix to explain how cross loadings are added 
      output$matrix <- renderUI({
        withMathJax(
          helpText(
            strong('Cross loadings are added in different orders as shown in the matrices below.'), 
            'Asterisks indicate main loadings, and the numbers indicate blocks of cross-loadings.
        Specific values of cross loadings with respect to the items and factors are printed in R console.',
            '\\begin{pmatrix}     *&1&2\\\\
                                3& *&4\\\\ 
                              5&6&* \\end{pmatrix}',
            HTML("<p>In <strong>“Same 1”</strong> order, cross-loadings are added column-wise, until each factor reflects all other indicators; blocks are completely filled in in the following order: 3, 5, 1, 6, 2, 4.</p>"),
            HTML("<p>In <strong>“Same 2”</strong> order, cross-loadings are added until a factor covers indicators of just one other factor before moving on to the next factor; blocks are completely filled in in the following order: 3, 6, 2, 5, 1, 4.</p>"),
            HTML("<p>In <strong>“Alternating 1”</strong> order, loadings are added row-wise, first alternating between blocks 1 and 2 until they are saturated, then between blocks 3 and 4 until they are saturated, then between blocks 5 and 6.</p>"),
            HTML("<p>In <strong>“Alternating 2”</strong> order, one cross-loading is added to each block, and then the cycle is repeated, so that no block gets saturated with cross-loadings until the very end of the series.</p>")
          )
        )
      })
      
      ### 3f figures ====
      ColorblindnessFriendlyValues4 <- c("Same1" = "#648FFF", "Same2" = "#D81B60", "Alt1" = "#FFB000", "Alt2" = "#004D40")
      
      output$plots <- renderPlotly({
        
        upperbound_rmsea = max(c(results$rmsea_same1_f,results$rmsea_dif1_f,results$rmsea_same2_f,results$rmsea_dif2_f,0.08)) + 0.005
        upperbound_srmr = max(c(results$srmr_same1_f,results$srmr_dif1_f,results$srmr_same2_f,results$srmr_dif2_f,0.08)) + 0.005
        lowerbound_cfi = min(c(results$cfi_same1_f,results$cfi_same2_f,results$cfi_dif1_f,results$cfi_dif2_f,0.9))-0.005
        
        plot1 <- ggplot(data=results,aes(x=number_crossloadings)) +
          scale_y_continuous(limits = c(NA,upperbound_rmsea), labels = function(x) sprintf("%.2f", x)) +
          # scale_x_continuous(breaks = seq(min(results$number_crossloadings), max(results$number_crossloadings), step_tick)) +
          geom_line(mapping = aes( y=rmsea_same1_f, color="Same1"), show.legend = FALSE)+
          suppressWarnings(geom_point(aes(y=rmsea_same1_f,
                                          color="Same1",
                                          # shape="Same1",
                                          text = paste0("# of cross Loadings: ", number_crossloadings, "<br>RMSEA: ", sprintf('%.3f', rmsea_same1_f)))))+
          
          geom_line(mapping = aes( y=rmsea_same2_f, color="Same2"), show.legend = FALSE)+
          suppressWarnings(geom_point(aes(y=rmsea_same2_f,
                                          color="Same2",
                                          # shape="Same2",
                                          text = paste0("# of cross Loadings: ", number_crossloadings, "<br>RMSEA: ", sprintf('%.3f', rmsea_same2_f)))))+
          
          geom_line(mapping = aes( y=rmsea_dif1_f, color="Alt1"), show.legend = FALSE)+
          suppressWarnings(geom_point(aes(y=rmsea_dif1_f,
                                          color="Alt1",
                                          # shape="Alt1",
                                          text = paste0("# of cross Loadings: ", number_crossloadings, "<br>RMSEA: ", sprintf('%.3f', rmsea_dif1_f)))))+
          
          geom_line(mapping = aes( y=rmsea_dif2_f, color="Alt2"), show.legend = FALSE)+
          suppressWarnings(geom_point(aes(y=rmsea_dif2_f,
                                          color="Alt2",
                                          # shape="Alt2",
                                          text = paste0("# of cross Loadings: ", number_crossloadings, "<br>RMSEA: ", sprintf('%.3f', rmsea_dif2_f)))))+
          geom_abline(color="grey",slope=0, intercept=0.08) +
          scale_color_manual(values = ColorblindnessFriendlyValues4, labels = c("Same1", "Same2", "Alt1", "Alt2")) +
          # labs(color = "Order") +
          ylab("RMSEA")
        p1 <- ggplotly(plot1,tooltip = c("text"))  %>% style(showlegend = FALSE)
        
        plot2 <- ggplot(data=results,aes(x=number_crossloadings)) +
          scale_y_continuous(limits = c(lowerbound_cfi, NA), labels = function(x) sprintf("%.2f", x)) +
          # scale_x_continuous(breaks = seq(min(results$number_crossloadings), max(results$number_crossloadings), step_tick)) +
          geom_line(mapping = aes( y=cfi_same1_f, color="Same1"), show.legend = FALSE)+
          suppressWarnings(geom_point(aes(y=cfi_same1_f,
                                          color="Same1",
                                          # shape="Same1",
                                          text = paste0("# of cross Loadings: ", number_crossloadings, "<br>CFI: ", sprintf('%.3f', cfi_same1_f)))))+
          
          geom_line(mapping = aes( y=cfi_same2_f, color="Same2"), show.legend = FALSE)+
          suppressWarnings(geom_point(aes(y=cfi_same2_f,
                                          color="Same2",
                                          # shape="Same2",
                                          text = paste0("# of cross Loadings: ", number_crossloadings, "<br>CFI: ", sprintf('%.3f', cfi_same2_f)))))+
          
          geom_line(mapping = aes( y=cfi_dif1_f, color="Alt1"), show.legend = FALSE)+
          suppressWarnings(geom_point(aes(y=cfi_dif1_f,
                                          color="Alt1",
                                          # shape="Alt1",
                                          text = paste0("# of cross Loadings: ", number_crossloadings, "<br>CFI: ", sprintf('%.3f', cfi_dif1_f)))))+
          
          geom_line(mapping = aes( y=cfi_dif2_f, color="Alt2"), show.legend = FALSE)+
          suppressWarnings(geom_point(aes(y=cfi_dif2_f,
                                          color="Alt2",
                                          # shape="Alt2",
                                          text = paste0("# of cross Loadings: ", number_crossloadings, "<br>CFI: ", sprintf('%.3f', cfi_dif2_f)))))+
          scale_color_manual(values = ColorblindnessFriendlyValues4, labels = c("Same1", "Same2", "Alt1", "Alt2")) +
          geom_abline(color="grey",slope=0, intercept=0.95) +
          # labs(color = "Order") +
          ylab("CFI")
        p2 <- ggplotly(plot2,tooltip = c("text"))  %>% style(showlegend = FALSE)
        
        plot3 <- ggplot(data=results,aes(x=number_crossloadings)) +
          scale_y_continuous(limits = c(NA,upperbound_srmr) , labels = function(x) sprintf("%.2f", x)) +
          # scale_x_continuous(breaks = seq(min(results$number_crossloadings), max(results$number_crossloadings), step_tick)) +
          geom_line(mapping = aes( y=srmr_same1_f, color="Same1"), show.legend = FALSE)+
          suppressWarnings(geom_point(aes(y=srmr_same1_f,
                                          color="Same1",
                                          # shape="Same1",
                                          text = paste0("# of cross Loadings: ", number_crossloadings, "<br>SRMR: ", sprintf('%.3f', srmr_same1_f)))))+
          
          geom_line(mapping = aes( y=srmr_same2_f, color="Same2"), show.legend = FALSE)+
          suppressWarnings(geom_point(aes(y=srmr_same2_f,
                                          color="Same2",
                                          # shape="Same2",
                                          text = paste0("# of cross Loadings: ", number_crossloadings, "<br>SRMR: ", sprintf('%.3f', srmr_same2_f)))))+
          
          geom_line(mapping = aes( y=srmr_dif1_f, color="Alt1"), show.legend = FALSE)+
          suppressWarnings(geom_point(aes(y=srmr_dif1_f,
                                          color="Alt1",
                                          # shape="Alt1",
                                          text = paste0("# of cross Loadings: ", number_crossloadings, "<br>SRMR: ", sprintf('%.3f', srmr_dif1_f)))))+
          
          geom_line(mapping = aes( y=srmr_dif2_f, color="Alt2"), show.legend = FALSE)+
          suppressWarnings(geom_point(aes(y=srmr_dif2_f,
                                          color="Alt2",
                                          # shape="Alt2",
                                          text = paste0("# of cross Loadings: ", number_crossloadings, "<br>SRMR: ", sprintf('%.3f', srmr_dif2_f)))))+
          geom_abline(color="grey",slope=0, intercept=0.08) +
          scale_color_manual(values = ColorblindnessFriendlyValues4, labels = c("Same1", "Same2", "Alt1", "Alt2")) +
          labs(color = "Order") +
          ylab("SRMR")
        
        p3 <- ggplotly(plot3,tooltip = c("text"))  #%>%
          # layout(legend = list(orientation="h",
          #                      entrywidth=70,
          #                      yanchor="bottom",
          #                      y=4,
          #                      xanchor="right",
          #                      x=1)) 
        
        subplot(
          p1,p2,p3,
          nrows = 1
        ) %>%
          layout(annotations = layout_3_figs)
      })
    }
    
    # Managing the layout of three figures 
    layout_3_figs <- list(
      list(
        x= -0.005,
        y =1.05,
        xanchor="left",
        text = "RMSEA",
        xref = "paper",
        yref = "paper",
        showarrow = FALSE
      ),
      list(
        x = 0.35,
        y = 1.05,
        text = "CFI",
        xref = "paper",
        yref = "paper",
        xanchor="left",
        showarrow = FALSE
      ),
      list(
        x = 0.68,
        y = 1.05,
        text = "SRMR",
        xref = "paper",
        yref = "paper",
        xanchor="left",
        showarrow = FALSE
      ) )
  }) 
}


shinyApp(ui = ui, server = server)
