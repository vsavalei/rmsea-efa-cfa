library(shiny)
library(plotly) 

source("functions_app.R")

#--------------------------------------------------------------------------------------------------------------------#
# UI code

# for adjusting the overlapping anchors on the slider 
# js <- paste("function MRdoesOverlap() {",
#             "   var $lastLabel = $('#sliderange .irs-grid-text:last');", # MR 
#             "   var $prevLastLabel = $lastLabel.prevAll('.irs-grid-text').first();",
#             "   return $lastLabel.offset().left < $prevLastLabel.offset().left + $prevLastLabel.width();",
#             "}\n",
#             "Shiny.addCustomMessageHandler('regrid', function(force) {",
#             "   if (MRdoesOverlap() | force) {",
#             "      console.log('Overlap detected - adjusting tick number');",
#             "      var $sld = $('#range').data('ionRangeSlider');", #range is the input name 
#             "      var ticks_n = $sld.options.grid_num;",
#             "      $sld.update({grid_num: Math.round(ticks_n)});",
#             "   }",
#             "});",
#             "function CLdoesOverlap() {",
#             "   var $lastLabel = $('#slidemax_cross .irs-grid-text:last');", # CL
#             "   var $prevLastLabel = $lastLabel.prevAll('.irs-grid-text').first();",
#             "   return $lastLabel.offset().left < $prevLastLabel.offset().left + $prevLastLabel.width();",
#             "}\n",
#             "Shiny.addCustomMessageHandler('regrid', function(force) {",
#             "   if (CLdoesOverlap() | force) {",
#             "      console.log('Overlap detected - adjusting tick number');",
#             "      var $sld = $('#aveloading_cross').data('ionRangeSlider');", #range is the input name
#             "      var ticks_n = $sld.options.grid_num;",
#             "      $sld.update({grid_num: Math.round(ticks_n)});",
#             "   }",
#             "});",
#             sep = "\n")

js1 <- paste("function MRdoesOverlap() {",
            "   var $lastLabel = $('#sliderange .irs-grid-text:last');", # MR 
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
			sliderInput("p", "Number of variables:", min=4, max=50, step=2, value=8), 
			
			sliderInput("fcor", "Factor correlation in the true model:", min=0, max=1, value=.2), 
			
			radioButtons("custom", "How do you want to define the factor loadings in the true model?", 
				c("Generate randomly"="rand","Specify manually"="spec")),
			
			conditionalPanel(
				condition = "input.custom == 'rand'",
				"Factor loadings are drawn from uniform distributions on [a,b] with mean L=(a+b)/2 and range R=(b-a). 
				 The user can set ranges (MR and CR) to zero to specify constant main loadings and constant cross-loadings. 
				 Variables' variances are assumed to be 1. For this reason, 
				 the maximum possible value of a crossloading is determined by the corresponding value of the main loading and the factor correlation. 
				 The upper limit on the choice of the average CL value below is computed from the average ML value specified above. This is done in 
				 an attempt to minimize impossible configurations (with negative error variances). However, 
				 impossible configurations can still be generated unless both MR and CR are set to zero. Check the printed residual 
				 variances for negative values. Plots will be omitted for all such configurations.",  
				br(),br(),
				sliderInput("aveloading", "Average Main Loading (ML)", min=0, max=1,value=.7),
				uiOutput("sliderange"), #MR
				uiOutput("slidemax_cross"), #CL
				uiOutput("sliderange_cross"), #CR
				actionButton("updateButton", "Compute fit index values!")
				),
			
		  conditionalPanel(
				condition = "input.custom == 'spec'",
				strong("Specify the main factor loading for each variable:"),br(),
				uiOutput("slideload"),actionButton("updateButtonSpec", "Plot!!")
				)
			),
		
		mainPanel(  
			conditionalPanel(
				condition = "input.custom == 'rand'", 
				h4(textOutput("printload")),h4(textOutput("printload_cross")),
				h4(textOutput("printres")),h4(textOutput("printres2")),h4(textOutput("plottext")), plotlyOutput("plot1"),
				plotlyOutput("plot2"),plotlyOutput("plot3")
				),
			
		conditionalPanel(
				condition = "input.custom == 'spec'",
				h4(textOutput("printloadspec")),
				br(),
				h4(textOutput("plottextspec")),
				plotlyOutput("plotspec"),
				h4(textOutput("plottextspecadd"))
				)
			
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
	
	#define input sliders   
	
  output$slidemax_cross <- renderUI({
    sliderInput("aveloading_cross", "Average Cross Loading (CL)", min = 0, 
                max = round((sqrt(1-(input$aveloading)^2+(input$fcor*input$aveloading)^2)-input$fcor*input$aveloading),2), 
                round = -2, step = 0.01, value = .2) 
  })
  
  #for randomly generated factor loadings:
	output$sliderange <- renderUI({
		sliderInput("range", "Main Loadings Range (MR)", min = 0, max = round((min(2*input$aveloading, 2*(1-input$aveloading))),2), 
		            value = min(.1,input$aveloading, (1-input$aveloading)), round = -3, step = 0.01) 
	})
	
	#for randomly generated factor cross-loadings:
	output$sliderange_cross <- renderUI({
	  sliderInput("range_cross", "Cross-Loadings Range (CR)", min = 0, max = round(min(1,2*input$aveloading_cross, 
	  2*(sqrt(1-(input$aveloading)^2+(input$fcor*input$aveloading)^2)-input$fcor*input$aveloading-input$aveloading_cross)),2), 
	  
	  value = min(0,input$slidermax_cross, (1-input$slidermax_cross)), round = -2, step = 0.01) 
	})
	
	#for user defined factor loadings: 
	numP <- reactive({input$p})  # get number of loadings from ui slider 
	
	output$slideload <- renderUI({   # create custom loading sliders
		lapply(1:(numP()), function(i) {
			sliderInput(inputId = paste0("load", i), # this creates input$load1, input$load2, etc.
				label = paste0("Loading ", i),
				min = 0, max = 1, value = .5)
		})
	})
	
	# button pressed to update randomly generated loadings
	# code needed
	
	observeEvent(input$updateButton,{      

	  loadtxt <- "The true model is a 2-factor model with the randomly generated main loadings of " 
	  genLoadingss <- runif(input$p, min=input$aveloading-.5*input$range, max=input$aveloading+.5*input$range)
	 
	  loadtxt_cross <- "The randomly generated values of crossloadings to be added to the true model, one by one, are " 
	  genLoadingss_cross <- runif(input$p, min=input$aveloading_cross-.5*input$range_cross, max=input$aveloading_cross+.5*input$range_cross)
	  
	  residualstxt <- "When these loadings are added first to one factor, then to the next, the residual variances in the final 
	  model with all the cross-loadings added are 
	  (i.e., 1 minus squared main loading, minus squared crossloading, and minus twice 
	  the factor correlation times the main loading times the crossloading)  " 
	  genResiduals <- 1-genLoadingss^2-genLoadingss_cross^2-2*genLoadingss*genLoadingss_cross*input$fcor 
	  
	  residualstxt2 <- "When these loadings are added to alternating factors, the residual variances in the final 
	  model with all the cross-loadings added are 
	  (i.e., 1 minus squared main loading, minus squared crossloading, and minus twice 
	  the factor correlation times the main loading times the crossloading)  " 
	  genLoadingss_cross_reordered<- c(genLoadingss_cross[c(TRUE, FALSE)], genLoadingss_cross[c(TRUE, FALSE)])
	  genResiduals2 <- 1-genLoadingss^2-genLoadingss_cross_reordered^2-2*genLoadingss*genLoadingss_cross_reordered*input$fcor 
	  
	  
	 #define output text (recap)
	  output$printload <- renderLoadingRecap(loadtxt, genLoadingss, isolate(input$p))
	  output$printload_cross <- renderLoadingRecap(loadtxt_cross, genLoadingss_cross, isolate(input$p))
	  output$printres <- renderLoadingRecap(residualstxt, genResiduals, isolate(input$p))
	  output$printres2 <- renderLoadingRecap(residualstxt2, genResiduals2, isolate(input$p))
	  
	  #replace rmseas with results
	  results <- as.data.frame(main.2f(isolate(input$p),isolate(input$fcor),genLoadingss,genLoadingss_cross))
	  
	  #define text
	  output$plottext <- renderText({ 
	    paste("In the plots below, the number of crossloadings in the true model is on the x-axis; this number 
	          varies from 0 to ", isolate(input$p), " (the number of variables). The cross-loadings in the true 
	          model are being added either a) to the first factor first and then to the second factor
	          or b) to the alternating factors. Their exact values are given by the list above. (If some 
	          values are missing from the plots, one of the residual variances is negative or the model 
	          failed to converge). The fitted model is a 2-factor model with no cross-loadings. 
	          You can hover over the curve to get specific fit index values.", sep="") 
	  })
	  
	  output$plot1 <- renderPlotly({
	    pt <- ggplot(data=results,aes(x=number_crossloadings))+ 
	      geom_line(data=results,aes(y=rmsea_same_f,
	                                 color="To 1st factor, then 2d",
	                                 text = paste0("# of cross Loadings: ", number_crossloadings,
	                                              "<br>RMSEA: ", sprintf('%.3f', rmsea_same_f))),group=1)+ 
	      geom_line(data=results,aes(y=rmsea_dif_f,color="To alternating factors",
	                                 text = paste0("# of cross Loadings: ", number_crossloadings,
	                                                                                                      "<br>RMSEA: ", sprintf('%.3f', rmsea_dif_f))),group=1)+ 
	      geom_abline(color="grey",slope=0, intercept=0.08) + labs(color = "How crossloadings are added") +
	      xlab("Number of crossloadings in the true model")+
	      ylab("RMSEA for the model with no crossloadings")
	    
	    pt <- ggplotly(pt, tooltip = "text")
	  })
	  
	  output$plot2 <- renderPlotly({
	    pt <- ggplot(data=results, aes(x=number_crossloadings))+ 
	      geom_line(data=results,
	                aes(y=cfi_same_f,
	                    color="To 1st factor, then 2d",
	                    text = paste0("# of cross Loadings: ", number_crossloadings,
	                                  "<br>CFI: ", sprintf('%.3f', cfi_same_f))),group=1)+ 
	      geom_line(data=results,aes(y=cfi_dif_f,
	                                 color="To alternating factors",
	                                 text = paste0("# of cross Loadings: ", number_crossloadings,
	                                               "<br>CFI: ", sprintf('%.3f', cfi_dif_f))),group=1)+ 
	      geom_abline(color="grey",slope=0, intercept=0.90) + labs(color = "How crossloadings are added") +
	      xlab("Number of crossloadings in the true model")+
	      ylab("CFI for the model with no crossloadings")
	    pt <- ggplotly(pt,tooltip = c("text"))
	  })
	  
	  output$plot3 <- renderPlotly({
	    pt <- ggplot(data=results, aes(x=number_crossloadings))+ 
	      geom_line(data=results,aes(y=srmr_same_f,color="To 1st factor, then 2d",
	                                 text = paste0("# of cross Loadings: ", number_crossloadings,
	                                               "<br>SRMR: ", sprintf('%.3f', srmr_same_f))),group=1)+ 
	      geom_line(data=results,aes(y=srmr_dif_f,color="To alternating factors",
	                                 text = paste0("# of cross Loadings: ", number_crossloadings,
	                                               "<br>SRMR: ", sprintf('%.3f', srmr_dif_f))),group=1)+ 
	      geom_abline(color="grey",slope=0, intercept=0.08) + labs(color = "How crossloadings are added") +
	      xlab("Number of crossloadings in the true model")+
	      ylab("SRMR for the model with no crossloadings")
	    pt <- ggplotly(pt,tooltip = c("text"))
	  })
	  
	  
	  
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
# 				tnoend=paste(as.character(round(specLoadingss[1:(isolate(input$p)-1)],3)), collapse=", ")
# 				tend=paste(as.character(round(specLoadingss[isolate(input$p)],3)), collapse=", ")
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
# 		if(input$p<=pm){
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
	# 	gs <- nrow(relspecx)/input$p #gridsize, as set in alpha.omega.pop.varyl
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
	# 		if (isolate(input$p)<=pm){ # show linetype legend
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