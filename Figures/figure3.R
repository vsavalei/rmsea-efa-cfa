library(ggplot2)
library(dplyr)
library("patchwork")
library(cowplot)
source("functions_app.R")

ColorblindnessFriendlyValues4 <- c("#648FFF", "#D81B60", "#FFB000", "#004D40")

ShapeValue2 <- c(15,19,17,3)

# Figures - Left 
set.seed(123)

p = 15
fcor = 0.2
ML = 0.7
MR = 0 
CL = 0.2
CR = 0 

genLoadingss <- runif(p, min=ML-.5*MR, max=ML+.5*MR) 

numCrossLoading <- runif(p*2, min=0.2 -.5*CR, max=0.2+.5*CR)

resultsAndOrder <- main.3f(p,fcor,genLoadingss,numCrossLoading)

results <- as.data.frame(resultsAndOrder$results)


upperbound_rmsea = max(c(results$rmsea_same1_f,results$rmsea_dif1_f,results$rmsea_same2_f,results$rmsea_dif2_f,0.08)) + 0.005
upperbound_srmr = max(c(results$srmr_same1_f,results$srmr_dif1_f,results$srmr_same2_f,results$srmr_dif2_f,0.08)) + 0.005
lowerbound_cfi = min(c(results$cfi_same1_f,results$cfi_same2_f,results$cfi_dif1_f,results$cfi_dif2_f,0.9))-0.005

# ensure the x-axis only has whole number 
if ((max(results$number_crossloadings) < 9)){
  step_tick <- 2
} else {
  step_tick <- floor(max(results$number_crossloadings)/ 4)
} 

plotL1 <- ggplot(data=results,aes(x=number_crossloadings))+ 
  scale_x_continuous(breaks = seq(min(results$number_crossloadings), max(results$number_crossloadings), step_tick)) +
  scale_y_continuous(limits = c(NA,upperbound_rmsea), labels = function(x) sprintf("%.2f", x)) +
  geom_line(mapping = aes( y=rmsea_same1_f, color="Same1"), show.legend = FALSE)+
  suppressWarnings(geom_point(aes(y=rmsea_same1_f,
                                  color="Same1", shape = "Same1")))+
  
  geom_line(mapping = aes( y=rmsea_same2_f, color="Same2"), show.legend = FALSE)+
  suppressWarnings(geom_point(aes(y=rmsea_same2_f,
                                  color="Same2", shape = "Same2")))+
  
  geom_line(mapping = aes( y=rmsea_dif1_f, color="Alt1"), show.legend = FALSE)+
  suppressWarnings(geom_point(aes(y=rmsea_dif1_f,
                                  color="Alt1", shape = "Alt1")))+
  
  geom_line(mapping = aes( y=rmsea_dif2_f, color="Alt2"), show.legend = FALSE)+
  suppressWarnings(geom_point(aes(y=rmsea_dif2_f,
                                  color="Alt2", shape = "Alt2")))+
  scale_color_manual(values = ColorblindnessFriendlyValues4, labels = c("Same1", "Same2", "Alt1", "Alt2")) +
  scale_shape_manual(values = ShapeValue2, labels = c("Same1", "Same2", "Alt1", "Alt2")) +
  geom_abline(color="grey",slope=0, intercept=0.08) +
  labs(color = "Order", shape = "Order")  + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(size = 10, hjust = 0.5))+
  ggtitle("RMSEA")


plotL2 <- ggplot(data=results,aes(x=number_crossloadings)) +
  scale_y_continuous(limits = c(lowerbound_cfi, NA), labels = function(x) sprintf("%.2f", x)) +
  scale_x_continuous(breaks = seq(min(results$number_crossloadings), max(results$number_crossloadings), step_tick)) +
  geom_line(mapping = aes( y=cfi_same1_f, color="Same1"), show.legend = FALSE)+
  suppressWarnings(geom_point(aes(y=cfi_same1_f,
                                  color="Same1",
                                  shape = "Same1")))+
  
  geom_line(mapping = aes( y=cfi_same2_f, color="Same2"), show.legend = FALSE)+
  suppressWarnings(geom_point(aes(y=cfi_same2_f,
                                  color="Same2",
                                  shape = "Same2")))+
  
  geom_line(mapping = aes( y=cfi_dif1_f, color="Alt1"), show.legend = FALSE)+
  suppressWarnings(geom_point(aes(y=cfi_dif1_f,
                                  color="Alt1",
                                  shape = "Alt1")))+
  
  geom_line(mapping = aes( y=cfi_dif2_f, color="Alt2"), show.legend = FALSE)+
  suppressWarnings(geom_point(aes(y=cfi_dif2_f,
                                  color="Alt2",
                                  shape = "Alt2")))+
  scale_color_manual(values = ColorblindnessFriendlyValues4, labels = c("Same1", "Same2", "Alt1", "Alt2")) +
  scale_shape_manual(values = ShapeValue2, labels = c("Same1", "Same2", "Alt1", "Alt2")) +
  geom_abline(color="grey",slope=0, intercept=0.08) +
  labs(color = "Order", shape = "Order")  + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(size = 10, hjust = 0.5))+
  ggtitle("CFI")


plotL3 <- ggplot(data=results,aes(x=number_crossloadings))+ 
  scale_x_continuous(breaks = seq(min(results$number_crossloadings), max(results$number_crossloadings), step_tick)) +
  scale_y_continuous(limits = c(NA,upperbound_srmr), labels = function(x) sprintf("%.2f", x)) +
  geom_line(mapping = aes( y=srmr_same1_f, color="Same1"), show.legend = FALSE)+
  suppressWarnings(geom_point(aes(y=srmr_same1_f,
                                  color="Same1", shape = "Same1")))+
  
  geom_line(mapping = aes( y=srmr_same2_f, color="Same2"), show.legend = FALSE)+
  suppressWarnings(geom_point(aes(y=srmr_same2_f,
                                  color="Same2", shape = "Same2")))+
  
  geom_line(mapping = aes( y=srmr_dif1_f, color="Alt1"), show.legend = FALSE)+
  suppressWarnings(geom_point(aes(y=srmr_dif1_f,
                                  color="Alt1", shape = "Alt1")))+
  
  geom_line(mapping = aes( y=srmr_dif2_f, color="Alt2"), show.legend = FALSE)+
  suppressWarnings(geom_point(aes(y=srmr_dif2_f,
                                  color="Alt2", shape = "Alt2")))+
  scale_color_manual(values = ColorblindnessFriendlyValues4, labels = c("Same1", "Same2", "Alt1", "Alt2")) +
  scale_shape_manual(values = ShapeValue2, labels = c("Same1", "Same2", "Alt1", "Alt2")) +
  geom_abline(color="grey",slope=0, intercept=0.08) +
  labs(color = "Order", shape = "Order")  + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(size = 10, hjust = 0.5))+
  ggtitle("SRMR")

# Figures - Right 
set.seed(123)

p = 24
fcor = 0.2
ML = 0.7
MR = 0 
CL = 0.2
CR = 0 

genLoadingss <- runif(p, min=ML-.5*MR, max=ML+.5*MR) 

numCrossLoading <- runif(p*2, min=0.2 -.5*CR, max=0.2+.5*CR)

resultsAndOrder <- main.3f(p,fcor,genLoadingss,numCrossLoading)

results <- as.data.frame(resultsAndOrder$results)


upperbound_rmsea = max(c(results$rmsea_same1_f,results$rmsea_dif1_f,results$rmsea_same2_f,results$rmsea_dif2_f,0.08)) + 0.005
upperbound_srmr = max(c(results$srmr_same1_f,results$srmr_dif1_f,results$srmr_same2_f,results$srmr_dif2_f,0.08)) + 0.005
lowerbound_cfi = min(c(results$cfi_same1_f,results$cfi_same2_f,results$cfi_dif1_f,results$cfi_dif2_f,0.9))-0.005

# ensure the x-axis only has whole number 
if ((max(results$number_crossloadings) < 9)){
  step_tick <- 2
} else {
  step_tick <- floor(max(results$number_crossloadings)/ 4)
} 

plotR1 <- ggplot(data=results,aes(x=number_crossloadings))+ 
  scale_x_continuous(breaks = seq(min(results$number_crossloadings), max(results$number_crossloadings), step_tick)) +
  scale_y_continuous(limits = c(NA,upperbound_rmsea), labels = function(x) sprintf("%.2f", x)) +
  geom_line(mapping = aes( y=rmsea_same1_f, color="Same1"), show.legend = FALSE)+
  suppressWarnings(geom_point(aes(y=rmsea_same1_f,
                                  color="Same1", shape = "Same1")))+
  
  geom_line(mapping = aes( y=rmsea_same2_f, color="Same2"), show.legend = FALSE)+
  suppressWarnings(geom_point(aes(y=rmsea_same2_f,
                                  color="Same2", shape = "Same2")))+
  
  geom_line(mapping = aes( y=rmsea_dif1_f, color="Alt1"), show.legend = FALSE)+
  suppressWarnings(geom_point(aes(y=rmsea_dif1_f,
                                  color="Alt1", shape = "Alt1")))+
  
  geom_line(mapping = aes( y=rmsea_dif2_f, color="Alt2"), show.legend = FALSE)+
  suppressWarnings(geom_point(aes(y=rmsea_dif2_f,
                                  color="Alt2", shape = "Alt2")))+
  scale_color_manual(values = ColorblindnessFriendlyValues4, labels = c("Same1", "Same2", "Alt1", "Alt2")) +
  scale_shape_manual(values = ShapeValue2, labels = c("Same1", "Same2", "Alt1", "Alt2")) +
  geom_abline(color="grey",slope=0, intercept=0.08) +
  labs(color = "Order", shape = "Order")  + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(size = 10, hjust = 0.5))+
  ggtitle("RMSEA")


plotR2 <- ggplot(data=results,aes(x=number_crossloadings)) +
  scale_y_continuous(limits = c(lowerbound_cfi, NA), labels = function(x) sprintf("%.2f", x)) +
  scale_x_continuous(breaks = seq(min(results$number_crossloadings), max(results$number_crossloadings), step_tick)) +
  geom_line(mapping = aes( y=cfi_same1_f, color="Same1"), show.legend = FALSE)+
  suppressWarnings(geom_point(aes(y=cfi_same1_f,
                                  color="Same1",
                                  shape = "Same1")))+
  
  geom_line(mapping = aes( y=cfi_same2_f, color="Same2"), show.legend = FALSE)+
  suppressWarnings(geom_point(aes(y=cfi_same2_f,
                                  color="Same2",
                                  shape = "Same2")))+
  
  geom_line(mapping = aes( y=cfi_dif1_f, color="Alt1"), show.legend = FALSE)+
  suppressWarnings(geom_point(aes(y=cfi_dif1_f,
                                  color="Alt1",
                                  shape = "Alt1")))+
  
  geom_line(mapping = aes( y=cfi_dif2_f, color="Alt2"), show.legend = FALSE)+
  suppressWarnings(geom_point(aes(y=cfi_dif2_f,
                                  color="Alt2",
                                  shape = "Alt2")))+
  scale_color_manual(values = ColorblindnessFriendlyValues4, labels = c("Same1", "Same2", "Alt1", "Alt2")) +
  scale_shape_manual(values = ShapeValue2, labels = c("Same1", "Same2", "Alt1", "Alt2")) +
  geom_abline(color="grey",slope=0, intercept=0.08) +
  labs(color = "Order", shape = "Order")  + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(size = 10, hjust = 0.5))+
  ggtitle("CFI")


plotR3 <- ggplot(data=results,aes(x=number_crossloadings))+ 
  scale_x_continuous(breaks = seq(min(results$number_crossloadings), max(results$number_crossloadings), step_tick)) +
  scale_y_continuous(limits = c(NA,upperbound_srmr), labels = function(x) sprintf("%.2f", x)) +
  geom_line(mapping = aes( y=srmr_same1_f, color="Same1"), show.legend = FALSE)+
  suppressWarnings(geom_point(aes(y=srmr_same1_f,
                                  color="Same1", shape = "Same1")))+
  
  geom_line(mapping = aes( y=srmr_same2_f, color="Same2"), show.legend = FALSE)+
  suppressWarnings(geom_point(aes(y=srmr_same2_f,
                                  color="Same2", shape = "Same2")))+
  
  geom_line(mapping = aes( y=srmr_dif1_f, color="Alt1"), show.legend = FALSE)+
  suppressWarnings(geom_point(aes(y=srmr_dif1_f,
                                  color="Alt1", shape = "Alt1")))+
  
  geom_line(mapping = aes( y=srmr_dif2_f, color="Alt2"), show.legend = FALSE)+
  suppressWarnings(geom_point(aes(y=srmr_dif2_f,
                                  color="Alt2", shape = "Alt2")))+
  scale_color_manual(values = ColorblindnessFriendlyValues4, labels = c("Same1", "Same2", "Alt1", "Alt2")) +
  scale_shape_manual(values = ShapeValue2, labels = c("Same1", "Same2", "Alt1", "Alt2")) +
  geom_abline(color="grey",slope=0, intercept=0.08) +
  labs(color = "Order", shape = "Order")  + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(size = 10, hjust = 0.5))+
  ggtitle("SRMR")


spacer <- plot_spacer()


plotL1 + plotL2 + plotL3 + spacer + plotR1 + plotR2 +plotR3+ plot_layout(ncol = 7, guides = "collect", widths  = c(1, 1,1, 0.1, 1, 1, 1))
