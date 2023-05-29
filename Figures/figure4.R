library(ggplot2)
library(dplyr)
library("patchwork")
library(cowplot)
source("/Users/mhhuang/Library/CloudStorage/OneDrive-UBC/rmsea-efa-cfa/functions_app.R")

ColorblindnessFriendlyValues2 <- c("Same" = "#648FFF", "Alt" = "#FFB000")
ShapeValue2 <- c("Same" = 15, "Alt" = 16)

# Figure 1 - Left 
set.seed(123)

p = 10
fcor = 0.2
ML = 0.7
MR = 0.3
CL = 0.2
CR = 0.2

genLoadingss <- runif(p, min=ML-.5*MR, max=ML+.5*MR) 

numCrossLoading <- runif(p*2, min=0.2 -.5*CR, max=0.2+.5*CR)

resultsAndOrder <- main.2f(p,fcor,genLoadingss,numCrossLoading)

results <- as.data.frame(resultsAndOrder$results)


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

plotL1 <- ggplot(data=results,aes(x=number_crossloadings))+ 
  scale_x_continuous(breaks = seq(min(results$number_crossloadings), max(results$number_crossloadings), step_tick)) +
  scale_y_continuous(limits = c(NA,upperbound_rmsea), labels = function(x) sprintf("%.2f", x)) +
  geom_line(aes(y=rmsea_same_f,
                color="Same"))+ 
  suppressWarnings(geom_point(aes(y=rmsea_same_f,
                                  color="Same",
                                  shape = "Same")))+
  geom_line(aes(y=rmsea_dif_f,color ="Alt"))+
  suppressWarnings(geom_point(aes(y=rmsea_dif_f,
                                  color ="Alt",
                                  shape = "Alt")))+
  scale_color_manual(values = ColorblindnessFriendlyValues2, labels = c("Same", "Alt")) +
  scale_shape_manual(values = ShapeValue2, labels = c("Same", "Alt")) +
  geom_abline(color="grey",slope=0, intercept=0.08) +
  labs(color = "Order", shape = "Order")  + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),legend.position = "none",plot.title = element_text(size = 10,hjust = 0.5))+
  ggtitle("RMSEA")


plotL2 <- ggplot(data=results, aes(x=number_crossloadings))+ 
  scale_x_continuous(breaks = seq(min(results$number_crossloadings), max(results$number_crossloadings), step_tick)) +
  scale_y_continuous(limits = c(lowerbound_cfi, NA), labels = function(x) sprintf("%.2f", x)) +
  geom_line(aes(y=cfi_same_f,
                color="Same"))+ 
  suppressWarnings(geom_point(aes(y=cfi_same_f,
                                  color="Same",
                                  shape = "Same")))+
  geom_line(aes(y=cfi_dif_f, 
                color="Alt"))+ 
  suppressWarnings(geom_point(aes(y=cfi_dif_f,
                                  color="Alt",
                                  shape = "Alt")))+
  scale_color_manual(values = ColorblindnessFriendlyValues2, labels = c("Same", "Alt")) +
  scale_shape_manual(values = ShapeValue2, labels = c("Same", "Alt")) +  
  geom_abline(color="grey",slope=0, intercept=0.95)+
  labs(color = "Order", shape = "Order")  + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),legend.position = "none",plot.title = element_text(size = 10,hjust = 0.5))+
  ggtitle("CFI")

plotL3 <- ggplot(data=results, aes(x=number_crossloadings))+ 
  scale_x_continuous(breaks = seq(min(results$number_crossloadings), max(results$number_crossloadings), step_tick)) +
  scale_y_continuous(limits = c(NA,upperbound_srmr) , labels = function(x) sprintf("%.2f", x)) +
  geom_line(aes(y=srmr_same_f,
                color="Same"))+ 
  geom_line(aes(y=srmr_dif_f,color="Alt"))+ 
  suppressWarnings(geom_point(aes(y=srmr_same_f,
                                  color="Same",
                                  shape = "Same")))+
  suppressWarnings(geom_point(aes(y=srmr_dif_f,
                                  color="Alt",
                                  shape = "Alt")))+
  scale_color_manual(values = ColorblindnessFriendlyValues2, labels = c("Same", "Alt")) +
  scale_shape_manual(values = ShapeValue2, labels = c("Same", "Alt")) + 
  geom_abline(color="grey",slope=0, intercept=0.08) + 
  labs(color = "Order", shape = "Order")  + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(size = 10, hjust = 0.5))+
  ggtitle("SRMR")

# Figure 1 - Right  
set.seed(123)

p = 16
fcor = 0.2
ML = 0.7
MR = 0.3
CL = 0.2
CR = 0.2

genLoadingss <- runif(p, min=ML-.5*MR, max=ML+.5*MR) 

numCrossLoading <- runif(p*2, min=0.2 -.5*CR, max=0.2+.5*CR)

resultsAndOrder <- main.2f(p,fcor,genLoadingss,numCrossLoading)

results <- as.data.frame(resultsAndOrder$results)

#send an alart for negative residual variances
residuals <- select(tail(results, n=1),matches("x[0-9]{1,2}"))


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

plotR1 <- ggplot(data=results,aes(x=number_crossloadings))+ 
  scale_x_continuous(breaks = seq(min(results$number_crossloadings), max(results$number_crossloadings), step_tick)) +
  scale_y_continuous(limits = c(NA,upperbound_rmsea), labels = function(x) sprintf("%.2f", x)) +
  geom_line(aes(y=rmsea_same_f,
                color="Same"))+ 
  suppressWarnings(geom_point(aes(y=rmsea_same_f,
                                  color="Same",
                                  shape = "Same")))+
  geom_line(aes(y=rmsea_dif_f,color ="Alt"))+
  suppressWarnings(geom_point(aes(y=rmsea_dif_f,
                                  color ="Alt",
                                  shape = "Alt")))+
  scale_color_manual(values = ColorblindnessFriendlyValues2, labels = c("Same", "Alt")) +
  scale_shape_manual(values = ShapeValue2, labels = c("Same", "Alt")) +
  geom_abline(color="grey",slope=0, intercept=0.08) +
  labs(color = "Order", shape = "Order")  + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),legend.position = "none",plot.title = element_text(size = 10,hjust = 0.5))+
  ggtitle("RMSEA")


plotR2 <- ggplot(data=results, aes(x=number_crossloadings))+ 
  scale_x_continuous(breaks = seq(min(results$number_crossloadings), max(results$number_crossloadings), step_tick)) +
  scale_y_continuous(limits = c(lowerbound_cfi, NA), labels = function(x) sprintf("%.2f", x)) +
  geom_line(aes(y=cfi_same_f,
                color="Same"))+ 
  suppressWarnings(geom_point(aes(y=cfi_same_f,
                                  color="Same",
                                  shape = "Same")))+
  geom_line(aes(y=cfi_dif_f, 
                color="Alt"))+ 
  suppressWarnings(geom_point(aes(y=cfi_dif_f,
                                  color="Alt",
                                  shape = "Alt")))+
  scale_color_manual(values = ColorblindnessFriendlyValues2, labels = c("Same", "Alt")) +
  scale_shape_manual(values = ShapeValue2, labels = c("Same", "Alt")) +  
  geom_abline(color="grey",slope=0, intercept=0.95)+
  labs(color = "Order", shape = "Order")  + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),legend.position = "none",plot.title = element_text(size = 10,hjust = 0.5))+
  ggtitle("CFI")

plotR3 <- ggplot(data=results, aes(x=number_crossloadings))+ 
  scale_x_continuous(breaks = seq(min(results$number_crossloadings), max(results$number_crossloadings), step_tick)) +
  scale_y_continuous(limits = c(NA,upperbound_srmr) , labels = function(x) sprintf("%.2f", x)) +
  geom_line(aes(y=srmr_same_f,
                color="Same"))+ 
  geom_line(aes(y=srmr_dif_f,color="Alt"))+ 
  suppressWarnings(geom_point(aes(y=srmr_same_f,
                                  color="Same",
                                  shape = "Same")))+
  suppressWarnings(geom_point(aes(y=srmr_dif_f,
                                  color="Alt",
                                  shape = "Alt")))+
  scale_color_manual(values = ColorblindnessFriendlyValues2, labels = c("Same", "Alt")) +
  scale_shape_manual(values = ShapeValue2, labels = c("Same", "Alt")) + 
  geom_abline(color="grey",slope=0, intercept=0.08) + 
  labs(color = "Order", shape = "Order")  + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(size = 10,hjust = 0.5))+
  ggtitle("SRMR")

spacer <- plot_spacer()

plotL1 + plotL2 + plotL3 + spacer + plotR1 + plotR2 +plotR3+ plot_layout(ncol = 7, guides = "collect", widths  = c(1, 1,1, 0.1, 1, 1, 1))

