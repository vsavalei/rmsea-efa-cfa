library(ggplot2)
library(dplyr)
library("patchwork")
library(cowplot)
source("functions_app.R")

ColorblindnessFriendlyValues4 <- c("#648FFF", "#D81B60", "#FFB000", "#004D40")

ShapeValue2 <- c(15,19,17,3)

# Figures - Left 
# set.seed(123)
# 
p = 15
fcor = 0.2
# ML = 0.7
# MR = 0 
# CL = 0.2
# CR = 0 
# 
# genLoadingss <- runif(p, min=ML-.5*MR, max=ML+.5*MR) 
# 
# numCrossLoading <- runif(p*2, min=0.2 -.5*CR, max=0.2+.5*CR)

genLoadingss <- c(0.7974023,0.5850553,0.6399342,0.6569682,0.8395585,0.8402682,0.7099277,0.6496467,
                  0.7461531,0.7343649,0.5542389,0.7790138,0.6577916,0.5804052,0.7170392)

numCrossLoading <- c(0.1886328,0.1878924,0.2650116,0.2258099,0.2388039,0.1810820,0.2342525,0.2275769,
                     0.1692982,0.2116642,0.1346637,0.2393963,0.1434427,0.1502740,0.1259607,0.2373200,
                     0.1332983,0.2676618,0.1945955,0.2589516,0.2152762,0.1338802,0.1381034,0.2572105,
                     0.1892597,0.2323606,0.1680277,0.2559470,0.1658718,0.1411056)

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
# set.seed(123)
# 
p = 24
fcor = 0.2
# ML = 0.7
# MR = 0 
# CL = 0.2
# CR = 0 
# 
# genLoadingss <- runif(p, min=ML-.5*MR, max=ML+.5*MR) 
# 
# numCrossLoading <- runif(p*2, min=0.2 -.5*CR, max=0.2+.5*CR)

genLoadingss <- c(0.5749046,0.6726384,0.7045855,0.6690647,0.6168168,0.6377049,0.7252198,0.6972734,
                  0.8268997,0.6339385,0.7815121,0.8070551,0.7774324,0.8050906,0.6727890,0.5665051,
                  0.7234369,0.7735156,0.8158941,0.5592366,0.7374656,0.7310614,0.5731757,0.6741507)

numCrossLoading <- c(0.1764058,0.1553328,0.2510339,0.1649118,0.2496762,0.1827019,0.2481070,0.1711574,
                     0.1910000,0.2548110,0.1415978,0.2704959,0.1538707,0.2659134,0.1472890,0.1726299,
                     0.1452313,0.2487741,0.2257318,0.1613649,0.1878208,0.1741678,0.2679642,0.1622507,
                     0.1921640,0.1474547,0.1881179,0.2504193,0.2032844,0.1576803,0.1621351,0.1799011,
                     0.1600096,0.2366973,0.1541555,0.2094425,0.1977551,0.1916701,0.1327101,0.1709655,
                     0.2105467,0.1250676,0.1992520,0.2097896,0.1770010,0.1532603,0.2206125,0.2579238)

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
