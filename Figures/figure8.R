library(ggplot2)
library(dplyr)
library("patchwork")
library(cowplot)
source("functions_app.R")

ColorblindnessFriendlyValues4 <- c("#648FFF", "#D81B60", "#FFB000", "#004D40")

ShapeValue2 <- c(15,19,17,3)

# Figures - Left 

p = 24
fcor = 0.2

genLoadingss <- c(0.6362733,0.7864915,0.6726931,0.8149052,0.8321402,0.5636669,0.7084316,0.8177257,0.7154305,0.6869844,0.8370500,0.6860002,0.7532712,0.7217900,0.5808774,0.8199475,0.6238263,0.5626179,0.6483762,0.8363511,0.8168618,0.7578410,0.7421520,0.8482809)

numCrossLoading <- c(0.062282320,0.083412187,0.017626410,0.037656808,-0.084336105,-0.141154541,0.185209693,0.160919618,0.076282111,0.118186967,-0.190154526,-0.008881612,0.103383815,-0.113436826,-0.072727597,-0.107349686,-0.142879991,-0.034181466,-0.034510269,-0.052461820,-0.139022101,-0.144477575,-0.106786360,-0.013615020,-0.093610944,0.143131086,-0.181667533,-0.023119970,0.119569938,-0.151240296,0.024379194,-0.117387444,-0.148987340,0.101323146,0.158018144,-0.050214890,0.066046078,-0.162063736,-0.046412145,-0.090246542,0.125856016,-0.020593463,0.124025741,0.124955804,0.117736928,-0.024067325,0.101790063,0.051688453)

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

genLoadingss <- c(0.5749046,0.6726384,0.7045855,0.6690647,0.6168168,0.6377049,0.7252198,0.6972734,0.8268997,0.6339385,0.7815121,0.8070551,0.7774324,0.8050906,0.6727890,0.5665051,0.7234369,0.7735156,0.8158941,0.5592366,0.7374656,0.7310614,0.5731757,0.6741507)

numCrossLoading <- c(0.1370821455,0.0808875153,0.3360903425,0.1064315191,0.3324698919,0.1538716243,0.3282854313,0.1230865054,0.1759999380,0.3461626373,0.0442607422,0.3879891191,0.0769886488,0.3757690293,0.0594373411,0.1270130562,0.0539500585,0.3300641527,0.2686182417,0.0969731729,0.1675221439,0.1311141402,0.3812379571,0.0993350820,0.1791040784,0.0598792013,0.1683144202,0.3344515469,0.2087583530,0.0871475966,0.0990268627,0.1464028953,0.0933589417,0.2978593600,0.0777480392,0.2251800742,0.1940135820,0.1777870168,0.0205602108,0.1225746907,0.2281244949,0.0001803759,0.1980052114,0.2261054884,0.1386693839,0.0753608147,0.2549667848,0.3544633645)

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
