Correlation analysis with bootstrapping confidence interval calculation
================
BS
21/08/2022

### load libraries

``` r
library(tidyverse)
library   (ggpubr)
library     (boot)
```

### data

``` r
Animal        <- c("MIDY737", "MIDY739",  "MIDY740", "MIDY744", "WT736", "WT738",  "WT741",  "WT743",  "WT745") 
IHC_without_air <- c(0.011437962, 0.017072698, 0.013136247, 0.011739773,                                          
                     0.01397391, 0.039885123, 0.029064475, 0.023870968, 0.020428169)
group           <- c(rep("MIDY", 4), rep("WT", 5))

# combined data
data <- data.frame(Animal, IHC_without_air, group)
rm(Animal, IHC_without_air, group)
```

### plot boxplot for ALOX15

``` r
boxplot_IHC <- ggplot(data, aes(x=reorder(group,-IHC_without_air), y=100*IHC_without_air))+
           geom_boxplot(fill = "#f0f8ff", color = "black", 
                        outlier.fill = NULL, lwd=0.3,
                        width=0.7/length(unique(data$group)))+
  scale_fill_manual(values=c(WT = "#0072B2", MIDY = "firebrick3"))+
  theme_bw()+
  theme(plot.margin = margin(-1,-1,-1,-1, "mm"))+
  geom_point(size = 2.2, shape = 21, aes(fill = group)) +
  scale_y_continuous(limits = c(1,4), breaks = c(1,2,3,4))+
  ylab(expression(V["v(ALOX15-positive cells/lung)"]))+
  xlab("")+
  theme(panel.border = element_rect(size = 0.5), 
        axis.text = element_text(size=10, color = "black"),     
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10), panel.grid =   element_blank())+
  stat_compare_means(label.x = 2, aes(label = paste0("p = ", ..p.format..)),
                     method = "wilcox.test", size = 3.5, vjust = 1)+
theme(legend.position = "none", 
        legend.box.spacing = unit(0.4, 'mm'), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 8))+
  theme(strip.text.x = element_text(size = 7), strip.background = element_blank(),legend.margin=margin(0,0,0,0), legend.box.margin=margin(-12,-12,-3,-12))
ggsave("IHC_ALOX15.svg", width = 2.5, height = 2.5)
```

\`
