Quantitative histomorphological analyses
================
BS
2022-11-05

### load libraries

``` r
library(tidyverse)
library(ggpubr)
library(grid)
library(cowplot)
```

### ALOX15 immunohistochemistry

``` r
Animal              <- c("MIDY737", "MIDY739",  "MIDY740", "MIDY744", "WT736", "WT738",  "WT741",  "WT743",  "WT745") 
IHC_without_air     <- c(0.011, 0.017, 0.013, 0.012, 0.014, 0.040, 0.030, 0.024, 0.021)
Condition           <- c(rep("MIDY", 4), rep("WT", 5))

# combined data
data_IHC <- data.frame(Animal, IHC_without_air, Condition)
```

### volume density of interstitial connective tissue in the lung (excluding air-filled alveolar spaces)

``` r
fibrotic_area <- c(0.122, 0.200,0.280,0.158, 0.153, 0.129,0.151,0.0559,0.150)

# combined data
data_fibrosis <- data.frame(Animal, fibrotic_area, Condition)
```

### plot data for ALOX15

#### prepare data

``` r
# prepare data for error bar calculation
ratios_error_bar <- data_IHC %>% 
  group_by(Condition) %>% 
  summarise(mean = mean(100*IHC_without_air), 
              sd = sd(100*IHC_without_air)) %>% 
  ungroup()

# prepare data for plotting
data_plot_IHC <- data_IHC %>%  
left_join(ratios_error_bar)
```

    ## Joining, by = "Condition"

``` r
# reorder data
data_plot_IHC$Condition <- factor(data_plot_IHC$Condition, levels = c("WT", "MIDY"))

# plot 
bar_plot <- ggplot(data_plot_IHC, aes(x=reorder(Condition,-IHC_without_air), y=100*IHC_without_air))+
    geom_bar(stat = "summary",fun = mean, fill = "white", alpha = 1, color = "black", lwd=0.5,
           width=0.6/length(unique(data_plot_IHC$Condition))) +
  geom_jitter(size = 2.3, width = 0.10, shape = 21, aes(fill = Condition), alpha = 0.8) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05)) +
  scale_fill_manual(values=c("#0072B2", "firebrick3"))+
  theme_bw()+
  ylab(expression(V["v(ALOX15-positive cells/lung)"]*" (%)"))+
  xlab("")+
  theme(panel.border = element_rect(size=0.8, color = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"), 
        axis.title.y = element_text(size = 10, colour = "black"), 
        panel.grid =   element_blank())+
   stat_compare_means(label.x = 1.95, aes(label = paste0("p = ", ..p.format..)),
                     method = "wilcox.test", size = 3.1, vjust = 1)+
   theme(legend.position = "bottom", 
        legend.box.spacing = unit(0.8, 'mm'), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 9))+
  theme(strip.background = element_blank(),legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-12,-12,-3,-12))
```

### add space for IHC pictures

``` r
p1 = rectGrob(width = 1, height = 1)


plot_grid(p1, p1, bar_plot, rel_widths = c(1,1,1), rel_heights = c(1,1,1),
          labels = c("A", "B", "C"),  label_size = 17, label_fontfamily = 'bold', 
          ncol = 3, nrow = 1,  legend = "bottom",
          common.legend = T)
```

    ## Warning in as_grob.default(plot): Cannot convert object of class character into
    ## a grob.

    ## Warning in as_grob.default(plot): Cannot convert object of class logical into a
    ## grob.

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## font family not found in Windows font database

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## font family not found in Windows font database

![](Quantitative-histomorphology_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
ggsave("IHC.svg", width = 7.1, height = 2.8)
```

### plot volume density of interstitial connective tissue in the lung (excluding air-filled alveolar spaces)

#### prepare data

``` r
# prepare data for error bar calculation
ratios_error_bar <- data_fibrosis %>% 
  group_by(Condition) %>% 
  summarise(mean = mean(100*fibrotic_area), 
          sd = sd(100*fibrotic_area)) %>% 
  ungroup()

# prepare data for plotting
data_plot_fibrotic <- data_fibrosis %>%  
left_join(ratios_error_bar)
```

    ## Joining, by = "Condition"

#### plot

``` r
# reorder data
data_plot_fibrotic$Condition <- factor(data_plot_fibrotic$Condition, levels = c("WT", "MIDY"))

# plot 
bar_plot_fibr <- ggplot(data_plot_fibrotic, aes(x=reorder(Condition,fibrotic_area), y=100*fibrotic_area))+
             geom_bar(stat = "summary",fun = mean, fill = "white", alpha = 1, color = "black", lwd=0.5,
           width=0.6/length(unique(data_plot_fibrotic$Condition))) +
  geom_jitter(size = 2.3, width = 0.10, shape = 21, aes(fill = Condition), alpha = 0.8) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05)) +
  scale_fill_manual(values=c("#0072B2", "firebrick3"))+
  theme_bw()+
  ylab(expression(V["v(interstitial connective tissue/lung)"]* " (%)"))+
  xlab("")+
  theme(panel.border = element_rect(size=0.8, color = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(size = 10, colour = "black"),
        axis.title.x = element_text(size = 10, colour = "black"), 
        axis.title.y = element_text(size = 10, colour = "black"), 
        panel.grid =   element_blank())+
  stat_compare_means(label.x = 0.9, aes(label = paste0("p = ", ..p.format..)),
                     method = "wilcox.test", size = 3.1, vjust = 1)+
  theme(legend.position = "bottom", 
        legend.box.spacing = unit(0.8, 'mm'), 
        legend.title = element_blank(), 
        legend.text = element_text(size = 8))+
  theme(strip.text.x = element_text(size = 7), strip.background = element_blank(),legend.margin=margin(0,0,0,0), legend.box.margin=margin(-12,-12,-3,-12))
```

### add space for staining pictures

``` r
plot_grid(p1, p1, bar_plot_fibr, rel_widths = c(1,1,1), rel_heights = c(1,1,1),
          labels = c("A", "B", "C"),  label_size = 17, label_fontfamily = 'bold',
          ncol = 3, nrow = 1,  legend = "bottom",
          common.legend = T)
```

    ## Warning in as_grob.default(plot): Cannot convert object of class character into
    ## a grob.

    ## Warning in as_grob.default(plot): Cannot convert object of class logical into a
    ## grob.

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## font family not found in Windows font database

    ## Warning in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y, :
    ## font family not found in Windows font database

![](Quantitative-histomorphology_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
ggsave("fibrosis.svg", width = 7.1, height = 2.8)
```
