correlation analysis of the small leucine rich proteoglycans
================
BS
22/08/2022

### load libraries

``` r
library(tidyverse)
library(GGally)
library(cowplot)
```

### load the data

``` r
data <- read.delim("SLRP.txt", header = T, sep = "\t", row.names = "Animal")
```

### plot and save

``` r
#https://stackoverflow.com/questions/30858337/how-to-customize-lines-in-ggpairs-ggally
lowerFn <- function(data, mapping, method = "lm", ...) {
    p <- ggplot(data = data, mapping = mapping) +
    geom_point(colour = "navy", size = 1.5)+
    geom_smooth(method = method, size = 0.5, color = "black", ...)
  p
}

scatter <- ggpairs(
  data[, 1:6], lower = list(continuous = wrap(lowerFn, method = "lm")),
  diag = "blank",
  upper = list(continuous = wrap("cor", method = "spearman", size = 3))
    
)+
  theme_bw()+
  theme(panel.grid =element_blank())+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(strip.text = element_text(size = 9))+
  theme(strip.background = element_blank())+
  theme(axis.text = element_text(size = 8))

plot_grid(
   ggmatrix_gtable(scatter),
   nrow = 1)
```

![](scatter-plot-of-SLRPs_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
ggsave("slrpgcorrplot.svg", width = 4.5, height = 4.5)
```
