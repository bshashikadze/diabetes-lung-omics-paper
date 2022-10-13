statistics and bioinformatics proteomics
================
BS
10/10/2022

### load libraries

``` r
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(cowplot)
library(grid)
```

### load the data

``` r
msempire_results  <- read.delim("msempire_results_tidy_volcano.txt")
matrisome         <- read.delim("Matrisome_manual.txt")       
string_output     <- read.delim("enrichment.Process_raw.tsv")  
```

## volcano plot

### clean-up

``` r
# do not include gene names starting with "LOC" in the plot
msempire_results_volcano <- msempire_results %>% 
  mutate(delabel = case_when(diff_abundant != "n.s." ~ accession))



msempire_results_volcano <- msempire_results_volcano %>%
               mutate(delabel = case_when(
               str_detect(delabel, "LOC") ~ "",
               TRUE ~ delabel))
```

### plot the volcano plot

``` r
volcanoplot <- ggplot(msempire_results_volcano %>%                       
                       arrange(desc(diff_abundant)), 
                       mapping = aes(x = l2fc, y = log10p, 
                                     fill=diff_abundant, label = delabel, 
                                     alpha = diff_abundant))+
         geom_point(aes(shape =diff_abundant, size = diff_abundant))+
         scale_shape_manual(values = c(n.s. = 16, downregulated_in_MIDY =21, upregulated_in_MIDY =21))+
         scale_size_manual(values=c(n.s. = 1.4, downregulated_in_MIDY =1.8, upregulated_in_MIDY =1.8))+
         scale_fill_manual(values=c("n.s." = "#4a4949", 
                                    "downregulated_in_MIDY"= "firebrick3", "upregulated_in_MIDY"="#0072B2"))+
         scale_alpha_manual(values= c("n.s." = 0.2, "downregulated_in_MIDY"= 1, "upregulated_in_MIDY"= 1))+
         geom_text_repel(data = subset(msempire_results_volcano, 
                                       adj_p_value < 0.05 & l2fc > 0.6 |adj_p_value < 0.01 & l2fc < -0.7),
                                       aes(label = delabel),
                                       size = 2,
                                       color = "black",
                                       box.padding = 0.3,
                                       alpha = 1)+
         theme_bw()+
         theme(plot.margin = margin(1,1,5,1, "mm"))+
         scale_x_continuous(limits = c(-1.9,1.9)) +
         scale_y_continuous(limits = c(0, 15.5)) +
         theme(panel.border = element_rect(size = 1), 
               panel.grid.major = element_line(), 
               panel.grid.minor = element_blank(),
               panel.background = element_blank(), 
               axis.line = element_blank(), 
               legend.position = "NONE")+
         theme(axis.title = element_text(size  = 9), 
               axis.text.x = element_text(size = 9, colour = "black", vjust = -0.1), 
               axis.text.y = element_text(size = 9, colour = "black"))+
         xlab("log2 fold change (MIDY/WT)")+
         ylab("-log10 P-value")
```

## dot plot (matrisome proteins)

### prepare the data

Significant proteins are classified according to doi:
10.1074/mcp.M111.014647

``` r
Matrisome_sig <- matrisome %>% 
  left_join(msempire_results_volcano, by = c("Genes" = "accession")) %>% 
  select("Category", "Genes", "l2fc", "log10p", "adj_p_value") %>% 
  drop_na() %>% 
  filter(adj_p_value < 0.05)
```

### plot the dot plot

``` r
#orderrows based on enrichment score for SKM
Matrisome_sig$Genes <- factor(Matrisome_sig$Genes, levels = Matrisome_sig$Genes
                             [order(Matrisome_sig$Category)])
dotplot <-  ggplot(Matrisome_sig,  mapping = aes(x = Genes, y= Category, fill = l2fc, size = log10p)) +
             geom_point(shape=21)+
             theme_bw() +
             theme(plot.margin = margin(1,1,1,-3, "mm"))+
             theme(panel.border = element_rect(size = 1),
                            axis.text.x = element_text(angle = 90, colour = "black", size = 8.5, vjust = 0.5, hjust = 1),
                            axis.title = element_text(size= 8.5),
                            axis.text.y = element_text(size = 8.5, colour = "black"))+
                      xlab("")+
                      ylab("")+
                      theme(plot.title = element_blank()) +
                      theme(panel.grid.major = element_line(), panel.grid.minor = element_blank())+
                      scale_size_continuous(name = "-log10(FDR)", range = c(2, 3.5), breaks = c(3,6,9,12))+
                      scale_fill_gradient(name = "log2 fold change", low = "firebrick3", high = "navy")+
             theme(legend.position = "bottom", 
                   legend.justification = "left",
                   legend.box.spacing = unit(-4, 'mm'), 
                   legend.title = element_text(size = 8.5),
                   legend.text = element_text(size  = 8.5))+
                   theme(legend.key.height= unit(5, 'mm'))+
                      theme(strip.background = element_blank(), strip.text = element_text(size = 8.5, face = "bold"))
```

## ORA plot (STRING pre-ranked)

### data preparation

<https://string-db.org/> Signed log-transformed p-values were used as
ranking metrics and false discovery rate was controlled at 1%. <GO:BP>
will be used for plotting

``` r
data_BP_revigo_input <- string_output %>% 
  select("TermID", "enrichment.score")
write.table(data_BP_revigo_input, "revigo_input.txt", sep = "\t", row.names = F, quote = F)
```

### perform redundancy removal with REVIGO and read the data back

parameters: Medium (0.7); would you like to remove obsolete GO terms -
Yes; What species… - Whole UniProt database; Semantic similarity -
SimRel

``` r
revigo_data <- read.csv("Revigo.csv", header = TRUE)
revigo_data_filtered <- revigo_data %>%
  filter(Eliminated == " False") %>% 
  select(TermID) %>% 
  left_join(string_output) %>% 
  filter(enrichment.score > 1) %>% 
  filter(false.discovery.rate < 0.01) %>% 
  mutate(log10FDR = -log10(false.discovery.rate))
```

    ## Joining, by = "TermID"

### plot

``` r
#orderrows based on enrichment score for SKM
revigo_data_filtered$term.description <- factor(revigo_data_filtered$term.description, levels = revigo_data_filtered$term.description
                             [order(revigo_data_filtered$enrichment.score)])
revigo_data_filtered$direction = factor(revigo_data_filtered$direction, levels=c('top','both ends','bottom'))
oraplot <- ggplot(revigo_data_filtered, aes(x = enrichment.score, y= term.description, fill = log10FDR, 
                                             size = genes.mapped), fill = NULL) +
            geom_point(shape = 21)+
            theme_bw() +
            theme(panel.border = element_rect(size = 1),
                            axis.text.x = element_text(angle = 90, colour = "black", size = 9),
                            axis.title = element_text(size= 9),
                            axis.text.y = element_text(size = 9, colour = "black"))+
                      xlab("Enrichment score")+
                      ylab("")+
                      theme(plot.title = element_text(size = 9, hjust=0.5,
                                                      face = "bold")) +
                      theme(plot.margin = margin(1,1,1,-3, "mm"))+
                      theme(panel.grid.major = element_line(), panel.grid.minor = element_blank())+
                      scale_size_continuous(name = "Genes mapped", range = c(2, 5))+
                      scale_fill_gradient(name = "-log10(FDR)", low = "#bf9deb", high = "#200b3b")+
                      theme(legend.position = "none", 
                            legend.box.spacing = unit(0.5, 'mm'), 
                            legend.title = element_text(size = 9), 
                            legend.text = 
                            element_text(size = 9))+
                      facet_grid(~direction, scales = "free")+
            theme(strip.background = element_blank(), strip.text = element_text(size = 9, face = "bold"))
```

### gets legend separately using cowplot to make it position more flexibly

``` r
# cluster plot
oraplot_for_legend <- ggplot(revigo_data_filtered, aes(x = enrichment.score, y= term.description, fill = log10FDR, 
                                             size = genes.mapped), fill = NULL) +
                  geom_point(shape = 21) +
                  theme_bw()+
                        scale_size_continuous(name = "Genes mapped", range = c(2, 5))+
                        scale_fill_gradient  (name = "-log10(FDR)", low = "#bf9deb", high = "#200b3b")+
                  theme(legend.position  = "bottom", 
                            legend.box.spacing = unit(-3, 'mm'), 
                            legend.title = element_text(size = 8.5), 
                            legend.text  = element_text(size = 8.5)) +
                  theme(legend.key.height= unit(5, 'mm'))

oraplot_legend <- cowplot::get_legend(oraplot_for_legend) 
```

### combine all plots

``` r
p0 =  rectGrob(width = 1, height = 1)
p2 <- ggarrange(p0, volcanoplot, widths = c(3.4,3.4), heights = c(3.4,3.4), labels = c("","B"))
p3 <- ggarrange(p2, dotplot, ncol = 1, widths = c(6.8,6.8), heights = c(3.4, 2.6), labels = c("", "C"))
p4 <- plot_grid(oraplot, oraplot_legend, ncol = 1, rel_heights = c(6,1), labels = c("D"))
p5 <- ggarrange(p3, p4, ncol = 1, widths = c(6.8,6.8), heights = c(6, 3.5))
ggsave("proteomics_fig1.svg", width = 6.8, height = 9.5)
```
