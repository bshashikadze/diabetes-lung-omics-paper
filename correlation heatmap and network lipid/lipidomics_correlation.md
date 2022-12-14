Correlation analysis, visualization as heatmap and network
================
BS
16/10/2022

### load libraries

``` r
library     (tidyverse)
library(ComplexHeatmap)
library      (circlize)
library        (ggpubr)
library         (psych)
library     (tidygraph)
library        (ggraph)
library        (igraph)
library       (ggforce)
library     (gridExtra)
library       (cowplot)
library  (RColorBrewer)
library    (ggnewscale)
library       (ggrepel)
library    (concaveman)
```

### load omics data

order of the compounds in annotation should correspond to the order of
compounds in the main data

``` r
data_raw  <- read.delim("lipidomics_raw.txt", sep = "\t", header = T, check.names = F) %>% 
  filter(if_all(-contains("Compound"), ~ (.x >= 0))) 
# class of each column
sapply(data_raw, class)
```

    ##    Compound     MIDY737     MIDY739     MIDY740     MIDY744       WT736 
    ## "character" "character" "character" "character" "character" "character" 
    ##       WT738       WT741       WT743       WT745 
    ## "character" "character" "character" "character"

``` r
# convert all but first column to numeric
data_raw[ , 2:10] <- apply(data_raw[ , 2:10], 2,            
                    function(x) as.numeric(as.character(x)))
# read annotation file
annotation_all <- read.delim("annotation.txt", sep = "\t", row.names = "Compound", header = T, check.names = F)
```

### make sure row order in annotation file and main file have the same order

``` r
order <- data_raw$Compound
annotation_all <- annotation_all %>% 
  arrange(factor(rownames(annotation_all), levels = order)) 
```

### define groups

``` r
Bioreplicate  <- names(data_raw)[-1]
Condition     <- Bioreplicate
Groups        <- data.frame(Bioreplicate, Condition)
write.table(Groups, "Groups.txt", row.names = F, quote = F, sep = "\t")
rm(Bioreplicate, Condition, Groups)
cat("Groups file was generated, rename the file as Groups_modified, 
    and modify the second column according to experimental condition")
```

    ## Groups file was generated, rename the file as Groups_modified, 
    ##     and modify the second column according to experimental condition

### prepare data for the correlation analysis

``` r
multivariate_prep_function <- function(data){
  data <- data %>% 
          column_to_rownames(var="Compound") %>% 
          log2() %>% 
          t() %>% 
          as.data.frame()
}
data_for_corr <- multivariate_prep_function(data_raw)
```

### calculation of correlation coefficients and significance

1)  matrix with all the pair-wise correlation coefficients
2)  matrix with a p-values, lower triangle contains raw p-values, upper
    contains p-values adjusted for all the pair-wise comparisons

``` r
corr_function <- function(data, corr_method = "spearman", padjmethod = "BH", adjusted = T) 
  
   {
   correlations <- psych::corr.test(data, method=corr_method, 
                                    adjust=padjmethod, 
                                    alpha=.05,ci=F,
                                    minlength=5,
                                    normal=TRUE)
    
   pvals        <- as.data.frame(correlations$p)
   diag(pvals)  <- NA #replaces MA in the diagonal
   
   #makes an empty list where all the outcome will be stored
   cor_data <- list() 
   #get correlation matrix
   corrs_r     <- as.data.frame(correlations$r)
   cor_data[[1]]   <- corrs_r
   cor_data[[2]]   <- pvals
  
   if (adjusted == TRUE) {
   
   pvals[lower.tri(pvals)] <- NA
   #copy values from upper triangle to lower triangle *mirror* 
   for(i in 1:nrow(pvals)) {for(j in 1:i) {pvals[i,j]=pvals[j,i]}}
   cor_data[[3]] <- pvals
   cat(paste("this function returned a list of 3. The first element contains ", corr_method, "correlation coeficients,
             the second element contains raw p-values in lower, and",
             padjmethod, "adjusted p-values in an upper triangle, and the third element contains adjusted p-values only"))
    
     return(cor_data)}
    
   else {pvals[upper.tri(pvals)] <- NA
   #copy values from lower triangle to upper triangle *mirror* 
   for(i in 1:nrow(pvals)) {for(j in 1:i) {pvals[j,i]=pvals[i,j]}}
   cor_data[[3]] <- pvals
   cat(paste("this function returned a list of 3. The first element contains ", corr_method, "correlation coeficients,
             the second element contains raw p-values in lower, and",
             padjmethod, "adjusted p-values in an upper triangle, and the third element contains raw p-values only"))
     
     return(cor_data)}
      
   }
#apply the function
data_corr <- corr_function(data_for_corr, corr_method = "spearman", adjusted = T, padjmethod = "BH")
```

    ## this function returned a list of 3. The first element contains  spearman correlation coeficients,
    ##              the second element contains raw p-values in lower, and BH adjusted p-values in an upper triangle, and the third element contains adjusted p-values only

## hierarchical clustering of the correlation matrix

### heatmap including all correlations

#### partition the data with the k-means clustering

``` r
set.seed(1234)
km_clust <- kmeans(data_corr[[1]], 4)$cluster
```

#### plot heatmap for all correlations

``` r
#makes a list with all necessary data for heatmap
hmap_all_data      <- list()
hmap_all_data[[1]] <- colorRamp2(c(-1, 0, 1), c("navy",  "white", "firebrick3"))
hmap_all_data[[2]] <- list('Substrate'  = c('EPA' = "#999999", 'AA' = "#E69F00", 'DHA' = "#56B4E9", 'LA' = "navy"),
                           'Enzyme'     = c('LOX' = "#F0E442", 'CYP' = "#008080", 'COX' = "#D55E00", 'NE' = "#CC79A7"),
                           'Regulation' = c('Up in MIDY' = "#0072B2", 'Down in MIDY' = "firebrick3"))
hmap_all_data[[3]] <- HeatmapAnnotation(df = annotation_all, 
                            show_legend    = F, 
                            which          = 'col',
                            annotation_name_gp   = gpar(fontsize = 8),
                            annotation_name_side = "right",
                            col                  = hmap_all_data[[2]],
                            height               = unit(0.8, "cm"),
                            simple_anno_size_adjust = TRUE)
hmap_all_data[[4]] <- HeatmapAnnotation(df          = annotation_all, 
                            show_legend             = F, 
                            which                   = 'row',
                            #annotation_name_gp= gpar(fontsize = 7),
                            #annotation_name_side = "top",
                            show_annotation_name = F,
                            col                     = hmap_all_data[[2]],
                            width                   = unit(0.8, "cm"),
                            simple_anno_size_adjust = TRUE)
# plotting
hmap <- Heatmap(as.matrix(data_corr[[1]]),
                show_heatmap_legend        = F,
                row_dend_width             = unit(0.4, "cm"),
                column_dend_height         = unit(0.4, "cm"),
                show_row_names = F,
                show_column_names = F,
                row_title = NULL,
                column_title = NULL,
                row_split = km_clust,
                column_split = km_clust,
                clustering_distance_columns = "euclidean",
                clustering_distance_rows    = "euclidean",
                clustering_method_rows      = "complete",
                clustering_method_columns   = "complete",
                right_annotation            = hmap_all_data[[4]],
                top_annotation              = hmap_all_data[[3]],
                width                       = unit(2.98, "in"),
                height                      = unit(3.35, "in"), 
                col                         = hmap_all_data[[1]], 
                border                      = TRUE,
                row_names_gp                = gpar(fontsize = 5)) 
ht <- draw(hmap)
```

![](lipidomics_correlation_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
#calculate actual plot size
w1 = ComplexHeatmap:::width(ht)
w1 = convertX(w1, "inch", valueOnly = TRUE)
h1 = ComplexHeatmap:::height(ht)
h1 = convertY(h1, "inch", valueOnly = TRUE)
c(w1, h1)
```

    ## [1] 3.997964 4.033703

``` r
#save
heatmap_file <- paste0("Correlation_lipids_spearman_full", ".svg")
svg(file = heatmap_file, width = w1, height = h1)
ht           <- draw(hmap)
#for cowplot
gb_Hmap = grid.grabExpr(draw(hmap))
```

### save correlation matrix as it is on the heatmap

``` r
# correlation matrix
corr_matrix <- data_corr[[1]] %>% 
  rownames_to_column("Compound") 
# get the row order (it is not necessary to get the column_order in this case as it is the same as row order)
roworder_ht <- row_order(ht) %>% 
unlist()
# order the lipids
compound_order <- corr_matrix %>% 
  select(Compound) %>% 
  mutate(order = seq_along(1:length(roworder_ht))) %>% 
  arrange(factor(order, levels = roworder_ht)) %>% 
  select(-order)
# arrange rows of the corr matrix
corr_matrix <- corr_matrix %>% 
arrange(factor(Compound, levels = compound_order$Compound))
# arrange columns
Compound       <- "Compound"
compound_order <- rbind(Compound, compound_order)
corr_matrix    <- corr_matrix[compound_order$Compound]
# save correlation matrix
write.table(corr_matrix, "correlations.txt", sep = "\t", row.names = F, quote = F)
```

\####plot legends separatelly this is a current work-around to get a
desired arrangement of the legends as in the figure.

``` r
cor_legend = grid.grabExpr(color_mapping_legend(hmap@matrix_color_mapping, plot = T, title = "Correlation", legend_direction = c("horizontal"), title_position = "topcenter", title_gp = gpar(fontsize = 7, fontface = "bold"),  labels_gp = gpar(fontsize = 7)))
Substrate = grid.grabExpr(color_mapping_legend(hmap@top_annotation@anno_list[["Substrate"]]@color_mapping, plot = T, title_gp = gpar(fontsize = 7, fontface = "bold"),  labels_gp = gpar(fontsize = 7)))
Enzyme = grid.grabExpr(scale = 1, color_mapping_legend(hmap@top_annotation@anno_list[["Enzyme"]]@color_mapping, plot = T, title_gp = gpar(fontsize = 7, fontface = "bold"),  labels_gp = gpar(fontsize = 7)))
Regulation = grid.grabExpr(color_mapping_legend(hmap@top_annotation@anno_list[["Regulation"]]@color_mapping, plot = T, title_gp = gpar(fontsize = 7, fontface = "bold"),  labels_gp = gpar(fontsize = 7)))
legend <- ggarrange(arrangeGrob(cor_legend, Regulation, ncol=1, heights = c(1,1)),
           Substrate, Enzyme, heights = c(1,1), widths =  c(1,1),
           ncol = 3) 
ggsave("legend.svg", width =2.5, height = 1)
rm(cor_legend, Substrate, Enzyme, Regulation, legend)
```

### heatmap including interesting correlations

n_clust, start_n and end_n is chosen after inspecting hierarchical
clustering and choosing specific cluster to be separately visualized.
For this row_names and column_names should be enabled to know the order
of the cluster. start_n and end_n corresponds to the order as it is on
the heatmap.

``` r
extract_features_function <- function(data, n_clust = 3, start_n = 1, end_n =12) 
  
  {
   #get an order of the rows in the main heatmap (row_order function)
   r_ord <- row_order(hmap)
   row_order_list = r_ord[[paste0(n_clust)]][start_n:end_n]
  
   #prepare annotation file for the subset
   annotation_subset <- data %>% 
   as.data.frame() %>% 
   rownames_to_column("Compound") %>% 
   rename(cluster_n = 2) %>% 
   left_join(annotation_all %>% 
              rownames_to_column("Compound")) %>% 
   rownames_to_column("row_order") %>% 
   filter(row_order %in% row_order_list) %>% 
    select(-row_order, -cluster_n)
   #extract subset compound names (will be used for the filtering)
   subs_names <- annotation_subset$Compound
  
   #get correlation coefficients
   subset_r <- annotation_subset %>% 
   select(Compound) %>% 
   left_join(data_corr[[1]] %>% 
              as.data.frame() %>% 
              rownames_to_column("Compound")) %>% 
   column_to_rownames("Compound") %>% 
   select(subs_names)
  
   #get an adjusted p-values
   subset_padj <- annotation_subset %>% 
   select(Compound) %>% 
   left_join(data_corr[[3]] %>% 
              as.data.frame() %>% 
              rownames_to_column("Compound")) %>% 
   select(subs_names, "Compound") %>% 
   column_to_rownames("Compound")
  
   #replace adjusted p values with stars (will be displayed on plot)
   subset_padj_stars <- subset_padj %>% 
   rownames_to_column("Compound") %>% 
   pivot_longer(names_to = "Compound2", values_to = "padj", -Compound) %>% 
   mutate(padj_star = case_when(padj < 0.001 &  padj >  0     ~ "***",
                         padj < 0.01         &  padj >  0.001 ~ "**",
                         padj < 0.05         &  padj >  0.01  ~ "*",
                         padj > 0.05                          ~ "",
                         is.na(padj)                          ~ "")) %>% 
   select(-padj) %>% 
   pivot_wider(names_from = "Compound", values_from = "padj_star") %>% 
   column_to_rownames("Compound2")
   
   #include quantitative values for the boxplot
   if (file.exists("Groups_modified.txt")) {
   Groups <- read.delim("Groups_modified.txt")
   data_quant <- annotation_subset %>% 
     select("Compound") %>% 
     left_join(data_raw) %>% 
     pivot_longer(names_to = "Bioreplicate", values_to = "Concentration",  -Compound) %>%
     mutate(log2conc = log2(Concentration)) 
   
   class(Groups$Bioreplicate) <- class(data_quant$Bioreplicate)
   
   data_quant <- data_quant %>% 
     left_join(Groups)
   }
   
   else{cat("Groups_modified file does not exist in the folder")}
   
   
   # store all data in the list
   subs_data      <- list()
   subs_data[[1]] <- annotation_subset
   subs_data[[2]] <- subset_r
   subs_data[[3]] <- subset_padj
   subs_data[[4]] <- subset_padj_stars
   subs_data[[5]] <- data_quant
   cat("this function returned a list of 5. the first element is annotation information. The second element is correlation coefficients, the 
        third element is adjusted p-values, the fourth element is adjusted p values replaced with starts, 
        the fifth element is quantitative data")
   
   return(subs_data)
   
}
#apply the function
subset_data  <- extract_features_function(km_clust)
```

    ## this function returned a list of 5. the first element is annotation information. The second element is correlation coefficients, the 
    ##         third element is adjusted p-values, the fourth element is adjusted p values replaced with starts, 
    ##         the fifth element is quantitative data

### plot heatmap for interesting correlations

``` r
hmap_subset_data      <- list()
hmap_subset_data[[1]] <- HeatmapAnnotation(df = column_to_rownames(subset_data[[1]], "Compound"), show_legend = F, 
                            which = 'col',
                            annotation_name_gp= gpar(fontsize = 8),
                            annotation_name_side = "right",
                            col = hmap_all_data[[2]],
                            height = unit(0.8, "cm"),
                            simple_anno_size_adjust = TRUE)
hmap_subset_data[[2]] <- HeatmapAnnotation(df = column_to_rownames(subset_data[[1]], "Compound"), 
                            show_legend = F, 
                            which = 'row',
                            show_annotation_name = F,
                            col = hmap_all_data[[2]],
                            width = unit(0.6, "cm"),
                            simple_anno_size_adjust = TRUE)
# plot
hmap_subset <- Heatmap(as.matrix(subset_data[[2]]), 
                show_heatmap_legend   = F,
                row_dend_width        = unit(0.25, "cm"),
                column_dend_height = unit(0.25, "cm"),
                show_row_names     = T,
                show_column_names  = T,
                row_title = NULL,
                column_title = NULL,
                clustering_distance_columns = "euclidean",
                clustering_distance_rows    = "euclidean",
                clustering_method_rows      = "complete",
                clustering_method_columns   = "complete",
                top_annotation   = hmap_subset_data[[1]],
                right_annotation = hmap_subset_data[[2]],
                width  = unit(1.5, "in"),
                height = unit(1.32, "in"),
                col    = hmap_all_data[[1]], 
                border = TRUE,
                row_names_gp    = gpar(fontsize = 8),
                column_names_gp = gpar(fontsize = 8),
                cell_fun = function(j, i, x, y, width, height, fill) {
                                   grid.text((subset_data[[4]][i, j]), x, y, 
                                             gp = gpar(fontsize = 5.5))})
ht_subs <- draw(hmap_subset)
w2 = ComplexHeatmap:::width(ht_subs)
w2 = convertX(w2, "inch", valueOnly = TRUE)
h2 = ComplexHeatmap:::height(ht_subs )
h2 = convertY(h2, "inch", valueOnly = TRUE)
c(w2, h2)
```

    ## [1] 3.103204 3.001944

``` r
# save
heatmap_file <- paste0("Correlation_lipids_spearman_subset", ".svg")
svg(file = heatmap_file, width = w2, height = h2)
ht_subs <- draw(hmap_subset)
# for cowplot
gb_hmap_subset = grid.grabExpr(draw(hmap_subset))
```

## univariate scatter plot of interesting features

### prepare data

``` r
# prepare data for error bar calculation
ratios_error_bar <- subset_data[[5]] %>% 
group_by(Compound, Condition) %>% 
summarise(mean = mean(Concentration), 
          sd = sd(Concentration)) %>% 
  ungroup()
```

    ## `summarise()` has grouped output by 'Compound'. You can override using the
    ## `.groups` argument.

``` r
# prepare data for plotting
data_plot <- subset_data[[5]] %>%  
left_join(ratios_error_bar)
```

    ## Joining, by = c("Compound", "Condition")

### plot bar chart with data

``` r
# reorder data
data_plot$Condition <- factor(data_plot$Condition, levels = c("WT", "MIDY"))

# plot 
bar_plot <- ggplot(data_plot, aes(x=Condition, y=Concentration))+
  geom_bar(stat = "summary",fun = mean, fill = "white", alpha = 1, color = "black", lwd=0.5,
           width=0.7/length(unique(data_plot$Condition))) +
  geom_jitter(size = 1.8, width = 0.15, shape = 21, aes(fill = Condition), alpha = 0.8) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05)) +
  scale_fill_manual(values=c("#0072B2", "firebrick3"))+
  theme_bw()+
  ylab("Lipid concentration (ng/g tissue)")+
  xlab("")+
  theme(panel.border = element_rect(size=0.8, color = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.text = element_text(size = 9, colour = "black"),
        axis.title.x = element_text(size = 9, colour = "black"), 
        axis.title.y = element_text(size = 9, colour = "black"), 
        panel.grid =   element_blank())+
  facet_wrap(~Compound, scales = "free_y", ncol = 6)+
  stat_compare_means(data =  data_plot %>% 
                       select(-log2conc) %>% 
                       mutate(Concentration = log2(Concentration)),
                     aes(label = sprintf("p = %5.2f", as.numeric(..p.format..))),
                     label.x = 2, 
                     hjust = 1, 
                     method = "t.test", 
                     size = 2.75, 
                     vjust = 2.2, 
                     angle =90)+
  facet_wrap(~Compound, scales = "free_y", ncol = 6)+
  theme(legend.position = "bottom", 
        legend.box.spacing = unit(0.8, 'mm'), 
        legend.key.height= unit(1, 'mm'),
        legend.title = element_blank(), 
        legend.text = element_text(size = 8))+
  theme(strip.text.x = element_text(size = 8), strip.background = element_blank(),legend.margin=margin(0,0,0,0), legend.box.margin=margin(-12,-12,-3,-12))
```

### correlation network

``` r
corr_netw_function <- function(r_threshold = 0.8, padj_threshold = 0.05) 
    {
     cor_netw_r <- data_corr[[1]] %>% 
     replace(., upper.tri(., diag = T), NA) %>% 
     rownames_to_column("x") %>% 
     pivot_longer(names_to  = "y", values_to = "r", values_drop_na = T, -x) %>% 
     filter(abs(r) >= r_threshold)
     
     x <- as.data.frame(cor_netw_r$x) 
     y <- as.data.frame(cor_netw_r$y)
     y_rename <- setNames(y, names(x)) 
     annot <- rbind(x,y_rename) %>% 
     unique() %>% 
     rename(Compound = 1) %>% 
     left_join(annotation_all %>% rownames_to_column("Compound"))
     
     corr_data_netw <- list()
     corr_data_netw[[1]] <- cor_netw_r
     corr_data_netw[[2]] <- annot 
  if (exists("subset_data")) {
    
     cor_netw_r_subs <- subset_data[[2]] %>% 
     replace(., upper.tri(., diag = T), NA) %>% 
     rownames_to_column("x") %>% 
     pivot_longer(names_to  = "y", values_to = "r", values_drop_na = T, -x) %>% 
      left_join(subset_data[[3]] %>% 
                  replace(., upper.tri(., diag = T), NA) %>% 
     rownames_to_column("x") %>% 
     pivot_longer(names_to  = "y", values_to = "padj", values_drop_na = T, -x)) %>% 
      filter(padj <= padj_threshold) %>% 
       select(-padj)
    
     annot_subs <- annot %>% 
     dplyr::filter(Compound %in% cor_netw_r_subs$x | Compound %in% cor_netw_r_subs$y)
  
     corr_data_netw[[3]] <- cor_netw_r_subs 
     corr_data_netw[[4]] <- annot_subs
     cat(paste("this function returned a list of 4. the first is pairwise correlations filtered for an absolute correlation coefficient more 
               than ",  r_threshold = 0.8, "the second element is annotation information for the compounds, the third and the fourth elements                are similar to the first and the second but fot the subset which was futher filtered to keep the correlations with a 
               significance less than ", padj_threshold))
       
       return(corr_data_netw)
  }
  else{cat(paste("this function returned a list of 2. the first is pairwise correlations filtered for an absolute correlation coefficient more
                 than ", r_threshold = 0.8, "the second element is annotation information for the compounds"))
               
      return(corr_data_netw)}
      
  }
#apply the function
data_corr_netw <- corr_netw_function()
```

    ## Joining, by = "Compound"
    ## Joining, by = c("x", "y")

    ## this function returned a list of 4. the first is pairwise correlations filtered for an absolute correlation coefficient more 
    ##                than  0.8 the second element is annotation information for the compounds, the third and the fourth elements                are similar to the first and the second but fot the subset which was futher filtered to keep the correlations with a 
    ##                significance less than  0.05

### prepare data to plot the network with ggplot

``` r
corr_netw_forggplot_function <- function(data, community_detection = TRUE, community_detection_alg = cluster_walktrap, netw_layout = "fr")  
  
  {
  #network with a tidygraph package
  corr.graph <- as_tbl_graph(data, directed = FALSE)
  corr.graph %>% activate(nodes) 
  
  #number of adjustent edges
  deg <- degree(corr.graph, mode="all") %>% 
  as.data.frame() %>% 
    rename(n_edge = 1) %>% 
    rownames_to_column("Compound")
  #plot layout with a ggraph (to get an x and y coordinates)
  #according to https://ggraph.data-imaginist.com/articles/Nodes.html 
  df_nodes   <- create_layout(corr.graph, layout = netw_layout) %>% 
    select(x,y,name) %>% 
    rename(Compound = name) %>% 
    left_join(annotation_all %>% rownames_to_column("Compound")) %>% 
    left_join(deg)
  
  #extract edges
  # according to the answer: https://stackoverflow.com/questions/63997659/get-edge-data-from-the-tidygraph-package
  df_edges <- corr.graph %>%
  activate(edges) %>%
  get.edgelist() %>%
  data.frame() %>% 
    rename(edge_from = X1, edge_to = X2)
  
  #add the correlation coefficients (can be used to manipulate the edge thickness)
  df_edges <- cbind(df_edges, data %>% select("r"))
  
  #add the start and end coordinates for the plotting
  df_edges <-  df_edges %>% 
    left_join(df_nodes  %>% select(x,y, Compound) %>%  rename(x_start = x, y_start = y), by = c("edge_from" = "Compound")) %>% 
    left_join(df_nodes  %>% select(x,y, Compound) %>%  rename(x_end = x, y_end = y), by = c("edge_to" = "Compound")) 

  
  corr_netw_ggplot_data <- list()
  
  if (community_detection == TRUE) {
    
    communityalg <- community_detection_alg
    community    <-  communityalg (corr.graph)$membership %>% 
    as.data.frame %>% 
    rename(Community = 1)
  
    community  <- cbind(df_nodes$Compound, community)
    names(community)[1] <- "Compound" 
    
    df_edges   <-  df_edges %>% 
    left_join(community %>% rename(comm_start = Community), by = c("edge_from" = "Compound")) %>% 
    left_join(community %>% rename(comm_end = Community), by = c("edge_to" = "Compound")) %>% 
    mutate(comm_col = case_when(comm_start == comm_end ~ "within",
                               TRUE ~ "between"))
    df_nodes   <- df_nodes %>% 
    left_join(community) 
    
    corr_netw_ggplot_data[[1]] <- df_edges
    corr_netw_ggplot_data[[2]] <- df_nodes
    return(corr_netw_ggplot_data)
    
  }
  
  else  { 
    corr_netw_ggplot_data[[1]] <- df_edges
    corr_netw_ggplot_data[[2]] <- df_nodes
    return(corr_netw_ggplot_data)}
  
}


set.seed(8889)
corr_plot_ggplot <- corr_netw_forggplot_function(data_corr_netw[[1]], community_detection = T, netw_layout = 'fr')
```

    ## Joining, by = "Compound"
    ## Joining, by = "Compound"
    ## Joining, by = "Compound"

### geom_mark_hull does not make circles when communitues contains only two entries, but works when I specify geom_mark_hull for specifically those communities that contains one entries

``` r
netwenzyme <- ggplot(mapping = aes(x=x,y=y)) +
              geom_mark_hull(data = corr_plot_ggplot[[2]], aes(group = as.factor(Community), #community
                                                           filter = Community == 12), 
                                                           fill = '#B15928',
                                                           col  = '#B15928',
                                                           concavity = 5,
                                                           expand = unit(2.5, "mm"),
                                                           radius = unit(2.5, "mm"),
                                                           alpha  = 0.5,
                                                           lwd = 0.5, 
                                                           show.legend = F) +
                geom_mark_hull(data = corr_plot_ggplot[[2]], aes(group = as.factor(Community), #community
                                              filter = Community == 4), 
                                                           fill = "#FFFF99",
                                                           col  = "#FFFF99",
                                                           concavity = 5,
                                                           expand = unit(2.5, "mm"),
                                                           radius = unit(2.5, "mm"),
                                                           alpha  = 0.5,
                                                           lwd = 0.5, 
                                                           show.legend = F) +
                geom_mark_hull(data = corr_plot_ggplot[[2]], aes(group = as.factor(Community), #community
                                              filter = Community == 5), 
                                                           fill = "#6A3D9A",
                                                           col  = "#6A3D9A",
                                                           concavity = 5,
                                                           expand = unit(2.5, "mm"),
                                                           radius = unit(2.5, "mm"),
                                                           alpha  = 0.5,
                                                           lwd = 0.5, 
                                                           show.legend = F) +
                geom_mark_hull(data = corr_plot_ggplot[[2]], aes(group = as.factor(Community), #community
                                                           fill      = as.factor(Community),  
                                                           colour    = as.factor(Community)), 
                                                           concavity = 6,
                                                           expand = unit(2.5, "mm"),
                                                           radius = unit(2.5, "mm"),
                                                           alpha  = 0.5,
                                                           lwd = 0.5, 
                                                           show.legend = F) +
                                        scale_fill_brewer(palette = "Paired") + 
                                        scale_colour_brewer(palette = "Paired") + 
  new_scale_color()+
  geom_segment(data=corr_plot_ggplot[[1]], aes(x=x_start, xend = x_end, y=y_start,yend = y_end, col = comm_col), 
               show.legend = F, size = 0.2) +
  scale_color_manual(values = c("between" = "firebrick3", "within" = "black"))+  
  new_scale_fill()+
  geom_point(data = corr_plot_ggplot[[2]], size= 2.8, shape = 21, lwd= 0.2, aes(fill = Enzyme))+
  scale_fill_manual(values = c('LOX' = "#F0E442", 'CYP' = "#008080", 'COX' = "#D55E00", 'NE' = "#CC79A7"))+
  scale_x_continuous(expand=c(0,1))+  # expand the x limits 
  scale_y_continuous(expand=c(0,1))+ # expand the y limits
  theme_void()+
  theme(plot.margin = margin(1,1,1,1, "mm"))+
  theme(legend.position = "bottom", legend.text = element_text(size = 8), legend.title = element_blank()) 
```

    ## Warning: Duplicated aesthetics after name standardisation: size

``` r
ggsave("netwenzyme.svg", width = 2.5, height = 2.5)
```

### geom_mark_hull does not make circles when communitues contains only two entries, but works when I specify geom_mark_hull for specifically those communities that contains one entries

``` r
netwsubstrate <- ggplot(mapping = aes(x=x,y=y)) +
                                                      geom_mark_hull(data = corr_plot_ggplot[[2]], aes(group = as.factor(Community), #community
                                                           filter = Community == 12), 
                                                           fill = '#B15928',
                                                           col  = '#B15928',
                                                           concavity = 5,
                                                           expand = unit(2.5, "mm"),
                                                           radius = unit(2.5, "mm"),
                                                           alpha  = 0.5,
                                                           lwd = 0.5, 
                                                           show.legend = F) +
                geom_mark_hull(data = corr_plot_ggplot[[2]], aes(group = as.factor(Community), #community
                                              filter = Community == 4), 
                                                           fill = "#FFFF99",
                                                           col  = "#FFFF99",
                                                           concavity = 5,
                                                           expand = unit(2.5, "mm"),
                                                           radius = unit(2.5, "mm"),
                                                           alpha  = 0.5,
                                                           lwd = 0.5, 
                                                           show.legend = F) +
                geom_mark_hull(data = corr_plot_ggplot[[2]], aes(group = as.factor(Community), #community
                                              filter = Community == 5), 
                                                           fill = "#6A3D9A",
                                                           col  = "#6A3D9A",
                                                           concavity = 5,
                                                           expand = unit(2.5, "mm"),
                                                           radius = unit(2.5, "mm"),
                                                           alpha  = 0.5,
                                                           lwd = 0.5, 
                                                           show.legend = F) +
                geom_mark_hull(data = corr_plot_ggplot[[2]], aes(group = as.factor(Community), #community
                                                           fill      = as.factor(Community),  
                                                           colour    = as.factor(Community)), 
                                                           concavity = 6,
                                                           expand = unit(2.5, "mm"),
                                                           radius = unit(2.5, "mm"),
                                                           alpha  = 0.5,
                                                           lwd = 0.5, 
                                                           show.legend = F) +
                                        scale_fill_brewer(palette = "Paired") + 
                                        scale_colour_brewer(palette = "Paired") + 
  new_scale_color()+
  geom_segment(data=corr_plot_ggplot[[1]], aes(x=x_start, xend = x_end, y=y_start,yend = y_end, col = comm_col), 
               show.legend = F, size = 0.2) +
  scale_color_manual(values = c("between" = "firebrick3", "within" = "black"))+  
  new_scale_fill()+
  geom_point(data = corr_plot_ggplot[[2]], size= 2.8, shape = 21, lwd= 0.2, aes(fill = Substrate))+
  scale_fill_manual(values = c('EPA' = "#999999", 'AA' = "#E69F00", 'DHA' = "#56B4E9", 'LA' = "navy"))+
  scale_x_continuous(expand=c(0,1))+  # expand the x limits 
  scale_y_continuous(expand=c(0,1))+ # expand the y limits
  theme_void()+
  theme(plot.margin = margin(1,1,1,1, "mm"))+
  theme(legend.position = "bottom", legend.text = element_text(size = 8), legend.title = element_blank()) 
```

    ## Warning: Duplicated aesthetics after name standardisation: size

``` r
ggsave("netwsubstrate.svg", width = 2.5, height = 2.5)
```

``` r
corr_plot_ggplot_subs <- corr_netw_forggplot_function(data_corr_netw[[3]], community_detection = F, netw_layout = "stress")
```

    ## Joining, by = "Compound"
    ## Joining, by = "Compound"

``` r
netwenzyme_subset <- ggplot(mapping = aes(x=x,y=y)) +
  geom_segment(data=corr_plot_ggplot_subs[[1]], aes(x=x_start, xend = x_end, y=y_start,yend = y_end, size = r), 
               color= "grey", show.legend = F) +
  scale_size_continuous(range = c(0.05, 0.8), guide = "none") +
  new_scale("size")+
  geom_point(data = corr_plot_ggplot_subs[[2]], show.legend = T, shape = 21, aes(fill = Enzyme, size = n_edge))+
  scale_fill_manual(values = c('LOX' = "#F0E442", 'CYP' = "#008080"))+
  scale_size_continuous(range = c(1, 4), guide = "none") +
  geom_text_repel(data = corr_plot_ggplot_subs[[2]], aes(label = Compound), size =3.2) +
  theme_void()+
  theme(plot.margin = margin(1,1,1,1, "mm"))+
  theme(legend.position = "bottom", legend.text = element_text(size = 8), legend.title = element_blank())
ggsave("netwenzyme_subset.svg", width = 1.5, height = 1.5)
```

## combine plots

``` r
p1 <- ggarrange(gb_Hmap, gb_hmap_subset, widths = c(w1, w2), heights = c(h1, h2), labels = c("A"), font.label = list(size = 17, face = 'bold'))
p2 <- ggarrange(netwenzyme, netwenzyme_subset, widths = c(2.4, 2.3), heights = c(2.5, 2.5), labels = c("C"), font.label = list(size = 17, face = 'bold'))
p3 <- ggarrange(netwsubstrate, p2, ncol = 2, widths = c(2.4, 4.7), heights = c(2.5,2.5), labels = c("B"), font.label = list(size = 17, face = 'bold'))
p4 <- ggarrange(p1, p3, ncol = 1, widths = c(w1+w2, w1+w2), heights = c(h1, 2.5))
p5 <- ggarrange(p4, bar_plot, ncol = 1, widths = c(w1+w2, w1+w2), heights = c(h1+2.5, 3), labels = c("", "D"), font.label = list(size = 17, face = 'bold'))
ggsave("correlationfigure.svg", height = h1+2.5+3, width = 7.1)
```
