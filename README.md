### Reproducibility repo for Shashikadze et al. "Multi-omics profiling of lungs from genetically diabetic pigs reveals reduced ALOX15 levels and impaired eicosanoid switching as mechanisms contributing to a proinflammatory milieu".

Proteomics dataset can be downloaded from [ProteomeExchange](http://www.proteomexchange.org/) repository (dataset identifier - PXD038014)

Lipidomics dataset can be found in supplementary materials of the paper

R scripts to reproduce the statistical analysis and visualization present in the manuscript
* [code](https://github.com/bshashikadze/diabetes-lung-omics-paper/blob/main/proteomics%20bioinformatics/proteomics%20bioinformatics.Rmd) to reproduce figure 1 (and supplementary figure 2)
* [code]([https://github.com/bshashikadze/diabetes-lung-omics-paper/blob/main/Quantitative%20histomorphological%20analyses/Quantitative-histomorphology.md](https://github.com/bshashikadze/diabetes-lung-omics-paper/blob/main/Quantitative%20histomorphological%20analyses/Quantitative%20histomorphology.Rmd)) to reproduce figure 2 (and supplementary figure 3)
* [code]([https://github.com/bshashikadze/diabetes-lung-omics-paper/blob/main/correlation%20heatmap%20and%20network%20lipid/lipidomics_correlation.md](https://github.com/bshashikadze/diabetes-lung-omics-paper/blob/main/correlation%20heatmap%20and%20network%20lipid/lipidomics_correlation.Rmd)) to reproduce figure 3
* [code]([https://github.com/bshashikadze/diabetes-lung-omics-paper/tree/main/lipidomics%20bioinformatics](https://github.com/bshashikadze/diabetes-lung-omics-paper/blob/main/lipidomics%20bioinformatics/lipidomics%20bioinformatics.Rmd))  to reproduce figure 4
* [code]([https://github.com/bshashikadze/diabetes-lung-omics-paper/tree/main/multiomics%20coinertia](https://github.com/bshashikadze/diabetes-lung-omics-paper/blob/main/multiomics%20coinertia/CIA.Rmd))  to reproduce figure 5
* [code]([https://github.com/bshashikadze/diabetes-lung-omics-paper/blob/main/correlation%20SLRPs/scatter-plot-of-SLRPs.md](https://github.com/bshashikadze/diabetes-lung-omics-paper/blob/main/correlation%20SLRPs/scatter%20plot%20of%20SLRPs.Rmd))  to reproduce supplementary figure 1
* [code]([https://github.com/bshashikadze/diabetes-lung-omics-paper/blob/main/prep%20suppl%20tables/suppl_tables.md](https://github.com/bshashikadze/diabetes-lung-omics-paper/blob/main/prep%20suppl%20tables/suppl_tables.Rmd))  to reproduce supplementary tables

Differential abundance analysis of proteomics data was performed using [MS-EmpiRe](https://github.com/zimmerlab/MS-EmpiRe) and [pepquantify](https://github.com/bshashikadze/pepquantify) ([script](https://github.com/bshashikadze/diabetes-lung-omics-paper/blob/main/DIA%20quant%20with%20MS-EmpiRe/DIA_NN_output_quant_with_MS_EmpiRe.Rmd))
