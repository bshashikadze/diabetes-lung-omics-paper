DIA proteomics (from DIA-NN) analysis with an MS-EmpiRe
================
BS
31/01/2023

initially directory should only contain the main output of DIA-NN and
contaminats fasta file from MaxQuant (latter only necessary if it was
used during DIA-NN search)

pepquantify - <https://github.com/bshashikadze/pepquantify> Ms-EmpiRe -
<https://github.com/zimmerlab/MS-EmpiRe>

### load libraries

``` r
library      (msEmpiRe)    # differential abundance analysis
library      (seqinr)      # filtering for contaminants
library      (tidyverse)   # data manipulation
library      (pepquantify) # prepare data for MS-EmpiRe analysis
```

### read the main output of DIA-NN

implicit reading is necessary only if any modifications other than
supported by pepquantify package is necessary, e.g. filtering for
contaminats

``` r
raw_diann <- read.delim("DIA-NN_output_precursors.tsv", sep = "\t", header = T) #can be downloaded from PRIDE repository
```

### remove contaminants (contaminants fasta file from MaxQuant common contaminants)

``` r
contamintants           <- read.fasta("contaminants.fasta")
contaminant.names       <- getName(contamintants) 
raw_diann_filtered      <- raw_diann %>% 
  filter(!str_detect(Protein.Group, str_c(contaminant.names, collapse="|")))

# move original file to the separate folder (this is because "pepquantify" will read automatically largest tsv file so it is necessary to leave only filtered data in the main directory)
dir.create("original")
```

    ## Warning in dir.create("original"): 'original' already exists

``` r
file.copy(from = paste0(getwd(), "/DIA-NN_output_precursors.tsv"),
          to   = paste0(getwd(), "/original/DIA-NN_output_precursors.tsv"))
```

    ## [1] TRUE

``` r
unlink("DIA-NN_output_precursors.tsv")

# save contaminants-removed data
write.table(raw_diann_filtered, "MIDY_Lung_DIA_nocontaminants.tsv", quote = F, sep = "\t", row.names = F)
```

## quantification

### function which performs MS-EmpiRe normalization and quanfications

``` r
msempire_calculation <- function(data, data2 = data_raw, seed=1234, fc_threshold = 1.5) {
  
  require(magrittr)
  # read the data in the expressionset format and perform msempire normalization and quantification  
  # (https://github.com/zimmerlab/MS-EmpiRe/blob/master/example.R)
  msempiredata  <- msEmpiRe::read.standard(data[[1]], data[[2]],
                                            prot.id.generator = function(pep) unlist(strsplit(pep, "\\.[0-9]*$"))[1],
                                            signal_pattern="Intensity.*")
  
  # msempire calculations
  set.seed(seed = seed)
  msempire_results <- msempiredata  %>%
    msEmpiRe::normalize() %>%
    msEmpiRe::de.ana() %>%
    write.table(paste0(data[[3]], "/msempire_results_raw.txt"), sep = "\t", row.names = F)
  
  
# tidy results (pepquantify package)
pepquantify::resultstidy(data, data2,  fc_threshold = fc_threshold)}
```

### load data and filter

``` r
data_raw <- pepquantify::read_diann(experimental_library = T, include_mod_in_pepreport = T)
```

    ## [1] "in the peptide output modifications will be included (so far (v.2.1.2) only works for Carbamidomethyl(C))"
    ## conditions file was generated. First rename file as: conditions_modified.txt. afterwards,
    ##          modify ONLY the second column according to the experimental conditions. Do not change the column headers

### normalization and quantification

``` r
msempire_data <- pepquantify::pepquantify_funs(data_raw, condition1 = "MIDY", condition2 = "WT", imputation = TRUE)
```

    ##  no peptides have met to the imputation criteria. you can relax the parameters or continue without the imputation

``` r
msempire_calculation(msempire_data, fc_threshold = 1.5)
```

    ## [1] "final error: 0.762124"
    ## [1] "final error: 0.76472"
    ## [1] 0.7621244
    ## [1] 0.7647205
    ## [1] "detecting mode"
    ## [1] "mode: -0.2"
    ## [1] "Final error: 0.760565"
    ## [1] "numeric" "numeric" "numeric" "numeric" "numeric" "numeric"
    ## three files were saved in the working directory:
    ##     1 - msempire_results_raw:     this is the raw results of MS-EmpiRe
    ##     2 - msempire_results_tidy:    this is the results that has been cleaned-up and can be used in suppl tables
    ##     3 - msempire_results_volcano: some columns were adjusted to make it suitable for the volcano plot
