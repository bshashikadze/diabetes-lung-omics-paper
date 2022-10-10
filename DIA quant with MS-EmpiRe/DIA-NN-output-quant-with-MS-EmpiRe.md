DIA proteomics (from DIA-NN) analysis with an MS-EmpiRe
================
B. Shashikadze
09/10/2022

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
raw_diann <- read.delim("MIDY_Lung_DIA.tsv", sep = "\t", header = T) 
```

### remove contaminants (contaminants fasta file from MaxQuant common contaminants)

``` r
contamintants           <- read.fasta("contaminants.fasta")
contaminant.names       <- getName(contamintants) 
raw_diann_filtered      <- raw_diann %>% 
  filter(!str_detect(Protein.Group, str_c(contaminant.names, collapse="|")))

# overwrite original file with contaminats filtered file
write.table(raw_diann_filtered, "MIDY_Lung_DIA.tsv", quote = F, sep = "\t", row.names = F)
```

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

``` r
data_raw <- pepquantify::read_diann(experimental_library = T)
```

    ## conditions file was generated. First rename file as: conditions_modified.txt. afterwards,
    ##          modify ONLY the second column according to the experimental conditions. Do not change the column headers

``` r
msempire_data <- pepquantify::pepquantify_funs(data_raw, condition1 = "MIDY", condition2 = "WT", imputation = TRUE)
```

    ##  no peptides have met to the imputation criteria. you can relax the parameters or continue without the imputation

``` r
msempire_calculation(msempire_data, fc_threshold = 1.5)
```
