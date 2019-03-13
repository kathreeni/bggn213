class 17
================

Section 1. Protein sequences from healthy and tumor tissue

``` r
library(bio3d)
```

Q1: Identify sequence regions that contain all 9-mer peptides that are only found in the tumor. Hint: You will need to first identify the sites of mutation in the above sequences and then extract the surrounding subsequence region. This subsequence should encompass all possible 9-mers in the tumor derived sequence. In other words extract the subsequence from 8 residues before and 8 residues after all point mutations in the tumor sequence.

``` r
seqs <- read.fasta("lecture18_sequences.fa")
```

Here we calculate the identity per aligned position, then use the information to find non-identical sites that do not contain gaps (indels)

The "mutant.sites" are now telling us which locations within our alignment sequence is there a misalignment

``` r
ide <- conserv(seqs$ali, method = "identity")
mutant.sites <- which(ide <1)
```

Exclude gap positions from analysis
===================================

``` r
gaps <- gap.inspect(seqs)
mutant.sites <- mutant.sites[mutant.sites %in% gaps$f.inds]
mutant.sites
```

    ## [1]  41  65 213 259

``` r
mutant.names <- paste0(seqs$ali["P53_wt", mutant.sites], mutant.sites, seqs$ali["P53_mutant", mutant.sites])
mutant.names
```

    ## [1] "D41L"  "R65W"  "R213V" "D259V"

``` r
## Sequence positions surounding each mutant site
start.position <- mutant.sites - 8
end.position <-  mutant.sites + 8

# Blank matrix to store sub-sequences
store.seqs <- matrix("-", nrow=length(mutant.sites), ncol=17)
rownames(store.seqs) <- mutant.names

## Extract each sub-sequence
for(i in 1:length(mutant.sites)) {
  store.seqs[i,] <- seqs$ali["P53_mutant",start.position[i]:end.position[i]]
}

store.seqs
```

    ##       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
    ## D41L  "S"  "P"  "L"  "P"  "S"  "Q"  "A"  "M"  "L"  "D"   "L"   "M"   "L"  
    ## R65W  "D"  "P"  "G"  "P"  "D"  "E"  "A"  "P"  "W"  "M"   "P"   "E"   "A"  
    ## R213V "Y"  "L"  "D"  "D"  "R"  "N"  "T"  "F"  "V"  "H"   "S"   "V"   "V"  
    ## D259V "I"  "L"  "T"  "I"  "I"  "T"  "L"  "E"  "V"  "-"   "-"   "-"   "-"  
    ##       [,14] [,15] [,16] [,17]
    ## D41L  "S"   "P"   "D"   "D"  
    ## R65W  "A"   "P"   "P"   "V"  
    ## R213V "V"   "P"   "Y"   "E"  
    ## D259V "-"   "-"   "-"   "-"

``` r
## First blank out the gap positions 
store.seqs[store.seqs == "-"] <- ""

## Output a FASTA file for further analysis
write.fasta(seqs=store.seqs, ids=mutant.names, file="subsequences.fa")
```

Section 2. Patient HLA typing results and HLA binding prediction:

``` r
results <- read.csv("result-3.csv")
```

``` r
#results$allele
#head(results$percentile_rank)
sort(results$allele)
```

    ##   [1] HLA-A*02:01 HLA-A*02:01 HLA-A*02:01 HLA-A*02:01 HLA-A*02:01
    ##   [6] HLA-A*02:01 HLA-A*02:01 HLA-A*02:01 HLA-A*02:01 HLA-A*02:01
    ##  [11] HLA-A*02:01 HLA-A*02:01 HLA-A*02:01 HLA-A*02:01 HLA-A*02:01
    ##  [16] HLA-A*02:01 HLA-A*02:01 HLA-A*02:01 HLA-A*02:01 HLA-A*02:01
    ##  [21] HLA-A*02:01 HLA-A*02:01 HLA-A*02:01 HLA-A*02:01 HLA-A*02:01
    ##  [26] HLA-A*02:01 HLA-A*02:01 HLA-A*02:01 HLA-A*68:01 HLA-A*68:01
    ##  [31] HLA-A*68:01 HLA-A*68:01 HLA-A*68:01 HLA-A*68:01 HLA-A*68:01
    ##  [36] HLA-A*68:01 HLA-A*68:01 HLA-A*68:01 HLA-A*68:01 HLA-A*68:01
    ##  [41] HLA-A*68:01 HLA-A*68:01 HLA-A*68:01 HLA-A*68:01 HLA-A*68:01
    ##  [46] HLA-A*68:01 HLA-A*68:01 HLA-A*68:01 HLA-A*68:01 HLA-A*68:01
    ##  [51] HLA-A*68:01 HLA-A*68:01 HLA-A*68:01 HLA-A*68:01 HLA-A*68:01
    ##  [56] HLA-A*68:01 HLA-B*07:02 HLA-B*07:02 HLA-B*07:02 HLA-B*07:02
    ##  [61] HLA-B*07:02 HLA-B*07:02 HLA-B*07:02 HLA-B*07:02 HLA-B*07:02
    ##  [66] HLA-B*07:02 HLA-B*07:02 HLA-B*07:02 HLA-B*07:02 HLA-B*07:02
    ##  [71] HLA-B*07:02 HLA-B*07:02 HLA-B*07:02 HLA-B*07:02 HLA-B*07:02
    ##  [76] HLA-B*07:02 HLA-B*07:02 HLA-B*07:02 HLA-B*07:02 HLA-B*07:02
    ##  [81] HLA-B*07:02 HLA-B*07:02 HLA-B*07:02 HLA-B*07:02 HLA-B*35:01
    ##  [86] HLA-B*35:01 HLA-B*35:01 HLA-B*35:01 HLA-B*35:01 HLA-B*35:01
    ##  [91] HLA-B*35:01 HLA-B*35:01 HLA-B*35:01 HLA-B*35:01 HLA-B*35:01
    ##  [96] HLA-B*35:01 HLA-B*35:01 HLA-B*35:01 HLA-B*35:01 HLA-B*35:01
    ## [101] HLA-B*35:01 HLA-B*35:01 HLA-B*35:01 HLA-B*35:01 HLA-B*35:01
    ## [106] HLA-B*35:01 HLA-B*35:01 HLA-B*35:01 HLA-B*35:01 HLA-B*35:01
    ## [111] HLA-B*35:01 HLA-B*35:01
    ## Levels: HLA-A*02:01 HLA-A*68:01 HLA-B*07:02 HLA-B*35:01

``` r
rank <- sort(results$percentile_rank)
```

make a subset, and then sort it based on that. \*review order and sort: order will give you the indices to sort by, so really you wat to do sort(order())
=========================================================================================================================================================

``` r
20^9
```

    ## [1] 5.12e+11
