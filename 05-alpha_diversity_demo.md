# Alpha diversity demo



## Alpha diversity estimation


First let`s load the required packages and data set



```r
library(mia)
```

```
## Loading required package: SummarizedExperiment
```

```
## Loading required package: MatrixGenerics
```

```
## Loading required package: matrixStats
```

```
## 
## Attaching package: 'MatrixGenerics'
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
##     colWeightedMeans, colWeightedMedians, colWeightedSds,
##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
##     rowWeightedSds, rowWeightedVars
```

```
## Loading required package: GenomicRanges
```

```
## Loading required package: stats4
```

```
## Loading required package: BiocGenerics
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
##     union, unique, unsplit, which.max, which.min
```

```
## Loading required package: S4Vectors
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following objects are masked from 'package:base':
## 
##     expand.grid, I, unname
```

```
## Loading required package: IRanges
```

```
## Loading required package: GenomeInfoDb
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## 
## Attaching package: 'Biobase'
```

```
## The following object is masked from 'package:MatrixGenerics':
## 
##     rowMedians
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     anyMissing, rowMedians
```

```
## Loading required package: SingleCellExperiment
```

```
## Loading required package: TreeSummarizedExperiment
```

```
## Loading required package: Biostrings
```

```
## Loading required package: XVector
```

```
## 
## Attaching package: 'Biostrings'
```

```
## The following object is masked from 'package:base':
## 
##     strsplit
```

```
## Loading required package: MultiAssayExperiment
```

```r
library(miaViz)
```

```
## Loading required package: ggplot2
```

```
## Loading required package: ggraph
```

```r
library(tidyverse)
```

```
## -- Attaching packages --------------------------------------- tidyverse 1.3.1 --
```

```
## v tibble  3.1.7     v dplyr   1.0.9
## v tidyr   1.2.0     v stringr 1.4.0
## v readr   2.1.2     v forcats 0.5.1
## v purrr   0.3.4
```

```
## -- Conflicts ------------------------------------------ tidyverse_conflicts() --
## x dplyr::collapse()   masks Biostrings::collapse(), IRanges::collapse()
## x dplyr::combine()    masks Biobase::combine(), BiocGenerics::combine()
## x purrr::compact()    masks XVector::compact()
## x dplyr::count()      masks matrixStats::count()
## x dplyr::desc()       masks IRanges::desc()
## x tidyr::expand()     masks S4Vectors::expand()
## x dplyr::filter()     masks stats::filter()
## x dplyr::first()      masks S4Vectors::first()
## x dplyr::lag()        masks stats::lag()
## x ggplot2::Position() masks BiocGenerics::Position(), base::Position()
## x purrr::reduce()     masks GenomicRanges::reduce(), IRanges::reduce()
## x dplyr::rename()     masks S4Vectors::rename()
## x dplyr::slice()      masks XVector::slice(), IRanges::slice()
```

```r
# library(vegan)

tse <- read_rds("data/Tengeler2020/tse.rds")

tse
```

```
## class: TreeSummarizedExperiment 
## dim: 151 27 
## metadata(0):
## assays(1): counts
## rownames(151): 1726470 1726471 ... 17264756 17264757
## rowData names(6): Kingdom Phylum ... Family Genus
## colnames(27): A110 A12 ... A35 A38
## colData names(4): patient_status cohort patient_status_vs_cohort
##   sample_name
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
## rowLinks: a LinkDataFrame (151 rows)
## rowTree: 1 phylo tree(s) (151 leaves)
## colLinks: NULL
## colTree: NULL
```
Then let's estimate multiple diversity indices.


```r
?estimateDiversity

tse <- estimateDiversity(tse, 
                              index = c("shannon","gini_simpson","faith"),
                              name = c("shannon","gini_simpson","faith"))
head(colData(tse))
```

```
## DataFrame with 6 rows and 7 columns
##      patient_status      cohort patient_status_vs_cohort sample_name   shannon
##         <character> <character>              <character> <character> <numeric>
## A110           ADHD    Cohort_1            ADHD_Cohort_1        A110   1.76541
## A12            ADHD    Cohort_1            ADHD_Cohort_1         A12   2.71644
## A15            ADHD    Cohort_1            ADHD_Cohort_1         A15   3.17810
## A19            ADHD    Cohort_1            ADHD_Cohort_1         A19   2.89199
## A21            ADHD    Cohort_2            ADHD_Cohort_2         A21   2.84198
## A23            ADHD    Cohort_2            ADHD_Cohort_2         A23   2.79794
##      gini_simpson     faith
##         <numeric> <numeric>
## A110     0.669537   7.39224
## A12      0.871176   6.29378
## A15      0.930561   6.60608
## A19      0.899210   6.79708
## A21      0.885042   6.65110
## A23      0.859813   5.96246
```

We can see that the variables are included in the data.
Similarly, let's calculate richness indices.


```r
tse <- estimateRichness(tse, 
                              index = c("chao1","observed"))
head(colData(tse))
```

```
## DataFrame with 6 rows and 10 columns
##      patient_status      cohort patient_status_vs_cohort sample_name   shannon
##         <character> <character>              <character> <character> <numeric>
## A110           ADHD    Cohort_1            ADHD_Cohort_1        A110   1.76541
## A12            ADHD    Cohort_1            ADHD_Cohort_1         A12   2.71644
## A15            ADHD    Cohort_1            ADHD_Cohort_1         A15   3.17810
## A19            ADHD    Cohort_1            ADHD_Cohort_1         A19   2.89199
## A21            ADHD    Cohort_2            ADHD_Cohort_2         A21   2.84198
## A23            ADHD    Cohort_2            ADHD_Cohort_2         A23   2.79794
##      gini_simpson     faith     chao1  chao1_se  observed
##         <numeric> <numeric> <numeric> <numeric> <numeric>
## A110     0.669537   7.39224        68  0.000000        68
## A12      0.871176   6.29378        51  0.000000        51
## A15      0.930561   6.60608        68  0.000000        68
## A19      0.899210   6.79708        62  0.000000        62
## A21      0.885042   6.65110        58  0.000000        58
## A23      0.859813   5.96246        61  0.247942        61
```

## Visualizing alpha diversity 

We can plot the distributions of individual indices:


```r
#individual plot
p <- as_tibble(colData(tse)) %>% 
  ggplot(aes(shannon)) +
  geom_histogram() 

print(p)
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](05-alpha_diversity_demo_files/figure-latex/distributions-1.pdf)<!-- --> 

```r
#multiple plots

p <- as_tibble(colData(tse)) %>% 
  pivot_longer(cols = c("shannon","gini_simpson","faith","chao1","observed"), names_to = "index", values_to = "alpha") %>% 
  ggplot(aes(alpha)) +
  geom_histogram() +
  facet_wrap(vars(index), scales = "free")


print(p)
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](05-alpha_diversity_demo_files/figure-latex/distributions-2.pdf)<!-- --> 

and the correlation between indices:


```r
p <- as_tibble(colData(tse)) %>% 
  pivot_longer(cols = c("shannon","gini_simpson","faith","chao1","observed"), names_to = "index", values_to = "alpha") %>% 
  full_join(.,., by = "sample_name") %>% 
  ggplot( aes(x = alpha.x, y = alpha.y)) + 
  geom_point() +
  geom_smooth() +
  facet_wrap(index.x ~ index.y, scales = "free")

print(p)
```

```
## `geom_smooth()` using method = 'loess' and formula 'y ~ x'
```

![](05-alpha_diversity_demo_files/figure-latex/scatterlots-1.pdf)<!-- --> 

## Comparing alpha diversity 

It is often interesting to look for any group differences:



```r
p <- as_tibble(colData(tse)) %>% 
  pivot_longer(cols = c("shannon","gini_simpson","faith","chao1","observed"), names_to = "index", values_to = "alpha") %>% 
  ggplot( aes(x = patient_status, y = alpha)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha =0.5) +
  facet_wrap(vars(index), scales = "free")

print(p)
```

![](05-alpha_diversity_demo_files/figure-latex/boxplots-1.pdf)<!-- --> 

Moreover, we can test the group differences by parametric or non-parametric tests:


```r
df1 <- as_tibble(colData(tse)) %>% 
  pivot_longer(cols = c("faith","chao1","observed"), names_to = "index", values_to = "alpha") %>% 
  group_by(index) %>% 
  nest() %>% 
  mutate(test_pval = map_dbl(data, ~ t.test(alpha ~ patient_status, data = .x)$p.value)) %>% 
  mutate(test = "ttest" ) 

df2 <- as_tibble(colData(tse)) %>% 
  pivot_longer(cols = c("shannon","gini_simpson"), names_to = "index", values_to = "alpha") %>% 
  group_by(index) %>% 
  nest() %>% 
  mutate(test_pval = map_dbl(data, ~ wilcox.test(alpha ~ patient_status, data = .x)$p.value))%>% 
  mutate(test = "wilcoxon" ) 

df <- rbind(df1,df2) %>% select(-data) %>% arrange(test_pval) %>% ungroup()

df
```

```
## # A tibble: 5 x 3
##   index        test_pval test    
##   <chr>            <dbl> <chr>   
## 1 shannon          0.488 wilcoxon
## 2 gini_simpson     0.685 wilcoxon
## 3 chao1            0.856 ttest   
## 4 observed         0.900 ttest   
## 5 faith            0.983 ttest
```
End of the demo.


## Exercises

Do "Alpha diversity basics" from the [exercises](https://microbiome.github.io/OMA/exercises.html).
