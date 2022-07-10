



# Importing microbiome data

This section demonstrates how to import taxonomic profiling data in R.


## Data access

*ADHD-associated changes in gut microbiota and brain in a mouse model*

Tengeler AC _et
al._ (2020) [**Gut microbiota from persons with
attention-deficit/hyperactivity disorder affects the brain in
mice**](https://doi.org/10.1186/s40168-020-00816-x). Microbiome
8:44.

In this study, mice are colonized with microbiota from participants
with ADHD (attention deficit hyperactivity disorder) and healthy
participants.  The aim of the study was to assess whether the mice
display ADHD behaviors after being inoculated with ADHD microbiota,
suggesting a role of the microbiome in ADHD pathology.

Download the data from
  [data](https://github.com/microbiome/course_2022_radboud/tree/main/data)
  subfolder.


## Importing microbiome data in R

**Import example data** by modifying the examples in the online book
section on [data exploration and
manipulation](https://microbiome.github.io/OMA/data-introduction.html#loading-experimental-microbiome-data).

The data files in our example are in _biom_ container, which is a
standard file format for microbiome data. Other file formats exist as
well, and import details vary by platform. Here, we import _biom_ data
files into a specific data container (structure) in R,
_TreeSummarizedExperiment_ (TreeSE) [Huang et
al. (2020)](https://f1000research.com/articles/9-1246).

The data container provides the basis for downstream data analysis in
the _miaverse_ data science framework.


<img src="https://raw.githubusercontent.com/FelixErnst/TreeSummarizedExperiment
/2293440c6e70ae4d6e978b6fdf2c42fdea7fb36a/vignettes/tse2.png" width="100%"/>

**Figure sources:** 

**Original article**
-   Huang R _et al_. (2021) [TreeSummarizedExperiment: a S4 class 
for data with hierarchical structure](https://doi.org/10.12688/
f1000research.26669.2). F1000Research 9:1246.

**Reference Sequence slot extension**
- Lahti L _et al_. (2020) [Upgrading the R/Bioconductor ecosystem for microbiome 
research](https://doi.org/10.7490/
f1000research.1118447.1) F1000Research 9:1464 (slides).




### Example solution

Let us first import the biom file into R / TreeSE container.


```r
# Load the mia R package
library(mia)

# Defining file paths
## Biom file (taxonomic profiles)
biom_file_path <- "data/Tengeler2020/Aggregated_humanization2.biom"

# Import the data into SummarizedExperiment container
se <- loadFromBiom(biom_file_path)

# Convert this data to TreeSE container (no direct importer exists)
tse <- as(se, "TreeSummarizedExperiment")
```


Check and clean up rowData (information on the taxonomic features).



```r
# Investigate the rowData of this data object
print(head(rowData(tse)))
```

```
## DataFrame with 6 rows and 6 columns
##             taxonomy1          taxonomy2           taxonomy3
##           <character>        <character>         <character>
## 1726470  "k__Bacteria   p__Bacteroidetes      c__Bacteroidia
## 1726471  "k__Bacteria   p__Bacteroidetes      c__Bacteroidia
## 17264731 "k__Bacteria   p__Bacteroidetes      c__Bacteroidia
## 17264726 "k__Bacteria   p__Bacteroidetes      c__Bacteroidia
## 1726472  "k__Bacteria p__Verrucomicrobia c__Verrucomicrobiae
## 17264724 "k__Bacteria   p__Bacteroidetes      c__Bacteroidia
##                      taxonomy4              taxonomy5           taxonomy6
##                    <character>            <character>         <character>
## 1726470       o__Bacteroidales      f__Bacteroidaceae     g__Bacteroides"
## 1726471       o__Bacteroidales      f__Bacteroidaceae     g__Bacteroides"
## 17264731      o__Bacteroidales  f__Porphyromonadaceae g__Parabacteroides"
## 17264726      o__Bacteroidales      f__Bacteroidaceae     g__Bacteroides"
## 1726472  o__Verrucomicrobiales f__Verrucomicrobiaceae     g__Akkermansia"
## 17264724      o__Bacteroidales      f__Bacteroidaceae     g__Bacteroides"
```

```r
# We notice that the rowData fields do not have descriptibve names.
# Hence, let us rename the columns in rowData
names(rowData(tse)) <- c("Kingdom", "Phylum", "Class", "Order", 
                        "Family", "Genus")

# We also notice that the taxa names are of form "c__Bacteroidia" etc.
# Goes through the whole DataFrame. Removes '.*[kpcofg]__' from strings, where [kpcofg] 
# is any character from listed ones, and .* any character.
rowdata_modified <- BiocParallel::bplapply(rowData(tse), 
                                           FUN = stringr::str_remove, 
                                           pattern = '.*[kpcofg]__')

# Genus level has additional '\"', so let's delete that also
rowdata_modified <- BiocParallel::bplapply(rowdata_modified, 
                                           FUN = stringr::str_remove, 
                                           pattern = '\"')

# rowdata_modified is a list, so convert this back to DataFrame format. 
# and assign the cleaned data back to the TSE rowData
rowData(tse) <- DataFrame(rowdata_modified)

# Recheck rowData after the modifications
print(head(rowData(tse)))
```

```
## DataFrame with 6 rows and 6 columns
##              Kingdom          Phylum            Class              Order
##          <character>     <character>      <character>        <character>
## 1726470     Bacteria   Bacteroidetes      Bacteroidia      Bacteroidales
## 1726471     Bacteria   Bacteroidetes      Bacteroidia      Bacteroidales
## 17264731    Bacteria   Bacteroidetes      Bacteroidia      Bacteroidales
## 17264726    Bacteria   Bacteroidetes      Bacteroidia      Bacteroidales
## 1726472     Bacteria Verrucomicrobia Verrucomicrobiae Verrucomicrobiales
## 17264724    Bacteria   Bacteroidetes      Bacteroidia      Bacteroidales
##                       Family           Genus
##                  <character>     <character>
## 1726470       Bacteroidaceae     Bacteroides
## 1726471       Bacteroidaceae     Bacteroides
## 17264731  Porphyromonadaceae Parabacteroides
## 17264726      Bacteroidaceae     Bacteroides
## 1726472  Verrucomicrobiaceae     Akkermansia
## 17264724      Bacteroidaceae     Bacteroides
```


Next let us add sample information (colData) to our TreeSE object.



```r
## Sample phenodata
sample_meta_file_path <- "data/Tengeler2020/Mapping_file_ADHD_aggregated.csv"

# Check what type of data it is
# read.table(sample_meta_file_path)

# It seems like a comma separated file and it does not include headers
# Let us read the file; note that sample names are in the first column
sample_meta <- read.table(sample_meta_file_path, sep=",", header=FALSE, row.names=1)

# Check the data
print(head(sample_meta))
```

```
##        V2       V3            V4   V5
## A110 ADHD Cohort_1 ADHD_Cohort_1 A110
## A12  ADHD Cohort_1 ADHD_Cohort_1  A12
## A15  ADHD Cohort_1 ADHD_Cohort_1  A15
## A19  ADHD Cohort_1 ADHD_Cohort_1  A19
## A21  ADHD Cohort_2 ADHD_Cohort_2  A21
## A23  ADHD Cohort_2 ADHD_Cohort_2  A23
```

```r
# Add headers for the columns (as they seem to be missing)
colnames(sample_meta) <- c("patient_status", "cohort",
                           "patient_status_vs_cohort", "sample_name")

# Add this sample data to colData of the taxonomic data object
# Note that the data must be given in a DataFrame format (required for our purposes)
colData(tse) <- DataFrame(sample_meta)

# Check the colData after modifications
print(head(colData(tse)))
```

```
## DataFrame with 6 rows and 4 columns
##      patient_status      cohort patient_status_vs_cohort sample_name
##         <character> <character>              <character> <character>
## A110           ADHD    Cohort_1            ADHD_Cohort_1        A110
## A12            ADHD    Cohort_1            ADHD_Cohort_1         A12
## A15            ADHD    Cohort_1            ADHD_Cohort_1         A15
## A19            ADHD    Cohort_1            ADHD_Cohort_1         A19
## A21            ADHD    Cohort_2            ADHD_Cohort_2         A21
## A23            ADHD    Cohort_2            ADHD_Cohort_2         A23
```


Add phylogenetic tree (rowTree).


```r
## Phylogenetic tree
tree_file_path <- "data/Tengeler2020/Data_humanization_phylo_aggregation.tre"

# Read the tree file
tree <- ape::read.tree(tree_file_path)

# Add tree to rowTree
rowTree(tse) <- tree
```

We can save the data object as follows.


```r
saveRDS(tse, file="tse.rds")
```


### Loading readily processed data

Alternatively, you can just load the already processed data set.

By using this readily processed data set we can skip the data import
step, and assume that the data has already been appropriately
preprocessed and available in the TreeSE container.



```r
tse <- readRDS("data/Tengeler2020/tse.rds")
```

In addition to our example data, further demonstration data sets are
available in the TreeSE data container: see [OMA demo data
sets](https://microbiome.github.io/OMA/containers.html#example-data).

