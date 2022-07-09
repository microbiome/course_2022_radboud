



# Importing microbiome data

This section demonstrates how to import microbiome profiling data in R.


## Data access

**Option 1**

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

**Option 2**

*Open data set of your own choice*, different options are listed in [OMA](https://microbiome.github.io/OMA/containers.html#example-data).


## Importing microbiome data in R

**Import example data** by modifying the examples in the online book
section on [data exploration and
manipulation](https://microbiome.github.io/OMA/data-introduction.html#loading-experimental-microbiome-data). The
data files in our example are in _biom_ format, which is a standard
file format for microbiome data. Other file formats exist as well, and
import details vary by platform.

Here, we import _biom_ data files into a specific data container (structure)
in R, _TreeSummarizedExperiment_ (TSE) [Huang et
al. (2020)](https://f1000research.com/articles/9-1246). This provides
the basis for downstream data analysis in the _miaverse_ data science
framework.

In this course, we focus on downstream analysis of taxonomic profiling
data, and assume that the data has already been appropriately
preprocessed and available in the TSE format. In addition to our
example data, further demonstration data sets are readily available in
the TSE format through
[microbiomeDataSets](https://bioconductor.org/packages/release/data/experiment/html/microbiomeDataSets.html).


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




## Importing data: example solutions



```r
# Defining file paths
biom_file_path <- "data/Tengeler2020/Aggregated_humanization2.biom"
sample_meta_file_path <- "data/Tengeler2020/Mapping_file_ADHD_aggregated.csv"
tree_file_path <- "data/Tengeler2020/Data_humanization_phylo_aggregation.tre"

library(mia)

# Imports the data
se <- loadFromBiom(biom_file_path)


names(rowData(se)) <- c("Kingdom", "Phylum", "Class", "Order", 
                        "Family", "Genus")

# Goes through the whole DataFrame. Removes '.*[kpcofg]__' from strings, where [kpcofg] 
# is any character from listed ones, and .* any character.
rowdata_modified <- BiocParallel::bplapply(rowData(se), 
                                           FUN = stringr::str_remove, 
                                           pattern = '.*[kpcofg]__')

# Genus level has additional '\"', so let's delete that also
rowdata_modified <- BiocParallel::bplapply(rowdata_modified, 
                                           FUN = stringr::str_remove, 
                                           pattern = '\"')

# rowdata_modified is a list, so it is converted back to DataFrame format. 
rowdata_modified <- DataFrame(rowdata_modified)

# And then assigned back to the SE object
rowData(se) <- rowdata_modified


# We use this to check what type of data it is
# read.table(sample_meta_file_path)

# It seems like a comma separated file and it does not include headers
# Let us read it and then convert from data.frame to DataFrame
# (required for our purposes)
sample_meta <- DataFrame(read.table(sample_meta_file_path, sep = ",", header = FALSE))

# Add sample names to rownames
rownames(sample_meta) <- sample_meta[,1]

# Delete column that included sample names
sample_meta[,1] <- NULL

# We can add headers
colnames(sample_meta) <- c("patient_status", "cohort", "patient_status_vs_cohort", "sample_name")

# Then it can be added to colData
colData(se) <- sample_meta

# Convert to tse format
tse <- as(se, "TreeSummarizedExperiment")

# Reads the tree file
tree <- ape::read.tree(tree_file_path)

# Add tree to rowTree
rowTree(tse) <- tree
```



### Loading the processed data

Alternatively, you can just load the already processed data set.


```r
tse <- readRDS("data/Tengeler2020/tse.rds")
```


