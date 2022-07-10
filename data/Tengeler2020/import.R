# Load the mia R package
library(mia)

# Defining file paths
## Biom file (taxonomic profiles)
biom_file_path <- "data/Tengeler2020/Aggregated_humanization2.biom"

# Import the data into SummarizedExperiment container
se <- loadFromBiom(biom_file_path)

# Convert this data to TreeSE container (no direct importer exists)
tse <- as(se, "TreeSummarizedExperiment")

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

## Sample phenodata
sample_meta_file_path <- "data/Tengeler2020/Mapping_file_ADHD_aggregated.csv"

# It seems like a comma separated file and it does not include headers
# Let us read the file; note that sample names are in the first column
sample_meta <- read.table(sample_meta_file_path, sep=",", header=FALSE, row.names=1)

# Add headers for the columns (as they seem to be missing)
colnames(sample_meta) <- c("patient_status", "cohort",
                           "patient_status_vs_cohort", "sample_name")

# Add this sample data to colData of the taxonomic data object
# Note that the data must be given in a DataFrame format (required for our purposes)
colData(tse) <- DataFrame(sample_meta)

## Phylogenetic tree
tree_file_path <- "data/Tengeler2020/Data_humanization_phylo_aggregation.tre"

# Read the tree file
tree <- ape::read.tree(tree_file_path)

# Add tree to rowTree
rowTree(tse) <- tree

saveRDS(tse, file="tse.rds")

