# Getting started {#start}

## Checklist (before the course)

### CSC Notebook

We will provide a temporary access to a cloud computing environment
that readily contains the available software packages. Instructions to access 
the environment will be sent to the registered participants. 

1. Read the [instructions](https://docs.csc.fi/cloud/csc_notebooks/guide_for_students/)
2. Go to the [CSC notebook frontpage](https://notebooks-beta.rahtiapp.fi/welcome)
3. Login
    a. Haka login
        * If you have a Finnish university account, you should be able to login with Haka
        1. Press **Login** button from the frontpage
        2. Press **Haka** button
        3. Select right organization
        4. Enter login information
    b. CSC login
        * You can create a CSC account by following the [instructions](https://research.csc.fi/accounts-and-projects)
        1. Press **Login** button from the frontpage
        2. Press **CSC** button
        3. Enter login information
    c. Special login
        * For those who cannot login with Haka or CSC account
        1. Contact Tuomas by email (first.v.last@utu.fi) if you are not able to login
        2. We give you a guest account
        3. Press **Special Login** button from the frontpage (below the **Login** button)
        4. Enter login information (username goes to **email** slot)
4. Join workspace
    a. Press **Join workspace** button (Top right corner)
    b. Enter the **Join Code** (Check your email)
5. Start session
    a. Press **ON** button
6. You can save files to **my-work** directory. They are kept stored even when the session is closed. **shared** folder is shared with all participants.

    
### (Your own computer)

Setting up the system on your own computer is not required for the
course but it can be useful for later use. The required software:

* [R (version >4.1.0)](https://www.r-project.org/) 

* [RStudio](https://www.rstudio.com/products/rstudio/download/);
  choose "Rstudio Desktop" to download the latest version. Optional
  but preferred. For further details, check the [Rstudio home
  page](https://www.rstudio.com/).

* Install and load the required R packages (see Section \@ref(packages))

* After a successful installation you can start with the
  case study examples in this training material


## Support and resources

 * We recommend to have a look at the additional reading tips and try out online material listed in Section \@ref(material).

**You can run the workflows by simply copy-pasting the examples.** For
further, advanced material, you can test and modify further examples
from the online book, and apply these techniques to your own data.

 * Online support on installation and other matters, join us at [Gitter](https://gitter.im/microbiome/miaverse?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)


## Installing and loading the required R packages {#packages}

Note that the CSC/RStudio environment has readily installed setup. You
may need the examples from this subsection if you are installing the
environment on your own computer. If you need to add new packages, you
can modify the examples below.

This section shows how to install and load all required packages into
the R session, if needed. Only uninstalled packages are installed.


```r
# List of packages that we need from cran and bioc 
cran_pkg <- c("BiocManager", "bookdown", "dplyr", "ecodist", "ggplot2", 
              "gridExtra", "kableExtra",  "knitr", "scales", "vegan", "matrixStats")
bioc_pkg <- c("yulab.utils","ggtree","ANCOMBC", "ape", "DESeq2", "DirichletMultinomial", "mia", "miaViz")

# Get those packages that are already installed
cran_pkg_already_installed <- cran_pkg[ cran_pkg %in% installed.packages() ]
bioc_pkg_already_installed <- bioc_pkg[ bioc_pkg %in% installed.packages() ]

# Get those packages that need to be installed
cran_pkg_to_be_installed <- setdiff(cran_pkg, cran_pkg_already_installed)
bioc_pkg_to_be_installed <- setdiff(bioc_pkg, bioc_pkg_already_installed)

# Reorders bioc packages, so that mia and miaViz are first
bioc_pkg <- c(bioc_pkg[ bioc_pkg %in% c("mia", "miaViz") ], 
              bioc_pkg[ !bioc_pkg %in% c("mia", "miaViz") ] ) 

# Combine to one vector
packages <- c(bioc_pkg, cran_pkg)
packages_to_install <- c( bioc_pkg_to_be_installed, cran_pkg_to_be_installed )
```


```r
# If there are packages that need to be installed, install them 
if( length(packages_to_install) ) {
   BiocManager::install(packages_to_install)
}
```

Now all required packages are installed, so let's load them into the session.
Some function names occur in multiple packages. That is why miaverse's packages
mia and miaViz are prioritized. Packages that are loaded first have higher priority.


```r
# Loading all packages into session. Returns true if package was successfully loaded.
loaded <- sapply(packages, require, character.only = TRUE)
as.data.frame(loaded)
```

```
##                      loaded
## mia                    TRUE
## miaViz                 TRUE
## yulab.utils            TRUE
## ggtree                 TRUE
## ANCOMBC                TRUE
## ape                    TRUE
## DESeq2                 TRUE
## DirichletMultinomial   TRUE
## BiocManager            TRUE
## bookdown               TRUE
## dplyr                  TRUE
## ecodist                TRUE
## ggplot2                TRUE
## gridExtra              TRUE
## kableExtra             TRUE
## knitr                  TRUE
## scales                 TRUE
## vegan                  TRUE
## matrixStats            TRUE
```

