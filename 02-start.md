# Getting started {#start}

## Checklist (before the course)

### Computer setup and installations

Setting up the system on your own computer is required to follow the full course and will be useful for later use if you intent to analyze microbiome data on your computer. The required software:

* [R (version >4.2.0)](https://www.r-project.org/) 

* [RStudio](https://www.rstudio.com/products/rstudio/download/);
  choose "Rstudio Desktop" to download the latest version. Optional
  but preferred. For further details, check the [Rstudio home
  page](https://www.rstudio.com/).

* Install and load the required R packages (see Section \@ref(packages))

* After a successful installation you can start with the
  case study examples in this training material


## Support and resources

 * We recommend to have a look at the additional reading tips and try out online material listed in Section \@ref(material).

 * **You can run the workflows by simply copy-pasting the examples.**
For further, advanced material, you can test and modify further
examples from the online book, and apply these techniques to your own
data.

 * Online support on installation and other matters, join us at [Gitter](https://gitter.im/microbiome/miaverse?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)


## Installing and loading the required R packages {#packages}

You may need the examples from this subsection if you are installing
the environment on your own computer. If you need to add new packages,
you can modify the examples below.

This section shows how to install and load all required packages into
the R session, if needed. Only uninstalled packages are installed.

Download the file [pkgs.csv](pkgs.csv). This contains the list of
packages that we recommend to preinstall. This can be done with the
following code.


```r
# List of packages that we need
pkg <- read.csv("pkgs.csv")[,1]

# List packages that are already installed
pkg_already_installed <- pkg[ pkg %in% installed.packages() ]

# List remaining packages that need to be installed
packages_to_install <- setdiff(pkg, pkg_already_installed)
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
loaded <- sapply(pkg, require, character.only = TRUE)
as.data.frame(loaded)
```



