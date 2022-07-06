# Microbiome data science with R/Bioconductor

**Welcome to Radboud Summer School, July 2022**

<img src="https://user-images.githubusercontent.com/60338854/121848694-1072a480-ccf3-11eb-9af2-7fdefd8d1794.png" alt="ML4microbiome" width="50%"/>

Figure source: Moreno-Indias _et al_. (2021) [Statistical and Machine Learning Techniques in Human Microbiome Studies: Contemporary Challenges and Solutions](https://doi.org/10.3389/fmicb.2021.635781). _Frontiers in Microbiology_ 12:11. 


## Rendering the book

You can render the book locally in R with:

```{r serve}
bookdown::serve_book()
``` 

## The miaverse framework

The [_miaverse_](https://microbiome.github.io) (mia = **MI**crobiome **A**nalysis) is an
R/Bioconductor framework for microbiome data science. It aims to
extend the capabilities of another popular framework,
[phyloseq](https://joey711.github.io/phyloseq/).

The miaverse framework consists of an efficient data structure, an
associated package ecosystem, demonstration data sets, and open
documentation. These are explained in more detail in the online book
[Orchestrating Microbiome Analysis](https://microbiome.github.io/OMA).

- Landing page (html): [miaverse teaching material](https://microbiome.github.io/course_2021_radboud/)
- Source code (github): [miaverse teaching material](https://github.com/microbiome/course_2021_radboud)

The source code of this repository is fully reproducible and contains
the Rmd files with executable code. All files can be rendered at one
go by running the file [main.R](main.R). You can check the file for
details on how to clone the repository and convert it into a gitbook,
although this is not necessary for the training.

To setup the website, activate the gh-pages branch and wait for a few minutes.