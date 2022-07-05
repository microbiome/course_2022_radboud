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

This training material walks you through an example workflow that
shows the standard steps of taxonomic data analysis covering data
access, exploration, analysis, visualization and reporoducible
reporting. **You can run the workflow by simply copy-pasting the
examples.** For advanced material, you can test and modify further
examples from the [OMA book](https://microbiome.github.io/OMA), or try
to apply the techniques to your own data.




# Learning goals

This course will teach the basics of biomedical data analysis with R/Bioconductor, a
popular open source environment for scientific data analysis. The participants get an
overview of the reproducible data analysis workflow in modern multi-omics, with a focus
on recent examples from published microbiome studies. After the course you will know
how to approach new tasks in biomedical data analysis by utilizing available
documentation and R tools. The teaching format allows adaptations according to the
student’s learning speed.

The teaching will follow open online documentation created by the course teachers,
extending the online book Orchestrating Microbiome Analysis
(https://microbiome.github.io/OMA). The training material walks you through the
standard steps of biomedical data analysis covering data access, exploration, analysis,
visualization, reproducible reporting, and best practices in open science. The openly
licensed teaching material will be available online during and after the course, following
national recommendations on open teaching materials.


**Target audience** The course is primarily designed for advanced MSc and PhD students, Postdocs, and
biomedical researchers who wish to learn new skills in scientific programming and
biomedical data analysis. Academic students and researchers from Finland and abroad are
welcome and encouraged to apply. Priority will be given for local students from Radboud.

**Venue** Radboud University. July 11-15, 2022.

## Acknowledgments

**Citation** "Microbiome data science with R/Bioconductor (2022). Tuomas Borman, Henrik Eckermann, Anna Aatsinki, Leo Lahti. URL: https://microbiome.github.io/course_2022_radboud/"

**Contact**
- [Leo Lahti](http://datascience.utu.fi), University of Turku
- [mia Collective](https://microbiome.github.io)

**License** All material is released under the open [CC BY-NC-SA 3.0 License](LICENSE).

- Landing page (html): [miaverse teaching material](https://microbiome.github.io/course_2021_radboud/)
- Source code (github): [miaverse teaching material](https://github.com/microbiome/course_2021_radboud)

The source code of this repository is fully reproducible and contains
the Rmd files with executable code. All files can be rendered at one
go by running the file [main.R](main.R). You can check the file for
details on how to clone the repository and convert it into a gitbook,
although this is not necessary for the training.
