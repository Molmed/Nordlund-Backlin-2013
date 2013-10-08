Genome-wide signatures of differential DNA methylation in pediatric acute lymphoblastic leukemia
======================

[Nordlund J](http://scholar.google.se/citations?user=ZztFeTEAAAAJ&hl=sv&oi=ao), [Bäcklin C](http://stackoverflow.com/users/840460/backlin), Wahlberg P, Busche S, Berglund EC, Eloranta M-L, Flaegstad T, Forestier E,  Frost B-M, Harila-Saari A, Heyman M, Jónsson OG, Larsson R, Palle J,  Rönnblom L, Schmiegelow K, Sinnett D, Söderhäll S, Pastinen T, Gustafsson MG, Lönnerholm G & Syvänen AC

Original publication: [Genome Biology 2013, 14:r105](http://genomebiology.com/2013/14/9/r105/abstract)

This code was produced by the groups of [Molecular Medicine](http://www.molmed.medsci.uu.se/) and [Cancer Pharmacology and Computational Medicine](http://www.medsci.uu.se/research/Cancer/Cancer+Pharmacology+and+Computational+Medicine/) at the [Department of Medical Sciences](http://www.medsci.uu.se) at [Uppsala University](http://www.uu.se).

System requirements
-------------------
R version 2.15.3 or later. Packages [`plyr`](http://plyr.had.co.nz/) and [`GEOquery`](http://www.bioconductor.org/packages/2.12/bioc/html/GEOquery.html) are required, but are installed automatically. 

Parallelization using the [`foreach`](http://cran.r-project.org/web/packages/foreach/index.html) and [`doSNOW`](http://cran.r-project.org/web/packages/doSNOW/index.html) packages is recommended, but not required.

Figures were annotated using the [`biomaRt`](http://www.bioconductor.org/packages/2.12/bioc/html/biomaRt.html) and [`GenomeGraphs`](http://www.bioconductor.org/packages/2.12/bioc/html/GenomeGraphs.html) packages. However, code for producing the figures is not included.

Instructions
------------
Download all scripts in this repo to a new directory. The most convenient way to do this on a linux/unix system is to clone the whole repo. Then run the files `setup.R` and `analyse.R`.

    git clone git@github.com:Molmed/Nordlund-Backlin-2013.git
    cd Nordlund-Backlin-2013
    R -f setup.R
    R -f analyse.R

`setup.R` will download all data from [GEO](www.ncbi.nlm.nih.gov/geo/), prepare it for use in R and store it in a new subfolder called `data`. `analyse.R` will run the analyses, produce the results and save them in a new subfolder called `results`. Plots are not produced.
