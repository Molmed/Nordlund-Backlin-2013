Genome-wide signatures of differential DNA methylation in pediatric acute lymphoblastic leukemia
======================

[Nordlund J](http://scholar.google.se/citations?user=ZztFeTEAAAAJ&hl=sv&oi=ao), [Bäcklin C](http://stackoverflow.com/users/840460/backlin), Wahlberg P, Busche S, Berglund EC, Eloranta M-L, Flaegstad T, Forestier E,  Frost B-M, Harila-Saari A, Heyman M, Jónsson OG, Larsson R, Palle J,  Rönnblom L, Schmiegelow K, Sinnett D, Söderhäll S, Pastinen T, Gustafsson MG, Lönnerholm G & Syvänen AC

Original publication: [Link to publication to be added here](#)

This code was produced by the groups of [Molecular Medicine](http://www.molmed.medsci.uu.se/) and [Cancer Pharmacology and Computational Medicine](http://www.medsci.uu.se/research/Cancer/Cancer+Pharmacology+and+Computational+Medicine/) at the [Department of Medical Sciences](http://www.medsci.uu.se) at [Uppsala University](http://www.uu.se).

System requirements
-------------------
R version 2.15.3 or later.

Parallelization using the [`foreach`](http://cran.r-project.org/web/packages/foreach/index.html) and [`doSNOW`](http://cran.r-project.org/web/packages/doSNOW/index.html) packages is recommended, but not required.

The package `analyse450k` called in [`03_analyse.R`](https://github.com/Molmed/Nordlund-Backlin-2013/blob/master/03_analyse.R) is only used for in house data management, and will be replaced once the dataset is published on [GEO](http://www.ncbi.nlm.nih.gov/geo/).

Figures were annotated using the [`biomaRt`](http://www.bioconductor.org/packages/2.12/bioc/html/biomaRt.html) and [`GenomeGraphs`](http://www.bioconductor.org/packages/2.12/bioc/html/GenomeGraphs.html) packages.

