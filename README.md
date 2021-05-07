
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-stable-ff69b4.svg)](https://www.tidyverse.org/lifecycle/#stable)
![last
commit](https://img.shields.io/github/last-commit/fjeanneret/multiSight.svg)
<!-- ![license](https://img.shields.io/github/license/fjeanneret/multiSight.svg) -->
<!-- badges: end -->

<img src="inst/app/www/multiSight_hexicon.png" align="right" alt="hexicon" />

# **multiSight**

:rocket: The purpose of this document is to help you become productive
as quickly as possible with the **multiSight** package.

  - The goal of **multiSight** is to handle multi-omics data and network
    inference in a easy-to-use R shiny package.

> You could use this tool with a graphical interface or only with script
> functions (see **Vignette** and *manual* for detailed examples).

# Installation

You can install the released version of **multiSight** from
[Bioconductor](https://www.bioconductor.org/) with:

``` r
#To install this package ensure you have BiocManager installed
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

#The following initializes usage of Bioc devel
BiocManager::install("multiSight")
```

# **multiSight** purpose

**multiSight** is an R package providing an user-friendly graphical
interface to analyze and explore your omic datasets in a **multi-omics**
manner by *DESeq2* (see **Biological Insights tab**), machine learning
methods with *biosigner* and *multi-block statistical analysis* (see
**Classification tab**) helped by *p-values pooling Stouffer’s* method.

Classification models are fitted to select few subsets of features,
using **biosigner** or **sPLS-DA methods**. *biosigner* provides one
model by omic block and one list of features named *biosignature*.
Nevertheless, sPLS-DA *biosignatures* are based on more features than
*biosigner*.

**Biosignatures** can be used:

  - To forecast phenotype (e.g. for diagnostic, histological subtyping);
  - To design ***Pathway*** and ***gene ontology*** **enrichment**
    (sPLS-DA biosignatures only);
  - To build ***Network inference***;
  - To find ***PubMed*** references to make assumptions easier and
    data-driven.

# :newspaper: App

**multiSight** enables you to get better biological insights for each
omic dataset helping by **four analytic modules** which content:

  - :memo: **Data input** & **results**;
  - :dart: **Classification** models building;
  - :books: **Biological databases** querying;
  - :seedling: **Network Inference** & **PubMed** querying.

> :point\_right: Run the application

``` r
run_app()
```

| :memo: Home                          | :dart: Classification                         | :books: Biological Insights                  | :seedling: Assumption                    |
| ------------------------------------ | --------------------------------------------- | -------------------------------------------- | ---------------------------------------- |
| ![home](inst/app/www/home/home1.png) | ![](inst/app/www/classification/classif1.png) | ![](inst/app/www/biologicalInsight/bio1.png) | ![](inst/app/www/networkInf/assump1.png) |

# What kind of data?

All types of omic data respecting input format is supported to build
**classification models**, **biosignatures** selection and **network
inference**.

  - Genomics;
  - Transcriptomics;
  - Proteomics;
  - Metabolomics;
  - Lipidomics;

> :point\_right: In fact all numeric matrices.

## Data inputs formats

You have to provide two types of data: **numeric matrices** and
**classes vector** as csv tables for all **same samples**.

### Omic data 1

|      | SIGIRR | MAOA | MANSC1 |   |
| ---- | ------ | ---- | ------ | - |
| AOFJ | 0      | 150  | 1004   | … |
| A13E | 34     | 0    | 0      |   |
|      |        | …    |        |   |

### Omic data 2

|      | ENSG00000139618 | ENSG00000226023 | ENSG00000198695 |   |
| ---- | --------------- | --------------- | --------------- | - |
| AOFJ | 25              | 42              | 423             | … |
| A13E | 0               | 154             | 4900            |   |
|      |                 | …               |                 |   |

… :point\_right: unlimited number of omic datasets.

### Omic data n

|      | 4292 | 5254 | 7432 |   |
| ---- | ---- | ---- | ---- | - |
| AOFJ | 25   | 42   | 423  | … |
| A13E | 0    | 154  | 4900 |   |
|      |      | …    |      |   |

### Omic classes

|      | Y     |
| ---- | ----- |
| AOFJ | condA |
| A13E | condB |
|      | …     |

# :dart: Classification tab

Two types of models have been implemented so far to answer different
questions: [**biosigner**](https://doi.org/10.3389/fmolb.2016.00026) &
[**sPLS-DA (DIABLO)**](https://doi.org/10.1093/bioinformatics/bty1054) .

  - To determine *small biosignatures* - biosigner.
  - To build *classification models* in a *multi-omics* way - DIABLO.
  - To select relevant biological *features* to *enrich* - DIABLO.

| Features selected                             | Performances                                  |
| --------------------------------------------- | --------------------------------------------- |
| ![](inst/app/www/classification/classif2.png) | ![](inst/app/www/classification/classif3.png) |

# :books: Biological insights tab

**Biological Insight** tab is dedicated to give biological sense to your
data.

  - You could process ***2 analysis in 2 clicks***: both *DESeq2* and
    *DIABLO features* ORAs for functional enrichment.

## Biological Annotation Databases

**multiSight** uses so far **several databases** to provide a large
panel of **enrichment analysis**, automatically after few clicks:

**Pathways** and **Gene Ontology** databases are implemented, helped by
**clusterProfiler** and **reactomePA** R Bioconductor packages.

  - Kegg;
  - Reactome;
  - wikiPathways;
  - Molecular Function (GO)
  - Cellular Component (GO)
  - Biological Process (GO)

## Visualizations

Two types of result visualization are given:

  - Classical **Enrichment tables** for each omic and each database
    (e.g.  Pathways id, p-value, padjust columns).
  - And, when more than one omic enriched: a *Multi-omics table* and a
    *multi-omics enrichment map* for **DESeq2** and *DIABLO selected
    features*.

| DESeq2 & DIABLO features                     | Enrichment tables                            | Enrichment Map                               |
| -------------------------------------------- | -------------------------------------------- | -------------------------------------------- |
| ![](inst/app/www/biologicalInsight/bio2.png) | ![](inst/app/www/biologicalInsight/bio3.png) | ![](inst/app/www/biologicalInsight/bio4.png) |

# :seedling: Assumption tab

> :point\_right: Some clicks (from 4 to number of PubMed queries)

**Assumption tab** aims to help biological hypothesis making by *network
inference* from feature relationship values (e.g correlation, partial
correlation) and by a *PubMed module*.

You can find both functions:

  - To compute *_network inference_* and to reveal feature
    relationships.
  - To get *_PubMed articles_* based on your personalized query without
    leaving app.

| Network Inference                        | PubMed query                             |
| ---------------------------------------- | ---------------------------------------- |
| ![](inst/app/www/networkInf/assump2.png) | ![](inst/app/www/networkInf/assump3.png) |

# :checkered\_flag: Results

You could retrieve different results computed by multiSight in Home tab
by:

  - Automatic report with all results in **HTML** and **.doc**
    documents.
  - **.RData** with all results obtained by the graphical application.

Note that tables could be downloaded in a separated way in relative
tabs.

> **MODELS**: classification models you can use on future data.

> **DESeq2**: differential expression analysis tables.

> **BIOSIGNATURES**: DESeq2 tables thresholding and DIABLO multi-omics
> features selection method

> **Functional ENRICHMENTS**: 6 databases functional enrichment for all
> omic datasets you provide enriched by Stouffer’s pooling p-value
> method giving a **multi-omics enrichmentt able** easily to discuss.

> **NETWORKS**: network inference analysis with all features selected
> from all omic datasets according to DESeq2 tables thresholding or
> multi-omics feature selection (correlation, partial correlation,
> mutual information).

> **BIBLIOGRAPHY** : a subset of PubMed articles relative to relations
> you choose in network inference tab.
