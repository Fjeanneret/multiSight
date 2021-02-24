
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
![last
commit](https://img.shields.io/github/last-commit/fjeanneret/MultiSight.svg)
![license](https://img.shields.io/github/license/fjeanneret/MultiSight.svg)
<!-- badges: end -->

<img src="inst/app/www/multiSight_hexicon.png" align="right" alt="hexicon" />

# **multiSight**

:rocket: This document is built to be efficient as quickly as possible
with **multiSight**.

  - The goal of **multiSight** is to handle multi-omic data and network
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

# To get up-to-date package
install.packages("remotes")
remotes::install_github("Fjeanneret/multiSight")
```

# What is **multiSight** ?

**multiSight** is a R package providing an user-friendly graphical
interface to analyze and explore your omic data sets in a **multi-omic**
manner by *DESEQ2* (see **Biological Insights tab**), machine learning
methods with *biosigner* and *multi-block statistical analysis* (see
**Classification tab**) helped by *p-values pooling Stouffer’s* method.

For each omic data set you furnish, it provides ***classification
models*** with *feature selection* you can use as biosignature:

  - To forecast phenotypes (e.g. to diagnostic tasks, histological
    subtyping)
  - To design ***Pathways*** and ***gene ontology*** **enrichments**
  - To build ***Network inference*** linked to ***PubMed*** querying to
    make assumptions easier and data-driven.

# :newspaper: App

**multiSight** enables you to get better biological insights for each
omic dataset helping by **four analytic modules** which content:

  - :memo: **Data input** & **results**;
  - :dart: **Classification** models building;
  - :books: **Biological databases** querying;
  - :seedling: **Network Inference** & **Pubmed** querying.

> :point\_right: Run the application

``` r
# run_app()
```

| :memo: Home                          | :dart: Classification                         | :books: Biological Insights                  | :seedling: Assumption                    |
| ------------------------------------ | --------------------------------------------- | -------------------------------------------- | ---------------------------------------- |
| ![home](inst/app/www/home/home1.png) | ![](inst/app/www/classification/classif1.png) | ![](inst/app/www/biologicalInsight/bio1.png) | ![](inst/app/www/networkInf/assump1.png) |

# What kind of data ?

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

| Omic data 1 |        |      |        |   |
| ----------- | ------ | ---- | ------ | - |
|             | SIGIRR | MAOA | MANSC1 |   |
| AOFJ        | 0      | 150  | 1004   | … |
| A13E        | 34     | 0    | 0      |   |
|             |        | …    |        |   |

| Omic data 2 |                 |                 |                 |   |
| ----------- | --------------- | --------------- | --------------- | - |
|             | ENSG00000139618 | ENSG00000226023 | ENSG00000198695 |   |
| AOFJ        | 25              | 42              | 423             | … |
| A13E        | 0               | 154             | 4900            |   |
|             |                 | …               |                 |   |

… :point\_right: unlimited number of omic data sets.

| Omic data n |      |      |      |   |
| ----------- | ---- | ---- | ---- | - |
|             | 4292 | 5254 | 7432 |   |
| AOFJ        | 25   | 42   | 423  | … |
| A13E        | 0    | 154  | 4900 |   |
|             |      | …    |      |   |

| Omic classes |       |
| ------------ | ----- |
|              | Y     |
| AOFJ         | condA |
| A13E         | condB |
|              | …     |

# :dart: Classification tab

Two types of models have been implemented so far to answer different
questions: **Biosigner** & **Diablo**.

  - To determine *small biosignatures* - Biosigner.
  - To build *classification models* in a *multi-omic* way - Diablo.
  - To select relevant biological *features* to *enrich* - Diablo.

| Features selected                             | Performances                                  |
| --------------------------------------------- | --------------------------------------------- |
| ![](inst/app/www/classification/classif2.png) | ![](inst/app/www/classification/classif3.png) |

# :books: Biological insights tab

**Biological Insight** tab is dedicated to give biological sense to your
data.

  - You could process ***2 analysis in 2 clicks***: *DESEQ2* and *Diablo
    features* enrichments.

## Biological Annotation Databases

**multiSight** uses so far **several databases** to provide large panel
of **enrichment analysis** automatically after few clicks:

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

  - Classical **Enrichment table** for each omic and each database
    (e.g.  Pathways id, p-value, padjust columns).
  - And, when more than one omic enriched: *Multi-omic table* and
    *multi-omic enrichment map* for **DESEQ2** and *Diablo selected
    features*

| DESEQ2 & Diablo features                     | Enrichment tables                            | Enrichment Map                               |
| -------------------------------------------- | -------------------------------------------- | -------------------------------------------- |
| ![](inst/app/www/biologicalInsight/bio2.png) | ![](inst/app/www/biologicalInsight/bio3.png) | ![](inst/app/www/biologicalInsight/bio4.png) |

# :seedling: Assumption tab

> :point\_right: Some clicks (from 4 to number of pubmed queries)

**Assumption tab** aims to help biological hypothesis making by *network
inference* with feature relation values (e.g correlation, partial
correlation) and *PubMed module* linked to relation tables.

Tools:

  - To compute *_network inference_* and reveal feature relationships.
  - To get *_Pubmed articles_* based on your personalized query without
    leaving app.

| Network Inference                        | PubMed query                             |
| ---------------------------------------- | ---------------------------------------- |
| ![](inst/app/www/networkInf/assump2.png) | ![](inst/app/www/networkInf/assump3.png) |

# :checkered\_flag: Results

You could retrieve different results computed by multiSight in Home tab
by:

  - Automatic report with all results in **HTML** and **.doc**
    documents.
  - **.RData** with all results obtained by all modules.

Note that tables could be download in a separated way in relative tabs.

> **MODELS**: Classification models you can use on future data.

> **DESEQ2** Differential expression analysis tables.

> **BIOSIGNATURES**: DESEQ2’s tables thresholding and diablo multi-omic
> features selection method

> **Functional ENRICHMENTS**: 6 databases functional enrichment for all
> omic data sets you provide enriched by Stouffer’s pooling p-value
> method giving a **multi-omic enrichmentt able** easily to discuss.

> **NETWORKS**: network inference with all features selected from all
> omic datasets according to DESEQ2 tables thresholding or multi-omic
> feature selection (correlation, partial correlation, mutual
> information).

> **BIBLIOGRAPHY** : Subset of pubMed articles relative to relations you
> choose in network inference tab.
