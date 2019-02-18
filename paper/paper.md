---
title: 'riskclustr: Functions to Study Etiologic Heterogeneity'
tags:
  - R
  - biostatistics
  - epidemiology
  - etiologic heterogeneity
  - clustering
  - disease subtypes
  - multinomial logistic regression
  - polytomous logistic regression
authors:
 - name: Emily C. Zabor
   orcid: 0000-0002-1402-4498
   affiliation: 1
affiliations:
 - name: Memorial Sloan Kettering Cancer Center
   index: 1
date: 19 February 2019
bibliography: paper.bib
---

# Summary

Etiologic heterogeneity refers to the concept of subtypes of disease that are influenced by different risk factors. In cancer epidemiology, it is well known that many types of cancer demonstrate etiologic heterogeneity with respect to subtypes formed by genetic markers and/or other tumor characteristics. Likely the best known example of this is in breast cancer research, where subtypes are often formed based on immunohistochemical staining of estrogen receptor (ER), progesterone receptor (PR) and human epidermal growth factor receptor 2 (HER2). Risk factors, including patient characteristics such as age and body mass index as well as hormonal risk factors such as age at menarche, parity, and menopausal status, have been shown to have different relative risks for these disease subtypes [@gaudet2011; @kwan2009; @yang2007; @phipps2008; @ma2010]. With the growing use of genetic testing, the search for disease subtypes will only become more common in cancer research, and across other disease areas. There may be three main goals in any study of etiologic heterogeneity: 1) quantifying the extent of etiologic heterogeneity for a given set of disease subtypes according to all known risk factors, 2) testing the etiologic heterogeneity of individual risk factors with respect to a given set of disease subtypes, and 3) identifying the most etiologically heterogeneous subtype solution from possibly high-dimensional disease marker data. The R [@rct2018] package `riskclustr` [@zabor2019] provides user-friendly functions to address all three areas of study.

In previous work we introduced a scalar measure that can be used to quantify the extent of etiologic heterogeneity of pre-defined disease subtypes based on known risk factors in the context of a case-control study [@begg2013]. Later work generalized this methodology to the context of case-only studies, as researchers are often working in hospital settings without access to control subjects [@begg2014]. The R package `riskclustr` introduces original functionality to implement these methods to estimate the extent of etiologic heterogeneity in case-control studies using the `d` function and in case-only studies using the `dstar` function. 

Separately, we compared available statistical methods for the study of etiologic heterogeneity in case-control studies [@zabor2017]. By unifying the notation of the different methods and employing simulations, we showed that the methods are all able to address two key research questions: 1) whether risk factor effects differ across subtypes of disease and 2) whether risk factor effects differ across levels of each individual disease marker of which the disease subtypes are comprised. Our research also showed that adaptation of polytomous logistic regression has statistical properties at least as good as the more sophisticated methods that have been proposed, provided that the number of disease markers is sufficiently small that the analysis is feasible [@zabor2017]. In R polytomous logistic regression can be implemented using the `multinom` function from the `nnet` package [@venables2002] or the `mlogit` function from the `mlogit` package [@croissant2018], but the additional calculations needed to perform an analysis of etiologic heterogeneity are cumbersome. To facilitate use of this method the R package `riskclustr` introduces functions `eh_test_subtype` and `eh_test_marker` that first fit a standard polytomous logistic regression model using the `mlogit` function from the `mlogit` package [@croissant2018] and then perform additional calculations to address the two preceding questions regarding etiologic heterogeneity. 

Finally, it is increasingly common for statistical or epidemiologic researchers to be confronted with high dimensional disease marker data and when disease subtypes are not pre-defined, this disease marker data can be clustered, and the optimally etiologically heterogeneous subtype solution can be identified. This methodology has been applied to breast cancer [@begg2013; @begg2015] and melanoma [@mauguen2017] and the statistical properties of the approach have been investigated using simulation studies (Zabor et al, under review). In R [@rct2018] unsupervised $k$-means clustering can be implemented using the `kmeans` function in the `stats` package [@rct2018], and the R package`riskclustr` includes a wrapper for the `kmeans` function that calculates the extent of etiologic heterogeneity using the `d` function for each cluster solution resulting from many random starts of the clustering algorithm, and returns the subtype solution that maximizes etiologic heterogeneity. 

The R package `riskclustr` was designed for use by researchers in epidemiology and biostatistics, and the open-source software package includes  user-friendly tutorials that include examples of how to use the various functions and cover details of all underlying statistical calculations.

# Acknowledgements

The research was supported by the National Cancer Institute, awards CA163251 and CA167237. Additional funding support was provided by the Core Grant (P30 CA008748).

# References
