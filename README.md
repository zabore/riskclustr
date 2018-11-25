
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Coverage status](https://codecov.io/gh/zabore/riskclustr/branch/master/graph/badge.svg)](https://codecov.io/github/zabore/riskclustr?branch=master) [![Build Status](https://travis-ci.com/zabore/riskclustr.svg?branch=master)](https://travis-ci.com/zabore/riskclustr) [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing) [![CRAN status](https://www.r-pkg.org/badges/version/riskclustr)](https://cran.r-project.org/package=riskclustr)

riskclustr
==========

Etiologic heterogeneity refers to the concept of subtypes of disease that are influenced by different risk factors. Likely the best known example of this is in breast cancer research, where subtypes are often formed based on immunohistochemical staining of estrogen receptor (ER), progesterone receptor (PR) and human epidermal growth factor receptor 2 (HER2). Risk factors including patient characteristics such as age and BMI as well as hormonal risk factors such as age at menarche, parity, and menopausal status have been shown to have different relative risks for these disease subtypes.

`riskclustr` is a collection of functions related to the study of etiologic heterogeneity both across disease subtypes and across individual disease markers. The included functions allow one to quantify the extent of etiologic heterogeneity in the context of a case-control study or case-only study, and provide p-values to test for etiologic heterogeneity with respect to either disease subtypes or individual disease markers in the context of a case-control study.

Installation
------------

You can install `riskclustr` by running:

``` r
devtools::install_github("zabore/riskclustr")
library(riskclustr)
```

This package is documented using [pkgdown](https://pkgdown.r-lib.org/articles/pkgdown.html), and the resulting website is available at <http://www.emilyzabor.com/riskclustr/>, where detailed Tutorials can be found covering all of the package functionality.

See <http://www.emilyzabor.com/riskclustr/reference/> for detailed function documentation.
