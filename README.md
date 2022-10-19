
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Coverage
status](https://codecov.io/gh/zabore/riskclustr/branch/master/graph/badge.svg)](https://codecov.io/github/zabore/riskclustr?branch=master)
[![R-CMD-check](https://github.com/zabore/riskclustr/workflows/R-CMD-check/badge.svg)](https://github.com/zabore/riskclustr/actions)
[![CRAN
status](https://www.r-pkg.org/badges/version/riskclustr)](https://cran.r-project.org/package=riskclustr)
[![status](https://joss.theoj.org/papers/10.21105/joss.01269/status.svg)](https://joss.theoj.org/papers/10.21105/joss.01269)

# riskclustr

Etiologic heterogeneity refers to the concept of subtypes of disease
that are influenced by different risk factors. Likely the best known
example of this is in breast cancer research, where subtypes are often
formed based on immunohistochemical staining of estrogen receptor (ER),
progesterone receptor (PR) and human epidermal growth factor receptor 2
(HER2). Risk factors including patient characteristics such as age and
BMI as well as hormonal risk factors such as age at menarche, parity,
and menopausal status have been shown to have different relative risks
for these disease subtypes.

`riskclustr` is a collection of functions related to the study of
etiologic heterogeneity both across disease subtypes and across
individual disease markers. The included functions allow one to quantify
the extent of etiologic heterogeneity in the context of a case-control
study or case-only study, and provide p-values to test for etiologic
heterogeneity with respect to either disease subtypes or individual
disease markers in the context of a case-control study.

## Installation

You can install the development version of `riskclustr` by running:

``` r
devtools::install_github("zabore/riskclustr")
library(riskclustr)
```

Or the prodcution version from CRAN by running:

``` r
install.packages("riskclustr")
library(riskclustr)
```

## Documentation

This package is documented using
[pkgdown](https://pkgdown.r-lib.org/articles/pkgdown.html), and the
resulting website is available at
<http://www.emilyzabor.me/riskclustr/>, where detailed Tutorials can be
found covering all of the package functionality.

See <http://www.emilyzabor.me/riskclustr/reference/> for detailed
function documentation.

## References

> Begg CB, Zabor EC, Bernstein JL, Bernstein L, Press MF, Seshan VE. A
> conceptual and methodological framework for investigating etiologic
> heterogeneity. *Stat Med.* 2013;32(29):5039-52. doi: 10.1002/sim.5902

> Begg CB, Seshan VE, Zabor EC, et al. Genomic investigation of
> etiologic heterogeneity: methodologic challenges. *BMC Med Res
> Methodol.* 2014;14:138. doi: 10.1186/1471-2288-14-138
