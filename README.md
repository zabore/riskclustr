# riskclustr

Functions related to the study of etiologic heterogeneity both across disease subtypes and across individual tumor markers.

You can install the current version using the command

```
devtools::install_github("zabore/riskclustr")
```

The primary functions in this package are 

1. `dest`, which estimates the D metric, a measure of the incremental explained risk variation for a set of pre-defined subtypes
2. `plrtm` and `plrsub`, which fit polytomous logistic regression models for subtypes defined by individual tumor markers or pre-defined subtypes, respectively, and result in tests of etiologic heterogeneity for each risk factor
