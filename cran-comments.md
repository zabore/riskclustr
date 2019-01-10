## Test environments
* local Windows 7 Professional install, R 3.5.1
* Ubuntu 14.04.5 LTS (on travis-ci.org), R 3.5.1 
* win-builder (devel and release)


## R CMD check results
There were no ERRORs or NOTEs

* On the local Windows 7 machine only, there was one WARNING:

    WARNING
   'qpdf' is needed for checks on size reduction of PDFs
   
   This warning is not applicable as no pdf vignettes are being produced, only html vignettes.
   
   This warning was not produced in the travis-ci or win-builder checks.
   
* One example is wrapped in \donttest due to longer run time, but is still made available for the user to run locally as needed. 
