# CLEAN

Installation:

```{r}
#In order to avoid unnecessary warnings, check for installation of devtools (>= 1.12.0):
packageVersion("devtools")
[1] ‘1.12.0’

#Companion data packages
devtools::install_github("uc-bd2k/CLEAN.Hs")
devtools::install_github("uc-bd2k/CLEAN.Mm")
devtools::install_github("uc-bd2k/CLEAN.Rn")

#Main package
devtools::install_github("uc-bd2k/CLEAN")

```
