Summary and Statistiks
========================================================

This page gives an overview, of the cohorts, methods and data types used by TCGA-WebTools.
For each data type, you'll find a summary below

## Variant Section
Overall there are ```{r echo=F} length(unique(variant.Table$Tumor_Sample_Barcode)) ``` patients within ```{r echo=F} length(getCancerTypes(unique(variant.Table$Tumor_Sample_Barcode))) ``` entities. The cohorts include the following tumor entities with given the number of patient samples

