# geuvadis-predixcan
Pipeline to train and analyze [PrediXcan](https://github.com/hakyim/PrediXcan) prediction weights using data from the Genetic European Variation in Health and Disease ([GEUVADIS](www.geuvadis.org/web/geuvadis)) project.

Genotypes and normalized gene expression measures for 373 EUR samples and 89 YRI samples from GEUVADIS were downloaded from the [EBI website](https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/analysis_results/).

New prediction weights for each population were computed using [elastic net regularized](https://en.wikipedia.org/wiki/Elastic_net_regularization) linear regression with the [`glmnet`](https://cran.r-project.org/web/packages/glmnet/index.html) package in R.
