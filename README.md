# BLUP
Calculations and statistical analyzes carried out as part of the master thesis, the topic of which was "Estimation of the breeding value for the production traits of the Polish Holstein - Friesian black and white variety using the BLUP method". The analysis included the following steps:
1. preliminary preparation of the data set
2. preparation of generalized linear models for cattle production traits
3. choosing the model best suited to the data, based on minimizing the AIC or BIC criterion
4. preparation of mixed models for cattle production traits based on previously constructed generalized linear models
5. estimation of variance components using REML and ML method
6. building a pedigree
7. building an inverse relationship matrix based on a constructed pedigree
8. estimation of breeding values of cattle by the BLUP method for each of the production traits (i.e. preparation of a data frame with individuals ranked from the highest to the lowest breeding value with the determination of the accuracy of their estimation)
Additionally, the script includes graphic diagnostics of mixed models and testing the significance of random effects using the likelihood ratio test (+ a few graphs presenting the analyzed data set).
