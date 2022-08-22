# ukbiobank-blood-trait-analysis

This is the GitHub repository associated with the Masters thesis of Andrew Elmore, and is broken into 3 parts:

1. The `biobank` directory is responsible for parsing the UK Biobank dataset and creating residuals per patient to be used for all subsequent GWAS and PGS calculations.

2. The `lassosum` directory used lassosum analysis, which used the lassosum repository outlined here: https://github.com/tshmak/lassosum

3. The `ml` direcotyr used Bespoke SNP recalculation, chosen by 2 separate processes, the first being a Bayesian regression uses original effect sizesas priors, the other being a basic linear regression with no priors.
