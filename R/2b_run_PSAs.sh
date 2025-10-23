#!/bin/bash

# # Define an array of parameters for the R script
# params=(
#     5
#     10
#     15
#     20
#     25
#     30
#     40
#     50
# )

# # Loop through each parameter and execute the R script
# for param in "${params[@]}"; do
#     echo $param;
#     R --slave --vanilla --args $param < ./utils/PSA_single_N.R
# done


# aggregate version (in parallel)
# 30, 35, 40, 45, 50
R --slave --vanilla --args 30 < ./utils/PSA_single_N.R & R --slave --vanilla --args 35 < ./utils/PSA_single_N.R & R --slave --vanilla --args 40 < ./utils/PSA_single_N.R & R --slave --vanilla --args 45 < ./utils/PSA_single_N.R & R --slave --vanilla --args 50 < ./utils/PSA_single_N.R

# 5, 10, 15, 20, 25
R --slave --vanilla --args 5 < ./utils/PSA_single_N.R & R --slave --vanilla --args 10 < ./utils/PSA_single_N.R & R --slave --vanilla --args 15 < ./utils/PSA_single_N.R & R --slave --vanilla --args 20 < ./utils/PSA_single_N.R  & R --slave --vanilla --args 25 < ./utils/PSA_single_N.R
