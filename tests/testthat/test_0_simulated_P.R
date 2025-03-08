source("setup_Pmatrix_TruncMVN.R")


source("test_funcs.R")
library(tidyverse)

small_test_convergence_control <- new_convergence_control(maxiters = 1000, MAP_over = 500)
large_test_convergence_control <- new_convergence_control(maxiters = 10000, MAP_over = 1000)



get_cosmic <- function() {
    P <- read.csv(
        "https://cog.sanger.ac.uk/cosmic-signatures-production/documents/COSMIC_v3.3.1_SBS_GRCh37.txt",
        sep = '\t'
    ) %>%
        column_to_rownames("Type") %>%
        as.matrix() %>%
        # as.numeric()
        return(P)
}

cosmic_matrix <- get_cosmic()

cor_matrix = cor(t(cosmic_matrix))

library(corpcor)
# Shrink the correlation matrix to ensure invertibility
shrunk_cor_matrix <- cov2cor(cov.shrink(cor_matrix))




# MVN
bayesNMF(
    M, rank = 4,
    likelihood = 'normal',
    prior = 'truncnormal',
    file = "log_files/simulated_P/mar8_modelNT_dataP_N4_g6_707neg05_startval_lowerupperinf",
    overwrite = TRUE,
    true_P = true_P,
    prior_parameters = list(Cor_p = shrunk_cor_matrix, Cor_p_inv = solve(shrunk_cor_matrix)),
    convergence_control = large_test_convergence_control
)


# NONMVN
bayesNMF(
    M, rank = 4,
    likelihood = 'normal',
    prior = 'truncnormal',
    file = "log_files/simulated_P/mar8_NONMVN_modelNT_dataP_N4_g6_707neg05_startval",
    overwrite = TRUE,
    true_P = true_P,
    convergence_control = large_test_convergence_control
)
