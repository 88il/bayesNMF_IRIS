source("setup_poisson_hierarchical.R")

source("test_funcs.R")
library(tidyverse)

small_test_convergence_control <- new_convergence_control(maxiters = 1000, MAP_over = 500)
large_test_convergence_control <- new_convergence_control(maxiters = 10000, MAP_over = 1000)

feb6_convg_ctrl <- new_convergence_control(maxiters = 1500, MAP_over = 500)
# # below is mutation6 large convergence control
# large_test_convergence_control <- new_convergence_control(maxiters = 10000, MAP_over = 10, MAP_every=1)


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


# # REORDER COSMIC rows to be grouped by center mutation
# true_P_reordered <- reorder_matrix_by_center_mutation(true_P)



# NONCHEAT
res <- bayesNMF(
    M, rank = 3,
    likelihood = 'normal',
    prior = 'truncnormal',
    file = "log_files/hierarchical_mvn/hierarchical_3param_N3_noncheat_feb12",
    overwrite = TRUE,
    true_P = true_P,
    prior_parameters = list(sigma2_prior = c(shape = 2, scale = 2),
                            rho_same_prior = c(a = 2, b = 2),
                            rho_diff_prior = c(a = 2, b = 2)),
    convergence_control = small_test_convergence_control
)
get_heatmap(res$MAP$P, true_P)


# CHEAT
res <- bayesNMF(
    M, rank = 3,
    likelihood = 'normal',
    prior = 'truncnormal',
    file = "log_files/hierarchical_mvn/hierarchical_3param_N3_cheat",
    overwrite = TRUE,
    true_P = true_P,
    inits = list(P= true_P[,1:3]), # starting at true value (cheating)
    prior_parameters = list(sigma2_prior = c(shape = 2, scale = 2),
                            rho_same_prior = c(a = 2, b = 2),
                            rho_diff_prior = c(a = 2, b = 2)),
    convergence_control = large_test_convergence_control
)
get_heatmap(res$MAP$P, true_P)
