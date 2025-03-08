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

M_pcawg <- read.csv(
    "Panc-Endocrine.csv", # Panc-AdenoCA.csv
    header = TRUE, row.names = 1, sep = ",")
M_pcawg <- as.matrix(M_pcawg)

library(corpcor)
# Shrink the correlation matrix to ensure invertibility
shrunk_cor_matrix <- cov2cor(cov.shrink(cor_matrix))


# Panc-Endocrine
panc_endocrine_columns <- c("SBS2", "SBS3", "SBS5", "SBS6", "SBS30", "SBS36")
true_P_panc_endocrine <- cosmic_matrix[, panc_endocrine_columns]

# MVN
res <- bayesNMF(
    M = M_pcawg, rank = 6,
    likelihood = 'normal',
    prior = 'truncnormal',
    file = "log_files/pcawg/mvn/mar6_panc_endocrine_N6_large",
    overwrite = FALSE,
    # true_P = true_P_panc_endocrine,
    prior_parameters = list(Cor_p = shrunk_cor_matrix, Cor_p_inv = solve(shrunk_cor_matrix)),
    convergence_control = large_test_convergence_control
)
get_heatmap(res$MAP$P, true_P_panc_endocrine) # cosmic_matrix

# HIERARCHICAL
res <- bayesNMF(
    M = M_pcawg, rank = 6,
    likelihood = 'normal',
    prior = 'truncnormal',
    file = "log_files/pcawg/hierarchical/mar6_panc_endocrine_N6_noTrueP",
    overwrite = FALSE,
    # true_P = true_P_panc_endocrine,
    prior_parameters = list(sigma2_prior = c(shape = 2, scale = 2),
                            rho_same_prior = c(a = 2, b = 2),
                            rho_diff_prior = c(a = 2, b = 2)),
    convergence_control = small_test_convergence_control
)
get_heatmap(res$MAP$P, true_P_panc_endocrine)

# NONMVN
res <- bayesNMF(
    M = M_pcawg, rank = 6,
    likelihood = 'normal',
    prior = 'truncnormal',
    file = "log_files/pcawg/nonmvn/mar6_panc_endocrine_N6_large",
    overwrite = FALSE,
    convergence_control = large_test_convergence_control
)
get_heatmap(res$MAP$P, true_P_panc_endocrine)
