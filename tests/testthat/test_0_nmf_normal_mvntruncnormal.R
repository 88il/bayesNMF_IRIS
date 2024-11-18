source("setup_poisson.R")
source("test_funcs.R")
library(tidyverse)

small_test_convergence_control <- new_convergence_control(maxiters = 1000, MAP_over = 500)
large_test_convergence_control <- new_convergence_control(maxiters = 2000, MAP_over = 500)

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

# make cor_matrix invertible
cor_matrix_1 <- cor_matrix
cutoff = quantile(abs(cor_matrix_1), 0.1)
cor_matrix_1[abs(cor_matrix_1) < cutoff] = 0

# example with diagonal correlation matrix
res <- bayesNMF(
    M, rank = 2,
    likelihood = 'normal',
    prior = 'truncnormal',
    file = "log_files/mvn/modelNT_dataP_N1",
    overwrite = TRUE,
    true_P = true_P,
    inits = list(P= true_P[,1:2]), # starting at true value (cheating)
    prior_parameters = list(Cor_p = diag(diag(cor_matrix_1)), Cor_p_inv = solve(diag(diag(cor_matrix_1)))),
    convergence_control = small_test_convergence_control
)
get_heatmap(res$MAP$P, true_P)

# should be the same as
res <- bayesNMF(
    M, rank = 5,
    likelihood = 'normal',
    prior = 'truncnormal',
    file = "log_files/mvn/modelNT_dataP_N1",
    overwrite = TRUE,
    true_P = true_P,
    prior_parameters = list(Sigmasq_p = rep(2.25, 96)),
    convergence_control = small_test_convergence_control
)
get_heatmap(res$MAP$P, true_P)

test_that("nmf_normal_mvntruncnormal works with 1 signature given N", {
    res <- bayesNMF(
        M, rank = 5,
        likelihood = 'normal',
        prior = 'truncnormal',
        file = "log_files/mvn/modelNT_dataP_N1",
        overwrite = TRUE,
        true_P = true_P,
        prior_parameters = list(Cor_p = diag(diag(cor_matrix_1)), Cor_p_inv = solve(diag(diag(cor_matrix_1)))),
        convergence_control = small_test_convergence_control
    )

    expect_equal(sum(is.na(res$MAP$P)), 0)
    expect_equal(sum(is.na(res$MAP$E)), 0)
    heatmap <- get_heatmap(res$MAP$P, true_P)
    expect_equal(class(heatmap), c('gg','ggplot'))
})

# test_that("nmf_normal_mvntruncnormal works with 2 signatures given N", {
#     res <- bayesNMF(
#         M, N = 2,
#         likelihood = 'normal',
#         prior = 'truncnormal',
#         file = "log_files/mvn/modelNT_dataP_N2",
#         overwrite = TRUE,
#         true_P = true_P
#     )
#
#     expect_equal(sum(is.na(res$MAP$P)), 0)
#     expect_equal(sum(is.na(res$MAP$E)), 0)
#
#     expect_equal(ncol(res$MAP$P), 2)
#     expect_equal(nrow(res$MAP$P), 96)
#     expect_equal(nrow(res$MAP$E), 2)
#
#     log_post <- get_proportional_log_posterior(
#         Theta = res$final_Theta,
#         M = M,
#         P = res$MAP$P,
#         E = res$MAP$E,
#         sigmasq = res$MAP$sigmasq,
#         likelihood = 'normal',
#         prior = 'truncnormal'
#     )
#     expect_true(!is.na(log_post))
# })
#
# test_that("nmf_normal_mvntruncnormal works with Poisson data generating function given N", {
#     res <- bayesNMF(
#         M, N = 5,
#         likelihood = 'normal',
#         prior = 'truncnormal',
#         file = "log_files/mvn/modelNT_dataP_N5",
#         overwrite = TRUE,
#         true_P = true_P
#     )
#
#     expect_equal(sum(is.na(res$MAP$P)), 0)
#     expect_equal(sum(is.na(res$MAP$E)), 0)
#
#     expect_equal(ncol(res$MAP$P), 5)
#     expect_equal(nrow(res$MAP$P), 96)
#     expect_equal(nrow(res$MAP$E), 5)
#
#     sig_sims <- diag(reassign_signatures(res$sim_mat))
#     sig_sims <- sig_sims[sig_sims != min(sig_sims)]
#     expect_gt(min(sig_sims), 0.8)
# })
#
# test_that("nmf_normal_mvntruncnormal works with Poisson data generating function given max N", {
#     res <- bayesNMF(
#         M, max_N = 7,
#         likelihood = 'normal',
#         prior = 'truncnormal',
#         file = "log_files/mvn/modelNT_dataP_maxN7",
#         overwrite = TRUE,
#         true_P = true_P
#     )
#
#     expect_equal(sum(is.na(res$MAP$P)), 0)
#     expect_equal(sum(is.na(res$MAP$E)), 0)
#
#     expect_lt(abs(sum(res$MAP$A) - 5), 2)
#
#     sig_sims <- diag(reassign_signatures(res$sim_mat))
#     sig_sims <- sig_sims[sig_sims != min(sig_sims)]
#     expect_gt(min(sig_sims), 0.8)
# })
#
# test_that("nmf_normal_mvntruncnormal works with Poisson data generating function given max N and custom a, b", {
#     res <- bayesNMF(
#         M, max_N = 7,
#         likelihood = 'normal',
#         prior = 'truncnormal',
#         file = "log_files/mvn/modelNT_dataP_maxN7_customab",
#         overwrite = TRUE,
#         true_P = true_P,
#         prior_parameters = list('a' = 0.8, 'b' = 0.4)
#     )
#
#     expect_equal(sum(is.na(res$MAP$P)), 0)
#     expect_equal(sum(is.na(res$MAP$E)), 0)
#
#     expect_lt(abs(sum(res$MAP$A) - 5), 2)
#
#     sig_sims <- diag(reassign_signatures(res$sim_mat))
#     sig_sims <- sig_sims[sig_sims != min(sig_sims)]
#     expect_gt(min(sig_sims), 0.8)
# })
#
#
# source("setup_poisson_sparse.R")
#
# test_that("nmf_normal_mvntruncnormal works with sparse Poisson data generating function given N", {
#     res <- bayesNMF(
#         M, N = 5,
#         likelihood = 'normal',
#         prior = 'truncnormal',
#         file = "log_files/mvn/modelNT_dataPS_N5",
#         overwrite = TRUE,
#         true_P = true_P
#     )
#
#     expect_equal(sum(is.na(res$MAP$P)), 0)
#     expect_equal(sum(is.na(res$MAP$E)), 0)
#
#     sig_sims <- diag(reassign_signatures(res$sim_mat))
#     sig_sims <- sig_sims[sig_sims != min(sig_sims)]
#     expect_gt(min(sig_sims), 0.8)
# })
#
# test_that("nmf_normal_mvntruncnormal works with sparse Poisson data generating function given max N", {
#     res <- bayesNMF(
#         M, max_N = 7,
#         likelihood = 'normal',
#         prior = 'truncnormal',
#         file = "log_files/mvn/modelNT_dataPS_maxN7",
#         overwrite = TRUE,
#         true_P = true_P
#     )
#
#     expect_equal(sum(is.na(res$MAP$P)), 0)
#     expect_equal(sum(is.na(res$MAP$E)), 0)
#
#     expect_lt(abs(sum(res$MAP$A) - 5), 2)
#
#     sig_sims <- diag(reassign_signatures(res$sim_mat))
#     sig_sims <- sig_sims[sig_sims != min(sig_sims)]
#     expect_gt(min(sig_sims), 0.8)
# })
#
# source("setup_normal.R")
#
# test_that("nmf_normal_mvntruncnormal works with Normal data generating function given N", {
#     res <- bayesNMF(
#         M, N = 5,
#         likelihood = 'normal',
#         prior = 'truncnormal',
#         file = "log_files/mvn/modelNT_dataN_N5",
#         overwrite = TRUE,
#         true_P = true_P
#     )
#
#     expect_equal(sum(is.na(res$MAP$P)), 0)
#     expect_equal(sum(is.na(res$MAP$E)), 0)
#
#     sig_sims <- diag(reassign_signatures(res$sim_mat))
#     sig_sims <- sig_sims[sig_sims != min(sig_sims)]
#     expect_gt(min(sig_sims), 0.8)
# })
#
# test_that("nmf_normal_mvntruncnormal works with Normal data generating function given max N", {
#     res <- bayesNMF(
#         M, max_N = 7,
#         likelihood = 'normal',
#         prior = 'truncnormal',
#         file = "log_files/mvn/modelNT_dataN_maxN7",
#         overwrite = TRUE,
#         true_P = true_P
#     )
#
#     expect_equal(sum(is.na(res$MAP$P)), 0)
#     expect_equal(sum(is.na(res$MAP$E)), 0)
#
#     expect_lt(abs(sum(res$MAP$A) - 5), 2)
#
#     sig_sims <- diag(reassign_signatures(res$sim_mat))
#     sig_sims <- sig_sims[sig_sims != min(sig_sims)]
#     expect_gt(min(sig_sims), 0.8)
# })
#
# source("setup_normal_sparse.R")
#
# test_that("nmf_normal_mvntruncnormal works with sparse Normal data generating function given N", {
#     res <- bayesNMF(
#         M, N = 5,
#         likelihood = 'normal',
#         prior = 'truncnormal',
#         file = "log_files/mvn/modelNT_dataNS_N5",
#         overwrite = TRUE,
#         true_P = true_P
#     )
#
#     expect_equal(sum(is.na(res$MAP$P)), 0)
#     expect_equal(sum(is.na(res$MAP$E)), 0)
#
#     sig_sims <- diag(reassign_signatures(res$sim_mat))
#     sig_sims <- sig_sims[sig_sims != min(sig_sims)]
#     expect_gt(min(sig_sims), 0.8)
# })
#
# test_that("nmf_normal_mvntruncnormal works with sparse Normal data generating function given max N", {
#     res <- bayesNMF(
#         M, max_N = 7,
#         likelihood = 'normal',
#         prior = 'truncnormal',
#         file = "log_files/mvn/modelNT_dataNS_maxN7",
#         overwrite = TRUE,
#         true_P = true_P
#     )
#
#     expect_equal(sum(is.na(res$MAP$P)), 0)
#     expect_equal(sum(is.na(res$MAP$E)), 0)
#
#     expect_lt(abs(sum(res$MAP$A) - 5), 2)
#
#     sig_sims <- diag(reassign_signatures(res$sim_mat))
#     sig_sims <- sig_sims[sig_sims != min(sig_sims)]
#     expect_gt(min(sig_sims), 0.8)
# })
