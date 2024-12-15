# source("setup_poisson.R")

source("setup_poisson_6.R")


source("test_funcs.R")
library(tidyverse)

small_test_convergence_control <- new_convergence_control(maxiters = 1000, MAP_over = 500)
large_test_convergence_control <- new_convergence_control(maxiters = 10000, MAP_over = 10, MAP_every=1)
# MAP_every changes ^ saves RDS object and can view PDF of results
# view pdf of results / chain
# look at what exactly is causing the NaN sampling
# can load in the covar and stuff and step by step sample


# 12: In rgamma(n, shape, rate) : NAs produced
# 13: In rnorm(dims$G, num/denom, 1/denom) : NAs produced
# 14: In rgamma(n, shape, rate) : NAs produced
# truncnorm sampling under the hood may use ^^ JL doesn't use any
# can try to find source code of truncnorm package

# can try to use MVN norm to sample, then set negatives to 0. this would be a mixture dist.



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

# change to 6x6 only center mutation. collapse the cosmic matrix
cosmic_matrix_6 <- matrix(nrow = 6, ncol = ncol(cosmic_matrix))
start <- 1
end <- 16
for (i in (1:6)){
    cosmic_matrix_6[i,] = as.matrix(colSums(cosmic_matrix[start:end,]))
    start <- end + 1
    end <- start + 15
}


cor_matrix = cor(t(cosmic_matrix_6))
# try covar matrix then scaling

# try without zeroing (to ensure positive definiteness)
# try with covar matrix (scale after)

library(Matrix)
sigma = cov(t(cosmic_matrix_6))
sigma = do.call(bdiag, lapply(1:16, function(i) {
    sigma[((i-1)*6+1):(i*6), ((i-1)*6+1):(i*6)]
})) %>% as.matrix()


block_cor = do.call(bdiag, lapply(1:16, function(i) {
    cor_matrix[((i-1)*6+1):(i*6), ((i-1)*6+1):(i*6)]
})) %>% as.matrix()


library(corpcor)
# Shrink the correlation matrix to ensure invertibility
shrunk_cor_matrix_6 <- cov2cor(cov.shrink(cor_matrix))

# make cor_matrix invertible
cor_matrix_1 <- cor_matrix
cutoff = quantile(abs(cor_matrix_1), 0.1)
cor_matrix_1[abs(cor_matrix_1) < cutoff] = 0


# example with diagonal correlation matrix
res <- bayesNMF(
    M, rank = 2,
    likelihood = 'normal',
    prior = 'truncnormal',
    file = "log_files/mvn/modelNT_dataP_N1_mutation6",
    overwrite = TRUE,
    true_P = true_P,
    inits = list(P= true_P[,1:2]), # starting at true value (cheating)
    prior_parameters = list(Cor_p = cor_matrix_1, Cor_p_inv = solve(cor_matrix_1)),
    convergence_control = large_test_convergence_control
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
