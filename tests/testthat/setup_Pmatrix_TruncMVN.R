source("test_funcs.R")

set.seed(123)
N = 8
G = 10 # 60

# Construct P (signatures) matrix.
# Instead of constructing P from sampling N signatures from COSMIC,
# we simulate P given a fixed Covar_P matrix.
covar_P <- build_covariance_matrix_P(7, 0.5, -0.1, 96) # reduce rho_same to reduce correlation b/w signatures

# could use cosmic covar and mean
cosmic_covar <- cov(t(cosmic_matrix))
cosmic_mean <- rowMeans(cosmic_matrix)


mean_vector <- sqrt(diag(covar_P)) # length 96
# other options: set to all 1, set to constant * dirichlet
# rowMeans of cosmic_matrix (* constant ?)

# Set lower, upper bounds for MVNtruncnorm sampler
lower <- rep(0, length(mean_vector))
upper <- rep(Inf, length(mean_vector))

lower_stricter <- pmax(lower, mean_vector - 10 * sqrt(diag(covar_P)))
upper_stricter <- pmin(upper, mean_vector + 10 * sqrt(diag(covar_P)))

# take the transpose of the sample so simulated_P is 96 x N
simulated_P <- t(tmvtnorm::rtmvnorm(
    n = N,
    mean = cosmic_mean, # mean_vector,
    # H = solve(covar_P),
    sigma = lqmm::make.positive.definite(cosmic_covar, tol=1e-3), # covar_P,
    # lower = lower_stricter, upper = upper_stricter,
    lower = lower, upper = upper,
    algorithm = "gibbs" #, burn.in.samples = 1 #
    , start.value = mean_vector # cosmic_mean # mean_vector
))

# normalize so that columns of P sum to 1 <-- NECESSARY ??
simulated_P <- apply(simulated_P, 2, function(col) col / sum(col))

# want low some < 0.7
get_heatmap(simulated_P, simulated_P)



K = nrow(simulated_P)

E <- matrix(rexp(N*G, 0.001), nrow = N, ncol = G) # M on avg 96k mutations (96 * 1k) can switch to 0.01 (more realistic for non melanoma)

M <- matrix(nrow = K, ncol = G)
for (k in 1:K) {
    M[k,] <- rpois(G, simulated_P[k,]%*%E)
}
M <- round(M)

true_P <- simulated_P
true_E <- E
rm(list = c('simulated_P','E','N','G'))
