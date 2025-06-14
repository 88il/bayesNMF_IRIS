#' Compute log target pdf for updating Pn with Poisson model
#'
#' @param M mutational catalog matrix, K x G
#' @param Pn vector length K, value of Pn to evaluate
#' @param n integer, signature index
#' @param Theta list of parameters
#' @param prior string, one of c('truncnormal','exponential')
#'
#' @return scalar
#' @noRd
log_target_Pn <- function(M, Pn, n, Theta, prior) {
    Theta$P[,n] <- Pn
    Mhat <- get_Mhat(Theta)
    if (prior == "truncnormal") {
        log_prior <- log(truncnorm::dtruncnorm(
            Pn, mean = Theta$Mu_p[,n], sd = sqrt(Theta$Sigmasq_p[,n]),
            a = 0, b = Inf
        ))
    } else if (prior == 'exponential') {
        log_prior <- dexp(
            Pn, Theta$Lambda_p[,n], log = TRUE
        )
    }
    log_likelihood <- rowSums(dpois(M, lambda = Mhat, log = TRUE))
    log_target <- log_prior + log_likelihood
    return(log_target)
}

#' Get mu and sigmasq for P[,n] in Normal model
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param prior string, one of c('truncnormal','exponential')
#' @param gamma double, tempering parameter
#'
#' @return list of two items, mu and sigmasq
#' @noRd
get_mu_sigmasq_Pn <- function(n, M, Theta, prior, dims, gamma = 1) {
    Mhat_no_n <- get_Mhat_no_n(Theta, n)

    # MVN case for P, when P has covariance matrix
    if (prior == 'truncnormal' & !is.null(Theta$Covar_p)) {
        cov_p = Theta$Covar_p

        sigma_M <- diag(Theta$sigmasq) # KxK
        sigma_M_inv <- solve(sigma_M) # KxK

        if (any(diag(sigma_M_inv) <= 0)) {
            stop("sigma_M_inv contains non-positive values.")
        }

        A_star <- sapply(1:dims$K, function(k) sum(Theta$E[n, ] * (M[k, ] - Mhat_no_n[k, ]))) # Kx1 vector
        A <- Theta$Covar_p_inv + sum(Theta$E[n, ] ** 2) * sigma_M_inv # KxK

        if (any(diag(Theta$Covar_p_inv) <= 0)) {
            stop("Theta$Covar_p_inv contains non-positive values.")
        }

        B <- Theta$Covar_p_inv %*% Theta$Mu_p[, n] + sigma_M_inv %*% A_star # Kx1
        A_inv <- solve(A)
        mu_P <- A_inv %*% B
        covar_P <- A_inv
        covar_P_inv <- A

        return(list(
            mu = mu_P,
            covar = covar_P,
            covar_inv = covar_P_inv
        ))
    }

    # non MVN case
    else {
        # broadcast En / sigmasq to residuals M - Mhat
        mu_num_term_1 <- (gamma * Theta$A[1,n] * sweep(
            (M - Mhat_no_n), # dim KxG
            2, # multiply each row by E[n,]
            Theta$E[n, ] , # length G
            "*"
        ) %>% # dim KxG
            rowSums()) / Theta$sigmasq # length K

        denom <- gamma * sum(Theta$A[1,n] * Theta$E[n, ] ** 2) / Theta$sigmasq

        if (prior == 'exponential') {
            mu_num_term_2 <- Theta$Lambda_p[, n] # length K
            mu_P <- (mu_num_term_1 - mu_num_term_2) / denom # length K
            sigmasq_P <- 1 / denom # length K
        } else if (prior == 'truncnormal') {
            mu_num_term_2 <- Theta$Mu_p[, n] / Theta$Sigmasq_p[,n] # length K
            denom <- denom + (1/Theta$Sigmasq_p[,n])
            mu_P <- (mu_num_term_1 + mu_num_term_2) / denom # length K
            sigmasq_P <- 1 / denom # length K
        }

        return(list(
            mu = mu_P,
            sigmasq = sigmasq_P
        ))
    }
}

#' Sample P[,n] for Normal likelihood
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param prior string, one of c('truncnormal','exponential')
#' @param gamma double, tempering parameter
#'
#' @return vector length K
#' @noRd
sample_Pn_normal <- function(n, M, Theta, dims, prior, gamma = 1) {
    # MVN case for P
    if (prior == 'truncnormal' & !is.null(Theta$Covar_p) ) {

        # Obtain posterior mu, sigma for P
        mu_sigmasq_P <- get_mu_sigmasq_Pn(n, M, Theta, prior, dims, gamma = gamma)
        mean_vector <- mu_sigmasq_P$mu[,1]

        # Perturb diagonal otherwise sigma is not positive definite
        epsilon <- 0 # 1e-5
        newsigma <- mu_sigmasq_P$covar + diag(epsilon, nrow(Theta$Covar_p))

        newsigma_inv <- mu_sigmasq_P$covar_inv
        newsigma_inv <- as.matrix(Matrix::forceSymmetric(newsigma_inv))

        eigenvalues <- eigen(newsigma_inv)$values
        if (any(eigenvalues <= 0)) {
            stop("inv Covariance matrix is not positive definite")
        }

        newsigma <- as.matrix(Matrix::forceSymmetric(newsigma))
        newsigma <- lqmm::make.positive.definite(newsigma, tol=1e-3)

        mean_vector[mean_vector < 0] <- 0

        # Set lower, upper bounds for MVNtruncnorm sampler
        lower <- rep(0, length(mean_vector))
        upper <- rep(Inf, length(mean_vector))

        # looser upper_stricter bound
        lower_stricter <- pmax(lower, mean_vector - 10 * sqrt(diag(newsigma)))
        upper_stricter <- pmin(upper, mean_vector + 10 * sqrt(diag(newsigma)))

        sample <- tmvtnorm::rtmvnorm(
            1, mean = mean_vector, H = newsigma_inv,
            lower = lower_stricter, upper = upper_stricter,
            algorithm = "gibbs",
            start.value = pmax(lower, pmin(upper, mean_vector)),
            burn.in.samples = 10
        )

        return(sample)
    }

    # non MVN case
    else {
        # Normal-exponential doesn't collapse like normal-truncated normal
        # sample from prior when Ann = 0 or gamma = 0
        if (prior == 'exponential' & (Theta$A[1,n] == 0 | gamma == 0)) {
            sampled <- stats::rexp(dims$K, Theta$Lambda_p[,n])
            return(sampled)
        }

        # Otherwise, compute mean and sd for sample from full conditional
        mu_sigmasq_P <- get_mu_sigmasq_Pn(n, M, Theta, prior, dims, gamma = gamma)

        # Sample from truncated normal
        sampled <- truncnorm::rtruncnorm(
            1, mean = mu_sigmasq_P$mu, sd = sqrt(mu_sigmasq_P$sigmasq),
            a = 0, b = Inf
        )
        return(sampled)
    }
}

#' sample P[,n] for Poisson likelihood
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param prior string, one of c('gamma','exponential')
#' @param gamma double, tempering parameter
#'
#' @return vector length K
#' @noRd
sample_Pn_poisson <- function(n, M, Theta, dims, prior, gamma = 1) {
    if (prior == 'gamma') {
        shape <- sapply(1:dims$K, function(k) {
            Theta$Alpha_p[k,n] + gamma * sum(Theta$Z[k,n,])
        })
        rate <- sapply(1:dims$K, function(k) {
            Theta$Beta_p[k,n] + gamma * Theta$A[1,n] * sum(Theta$E[n,])
        })
    } else if (prior == 'exponential') {
        shape <- sapply(1:dims$K, function(k) {
            1 + gamma * sum(Theta$Z[k,n,])
        })
        rate <- sapply(1:dims$K, function(k) {
            Theta$Lambda_p[k,n] + gamma * Theta$A[1,n] * sum(Theta$E[n,])
        })
    }
    sampled <- sapply(1:dims$K, function(k) { rgamma(1, shape[k], rate[k]) })
    return(sampled)
}

#' Sample P[,n] wrapper function
#'
#' @param n signature index
#' @param M mutational catalog matrix, K x G
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param likelihood string, one of c('normal','poisson')
#' @param prior string, one of c('gamma','exponential','truncnormal')
#' @param gamma double, tempering parameter. DO NOT CHANGE, KEEP gamma = 1
#'
#' @return vector length K
#' @noRd
sample_Pn <- function(n, M, Theta, dims, likelihood, prior, fast, gamma = 1) {
    if (likelihood == 'normal' | (likelihood == 'poisson' & fast)) {
        proposal <- sample_Pn_normal(n, M, Theta, dims, prior, gamma)
        if (likelihood == 'poisson' & fast) {
            acceptance <- exp(
                log_target_Pn(M, proposal, n, Theta, prior) -
                log_target_Pn(M, Theta$P[,n], n, Theta, prior)
            )

            acceptance[is.na(acceptance)] <- 0.5 # NaN if 0/0, give 50-50 chance
            acceptance[acceptance > 1] <- 1 # cap acceptance probability at 1

            # accept with probability acceptance
            accepted <- runif(dims$K) < acceptance
            sampled <- Theta$P[,n]
            sampled[accepted] <- proposal[accepted]
        } else {
            # likelihood == 'normal'
            sampled <- proposal
            acceptance <- 1
        }
    } else {
        # likelihood == 'poisson' & !fast
        sampled <- sample_Pn_poisson(n, M, Theta, dims, prior, gamma)
        acceptance <- 1
    }

    # print("in samplePn")
    # print(sampled)

    return(list(
        sampled = sampled,
        acceptance = acceptance
    ))
}
