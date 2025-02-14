#### Truncated Normal Prior

#' Sample Prior Parameter: Mu for Pn
#'
#' @param n signature index
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return vector of length K
#' @noRd
sample_Mu_Pn <- function(n, Theta, dims, gamma) {
    num <- Theta$M_p[,n]/Theta$S_p[,n] + gamma*Theta$P[,n]/Theta$Sigmasq_p[,n]
    denom <- 1/Theta$S_p[,n] + gamma*1/Theta$Sigmasq_p[,n]
    rnorm(dims$K, num/denom, 1/denom)
}

#' Sample Prior Parameter: Mu for En
#'
#' @param n signature index
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return vector of length G
#' @noRd
sample_Mu_En <- function(n, Theta, dims, gamma) {
    num <- Theta$M_e[n,]/Theta$S_e[n,] + gamma*Theta$E[n,]/Theta$Sigmasq_e[n,]
    denom <- 1/Theta$S_e[n,] + gamma*1/Theta$Sigmasq_e[n,]
    rnorm(dims$G, num/denom, 1/denom)
}

#' Sample Prior Parameter: Sigmasq for Pn
#'
#' @param n signature index
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return vector of length K
#' @noRd
sample_Sigmasq_Pn <- function(n, Theta, dims, gamma) {
    invgamma::rinvgamma(
        n = dims$K,
        shape = Theta$A_p[,n] + gamma*1/2,
        rate = Theta$B_p[,n] + gamma*(Theta$P[,n] - Theta$Mu_p[,n])**2/2
    )
}

#' Sample Prior Parameter: Sigmasq for En
#'
#' @param n signature index
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return vector of length G
#' @noRd
sample_Sigmasq_En <- function(n, Theta, dims, gamma) {
    invgamma::rinvgamma(
        n = dims$G,
        shape = Theta$A_e[n,] + gamma*1/2,
        rate = Theta$A_e[n,] + gamma*(Theta$E[n,] - Theta$Mu_e[n,])**2/2
    )
}

#### Exponential Prior

#' Sample Prior Parameter: Lambda for Pn
#'
#' @param n signature index
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return vector of length K
#' @noRd
sample_Lambda_Pn <- function(n, Theta, dims, gamma) {
    rgamma(dims$K, Theta$A_p[,n] + gamma * 1, Theta$B_p[,n] + gamma * Theta$P[,n])
}

#' Sample Prior Parameter: Lambda for En
#'
#' @param n signature index
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return vector of length G
#' @noRd
sample_Lambda_En <- function(n, Theta, dims, gamma) {
    rgamma(dims$G, Theta$A_e[n,] + gamma * 1, Theta$B_e[n,] + gamma * Theta$E[n,])
}


#### Gamma Prior

#' Sample Prior Parameter: Beta for Pn
#'
#' @param n signature index
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return vector of length K
#' @noRd
sample_Beta_Pn <- function(n, Theta, dims, gamma) {
    rgamma(
        dims$K,
        Theta$A_p[,n] + gamma*Theta$Alpha_p[,n],
        Theta$B_p[,n] + gamma*Theta$P[,n]
    )
}

#' Sample Prior Parameter: Beta for En
#'
#' @param n signature index
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return vector of length G
#' @noRd
sample_Beta_En <- function(n, Theta, dims, gamma) {
    rgamma(
        dims$G,
        Theta$A_e[n,] + gamma*Theta$Alpha_e[n,],
        Theta$B_e[n,] + gamma*Theta$E[n,]
    )
}

#' Sample Prior Parameter: Alpha for Pkn
#'
#' @param k mutation type index
#' @param n signature index
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return scalar
#' @noRd
sample_Alpha_Pkn <- function(k, n, Theta, dims, gamma) {
    logpdf_prop <- function(x) {
        (gamma * Theta$Beta_p[k,n] + Theta$C_p[k,n] - 1) * log(x) -
        Theta$D_p[k,n] * x +
        gamma * (x - 1) * log(Theta$P[k,n]) -
        gamma * lgamma(x)
    }

    armspp::arms(
        n_samples = 1,
        log_pdf = logpdf_prop,
        lower = 1e-3,
        upper = 10000
    )
}

#' Sample Prior Parameter: Alpha for Eng
#'
#' @param n signature index
#' @param g tumor genome index
#' @param Theta list of parameters
#' @param dims list of dimensions
#' @param gamma double, tempering parameter
#'
#' @return scalar
#' @noRd
sample_Alpha_Eng <- function(n, g, Theta, dims, gamma) {
    logpdf_prop <- function(x) {
        (gamma * Theta$Beta_e[n,g] + Theta$C_e[n,g] - 1) * log(x) -
        Theta$D_e[n,g] * x +
        gamma * (x - 1) * log(Theta$E[n,g]) -
        gamma * lgamma(x)
    }

    armspp::arms(
        n_samples = 1,
        log_pdf = logpdf_prop,
        lower = 1e-3,
        upper = 10000
    )
}


# TODO: hierarchical bayesian model for Covar_P

#' Compute the unnormalized posterior for a single parameter
#'
#' @param param current value of the parameter being resampled
#' @param param_name name of the parameter ("sigma2", "rho_same", or "rho_diff")
#' @param P current signature matrix (K x N)
#' @param Theta list of parameters
#' @param dims list of dimensions
#'
#' @return log posterior value (unnormalized)
unnormalized_posterior <- function(param, param_name, P, Theta, dims) {
    print('inside unnormalized_posterior')

    # Update the relevant parameter in Theta
    Theta[[param_name]] <- param

    # Build the covariance matrix for P
    Sigma_P <- build_covariance_matrix_P(
        Theta$sigma2,
        Theta$rho_same,
        Theta$rho_diff,
        dims$K
    )

    # Log-likelihood of P (trunc MVN) under the hierarchical prior
    log_likelihood_P <- sum(sapply(1:dims$N, function(n) {

        # NOT dropping normalizing constants
        # 2/13/25 flat sigs --> maybe see  more rho_diff
        # SBS5 is flat sig
        # SBS1 C>T sig
        tmvtnorm::dtmvnorm(
            x = P[, n],  # nth column of P
            mean = Theta$Mu_p[, n],  # current mean for the nth signature
            sigma = Sigma_P,  # covar_P
            lower = rep(0, dims$K),
            upper = rep(Inf, dims$K),
            log = TRUE  # log density
        )
    }))

    print(param)



    # Get the prior distribution for the current parameter
    if (param_name == "sigma2") {
        log_prior <- invgamma::dinvgamma(
            param,
            shape = Theta$sigma2_prior["shape"],
            rate = Theta$sigma2_prior["scale"],
            log = TRUE
        )
    } else if (param_name == "rho_same") {
        log_prior <- dbeta(
            param,
            Theta$rho_same_prior["a"],
            Theta$rho_same_prior["b"],
            log = TRUE
        )
    } else if (param_name == "rho_diff") { # f(x) = 2x - 1 => f^-1(x) = (x+1)/2
        log_prior <- dbeta(
            (param+1)/2,
            Theta$rho_diff_prior["a"],
            Theta$rho_diff_prior["b"],
            log = TRUE
        )
    } else {
        stop("Invalid parameter name.")
    }

    print('returning from unnormalized_posterior')

    print(log_likelihood_P)

    print(log_prior)

    # Return the sum of the log-likelihood and the log-prior
    return(log_likelihood_P + log_prior)
}



#' #' Resample hyperparameters using ARMS (Adaptive Rejection Metropolis Sampling)
#' #'
#' #' @param Theta list of parameters
#' #' @param P current signature matrix (K x N)
#' #' @param dims list of dimensions
#' #'
#' #' @return list with updated values of sigma2, rho_same, and rho_diff
#' resample_hyperparameters <- function(Theta, P, dims) {
#'     # Helper function for resampling a single parameter
#'     resample_param <- function(param_name, lower, upper) {
#'         log_pdf <- function(x) {
#'             if (x < lower || x > upper) return(-Inf)  # Enforce bounds
#'             return(unnormalized_posterior(x, param_name, P, Theta, dims))
#'         }
#'         armspp::arms(
#'             n_samples = 1,
#'             log_pdf = log_pdf,
#'             lower = lower,
#'             upper = upper
#'         )
#'     }
#'
#'     # Resample each parameter
#'     sigma2_new <- resample_param("sigma2", lower = 1e-6, upper = 100)
#'     rho_same_new <- resample_param("rho_same", lower = -1, upper = 1)
#'     rho_diff_new <- resample_param("rho_diff", lower = -1, upper = 1)
#'
#'     # Update Theta and return
#'     Theta$sigma2 <- sigma2_new
#'     Theta$rho_same <- rho_same_new
#'     Theta$rho_diff <- rho_diff_new
#'     return(Theta)
#' }


# Using MH algorithm
resample_param_rw <- function(current, log_pdf, proposal_sd, lower, upper) {
    print('inside resample_param_rw')

    proposed <- rnorm(1, mean = current, sd = proposal_sd)
    if (proposed < lower || proposed > upper) {
        return(current)
    }
    log_p_current <- log_pdf(current) # unnormalized log posterior probability
    log_p_proposed <- log_pdf(proposed)
    acceptance_ratio <- exp(log_p_proposed - log_p_current)
    if (runif(1) < acceptance_ratio) {
        return(proposed)
    }
    return(current)
}


resample_hyperparameters <- function(Theta, P, dims) {
    print('inside resample_hyperparameters')

    # resampling a single parameter
    resample_param <- function(param_name, current_value, lower, upper) {
        log_pdf <- function(x) {
            if (is.na(x) || is.nan(x) || x < lower || x > upper) {
                return(-Inf)
            }
            value <- unnormalized_posterior(x, param_name, P, Theta, dims)
            cat("Param:", param_name, "Value:", x, "Log Posterior:", value, "\n")
            return(value)
        }
        # Attempt ARMS sampling
        tryCatch({
            armspp::arms(
                n_samples = 1,
                log_pdf = log_pdf,
                lower = lower,
                upper = upper,
                initial = current_value,
                max_iter = 10000  # Allow more iterations
            )
        }, error = function(e) {
            cat("ARMS failed for", param_name, "- falling back to random-walk.\n")

            # Fallback to random-walk Metropolis-Hastings
            # smaller sd for rho_diff
            proposal_sd_value <- if (param_name == "rho_diff") 0.1 else 0.1
            resample_param_rw(current_value, log_pdf, proposal_sd = proposal_sd_value, lower, upper)
        })
    }

    # Resample each parameter
    Theta$sigma2 <- resample_param("sigma2", Theta$sigma2, lower = 1e-6, upper = 100)
    Theta$rho_same <- resample_param("rho_same", Theta$rho_same, lower = 0, upper = 1)
    Theta$rho_diff <- resample_param("rho_diff", Theta$rho_diff, lower = -1, upper = 1)

    print('done resampling in resample_hyperparameters')

    return(Theta)
}

