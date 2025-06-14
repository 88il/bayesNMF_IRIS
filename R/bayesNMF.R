#' Bayesian Non-Negative Matrix Factorization
#' @description Bayesian Non-Negative Matrix Factorization.
#'
#' @param M matrix, data with samples as columns and features as rows
#' @param rank integer or vector, number of latent factors, integers
#' (e.g., \code{rank = 5}) or vector (e.g., \code{rank = 1:5}) are accepted.
#' If an integer is provided, rank is constant and will not be learned.
#' If a vector is provided, values must be sequential and start at 0 or 1
#' (e.g., \code{1:5} or \code{0:5}) for \code{learn_rank_method = "BFI"} or
#' \code{"SBFI"}, but can be non-sequential (e.g., \code{c(2,4,5)}) for
#' \code{learn_rank_method = "heuristic"}.
#' @param learn_rank_method string, method used to learn latent rank, one of
#' c("SBFI","BFI","heuristic"), used if a vector is provided for \code{rank}.
#' Bayesian Factor Inclusion (BFI) learns rank automatically as a part of the
#' Bayesian model with a binary Factor inclusion matrix. Sparse Bayesian
#' Factor inclusion (SBFI) is a variant of BFI with a sparse prior based on
#' the regularization of BIC. The heuristic approach
#' fits a fixed-rank Bayesian NMF for each value in \code{rank}, and selects
#' the model that minimizes BIC.
#' @param likelihood string, one of \code{c('normal','poisson')}, represents the
#' distribution used for likelihood f(M|P, E).
#' @param prior string, one of \code{c('truncnormal','exponential','gamma')},
#' represents the distribution used for priors on P and E, f(P) and f(E).
#' @param fast boolean, if \code{likelihood = 'poisson'} and \code{fast = TRUE},
#' then updates from the corresponding \code{likelihood = 'normal'} model
#' are used as proposals in an efficient Gibb's sampler. Only available for
#' \code{likelihood = 'poisson'} and \code{prior = c('truncnormal', 'exponential')}.
#' Defaults \code{TRUE} when possible.
#' @param inits (optional) named list, initial values for parameters, often P and E.
#' May also provide sigmasq if \code{likelihood = "normal"},
#' if \code{likelihood = "poisson"} and  \code{fast = TRUE}.
#' May provide initial values for A if \code{rank} is a vector.
#' @param fixed (optional) named list, parameters values to fix at constant
#' values, rather than learning them through Gibbs updates.
#' @param clip numeric, if \code{rank} is a vector and
#' \code{learn_rank_method = "BFI"} or \code{"SBFI"}, prior probabilities of
#' factor inclusion will be clipped by \code{clip}/N away from 0 and 1.
#' @param prior_parameters (optional) list, specification of prior parameters.
#' @param recovery boolean, if TRUE, allows for the addition of priors
#' based on previously discovered factors. In this case, \code{rank} denotes
#' the number of additional latent factors on top of those with priors in
#' \code{recovery_priors}.
#' @param recovery_priors "cosmic" or list, prior parameters for recovered
#' latent factors. If \code{recovery_priors = "cosmic"}, pre-computed priors
#' based on the 79 COSMIC cancer mutational signatures are used.
#' The function \code{get_recovery_priors} can be used to create recovery priors
#' given a reference factor matrix P and a likelihood-prior specification.
#' @param file string, file name (without extension) used for the log, rds,
#' and pdf files created by this function.
#' @param true_P (optional) matrix, reference latent factors matrix P to
#' compare estimated factors to with a heatmap. Not used in model, only
#' in evaluation.
#' @param convergence_control list, specification of convergence parameters.
#' See documentation for \code{new_convergence_control}.
#' @param store_logs boolean, if \code{store_logs = TRUE}, each iteration of the
#' Gibb's sampler is stored in resulting \code{.rds} file. Otherwise, only the
#' samples used to compute MAP are saved.
#' @param overwrite boolean, if \code{overwrite = TRUE}, the .log, .rds, and
#' .pdf files of previous runs with the same \code{file} will be overwritten.
#' This will permanently delete previous results, so be careful.
#'
#' @return list, also stored in the \code{.rds} file named by \code{file}
#' \itemize{
#'  \item \code{MAP}: list, maximum a-posteriori estimates of all parameters
#'  \item \code{credible_intervals}: list, lower and upper bounds for 95% credible intervals of all parameters
#'  \item \code{posterior_samples}: list, Gibb's samples used to compute MAP and credible intervals
#'  \item \code{metrics}: data.frame, contains MAP metrics at each convergence checkpoint, also plotted in the \code{.pdf} file named by \code{file}
#'  \item \code{model}: list, parameters of initial \code{bayesNMF} call
#'  \item \code{totaliters}: integer
#'  \item \code{converged_at}: integer, iteration determined as the convergence point
#'  \item \code{final_Theta}: list, parameter values of final iteration
#'  \item \code{time}: list, holds averge seconds per iteration and total seconds
#'  \item \code{logs}: list, only included if \code{store_logs = TRUE}, all Gibb's samples
#'  \item \code{sim_mat}: matrix, only included if \code{true_P} is provided, pairwise similarities between estimated and true factors
#'  \item \code{heatmap}: ggplot2 object, only included if \code{true_P} is provided, heatmap of pairwise similarities between estimated and true factors
#' }
#' @export
bayesNMF <- function(
    M, rank,
    range_N = NULL,
    learn_rank_method = "SBFI",
    likelihood = "poisson",
    prior = "truncnormal",
    fast = likelihood == "poisson" & prior %in% c('truncnormal','exponential'),
    inits = NULL,
    fixed = NULL,
    clip = 0.4,
    prior_parameters = NULL,
    recovery = FALSE,
    recovery_priors = "cosmic",
    file = paste0('nmf_', likelihood, '_', prior),
    true_P = NULL,
    convergence_control = new_convergence_control(),
    store_logs = TRUE,
    overwrite = FALSE
) {
    # set up file names
    final_file = file
    savefile = paste0(file, '.rds')
    logfile = paste0(file, '.log')
    tail = 0
    while (!overwrite & (file.exists(savefile) | file.exists(logfile))) {
        tail = tail + 1
        final_file = paste0(file, '_', tail)
        savefile = paste0(file, '_', tail, '.rds')
        logfile = paste0(file, '_', tail, '.log')
    }

    # if learning rank, consider "learn_rank_method"
    learn_A <- length(rank) > 1
    if (learn_A) {
        ################################################################
        # --------- Sparse Bayesian Factor Inclusion (SBFI) -----------#
        ################################################################
        if (learn_rank_method == "SBFI") {
            max_N = max(rank)
            min_N = min(rank)

            # rank must be 0/1:max for SBFI
            if (!setequal(min_N:max_N, rank)) {
                stop(paste(
                    "rank =", paste0("c(", paste(rank, collapse = ','), ")"),
                    "is not sequential. rank must be sequential",
                    "to use learn_rank_method = 'BFI'. Try a sequential rank",
                    "or learn_rank_method = 'heuristic'."
                ))
            }
            if (min_N > 1) {
                stop(paste(
                    "Minimum rank > 1 is not permitted with learn_rank_method =",
                    "'BFI'. Try rank =", paste0(1, ":", max_N),
                    "or use learn_rank_method = 'heuristic'."
                ))
            }

            # run bayesNMF
            tryCatch({
                sink(file = logfile)
                res <- inner_bayesNMF(
                    M = M,
                    N = NULL,
                    max_N = max_N,
                    range_N = range_N,
                    sparse_rank = TRUE, # updates A with BIC
                    likelihood = likelihood,
                    prior = prior,
                    fast = fast,
                    inits = inits,
                    fixed = fixed,
                    clip = clip,
                    prior_parameters = prior_parameters,
                    recovery = recovery,
                    recovery_priors = recovery_priors,
                    file = final_file,
                    true_P = true_P,
                    convergence_control = convergence_control,
                    store_logs = store_logs,
                    overwrite = overwrite
                )
                sink()
            }, error = function(e) {
                sink()
                traceback()
                stop(e)
            }, interrupt = function(e) {
                sink()
                stop("Interrupted", call. = FALSE)
            })

        ################################################################
        # ------------- Bayesian Factor Inclusion (BFI) ---------------#
        ################################################################
        } else if (learn_rank_method == "BFI") {
            max_N = max(rank)
            min_N = min(rank)

            # rank must be 0/1:max for BFI
            if (!setequal(min_N:max_N, rank)) {
                stop(paste(
                    "rank =", paste0("c(", paste(rank, collapse = ','), ")"),
                    "is not sequential. rank must be sequential",
                    "to use learn_rank_method = 'BFI'. Try a sequential rank",
                    "or learn_rank_method = 'heuristic'."
                ))
            }
            if (min_N > 1) {
                stop(paste(
                    "Minimum rank > 1 is not permitted with learn_rank_method =",
                    "'BFI'. Try rank =", paste0(1, ":", max_N),
                    "or use learn_rank_method = 'heuristic'."
                ))
            }

            # run bayesNMF
            tryCatch({
                sink(file = logfile)
                res <- inner_bayesNMF(
                    M = M,
                    N = NULL,
                    max_N = max_N,
                    range_N = range_N,
                    sparse_rank = FALSE, # A updated based on log likelihood
                    likelihood = likelihood,
                    prior = prior,
                    fast = fast,
                    inits = inits,
                    fixed = fixed,
                    clip = clip,
                    prior_parameters = prior_parameters,
                    recovery = recovery,
                    recovery_priors = recovery_priors,
                    file = final_file,
                    true_P = true_P,
                    convergence_control = convergence_control,
                    store_logs = store_logs,
                    overwrite = overwrite
                )
                sink()
            }, error = function(e) {
                sink()
                traceback()
                stop(e)
            }, interrupt = function(e) {
                sink()
                stop("Interrupted", call. = FALSE)
            })

        ################################################################
        # -------------- Heuristic Approach (min BIC) -----------------#
        ################################################################
        } else if (learn_rank_method == "heuristic") {

            # set up storage for results
            BICs <- data.frame(
                rank = rank,
                BIC = rep(NA, length(rank))
            )
            all_models <- list()

            for (r in rank) {
                # fit a model for each rank
                tryCatch({
                    this_file = paste0(final_file, "_rank", r)
                    sink(file = paste0(this_file, '.log'))
                    res_N <- inner_bayesNMF(
                        M = M, N = r, likelihood = likelihood, prior = prior,
                        fast = fast, inits = inits, fixed = fixed, clip = clip,
                        prior_parameters = prior_parameters,
                        recovery = recovery, recovery_priors = recovery_priors,
                        file = this_file, true_P = true_P,
                        convergence_control = convergence_control,
                        store_logs = store_logs, overwrite = overwrite
                    )
                    sink()
                }, error = function(e) {
                    sink()
                    traceback()
                    stop(e)
                }, interrupt = function(e) {
                    sink()
                    stop("Interrupted", call. = FALSE)
                })

                # log results of this run
                BICs$BIC[BICs$rank == r] <- res_N$metrics$BIC[
                    res_N$metrics$sample_idx == res_N$converged_at
                ]
                all_models[[r]] <- res_N

                # update and save all results
                best_rank <- BICs %>%
                    dplyr::filter(BIC == min(BIC, na.rm = TRUE)) %>%
                    dplyr::pull(rank)
                plot <- BICs %>%
                    ggplot2::ggplot(ggplot2::aes(x = rank, y = BIC)) +
                    ggplot2::geom_point() +
                    ggplot2::geom_line()
                res <- list(
                    best_rank = best_rank,
                    best_model = all_models[[best_rank]],
                    BICs = BICs,
                    all_models = all_models,
                    plot = plot
                )
                saveRDS(res, file = savefile)
            }
        } else {
            stop(paste("learn_rank_method =", paste0("'", learn_rank_method, "'"),
                       "not defined,", "must be one of c('BFI','heuristic')."))
        }
    } else {
        ################################################################
        # ---------------------- Fixed Rank ---------------------------#
        ################################################################
        tryCatch({
            sink(file = logfile)
            res <- inner_bayesNMF(
                M = M,
                N = rank,
                max_N = NULL,
                likelihood = likelihood,
                prior = prior,
                fast = fast,
                inits = inits,
                fixed = fixed,
                clip = clip,
                prior_parameters = prior_parameters,
                recovery = recovery,
                recovery_priors = recovery_priors,
                file = final_file,
                true_P = true_P,
                convergence_control = convergence_control,
                store_logs = store_logs,
                overwrite = overwrite
            )
            sink()
        }, error = function(e) {
            sink()
            traceback()
            stop(e)
        })
    }
    return(res)
}


#' Bayesian Non-Negative Matrix Factorization
#' @description Perform single-study Bayesian NMF with the provided likelihood and prior
#' combination. Exact rank `N` or maximum rank `max_N` must be provided.
#'
#' @param M mutational catalog matrix, K x G
#' @param N fixed number of latent factors
#' @param max_N maximum number of latent factors if learning rank
#' @param likelihood string, one of c('normal','poisson')
#' @param prior string, one of c('truncnormal','exponential','gamma')
#' @param fast boolean, if `likelihood == 'poisson'` and `fast = TRUE`, updates
#' from the corresponding `likelihood == 'normal'` model are used as proposals
#' in an efficient Gibb's sampler. Defaults TRUE when possible.
#' @param inits (optional) list of initial values for P and E (and sigmasq
#' if `likelihood = "normal"`)
#' @param fixed (optional) list of parameters to fix and not include in Gibbs
#' updates.
#' @param clip numeric, prior probabilities of inclusion will be clipped by
#' `clip`/N away from 0 and 1
#' @param prior_parameters list, optional specification of prior parameters
#' @param recovery boolean, whether to set priors of a subset of factors at
#' previously discovered factors
#' @param recovery_priors "cosmic" or list of prior parameters
#' @param file file name without extension of log, save, and plot files
#' @param true_P (optional) true latent factors matrix P to compare to in a heatmap
#' @param convergence_control list, specification of convergence parameters.
#' See documentation for `new_convergence_control`.
#' @param store_logs boolean, whether to store each sample in resulting .rds file
#' @param overwrite if `overwrite = TRUE`, the log, safe, and plot files of
#' previous runs with the same `file` will be overwritten
#'
#' @return list
#' @noRd
inner_bayesNMF <- function(
        M,
        N = NULL,
        max_N = NULL,
        range_N = NULL,
        sparse_rank = FALSE,
        likelihood = "poisson",
        prior = "truncnormal",
        fast = (likelihood == "poisson" &
                prior %in% c('truncnormal','exponential')),
        inits = NULL,
        fixed = NULL,
        clip = 0.4,
        prior_parameters = NULL,
        recovery = FALSE,
        recovery_priors = "cosmic",
        file = paste0('nmf_', likelihood, '_', prior),
        true_P = NULL,
        convergence_control = new_convergence_control(
            maxiters = ifelse(recovery, 5000, 2000)
        ),
        store_logs = TRUE,
        overwrite = FALSE
) {
    savefile = paste0(file, '.rds')
    plotfile = paste0(file, '.pdf')

    START = Sys.time()
    print(START)
    print(paste("maxiters =", convergence_control$maxiters))
    print(paste("fast", fast))

    # check recovery/discovery
    if (recovery) {
        if (is.character(recovery_priors)) {
            if (recovery_priors == "cosmic") {
                if (prior == 'truncnormal') {
                    recovery_priors <- normal_truncnormal_recovery_priors
                } else {
                    stop("Recovery priors not defined for prior combination")
                }
            }
        }
    } else {
        recovery_priors = list()
    }

    # check N/max_N combination is valid
    N <- validate_N(N, max_N, recovery_priors)
    if (is.null(range_N)) {
        range_N = 0:N
    }
    if (recovery) {
        scale_by <- sqrt(mean(M)/N) / mean(recovery_priors$Mu_p)
        recovery_priors$Mu_p <- recovery_priors$Mu_p * scale_by
        recovery_priors$Sigmasq_p <- recovery_priors$Sigmasq_p * (scale_by**2)
    }

    # check prior and likelihood are valid
    validate_model(likelihood, prior, fast)

    # set up tempering schedule
    learn_A <- !is.null(max_N)
    if (learn_A) {
        gamma_sched <- get_gamma_sched(
            len = convergence_control$maxiters,
            n_temp = round(0.5 * convergence_control$maxiters)
        )
    } else {
        gamma_sched <- rep(1, convergence_control$maxiters)
    }
    print(paste("learn_A",learn_A))

    # set up dimensions
    dims = list(
        K = dim(M)[1],
        G = dim(M)[2],
        N = N,
        S = 1
    )

    # set up Theta
    print("initializing Theta")
    Theta <- initialize_Theta(
        M = M,
        likelihood = likelihood,
        prior = prior,
        fast = fast,
        learn_A = learn_A,
        dims = dims,
        inits = inits, fixed = fixed,
        prior_parameters = prior_parameters,
        recovery = recovery,
        recovery_priors = recovery_priors,
        clip = clip,
        range_N = range_N
    )
    Theta$prob_inclusion <- Theta$A
    Theta$P_acceptance <- Theta$P
    Theta$E_acceptance <- Theta$E

    print("done initializing Theta")

    # set up metrics and logs
    metrics <- list(
        sample_idx = list(),
        RMSE = list(),
        KL = list(),
        loglik = list(),
        logpost = list(),
        N = list(),
        n_params = list(),
        MAP_A_counts = list(),
        BIC = list()
    )

    logs <- list(
        P = list(),
        E = list(),
        P_acceptance = list(),
        E_acceptance = list(),
        A = list(),
        prob_inclusion = list(),
        n = list(),

        logprior_every = list(), # IL added 4/25
        loglik_every = list(),
        logpost_every = list()
    )
    if (likelihood == "normal" | (likelihood == "poisson" & fast)) {
        logs$sigmasq <- list()
    } else {
        # likelihood == "poisson" & !fast
        logs$Z <- list()
    }

    if(!is.null(Theta$sigma2_prior)) {
        logs$sigma2 <- list()
        logs$rho_same <- list()
        logs$rho_diff <- list()
    }

    # initialize convergence status
    convergence_status <- list(
        prev_MAP_metric = Inf,
        best_MAP_metric = Inf,
        inarow_no_change = 0,
        inarow_no_best = 0,
        converged = FALSE,
        best_iter = NULL
    )

    PREV = Sys.time()
    print(paste("starting iterations,", PREV))
    avg_time = 0

    # avoid undesirable effect of sample(x) when x is a single value
    update_P_columns = which(!Theta$is_fixed$P)
    if (length(update_P_columns) > 1) {
        update_P_columns = sample(update_P_columns)
    }
    update_A_cols = which(!Theta$is_fixed$A)
    if (length(update_A_cols) > 1) {
        update_A_cols = sample(update_A_cols)
    }

    # print("before Gibbs starts printing Theta$P")
    # print(Theta$P)

    # Gibbs sampler
    iter = 1
    logiter = 1
    done = FALSE
    first_MAP = TRUE
    stop = NULL
    START_ITER <- Sys.time()
    while (iter <= convergence_control$maxiters & !done) {
        print(iter)

        # update non-fixed columns of P
        if (length(update_P_columns) > 1) {
            update_P_columns = sample(update_P_columns)
        }

        for (n in update_P_columns) {
            sample_Pn_out <- sample_Pn(
                n, M, Theta, dims,
                likelihood, prior, fast
            )
            Theta$P[, n] <- sample_Pn_out$sampled
            Theta$P_acceptance[, n] <- sample_Pn_out$acceptance
        }

        # update E
        if (!Theta$is_fixed$E) {
            for (n in sample(1:dims$N)) {
                sample_En_out <- sample_En(
                    n, M, Theta, dims,
                    likelihood, prior, fast
                )
                Theta$E[n, ] <- sample_En_out$sampled
                Theta$E_acceptance[n, ] <- sample_En_out$acceptance
            }
        }

        # if Normal likelihood, update sigmasq
        if (likelihood == 'normal' | (likelihood == 'poisson' & fast)) {
            if (!Theta$is_fixed$sigmasq) {
                Theta$sigmasq <- sample_sigmasq_normal(
                    M, Theta, dims, gamma = 1
                )
            }
        }

        # update A and n
        Theta$n <- sample_n(Theta, dims, clip, gamma = gamma_sched[iter])
        Theta$q <- update_q(Theta, dims, clip)
        if (length(update_A_cols) > 1) {
            update_A_cols = sample(update_A_cols)
        }
        for (n in update_A_cols) {
            sample_An_out <- sample_An(
                n, M, Theta, dims,
                likelihood, prior,
                sparse_rank = sparse_rank,
                gamma = gamma_sched[iter]
            )
            Theta$A[1, n] <- sample_An_out$sampled
            Theta$prob_inclusion[1,n] <- sample_An_out$prob_inclusion
        }

        # if Poisson likelihood, update latent counts Z
        if (likelihood == 'poisson' & !fast) {
            for (k in sample(1:dims$K)) {
                for (g in sample(1:dims$G)) {
                    Theta$Z[k,,g] <- sample_Zkg_poisson(
                        k, g, M, Theta, dims,
                        gamma = 1
                    )
                }
            }
        }

        # update priors
        for (n in 1:dims$N) {
            if (prior == "truncnormal") {
                if (!Theta$is_fixed$prior_P[n]) {
                    # Update mean and variance for P
                    Theta$Mu_p[,n] <- sample_Mu_Pn(
                        n, Theta, dims, gamma = 1
                    )
                    Theta$Sigmasq_p[,n] <- sample_Sigmasq_Pn(
                        n, Theta, dims, gamma = 1
                    )
                }

                # Update mean and variance for E
                Theta$Mu_e[n,] <- sample_Mu_En(
                    n, Theta, dims, gamma = 1
                )
                Theta$Sigmasq_e[n,] <- sample_Sigmasq_En(
                    n, Theta, dims, gamma = 1
                )
            } else if (prior == "exponential") {
                if (!Theta$is_fixed$prior_P[n]) {
                    Theta$Lambda_p[,n] <- sample_Lambda_Pn(
                        n, Theta, dims, gamma = 1
                    )
                }
                Theta$Lambda_e[n,] <- sample_Lambda_En(
                    n, Theta, dims, gamma = 1
                )
            } else if (prior == "gamma") {
                if (!Theta$is_fixed$prior_P[n]) {
                    Theta$Beta_p[,n] <- sample_Beta_Pn(
                        n, Theta, dims, gamma = 1
                    )
                    for (k in 1:dims$K) {
                        Theta$Alpha_p[k,n] <- sample_Alpha_Pkn(
                            k, n, Theta, dims, gamma = 1
                        )
                    }
                }
                Theta$Beta_e[n,] <- sample_Beta_En(
                    n, Theta, dims, gamma = 1
                )
                for (g in 1:dims$G) {
                    Theta$Alpha_e[n,g] <- sample_Alpha_Eng(
                        n, g, Theta, dims, gamma = 1
                    )
                }
            }
        }


        # update covar_P hyperpriors for Hierarchical method
        if (!is.null(Theta$sigma2_prior) & !is.null(Theta$rho_same_prior) & !is.null(Theta$rho_diff_prior)) {

            hyperparam_resample <- resample_hyperparameters(
                Theta, Theta$P, dims
            )
            Theta$sigma2 <- hyperparam_resample$sigma2
            Theta$rho_same <- hyperparam_resample$rho_same
            Theta$rho_diff <- hyperparam_resample$rho_diff

            # recompute Covar_p based on updated hyperparameters
            Theta$Covar_p <- build_covariance_matrix_P(
                Theta$sigma2,
                Theta$rho_same,
                Theta$rho_diff,
                dims$K
            )
        }


        # log on original scale
        # only if storing logs or if we will use it for MAP
        if (store_logs | iter >= convergence_control$MAP_every + 1) {
            logs$P[[logiter]] <- Theta$P
            logs$E[[logiter]] <- Theta$E
            logs$P_acceptance[[logiter]] <- Theta$P_acceptance
            logs$E_acceptance[[logiter]] <- Theta$E_acceptance
            logs$A[[logiter]] <- Theta$A
            logs$n[[logiter]] <- Theta$n
            logs$prob_inclusion[[logiter]] <- Theta$prob_inclusion
            if (likelihood == "normal" | (likelihood == "poisson" & fast)) {
                logs$sigmasq[[logiter]] <- Theta$sigmasq
            } else {
                # likelihood == "poisson" & !fast
                logs$Z[[logiter]] <- Theta$Z
            }

            # store hierarchical parameters logs
            if(!is.null(Theta$sigma2_prior)) {
                logs$sigma2[[logiter]] <- Theta$sigma2
                logs$rho_same[[logiter]] <- Theta$rho_same
                logs$rho_diff[[logiter]] <- Theta$rho_diff
            }

            # track logposterior and logprior
            logs$logprior_every[[logiter]] <- get_logprior(
                Theta, likelihood, prior, dims
            )
            logs$loglik_every[[logiter]] <- get_loglik_normal(M, Theta, dims)
            logs$logpost_every[[logiter]] <- logs$loglik_every[[logiter]] + logs$logprior_every[[logiter]]

            logiter = logiter + 1
        }

        # periodically check convergence and log progress
        if (
            (iter %% convergence_control$MAP_every == 0 &
             iter >= convergence_control$MAP_over + convergence_control$MAP_every)
            | iter == convergence_control$maxiters
        ) {
            # get MAP over past convergence_control$MAP_over iterations
            burn_in <- logiter - convergence_control$MAP_over
            keep <- burn_in:(logiter - 1)
            MAP <- get_MAP(logs, keep, dims)

            # log metrics
            out <- update_metrics(
                metrics, MAP, iter, Theta, M,
                likelihood, prior, dims
            )
            metrics <- out$metrics
            Theta_MAP <- out$Theta_MAP
            Mhat_MAP <- out$Mhat_MAP

            # check convergence
            convergence_status <- check_converged(
                iter, gamma_sched[iter],
                Mhat_MAP, M,
                convergence_status,
                convergence_control,
                first_MAP,
                Theta = Theta_MAP,
                likelihood = likelihood,
                prior = prior,
                dims = dims
            )

            # check whether convergence is allowed (i.e., whether tempering is over)
            if (gamma_sched[iter] == 1 & first_MAP & MAP$top_counts[1] >= convergence_control$minA) {
                first_MAP = FALSE
                # forces convergence after gamma == 1
                convergence_status$best_MAP_metric = Inf
            }

            # if not storing logs, store current best MAP to return if needed
            # and remove logs to save memory
            if (
                !store_logs &
                !is.null(convergence_status$best_iter)
            ) {
                # update stored MAP if this is the best_iter
                if (convergence_status$best_iter == iter) {
                    store_MAP <- MAP
                    store_MAP_iter <- iter
                    store_credible_intervals <- get_credible_intervals(logs, store_MAP$idx)
                }

                # remove logs to save memory
                drop <- 1:convergence_control$MAP_every
                for (name in names(logs)) {
                    logs[[name]][drop] <- NULL
                }
                logiter = logiter - length(drop)
            }

            # log progress
            NOW = Sys.time()
            diff = as.numeric(difftime(NOW, PREV, units = "secs"))
            PREV = NOW
            log_MAP(
                iter, done, diff, convergence_control,
                convergence_status, gamma_sched, MAP, learn_A
            )

            # if converged, stop
            if (convergence_status$converged){
                stop = convergence_status$best_iter
                log_converged(convergence_control, convergence_status)
                done = TRUE

                # re-compute MAP at stop, compute 95% credible intervals
                if (store_logs) {
                    burn_in <- stop - convergence_control$MAP_over
                    keep <- burn_in:stop
                    MAP <- get_MAP(logs, keep, dims, final = TRUE)
                    credible_intervals <- get_credible_intervals(logs, MAP$idx)
                } else {
                    MAP <- store_MAP
                    keep_sigs <- which(MAP$A[1,] == 1)
                    if (dims$N > 1) {
                        MAP$P <- MAP$P[, keep_sigs]
                        MAP$E <- MAP$E[keep_sigs, ]
                    } else {
                        if (length(keep_sigs) == 0) {
                            MAP$P <- matrix(0, nrow = dims$K, ncol = 1)
                            MAP$E <- matrix(0, nrow = 1, ncol = dims$G)
                        }
                    }
                    credible_intervals <- store_credible_intervals
                }
            }

            # plot metrics
            plot_metrics(metrics, plotfile, stop, learn_A, gamma_sched, iter, true_P)

            # save results
            metrics_df = data.frame(matrix(nrow = length(unlist(metrics$BIC)), ncol = 0))
            for (name in names(metrics)) {
                metrics_df[,name] <- unlist(metrics[[name]])
            }
            keep_sigs = as.logical(MAP$A[1, ])
            res <- list(
                MAP = MAP,
                metrics = metrics_df,
                model = list(
                    M = M,
                    true_P = true_P,
                    likelihood = likelihood,
                    prior = prior,
                    prior_parameters = prior_parameters,
                    fixed = fixed,
                    inits = inits,
                    convergence_control = convergence_control,
                    dims = dims,
                    gamma_sched = gamma_sched
                ),
                totaliters = iter,
                converged_at = stop,
                final_Theta = Theta,
                time = list(
                    "avg_secs_per_iter" = as.numeric(difftime(Sys.time(), START_ITER, units = "secs"))/iter,
                    "total_secs" = as.numeric(difftime(Sys.time(), START, units = "secs"))
                )
            )
            if (done) {
                credible_intervals$P[[1]] <- credible_intervals$P[[1]][,res$MAP$A[1,]==1]
                credible_intervals$P[[2]] <- credible_intervals$P[[2]][,res$MAP$A[1,]==1]
                credible_intervals$E[[1]] <- credible_intervals$E[[1]][res$MAP$A[1,]==1,]
                credible_intervals$E[[2]] <- credible_intervals$E[[2]][res$MAP$A[1,]==1,]

                res$credible_intervals <- credible_intervals
                # IL 2/6 returns last 1000 iterations
                posterior_samples <- list()
                for (name in names(logs)) {
                    posterior_samples[[name]] <- logs[[name]][MAP$idx]
                }
                res$posterior_samples <- posterior_samples
            }
            if (store_logs) {
                res$logs = logs
            }
            class(res) = 'bayesNMF'
            saveRDS(res, file = savefile)
        }

        iter = iter + 1
    }

    # similarity matrix with true_P, if provided
    if (!is.null(true_P)) {
        sim_mat <- pairwise_sim(
            res$MAP$P, true_P,
            name1 = "estimated", name2 = "true",
            which = "cols"
        )
        heatmap <- get_heatmap(res$MAP$P, true_P)

        res$sim_mat <- sim_mat
        res$heatmap <- heatmap
        class(res) = 'bayesNMF'
        saveRDS(res, file = savefile)
    }

    return(res)
}
