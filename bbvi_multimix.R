require(mvtnorm)
require(randtoolbox)
require(distr)
require(LaplacesDemon)

#' Random Sample Generation
#'
#' @param n Size of mixture sample to be generated
#' @param props An array of K mixing proportions per component, can be unnormalized
#' @param means a list of p-dimensional vectors of the means
#' @param covars a list of pxp covariance matrices
#' @param seed Random seed (defaults to 923)
#' 
#' @return An array of n observations corresponding to the Gaussian mixture distribution
generate_mixture <- function(n, props, means, covars, seed=923){
    set.seed(seed)
    props <- props/sum(props)
    sizes <- round(n * props, 0)
    
    samp <- NULL
    cls <- NULL
    for(i in seq_along(sizes)){
        cls <- c(cls,rep(i, times = sizes[i]))
        samp <- rbind(samp, rmvnorm(sizes[i], mean=means[[i]], sigma=covars[[i]]))
    }

    return(list(samp=samp, sizes=sizes, Cls = cls))
}

#' Cartesian Product of Indices
#' 
#' For use in naming parameters when initializing the method.
#'
#' @param seq1 A sequence of indices
#' @param seq2 A sequence of indices
#' 
#' @return An array of strings containing cartesian products of seq1 and seq2
cartesian_index <- function(seq1, seq2){
    combs <- expand.grid(seq1,seq2)
    sapply(1:nrow(combs), function(x){
        return(paste0(combs[x,1], "_" ,combs[x,2]))
    })
}

#' Unpack List of Prior Distributions
#' 
#' For use in initializing the method.
#'
#' @param priors a list of prior specifications
#' @param K the number of mixture components
#' @param p the dimension of the Gaussian distributions
#' @param N the number of observations
#' 
#' @return An updated list of prior specifications, completed in case of missing priors
unpack_priors <- function(priors, K, p, N){
    
    if(is.null(priors[["m0"]])){
        priors[["m0"]] <- rep(0, times = p)
    }

    if(is.null(priors[["S0"]])){
        priors[["S0"]] <- diag(1, nrow = p, ncol = p)
    }

    if(is.null(priors[["Sigma"]])){
        priors[["Sigma"]] <- diag(1, nrow = p, ncol = p)
    }

    if(is.null(priors[["Phi0"]])){
        priors[["Phi0"]] <- rep(1/K, times = K)
    }

    return(priors)
}

#' Fit a Multivariate Mixture of Gaussian Using BBVI
#' 
#' Fits a set of K (pre-specified) components each a p-dimensional Gaussian distribution to a dataset.
#' The method used Black Box Variational Inference using either its Naive, Rao-Blackwellized, or James-Stein
#' flavor.
#'
#' @param data an N x p matrix of observations
#' @param clusters the number K of components
#' @param method the flavor of BBVI to use in model fitting. Can be one of method = c("Naive","RB","JS+") for
#' the Naive, Rao-Blackwellized, and James-Stein flavor, respectively. Defaults to method = "Naive"
#' @param learn_rate the learning rate to use for the stochastic optimization. Can be rate_constant(rho) for a
#' certain value of rho, rate_adagrad(eta) for a value of eta to use the AdaGrad method, and rate_rmsprop(eta, beta)
#' for the RMSProp modification of AdaGrad. Defaults to rate_constant(1e-3)
#' @param max_iter the total number of iterations to stop the method in case of non-convergence
#' @param min_iter the minimum number of iterations to perform before convergence is assessed
#' @param mc_size the number S of monte carlo samples to draw for generating the BBVI Update steps
#' @param priors a list of prior specifications, defaults to standard priors if priors = NULL.
#' @param converge the convergence threshold delta
#' @param criterion the criterion in which to test for convergence, can be one of criterion = c("param","elbo").
#' @param var_threshold the lower and upper threshold on the variance, defaults to var_threshold = c(0.01, 1e3).
#' @param verbose logical, a switch for whether a verbose execution should be done.
#' @param seed the seed to use for starting the random samples
#' 
#' @return An updated list of prior specifications, completed in case of missing priors
bbvi_multimix <- function(data, clusters = 2, method = "Naive", learn_rate = rate_constant(), mc_size = 1000, max_iter = 1000, min_iter = 1, priors = NULL, seed = 923,
                            converge = 1e-4, crit = "param", var_threshold = c(0.01, 1e3), verbose = FALSE){
    p <- ncol(data) # number of components in the data
    N <- nrow(data) # number of rows in the data
    K <- clusters

    priors <- unpack_priors(priors, K, p, N)

    seq_K <- 1:K
    seq_p <- 1:p
    seq_N <- 1:N

    # the total number of parameters given by the following:
    # - K x p mean parameters m_k
    # - K x p (diagonal) variance matrix parameters S_k
    # - N x K mixture probabilities per observation of data
    mean_params <- paste0("m_", cartesian_index(seq_K, seq_p))
    var_params  <- paste0("ssq_", cartesian_index(seq_K, seq_p))
    prob_params <- paste0("phi_", cartesian_index(seq_N, seq_K))

    # initialize these variables
    # - setting m to a random location on the dataset
    # - setting ssq to 1 for all
    # - setting 1/K for all rows of data
    set.seed(923)
    init_locs <- as.matrix(data[sample(1:N, size = K), ])

    params <- list(max_iter)
    elbo <- list(max_iter - 1)
    steps <- list(max_iter - 1)

    params[[1]] <- data.frame(t(c(
        as.vector(init_locs),
        rep(1, times = length(var_params)),
        rep(1/K, times = length(prob_params))
    )))

    param_names <- c(mean_params, var_params, prob_params)
    colnames(params[[1]]) <- param_names
    G <- 0
    converged <- FALSE
    start.time <- Sys.time()

    set.seed(seed)
    for(t in 2:max_iter){
        old.params <- params[[t - 1]]

        samps <- generate_samples(mc_size, data, K, old.params, priors, method)
        updates <- generate_bbvi(samps)
        
        #G <- G + elbo %*% t(elbo) # for AdaGrad
        #G <- beta * G + (1 - beta) * (elbo %*% t(elbo)) # for AdaGrad
        updater <- learn_rate(t, old = old.params, step = updates, G)
        new.params <- updater$new
        G <- updater$G

        # hard thresholding on ssq
        for(j in which(grepl("ssq_", param_names, fixed = TRUE))){
            if(new.params[1,j] < var_threshold[1]) new.params[1,j] <- var_threshold[1]
            if(new.params[1,j] > var_threshold[2]) new.params[1,j] <- var_threshold[2]
        }
        

        # normalizing phi
        mu.new <- as.numeric(new.params[1, grepl("m_", param_names, fixed = TRUE)])
        mu.new <- matrix(mu.new, nrow = K, ncol = p, byrow = FALSE)

        phi.new <- as.numeric(new.params[1,grepl("phi_", param_names, fixed = TRUE)])
        phi.new <- matrix(phi.new, nrow = N, ncol = K, byrow = FALSE)
        
        for(i in 1:N){
            probs <- sapply(as.numeric(phi.new[i,]), function(x){
                if(x < 0.01) return(0.01)
                if(is.na(x) || is.infinite(x) || is.nan(x)) return(0.01)
                return(x)
            })
            total_prob <- sum(probs, na.rm = TRUE)
            phi.new[i,] <- probs/total_prob
            
        }

        new.params[1, grepl("phi_", param_names, fixed = TRUE)] <- as.vector(phi.new)

        bayesLik <- bayes.logLik(data, mu.new, phi.new, priors)
        pDIC <- 2 * (bayesLik - samps$simLikelihood)
        elpd <- bayesLik - pDIC
        DIC <- -2 * bayesLik + pDIC

        params[[t]] <- new.params
        elbo[[t-1]] <- data.frame(iter = t, elapsed = difftime(Sys.time(), start.time, units = "secs"), elbo = samps$ELBO,
                                    bayesLikelihood = bayesLik, simLikelihood = samps$simLikelihood, elpd = elpd, DIC = DIC)
        steps[[t-1]] <- updates

        # check for convergence
        if (t > 2){
            elbo.chg <- log(abs(samps$ELBO/elbo[[t-1]]$elbo[1]))
        }else{
            elbo.chg <- 1.0000
        }
        
        norm_old <- sum(old.params[1,]^2)
        norm_chg <- sum((new.params[1,] - old.params[1,])^2)
        lambda <- sqrt(norm_chg)/sqrt(norm_old)

        if((t > min_iter & crit == "param" & lambda < converge) | (t > min_iter & crit == "elbo" & elbo.chg < converge)){
            message(paste0("Algorithm converged after ", t, " iterations."))
            converged <- TRUE
            params[(t+1):max_iter] <- NULL
            elbo[t:(max_iter - 1)] <- NULL
            steps[t:(max_iter - 1)] <- NULL
            break
        }

        if(t%%100 == 0 & verbose) message(paste0("BBVI-",method,": Iteration ", t,
                                                    " | lambda: ", round(lambda, 4),
                                                    " | ELBO: ", round(samps$ELBO,2),
                                                    " | ELBO Change: ", round(elbo.chg, 4)))
    }

    if(! converged) message(paste0("Algorithm did not converge after ", t, " steps. Results may not be reliable."))
    
    paramlist <- do.call('rbind', params)
    elbolist <- do.call('rbind', elbo)
    steplist <- do.call('rbind', steps)

    return(list(
        trace = paramlist,
        elbo = elbolist,
        updates = steplist,
        iterations = t,
        converged = converged,
        elapsed.time = difftime(Sys.time(), start.time, units = "secs"),
        N = N,
        K = K,
        p = p,
        data = data,
        method = method
    ))
}

#' Fit a Multivariate Mixture of Gaussians Using QMC-BBVI
#' 
#' Fits a set of K (pre-specified) components each a p-dimensional Gaussian distribution to a dataset.
#' The method used Black Box Variational Inference using either its Naive, Rao-Blackwellized, or James-Stein
#' flavor.
#'
#' @param data an N x p matrix of observations
#' @param clusters the number K of components
#' @param method the flavor of BBVI to use in model fitting. Can be one of method = c("Naive","RB","JS+") for
#' the Naive, Rao-Blackwellized, and James-Stein flavor, respectively. Defaults to method = "Naive"
#' @param learn_rate the learning rate to use for the stochastic optimization. Can be rate_constant(rho) for a
#' certain value of rho, rate_adagrad(eta) for a value of eta to use the AdaGrad method, and rate_rmsprop(eta, beta)
#' for the RMSProp modification of AdaGrad. Defaults to rate_constant(1e-3)
#' @param max_iter the total number of iterations to stop the method in case of non-convergence
#' @param min_iter the minimum number of iterations to perform before convergence is assessed
#' @param mc_size the number S of monte carlo samples to draw for generating the BBVI Update steps
#' @param priors a list of prior specifications, defaults to standard priors if priors = NULL.
#' @param converge the convergence threshold delta
#' @param criterion the criterion in which to test for convergence, can be one of criterion = c("param","elbo").
#' @param var_threshold the lower and upper threshold on the variance, defaults to var_threshold = c(0.01, 1e3).
#' @param verbose logical, a switch for whether a verbose execution should be done.
#' @param seed the seed to use for starting the random samples
#' 
#' @return An updated list of prior specifications, completed in case of missing priors
qmcvi_multimix <- function(data, clusters = 2, method = "Naive", learn_rate = rate_constant(), mc_size = 1000, max_iter = 1000, min_iter = 1, priors = NULL, seed = 923,
                            converge = 1e-4, crit = "param", var_threshold = c(0.01, 1e3), verbose = FALSE){
    p <- ncol(data) # number of components in the data
    N <- nrow(data) # number of rows in the data
    K <- clusters

    priors <- unpack_priors(priors, K, p, N)

    seq_K <- 1:K
    seq_p <- 1:p
    seq_N <- 1:N

    # the total number of parameters given by the following:
    # - K x p mean parameters m_k
    # - K x p (diagonal) variance matrix parameters S_k
    # - N x K mixture probabilities per observation of data
    mean_params <- paste0("m_", cartesian_index(seq_K, seq_p))
    var_params  <- paste0("ssq_", cartesian_index(seq_K, seq_p))
    prob_params <- paste0("phi_", cartesian_index(seq_N, seq_K))

    # initialize these variables
    # - setting m to a random location on the dataset
    # - setting ssq to 1 for all
    # - setting 1/K for all rows of data
    set.seed(923)
    init_locs <- as.matrix(data[sample(1:N, size = K), ])

    params <- list(max_iter)
    elbo <- list(max_iter - 1)
    steps <- list(max_iter - 1)

    params[[1]] <- data.frame(t(c(
        as.vector(init_locs),
        rep(1, times = length(var_params)),
        rep(1/K, times = length(prob_params))
    )))

    param_names <- c(mean_params, var_params, prob_params)
    colnames(params[[1]]) <- param_names
    G <- 0
    converged <- FALSE
    start.time <- Sys.time()

    set.seed(seed)
    for(t in 2:max_iter){
        old.params <- params[[t - 1]]

        samps <- generate_qmc_samples(mc_size, data, K, old.params, priors, method)
        updates <- generate_bbvi(samps)
        
        updater <- learn_rate(t, old = old.params, step = updates, G)
        new.params <- updater$new
        G <- updater$G

        # hard thresholding on ssq
        for(j in which(grepl("ssq_", param_names, fixed = TRUE))){
            if(new.params[1,j] < var_threshold[1]) new.params[1,j] <- var_threshold[1]
            if(new.params[1,j] > var_threshold[2]) new.params[1,j] <- var_threshold[2]
        }
        

        # normalizing phi
        mu.new <- as.numeric(new.params[1, grepl("m_", param_names, fixed = TRUE)])
        mu.new <- matrix(mu.new, nrow = K, ncol = p, byrow = FALSE)

        phi.new <- as.numeric(new.params[1,grepl("phi_", param_names, fixed = TRUE)])
        phi.new <- matrix(phi.new, nrow = N, ncol = K, byrow = FALSE)
        
        for(i in 1:N){
            probs <- sapply(as.numeric(phi.new[i,]), function(x){
                if(x < 0.01) return(0.01)
                if(is.na(x) || is.infinite(x) || is.nan(x)) return(0.01)
                return(x)
            })
            total_prob <- sum(probs, na.rm = TRUE)
            phi.new[i,] <- probs/total_prob
            
        }

        new.params[1, grepl("phi_", param_names, fixed = TRUE)] <- as.vector(phi.new)

        bayesLik <- bayes.logLik(data, mu.new, phi.new, priors)
        pDIC <- 2 * (bayesLik - samps$simLikelihood)
        elpd <- bayesLik - pDIC
        DIC <- -2 * bayesLik + pDIC

        params[[t]] <- new.params
        elbo[[t-1]] <- data.frame(iter = t, elapsed = difftime(Sys.time(), start.time, units = "secs"), elbo = samps$ELBO,
                                    bayesLikelihood = bayesLik, simLikelihood = samps$simLikelihood, elpd = elpd, DIC = DIC)
        steps[[t-1]] <- updates

        # check for convergence
        if (t > 2){
            elbo.chg <- log(abs(samps$ELBO/elbo[[t-1]]$elbo[1]))
        }else{
            elbo.chg <- 1.0000
        }
        
        norm_old <- sum(old.params[1,]^2)
        norm_chg <- sum((new.params[1,] - old.params[1,])^2)
        lambda <- sqrt(norm_chg)/sqrt(norm_old)

        if((t > min_iter & crit == "param" & lambda < converge) | (t > min_iter & crit == "elbo" & elbo.chg < converge)){
            message(paste0("Algorithm converged after ", t, " iterations."))
            converged <- TRUE
            params[(t+1):max_iter] <- NULL
            elbo[t:(max_iter - 1)] <- NULL
            steps[t:(max_iter - 1)] <- NULL
            break
        }

        if(t%%100 == 0 & verbose) message(paste0("QMCVI-",method,": Iteration ", t,
                                                    " | lambda: ", round(lambda, 4),
                                                    " | ELBO: ", round(samps$ELBO,2),
                                                    " | ELBO Change: ", round(elbo.chg, 4)))
    }

    if(! converged) message(paste0("Algorithm did not converge after ", t, " steps. Results may not be reliable."))
    
    paramlist <- do.call('rbind', params)
    elbolist <- do.call('rbind', elbo)
    steplist <- do.call('rbind', steps)

    return(list(
        trace = paramlist,
        elbo = elbolist,
        updates = steplist,
        iterations = t,
        converged = converged,
        elapsed.time = difftime(Sys.time(), start.time, units = "secs"),
        N = N,
        K = K,
        p = p,
        data = data,
        method = method
    ))
}

#' Fit a Multivariate Mixture of Gaussian Using YOASOVI
#' 
#' Fits a set of K (pre-specified) components each a p-dimensional Gaussian distribution to a dataset.
#' The method uses the Naive and Rao-Blackwellized methods only.
#'
#' @param data an N x p matrix of observations
#' @param clusters the number K of components
#' @param method the flavor of BBVI to use in model fitting. Can be one of method = c("Naive","RB","JS+") for
#' the Naive, Rao-Blackwellized, and James-Stein flavor, respectively. Defaults to method = "Naive"
#' @param learn_rate the learning rate to use for the stochastic optimization. Can be rate_constant(rho) for a
#' certain value of rho, rate_adagrad(eta) for a value of eta to use the AdaGrad method, and rate_rmsprop(eta, beta)
#' for the RMSProp modification of AdaGrad. Defaults to rate_constant(1e-3)
#' @param max_iter the total number of iterations to stop the method in case of non-convergence
#' @param min_iter the minimum number of iterations to perform before convergence is assessed
#' @param mc_size the number S of monte carlo samples to draw for generating the BBVI Update steps
#' @param priors a list of prior specifications, defaults to standard priors if priors = NULL.
#' @param converge the convergence threshold delta
#' @param criterion the criterion in which to test for convergence, can be one of criterion = c("param","elbo").
#' @param var_threshold the lower and upper threshold on the variance, defaults to var_threshold = c(0.01, 1e3).
#' @param verbose logical, a switch for whether a verbose execution should be done.
#' @param seed the seed to use for starting the random samples
#' 
#' @return An updated list of prior specifications, completed in case of missing priors
yoasovi_multimix <- function(data, clusters = 2, method = "Naive", learn_rate = rate_constant(), mc_size = 1, max_iter = 1000, min_iter = 1, priors = NULL, seed = 923,
                            patience = 100, crit = "param", var_threshold = c(0.01, 1e3), seed_init = 923, verbose = FALSE){
    p <- ncol(data) # number of components in the data
    N <- nrow(data) # number of rows in the data
    K <- clusters

    priors <- unpack_priors(priors, K, p, N)

    seq_K <- 1:K
    seq_p <- 1:p
    seq_N <- 1:N

    # the total number of parameters given by the following:
    # - K x p mean parameters m_k
    # - K x p (diagonal) variance matrix parameters S_k
    # - N x K mixture probabilities per observation of data
    mean_params <- paste0("m_", cartesian_index(seq_K, seq_p))
    var_params  <- paste0("ssq_", cartesian_index(seq_K, seq_p))
    prob_params <- paste0("phi_", cartesian_index(seq_N, seq_K))

    # initialize these variables
    # - setting m to a random location on the dataset
    # - setting ssq to 1 for all
    # - setting 1/K for all rows of data
    set.seed(seed_init)
    init_locs <- as.matrix(data[sample(1:N, size = K), ])

    params <- list(max_iter)
    elbo <- list(max_iter - 1)
    steps <- list(max_iter - 1)

    params[[1]] <- data.frame(t(c(
        as.vector(init_locs),
        rep(1, times = length(var_params)),
        rep(1/K, times = length(prob_params))
    )))

    param_names <- c(mean_params, var_params, prob_params)
    colnames(params[[1]]) <- param_names
    G <- 0
    converged <- FALSE
    start.time <- Sys.time()
    num_rejected <- 0

    set.seed(seed)
    for(t in 2:max_iter){
        old.params <- params[[t - 1]]
        ELBO.Prev <- ifelse(t > 3, elbo[[t - 2]]$elbo, -9e6)

        samps <- generate_reject_samples(mc_size, ELBO.Prev, tempered = t, data, K, old.params, priors, method)

        if(t <= 2 | samps$accept){
            updates <- generate_bbvi(samps)
            updater <- learn_rate(t, old = old.params, step = updates, G)
            new.params <- updater$new
            G <- updater$G

            # hard thresholding on ssq
            for(j in which(grepl("ssq_", param_names, fixed = TRUE))){
                if(new.params[1,j] < var_threshold[1]) new.params[1,j] <- var_threshold[1]
                if(new.params[1,j] > var_threshold[2]) new.params[1,j] <- var_threshold[2]
            }
            
            # normalizing phi
            mu.new <- as.numeric(new.params[1, grepl("m_", param_names, fixed = TRUE)])
            mu.new <- matrix(mu.new, nrow = K, ncol = p, byrow = FALSE)

            phi.new <- as.numeric(new.params[1,grepl("phi_", param_names, fixed = TRUE)])
            phi.new <- matrix(phi.new, nrow = N, ncol = K, byrow = FALSE)
            
            for(i in 1:N){
                probs <- sapply(as.numeric(phi.new[i,]), function(x){
                    if(x < 0.01) return(0.01)
                    if(is.na(x) || is.infinite(x) || is.nan(x)) return(0.01)
                    return(x)
                })
                total_prob <- sum(probs, na.rm = TRUE)
                phi.new[i,] <- probs/total_prob
                
            }

            new.params[1, grepl("phi_", param_names, fixed = TRUE)] <- as.vector(phi.new)

            bayesLik <- bayes.logLik(data, mu.new, phi.new, priors)
            pDIC <- 2 * (bayesLik - samps$simLikelihood)
            elpd <- bayesLik - pDIC
            DIC <- -2 * bayesLik + pDIC

            params[[t]] <- new.params
            elbo[[t-1]] <- data.frame(iter = t, elapsed = difftime(Sys.time(), start.time, units = "secs"), probaccept = samps$probaccept, elbo = samps$ELBO,
                                        bayesLikelihood = bayesLik, simLikelihood = samps$simLikelihood, elpd = elpd, DIC = DIC)
            steps[[t-1]] <- updates
            num_rejected <- 0

        }else{
            retain.params <- params[[t - 1]]
            retain.elbo <- elbo[[t - 2]]
            retain.steps <- steps[[t - 2]]

            retain.elbo$iter = t

            params[[t]] <- retain.params
            elbo[[t-1]] <- retain.elbo
            steps[[t-1]] <- retain.steps
            num_rejected <- num_rejected + 1

        }

        if(t > min_iter & num_rejected >= patience){
            message(paste0("Algorithm converged after ", t, " iterations."))
            converged <- TRUE
            params[(t+1):max_iter] <- NULL
            elbo[t:(max_iter - 1)] <- NULL
            steps[t:(max_iter - 1)] <- NULL
            break
        }
        

        if(t%%100 == 0 & verbose) message(paste0("YOASOVI-",method,": Iteration ", t, " | ELBO: ", round(elbo[[t-1]]$elbo,2)))
    }

    if(! converged) message(paste0("Algorithm did not converge after ", t, " steps. Results may not be reliable."))
    
    paramlist <- do.call('rbind', params)
    elbolist <- do.call('rbind', elbo)
    steplist <- do.call('rbind', steps)

    return(list(
        trace = paramlist,
        elbo = elbolist,
        updates = steplist,
        iterations = t,
        converged = converged,
        elapsed.time = difftime(Sys.time(), start.time, units = "secs"),
        N = N,
        K = K,
        p = p,
        data = data,
        method = method
    ))
}

#' Complete Joint Log-Likelihood For Naive and JS+ Estimator
#'
#' @param y An array of N x p observations containing the observed data
#' @param z An N x K matrix of one-hot vector encodings for the K components
#' @param mu An array of K x p cluster means
#' @param mu0 The prior mean for each mu \sim N(mu0, tausq), defaults to 0
#' @param tausq The prior variance for each mu \sim N(mu0, tausq), defaults to 1
#' @param sigsq The data variance for each y \sim N(zi^T mu, sigsq), defaults to 1
#' @param phi0 The prior probabilities for zi, zi \sim Categorical(phi0), defaults to 1/K
#' 
#' @return The complete joint likelihood p(y, zu, mu | mu0, tausq, sigsq, phi0)
logp.all <- function(data, z, mu, priors){
    p <- ncol(data) # number of components in the data
    N <- nrow(data) # number of rows in the data
    K <- ncol(z)

    m0 <- priors[["m0"]]
    S0 <- priors[["S0"]]
    Sig <- priors[["Sigma"]]
    Phi0 <- priors[["Phi0"]]

    term1 <- 0
    term2 <- 0
    for(i in 1:N){
        this.z <- as.numeric(z[i, ])
        k <- which(this.z == 1)
        x <- as.numeric(data[i,])
        mk <- as.numeric(mu[k, ])
        lik.norm <- -0.5 * p * log(2 * pi) - 0.5 * log(det(Sig)) - 0.5 * t(x - mk) %*% solve(Sig) %*% (x - mk)
        lik.catg <- sum(this.z * log(Phi0))

        term1 <- term1 + lik.norm
        term2 <- term2 + lik.catg
    }

    term3 <- 0
    for(k in 1:K){
        mk <- as.numeric(mu[k, ])
        lik.norm <- -0.5 * p * log(2 * pi) - 0.5 * log(det(S0)) - 0.5 * t(mk - m0) %*% solve(S0) %*% (mk - m0)
        term3 <- term3 + lik.norm
    }

    logp <- term1 + term2 + term3
    return(logp)
}

#' Rao-Blackwellized Joint Log-Likelihood For mu
#'
#' @param y An array of N observations containing the observed data
#' @param z An N x K matrix of one-hot vector encodings for the K components
#' @param mu An array of K cluster means
#' @param mu0 The prior mean for each mu \sim N(mu0, tausq), defaults to 0
#' @param tausq The prior variance for each mu \sim N(mu0, tausq), defaults to 1
#' @param sigsq The data variance for each y \sim N(zi^T mu, sigsq), defaults to 1
#' @param phi0 The prior probabilities for zi, zi \sim Categorical(phi0), defaults to 1/K
#' 
#' @return The joint likelihood p(y, zu, mu | mu0, tausq, sigsq, phi0) with terms only for mu
logp.mu <- function(data, z, mu, priors){
    p <- ncol(data) # number of components in the data
    N <- nrow(data) # number of rows in the data
    K <- ncol(z)

    m0 <- priors[["m0"]]
    S0 <- priors[["S0"]]
    Sig <- priors[["Sigma"]]
    Phi0 <- priors[["Phi0"]]

    logp <- numeric(K)
    for(k in 1:K){
        mk <- as.numeric(mu[k, ])

        term1 <- 0
        for(i in 1:N){
            this.z <- as.numeric(z[i, ])
            if(which(this.z == 1) != k) next

            x <- as.numeric(data[i,])
            lik.norm <- -0.5 * p * log(2 * pi) - 0.5 * log(det(Sig)) - 0.5 * t(x - mk) %*% solve(Sig) %*% (x - mk)

            term1 <- term1 + lik.norm
        }
        
        term2 <- -0.5 * p * log(2 * pi) - 0.5 * log(det(S0)) - 0.5 * t(mk - m0) %*% solve(S0) %*% (mk - m0)
        logp[k] <- term1 + term2
    }

    return(logp)
}

#' Rao-Blackwellized Joint Log-Likelihood For z
#'
#' @param y An array of N observations containing the observed data
#' @param z An N x K matrix of one-hot vector encodings for the K components
#' @param mu An array of K cluster means
#' @param mu0 The prior mean for each mu \sim N(mu0, tausq), defaults to 0
#' @param tausq The prior variance for each mu \sim N(mu0, tausq), defaults to 1
#' @param sigsq The data variance for each y \sim N(zi^T mu, sigsq), defaults to 1
#' @param phi0 The prior probabilities for zi, zi \sim Categorical(phi0), defaults to 1/K
#' 
#' @return The joint likelihood p(y, zu, mu | mu0, tausq, sigsq, phi0) with terms only for z
logp.z <- function(data, z, mu, priors){
    p <- ncol(data) # number of components in the data
    N <- nrow(data) # number of rows in the data
    K <- ncol(z)

    m0 <- priors[["m0"]]
    S0 <- priors[["S0"]]
    Sig <- priors[["Sigma"]]
    Phi0 <- priors[["Phi0"]]

    logp <- numeric(N)
    for(i in 1:N){
        this.z <- as.numeric(z[i, ])
        mk <- as.numeric(mu[which(this.z == 1), ])
        x <- as.numeric(data[i, ])
        logp[i] <- -0.5 * p * log(2 * pi) - 0.5 * log(det(Sig)) - 0.5 * t(x - mk) %*% solve(Sig) %*% (x - mk)
    }

    return(logp)
}

#' Complete Mean-Field Variational Approximation to the Posterior
#'
#' @param z An N x K matrix of one-hot vector encodings for the K components
#' @param mu An array of K cluster means
#' @param phi An N x K matrix of cluster probabilities such that z \sim q(phi)
#' @param m An array of K posterior means such that mu \sim q(m, ssq)
#' @param ssq An array of K posterior variances such that mu \sim q(m, ssq)
#' 
#' @return The complete mean-field approximation q(z, mu | phi, m, s)
logq.all <- function(z, mu, phi, m, ssq){
    p <- ncol(mu) # number of components in the data
    N <- nrow(z) # number of rows in the data
    K <- ncol(z)

    term1 <- 0
    for(i in 1:nrow(z)){
        term1 <- term1 + sum(z[i,] * log(phi[i, ]))
    }

    term2 <- sapply(1:K, function(k){
        Sk <- ssq[[k]]
        mk <- as.numeric(m[k,])
        muk <- as.numeric(mu[k,])
        lik.norm <- -0.5 * p * log(2 * pi) - 0.5 * log(det(Sk)) - 0.5 * t(muk - mk) %*% solve(Sk) %*% (muk - mk)
        return(lik.norm)
    })

    logq <- term1 + sum(term2)
    return(logq)
}

#' Rao-Blackwellized Mean-Field Variational Approximation to the Posterior for mu
#'
#' @param z An N x K matrix of one-hot vector encodings for the K components
#' @param mu An array of K cluster means
#' @param phi An N x K matrix of cluster probabilities such that z \sim q(phi)
#' @param m An array of K posterior means such that mu \sim q(m, ssq)
#' @param ssq An array of K posterior variances such that mu \sim q(m, ssq)
#' 
#' @return The mean-field approximation q(z, mu | phi, m, s) with relevant terms for mu
logq.mu <- function(z, mu, phi, m, ssq){
    p <- ncol(mu) # number of components in the data
    N <- nrow(z) # number of rows in the data
    K <- ncol(z)

    logq <- sapply(1:K, function(k){
        Sk <- ssq[[k]]
        mk <- as.numeric(m[k,])
        muk <- as.numeric(mu[k,])
        lik.norm <- -0.5 * p * log(2 * pi) - 0.5 * log(det(Sk)) - 0.5 * t(muk - mk) %*% solve(Sk) %*% (muk - mk)
        return(lik.norm)
    })

    return(logq)
}

#' Rao-Blackwellized Mean-Field Variational Approximation to the Posterior for z
#'
#' @param z An N x K matrix of one-hot vector encodings for the K components
#' @param mu An array of K cluster means
#' @param phi An N x K matrix of cluster probabilities such that z \sim q(phi)
#' @param m An array of K posterior means such that mu \sim q(m, ssq)
#' @param ssq An array of K posterior variances such that mu \sim q(m, ssq)
#' 
#' @return The mean-field approximation q(z, mu | phi, m, s) with relevant terms for z
logq.z <- function(z, mu, phi, m, ssq){

    logq <- sapply(1:nrow(z), function(i){
        this.phi <- phi[i, ] + 0.0000001 * (phi[i, ] == 0)
        sum(z[i,] * log(this.phi))
    })

    return(logq)
}

#' Gradient Function for Mu, given m
#'
#' @param mu An array of K cluster means
#' @param m An array of K posterior means such that mu \sim q(m, ssq)
#' @param ssq An array of K posterior variances such that mu \sim q(m, ssq)
#' 
#' @return The gradient function for mu, given m
grad.mu_m <- function(mu, m, ssq){
    K <- nrow(mu)
    p <- ncol(mu)

    grads <- t(sapply(1:K, function(k){
        Sk <- ssq[[k]]
        mk <- as.numeric(m[k,])
        muk <- as.numeric(mu[k,])

        return(solve(Sk) %*% (muk - mk))
    }))

    return(grads)
}

#' Gradient Function for Mu, given ssq
#'
#' @param mu An array of K cluster means
#' @param m An array of K posterior means such that mu \sim q(m, ssq)
#' @param ssq An array of K posterior variances such that mu \sim q(m, ssq)
#' 
#' @return The gradient function for mu, given ssq
grad.mu_ssq <- function(mu, m, ssq){
    K <- nrow(mu)
    p <- ncol(mu)

    grads <- t(sapply(1:K, function(k){
        Sk <- ssq[[k]]
        mk <- as.numeric(m[k,])
        muk <- as.numeric(mu[k,])

        grad <- -0.5 * (solve(Sk) - solve(Sk) %*% (muk - mk) %*% t(muk - mk) %*% solve(Sk))
        return(diag(grad))
    }))

    return(grads)
}

#' Gradient Function for z
#'
#' @param phi An N x K matrix of cluster probabilities such that z \sim q(phi)
#' @param phi An N x K matrix of cluster probabilities such that z \sim q(phi)
#' 
#' @return The gradient function for z
grad.z <- function(z, phi){
    N <- nrow(z)
    grad <- sapply(1:N, function(i){
        matrix(z[i,]/phi[i,], nrow=1, byrow=TRUE)
    })

    return(t(grad))
}

#' Positive Part Clipping Function for James-Stein Estimator
#'
#' @param arr A single-dimensional array
#' 
#' @return The same length array, with negative values coerced to zero
positive <- function(arr){
    sapply(arr, function(x){
        if(is.na(x) || is.infinite(x) || is.nan(x)) return(0)
        if(x > 0) return(x)
        return(0)
    })
}


#' Unpack All Parameters From Array Form
#' 
#' In general the algorithm represents all parameters in the algorithm in a single array.
#' The following function unpacks them into their correct form, i.e. a K x p matrix for
#' the cluster centroids, an N x K matrix for the component probabilities per observation,
#' and a list of p x p (diagonal) matrices for the variance-covariance matrix of the
#' multivariate normal distributions on mu.
#'
#' @param params A named horizontal array of parameter values
#' @param N the number of observations
#' @param K the number of components in the mixture distribution
#' @param p the number of dimensions of the multivariate normal distributions
#' 
#' @return The same length array, with negative values coerced to zero
unpack_parameters <- function(params, N, K, p){
    phi <- as.numeric(params[,grepl("phi_",colnames(params))])
    phi <- matrix(phi, nrow = N, ncol = K, byrow = FALSE)

    m <- as.numeric(params[,grepl("m_",colnames(params))])
    m <- matrix(m, nrow = K, ncol = p, byrow = FALSE)

    ssq <- list(K)
    for(k in 1:K){
        this.ssq <- as.numeric(params[,grepl(paste0("ssq_",k,"_"),colnames(params))])
        this.ssq <- diag(this.ssq, nrow = p, ncol = p)
        ssq[[k]] <- this.ssq
    }

    return(list(
        phi = phi,
        m = m,
        ssq = ssq
    ))
}

#' Draw Monte Carlo Samples of the ELBO Gradient
#'
#' @param S Number of Monte Carlo samples to draw, defaults to S = 1000
#' @param y The set of observations
#' @param phi An N x K matrix of cluster probabilities such that z \sim q(phi)
#' @param m An array of K posterior means such that mu \sim q(m, ssq)
#' @param ssq An array of K posterior variances such that mu \sim q(m, ssq)
#' @param mu0 The prior mean for each mu \sim N(mu0, tausq), defaults to 0
#' @param tausq The prior variance for each mu \sim N(mu0, tausq), defaults to 1
#' @param sigsq The data variance for each y \sim N(zi^T mu, sigsq), defaults to 1
#' @param phi0 The prior probabilities for zi, zi \sim Categorical(phi0), defaults to 1/K
#' @param method Which type of ELBO estimator to use, method = c("Naive","JS+","RB","RB+")
#' 
#' @return S samples from the corresponding mean-field variational approximations
generate_samples <- function(S = 1000, data, K, params, priors, method="Naive"){
    N <- nrow(data)
    p <- ncol(data)

    this.params <- unpack_parameters(params, N, K, p)
    phi <- this.params$phi
    m <- this.params$m
    ssq <- this.params$ssq

    elbograd <- list(S)
    ELBO <- list(S)
    likelihoods <- list(S)

    for(s in 1:S){
        # draw samples
        z <- t(sapply(1:nrow(phi), function(i){
            rmultinom(1,1,prob = phi[i,])
        }))

        mu <- t(sapply(1:K, function(k){
            rmvnorm(1, mean = m[k,], sigma = ssq[[k]])
        }))

        # compute current ELBO
        ELBO[[s]] <- logp.all(data, z, mu, priors) - logq.all(z, mu, phi, m, ssq)

        # compute the simulated likelihood
        likelihoods[[s]] <- sim.logLik(data, mu, z, priors)

        # compute ELBO gradient
        if(method %in% c("Naive","JS","JS+")){
            grad <- c(
                as.vector(grad.mu_m(mu, m, ssq)),
                as.vector(grad.mu_ssq(mu, m, ssq)),
                as.vector(grad.z(z, phi))
            )

            diff <- rep(logp.all(data, z, mu, priors) - logq.all(z, mu, phi, m, ssq),
                times = length(grad))

            elbograd[[s]] <- grad * diff

        }else{
            elbo.mu_m <- as.vector(grad.mu_m(mu, m, ssq)) * (as.vector(logp.mu(data, z, mu, priors)) - as.vector(logq.mu(z, mu, phi, m, ssq)))
            elbo.mu_ssq <- as.vector(grad.mu_ssq(mu, m, ssq)) *
                (rep(as.vector(logp.mu(data, z, mu, priors)), times = p) - rep(as.vector(logq.mu(z, mu, phi, m, ssq)), times = p))
            elbo.z <- as.vector(grad.z(z, phi)) *
                (rep(logp.z(data, z, mu, priors), times = K) - rep(logq.z(z, mu, phi, m, ssq), times = K))
            elbograd[[s]] <- c(elbo.mu_m, elbo.mu_ssq, elbo.z)

        }

    }

    ELBO <- do.call('rbind', ELBO)
    elbograd <- do.call('rbind', elbograd)
    likelihoods <- do.call('rbind', likelihoods)
    
    list_data <- list(Z = elbograd, ELBO = mean(ELBO), simLikelihood = mean(likelihoods), S = S, method=method)
    return(list_data)
}

owen_scramble <- function(seq) {
  n <- length(seq)
  for (i in 1:n) {
    seq[i] <- seq[i] ^ runif(1)
  }
  return(seq)
}

#' Draw Quasi-Monte Carlo Samples of the ELBO Gradient
#'
#' @param S Number of Monte Carlo samples to draw, defaults to S = 1000
#' @param y The set of observations
#' @param phi An N x K matrix of cluster probabilities such that z \sim q(phi)
#' @param m An array of K posterior means such that mu \sim q(m, ssq)
#' @param ssq An array of K posterior variances such that mu \sim q(m, ssq)
#' @param mu0 The prior mean for each mu \sim N(mu0, tausq), defaults to 0
#' @param tausq The prior variance for each mu \sim N(mu0, tausq), defaults to 1
#' @param sigsq The data variance for each y \sim N(zi^T mu, sigsq), defaults to 1
#' @param phi0 The prior probabilities for zi, zi \sim Categorical(phi0), defaults to 1/K
#' @param method Which type of ELBO estimator to use, method = c("Naive","JS+","RB","RB+")
#' 
#' @return S samples from the corresponding mean-field variational approximations
generate_qmc_samples <- function(S = 1000, data, K, params, priors, method="Naive"){
    N <- nrow(data)
    p <- ncol(data)

    this.params <- unpack_parameters(params, N, K, p)
    phi <- this.params$phi
    m <- this.params$m
    ssq <- this.params$ssq

    elbograd <- list(S)
    ELBO <- list(S)
    likelihoods <- list(S)

    # the maximum dimension needed for a Sobol sequence is length(mean)
    sobol_norms <- sobol(S, dim = ncol(m))
    sobol_norms <- apply(sobol_norms, 2, owen_scramble)

    sobol_multinoms <- sobol(S, dim = ncol(phi))
    sobol_multinoms <- apply(sobol_multinoms, 2, owen_scramble)
   

    for(s in 1:S){
        # draw samples
        z <- t(sapply(1:nrow(phi), function(i){
            qmc_rmultinom(1,1,prob = phi[i,], sobol_samples = sobol_multinoms[s,])
        }))

        mu <- t(sapply(1:K, function(k){
            qmc_rmvnorm(1, mean = m[k,], sigma = ssq[[k]], sobol_samples = matrix(sobol_norms[s,], nrow = 1))
        }))

        # compute current ELBO
        ELBO[[s]] <- logp.all(data, z, mu, priors) - logq.all(z, mu, phi, m, ssq)

        # compute the simulated likelihood
        likelihoods[[s]] <- sim.logLik(data, mu, z, priors)

        # compute ELBO gradient
        if(method %in% c("Naive","JS","JS+")){
            grad <- c(
                as.vector(grad.mu_m(mu, m, ssq)),
                as.vector(grad.mu_ssq(mu, m, ssq)),
                as.vector(grad.z(z, phi))
            )

            diff <- rep(logp.all(data, z, mu, priors) - logq.all(z, mu, phi, m, ssq),
                times = length(grad))

            elbograd[[s]] <- grad * diff

        }else{
            elbo.mu_m <- as.vector(grad.mu_m(mu, m, ssq)) * (as.vector(logp.mu(data, z, mu, priors)) - as.vector(logq.mu(z, mu, phi, m, ssq)))
            elbo.mu_ssq <- as.vector(grad.mu_ssq(mu, m, ssq)) *
                (rep(as.vector(logp.mu(data, z, mu, priors)), times = p) - rep(as.vector(logq.mu(z, mu, phi, m, ssq)), times = p))
            elbo.z <- as.vector(grad.z(z, phi)) *
                (rep(logp.z(data, z, mu, priors), times = K) - rep(logq.z(z, mu, phi, m, ssq), times = K))
            elbograd[[s]] <- c(elbo.mu_m, elbo.mu_ssq, elbo.z)

        }

    }

    ELBO <- do.call('rbind', ELBO)
    elbograd <- do.call('rbind', elbograd)
    likelihoods <- do.call('rbind', likelihoods)
    
    list_data <- list(Z = elbograd, ELBO = mean(ELBO), simLikelihood = mean(likelihoods), S = S, method=method)
    return(list_data)
}

#' Perform Rejection Sampling With A Single Sample
#'
#' @param ELBO.Prev Previous ELBO Value
#' @param data The set of observations
#' @param K The number of clusters in the multivariate normal mixture
#' @param params The list of parameters drawn
#' @param priors The list of priors set for the mixture
#' @param method Which type of ELBO estimator to use, method = c("Naive","JS+","RB","RB+")
#' 
#' @return S samples from the corresponding mean-field variational approximations
generate_reject_samples <- function(S = 10, ELBO.Prev, tempered = 10, data, K, params, priors, method="Naive"){
    N <- nrow(data)
    p <- ncol(data)

    this.params <- unpack_parameters(params, N, K, p)
    phi <- this.params$phi
    m <- this.params$m
    ssq <- this.params$ssq

    elbograd <- list()
    ELBO <- list()
    likelihoods <- list()
    num_accept <- 0

    for(s in 1:S){
        # draw samples
        z <- t(sapply(1:nrow(phi), function(i){
            rmultinom(1,1,prob = phi[i,])
        }))

        mu <- t(sapply(1:K, function(k){
            rmvnorm(1, mean = m[k,], sigma = ssq[[k]])
        }))

        # compute current ELBO
        ELBO.New <- logp.all(data, z, mu, priors) - logq.all(z, mu, phi, m, ssq)

        # compute the simulated likelihood
        Lik.New <- sim.logLik(data, mu, z, priors)

        # compute ELBO gradient
        if(method %in% c("Naive","JS","JS+")){
            grad <- c(
                as.vector(grad.mu_m(mu, m, ssq)),
                as.vector(grad.mu_ssq(mu, m, ssq)),
                as.vector(grad.z(z, phi))
            )

            diff <- rep(logp.all(data, z, mu, priors) - logq.all(z, mu, phi, m, ssq),
                times = length(grad))

            elbograd.new <- grad * diff

        }else{
            elbo.mu_m <- as.vector(grad.mu_m(mu, m, ssq)) * (as.vector(logp.mu(data, z, mu, priors)) - as.vector(logq.mu(z, mu, phi, m, ssq)))
            elbo.mu_ssq <- as.vector(grad.mu_ssq(mu, m, ssq)) *
                (rep(as.vector(logp.mu(data, z, mu, priors)), times = p) - rep(as.vector(logq.mu(z, mu, phi, m, ssq)), times = p))
            elbo.z <- as.vector(grad.z(z, phi)) *
                (rep(logp.z(data, z, mu, priors), times = K) - rep(logq.z(z, mu, phi, m, ssq), times = K))
            elbograd.new <- c(elbo.mu_m, elbo.mu_ssq, elbo.z)

        }

        elbograd.new <- matrix(elbograd.new, nrow = 1, ncol = length(elbograd.new))

        probaccept <- min(1, 1 + (abs(ELBO.Prev) - abs(ELBO.New)) * 3 * log(tempered) /(abs(ELBO.Prev)))
        if(runif(1) < probaccept){
            num_accept <- num_accept + 1
            ELBO[[length(ELBO) + 1]] <- ELBO.New
            likelihoods[[length(likelihoods) + 1]] <- Lik.New
            elbograd[[length(elbograd) + 1]] <- elbograd.new
        }else if(s == 1){
            # guarantees there is at least one element even if not accepted
            # for error catching
            ELBO[[length(ELBO) + 1]] <- ELBO.New
            likelihoods[[length(likelihoods) + 1]] <- Lik.New
            elbograd[[length(elbograd) + 1]] <- elbograd.new
        }
    }

    accept = ifelse(num_accept > 0, TRUE, FALSE)
    ELBO <- do.call('rbind', ELBO)
    elbograd <- do.call('rbind', elbograd)
    likelihoods <- do.call('rbind', likelihoods)
    
    list_data <- list(accept = accept, probaccept = probaccept, Z = elbograd, ELBO = mean(ELBO), simLikelihood = mean(likelihoods), S = S, method=method)
    return(list_data)
}

#' Generate the BBVI Update Steps
#'
#' @param dat The generated samples drawn from `generate_samples`
#' 
#' @return A (N + 2) x (K) matrix of update steps, plus a column of parameter labels
generate_bbvi <- function(dat, method = NULL){
    if(is.null(method)) method = dat$method
    Z = dat$Z
    S = dat$S

    # Now estimate the mean
    means <- apply(Z, 2, mean)
    norm_est = sum(means^2)
    m = length(means) - 1

    if(method %in% c("JS+","RB+")){
        # Estimate the variance
        var_index <- sample(1:S, size=round(S/3), replace=TRUE)
        var_samples <- Z[var_index, ]
        variances <- apply(var_samples, 2, var)
        
        means <- positive(1 - ((m - 3) * variances)/norm_est) * means
    }

    return(means)
}

#' Constant Rate Function Constructor
#' 
#' Returns a function to use as a learning rate for the BBVI Update Steps
#'
#' @param rho The constant learning rate, at each step rho(t) = rho
#' 
#' @return Returns a function that the BBVI Update step uses for updating parameters
rate_constant <- function(rho = 1e-6, ...){
    rate_function <- function(t, old, step, ...){
        new <- old + rep(rho, times = length(step)) * step
        return(list(new = new, G = 0))
    }

    return(rate_function)
}

#' AdaGrad Rate Function Constructor
#' 
#' Returns a function to use as a learning rate for the BBVI Update Steps
#'
#' @param eta The tuning parameter eta for performing AdaGrad
#' 
#' @return Returns a function that the BBVI Update step uses for updating parameters
rate_adagrad <- function(eta = 1, ...){
    rate_function <- function(t, old, step, G, ...){
        G <- G + step %*% t(step)
        perturb <- runif(nrow(G), 0.001, 0.01)
        grads <- eta * 1/sqrt(diag(G) + perturb)

        new <- old + grads * step
        return(list(new = new, G = G))
    }

    return(rate_function)
}

#' RMSProp Rate Function Constructor
#' 
#' Returns a function to use as a learning rate for the BBVI Update Steps
#'
#' @param eta The tuning parameter eta for performing RMSProp
#' @param beta The tuning parameter for decay in performing RMSProp
#' 
#' @return Returns a function that the BBVI Update step uses for updating parameters
rate_rmsprop <- function(eta = 1, beta = 0.9, ...){
    rate_function <- function(t, old, step, G, ...){
        G <- beta * G + (1 - beta) * (step %*% t(step))
        perturb <- runif(nrow(G), 0.001, 0.01)
        grads <- eta * 1/sqrt(diag(G) + perturb)

        new <- old + grads * step
        return(list(new = new, G = G))
    }

    return(rate_function)
}

#' Plot A Bivariate Mixture of Gaussian
#' 
#' Plots the data used to train the multivariate gaussian mixture using BBVI.
#' Note that this only works for bivariate normal distributions.
#'
#' @param bbvi The object outputted by multimix_bbvi
#' 
#' @return Plots the scatter of observations, with a super-imposed contour plot of the
#' resulting densities of the mixing distributions.
plot_multimix <- function(bbvi, contours = TRUE, colors=c("black","red"), add = FALSE){
    N <- bbvi$N 
    K <- bbvi$K 
    p <- bbvi$p

    last_iter <- bbvi$iterations
    params <- bbvi$trace[last_iter, ]
    rownames(params) <- NULL
    this.params <- unpack_parameters(as.matrix(params), N, K, p)
    
    x <- seq(from = min(bbvi$data[,1]), to = max(bbvi$data[,1]), length = 100)
    y <- seq(from = min(bbvi$data[,2]), to = max(bbvi$data[,2]), length = 100)

    if(! add){
        par(bg = "white")
        plot(bbvi$data, col = colors[1])
    }

    if(contours){
        for(k in 1:K){
            mu <- as.numeric(this.params$m[k,])
            sigma <- this.params$ssq[[k]]
            f <- function(x, y) dmvnorm(cbind(x, y), mu, sigma)
            z <- outer(x, y, f)
            contour(x, y, z, add = TRUE, col = colors[2])
        }
    }

    points(this.params$m, pch="+", cex=2, col = colors[2])
    return(this.params$m)
}

#' Summary Mixtures of the BBVI Results
#' 
#' Returns a table of performance measures on the BBVI runs used for the update step.
#'
#' @param bbvi The object outputted by multimix_bbvi
#' 
#' @return A dataframe containing the number of iterations, time (in minutes),
#' expected predictive log-likelihood (elpd), and DIC of the resulting multivariate
#' mixture model.
summary_multimix <- function(bbvi){
    N <- bbvi$N 
    K <- bbvi$K 
    p <- bbvi$p

    last_iter <- bbvi$iterations
    params <- bbvi$trace[last_iter, ]
    rownames(params) <- NULL
    this.params <- unpack_parameters(as.matrix(params), N, K, p)

    performance_data <- data.frame(
        method = bbvi$method,
        iterations = bbvi$iterations,
        time = abs(as.numeric(bbvi$elapsed.time/60)),
        ELBO = as.numeric(bbvi$elbo[bbvi$iterations - 1, "elbo"]),
        elpd = as.numeric(bbvi$elbo[bbvi$iterations - 1, "elpd"]),
        DIC = as.numeric(bbvi$elbo[bbvi$iterations - 1, "DIC"])
    )

    return(performance_data)
}

#' Log-Likelihood of Simulated Parameters
#' 
#' For use during sampling. Takes in a random sample of parameters mu and z as computes for
#' the log-likelihood of the data.
#'
#' @param data An N x p matrix of observations for fitting the multivariate mixture model
#' @param mu A K x p matrix of means for the K components and their locations in p dimensions
#' @param z An N x K matrix of assignments for each observation to the K mixing components
#' @param priors A list of prior distributions (and assumptions)
#' 
#' @return The data log-likelihood for the given simulated parameters.
sim.logLik <- function(data, mu, z, priors){
    N = nrow(data)
    p = ncol(data)
    K = ncol(z)

    Sig <- priors[["Sigma"]]

    logp <- 0
    for(i in 1:N){
        this.z <- as.numeric(z[i, ])
        k <- which(this.z == 1)
        x <- as.numeric(data[i,])
        mk <- as.numeric(mu[k, ])
        lik.norm <- -0.5 * p * log(2 * pi) - 0.5 * log(det(Sig)) - 0.5 * t(x - mk) %*% solve(Sig) %*% (x - mk)
        logp <- logp + lik.norm
    }

    return(logp)
}

#' Log-Likelihood of Average Parameters
#' 
#' For use during sampling. Takes in the posterior mean values of mu and phi and computes
#' for the data log-likelihood.
#'
#' @param data An N x p matrix of observations for fitting the multivariate mixture model
#' @param mu A K x p matrix of means for the K components and their locations in p dimensions
#' @param phi An N x K matrix of probabilities for each observation to the K mixing components
#' @param priors A list of prior distributions (and assumptions)
#' 
#' @return The data log-likelihood for the given posterior average parameters.
bayes.logLik <- function(data, mu, phi, priors){
    N = nrow(data)
    p = ncol(data)
    K = ncol(phi)

    Sig <- priors[["Sigma"]]

    logp <- 0
    for(i in 1:N){
        x <- as.numeric(data[i,])

        for(k in 1:K){
            mk <- as.numeric(mu[k, ])
            lik.norm <- -0.5 * p * log(2 * pi) - 0.5 * log(det(Sig)) - 0.5 * t(x - mk) %*% solve(Sig) %*% (x - mk)
            logp <- logp + as.numeric(phi[i,k]) * lik.norm
        }
    }

    return(logp)
}

#' Quasi-Monte Carlo Random Samples From the (Univariate) Normal Distribution
#'
#' @param n The number of QMC samples to generate
#' @param mean The mean of the normal distribution
#' @param sd The standard deviation of the normal distribution
#' 
#' @return A vector of length n samples from the normal distribution
qmc_rnorm <- function(n, mean, sd){
 sobol_samples <- sobol(n, dim = 1)
 normal_samples <- qnorm(sobol_samples, mean = mean, sd = sd)
 return(normal_samples)
}

#' Quasi-Monte Carlo Random Samples From the (Multivariate) Normal Distribution
#'
#' @param n The number of QMC samples to generate
#' @param mean A p-size vector of means of the normal distribution
#' @param sigma A pxp covariance matrix for the normal distribution
#' 
#' @return A vector of length n samples from the normal distribution
qmc_rmvnorm <- function(n, mean, sigma, sobol_samples = NULL){
  Sig <- chol(sigma)

  if(is.null(sobol_samples)){
    sobol_samples <- sobol(n, dim = length(mean))
  }
  
  mvtsamp <- t(sapply(1:n, function(i){
    as.numeric(qnorm(sobol_samples[i,]) %*% Sig) + mean
  }))
  
  return(mvtsamp)
}

#' Quasi-Monte Carlo Random Samples From The Multinomial Distribution
#'
#' @param n The number of QMC samples to generate
#' @param size The number of N objects to place into K boxes
#' @param prob The probability for each of the K boxes
#' 
#' @return A vector of length n samples from the multinomial distribution
qmc_rmultinom <- function(n, size, prob, sobol_samples = NULL) {
  k <- length(prob)
  prob <- prob/sum(prob)
  multinom_samples <- matrix(0, nrow = n, ncol = k)
  
  if(is.null(sobol_samples)){
    sobol_samples <- sobol(n * length(prob), dim = 1)
  }
  
  for (i in 1:n) {
    cum_prob <- cumsum(prob)
    u <- sobol_samples[(i - 1) * k + 1:i * k]
    counts <- numeric(k)
    
    for (j in 1:size) {
      sample_index <- min(which(u[j] <= cum_prob))
      counts[sample_index] <- counts[sample_index] + 1
    }
    
    multinom_samples[i, ] <- counts
  }
  
  return(multinom_samples)
}