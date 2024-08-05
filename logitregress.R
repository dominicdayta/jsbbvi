require(mvtnorm)
require(randtoolbox)
require(distr)
require(pROC)

#' Random Sample Generation
#'
#' @param n Size of mixture sample to be generated
#' @param beta_coeff An array of p+1 beta coefficients, including the intercept and slopes for p variables
#' @param means An array of p means for the X variables
#' @param seed Random seed (defaults to 923)
#' 
#' @return An array of n observations corresponding to the Gaussian mixture distribution
generate_mixture <- function(n, beta_coeff, means = NULL, seed=923){
    set.seed(seed)
    p <- length(beta_coeff) - 1

    if(is.null(means)) means <- rnorm(p, 0, sd = 1)
    vars <- diag(1/rgamma(p, 3, 4), nrow = p, ncol = p)
    X <- rmvnorm(n, mean = means, sigma = vars)
    X <- cbind(X0 = rep(1, times = n), X)
    
    y <- 1/(1 + exp(- X %*% beta_coeff))
    t <- sapply(y, function(x){
        rbinom(1, size = 1, prob = x)
    })

    
    samp <- matrix(as.numeric(X[,-1]), nrow = n, ncol = p)
    colnames(samp) <- paste0("X", 1:p)

    return(list(samp=samp, categories=t, coeffs = beta_coeff))
}

#' Fit a Bayesian Logistic Regression Model Using BBVI
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
bbvi_logitreg <- function(categories, X, method = "Naive", learn_rate = rate_constant(), mc_size = 1000, max_iter = 1000, min_iter = 1, priors = NULL, seed = 923,
                            converge = 1e-4, crit = "param", gamma_threshold = c(1e-1, 50), verbose = FALSE){
    p <- ncol(X)
    n <- nrow(X)

    priors <- unpack_priors(priors)

    seq_p <- 1:(p+1)

    # the total number of parameters given by the following:
    # - (p + 1) x 1 beta_1, alpha_star, and beta_star parameters
    coef_params <- paste0("beta_1_", seq_p)
    shape_params  <- paste0("alpha_star_", seq_p)
    rate_params <- paste0("beta_star_", seq_p)

    # initialize these variables by a random draw in the prior distributions4
    set.seed(923)

    params <- list(max_iter)
    elbo <- list(max_iter - 1)
    steps <- list(max_iter - 1)

    params[[1]] <- data.frame(t(c(
        rep(0, times = p + 1), # beta_1
        rep(2, times = p + 1), # alpha_star
        rep(3, times = p + 1)  # beta_star
    )))

    param_names <- c(coef_params, shape_params, rate_params)
    colnames(params[[1]]) <- param_names
    G <- 0
    converged <- FALSE
    start.time <- Sys.time()

    set.seed(seed)
    for(t in 2:max_iter){
        old.params <- params[[t - 1]]

        samps <- generate_samples(S = mc_size, categories, X, old.params, priors, method=method)
        updates <- generate_bbvi(samps)
        
        updater <- learn_rate(t, old = old.params, step = updates, G)
        new.params <- updater$new
        G <- updater$G

        # hard thresholding on ssq
        for(j in which(grepl("alpha_star_", param_names, fixed = TRUE))){
            if(new.params[1,j] < gamma_threshold[1]) new.params[1,j] <- gamma_threshold[1]
            if(new.params[1,j] > gamma_threshold[2]) new.params[1,j] <- gamma_threshold[2]
        }

        for(j in which(grepl("beta_star_", param_names, fixed = TRUE))){
            if(new.params[1,j] < gamma_threshold[1]) new.params[1,j] <- gamma_threshold[1]
            if(new.params[1,j] > gamma_threshold[2]) new.params[1,j] <- gamma_threshold[2]
        }

        alpha_star <- as.numeric(new.params[1, which(grepl("alpha_star_", param_names, fixed = TRUE))])
        beta_star <- as.numeric(new.params[1, which(grepl("alpha_star_", param_names, fixed = TRUE))])
        beta_1 <- as.numeric(new.params[1, which(grepl("beta_1", param_names, fixed = TRUE))])

        bayesLik <- bayes.logLik(categories, X, beta_1, alpha_star, beta_star, priors)
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
        n = n,
        p = p,
        data = list(categories = categories, X = X),
        method = method
    ))
}

#' Fit a Bayesian Logistic Regression Model Using YOASOVI
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
yoasovi_logitreg <- function(categories, X, method = "Naive", learn_rate = rate_constant(), mc_size = 1, max_iter = 1000, min_iter = 1, priors = NULL, seed = 923,
                            patience = 100, crit = "param", gamma_threshold = c(1e-1, 50), verbose = FALSE){
    p <- ncol(X)
    n <- nrow(X)

    priors <- unpack_priors(priors)

    seq_p <- 1:(p+1)

    # the total number of parameters given by the following:
    # - (p + 1) x 1 beta_1, alpha_star, and beta_star parameters
    coef_params <- paste0("beta_1_", seq_p)
    shape_params  <- paste0("alpha_star_", seq_p)
    rate_params <- paste0("beta_star_", seq_p)

    # initialize these variables by a random draw in the prior distributions4
    set.seed(923)

    params <- list(max_iter)
    elbo <- list(max_iter - 1)
    steps <- list(max_iter - 1)

    params[[1]] <- data.frame(t(c(
        rep(0, times = p + 1), # beta_1
        rep(2, times = p + 1), # alpha_star
        rep(3, times = p + 1)  # beta_star
    )))

    param_names <- c(coef_params, shape_params, rate_params)
    colnames(params[[1]]) <- param_names
    G <- 0
    converged <- FALSE
    start.time <- Sys.time()
    num_rejected <- 0

    set.seed(seed)
    for(t in 2:max_iter){
        old.params <- params[[t - 1]]
        ELBO.Prev <- ifelse(t > 3, elbo[[t - 2]]$elbo, -9e6)

        samps <- generate_reject_samples(mc_size, ELBO.Prev, t, categories, X, old.params, priors, method=method)
        updates <- generate_bbvi(samps)

        if(t <= 2 | samps$accept){
            updater <- learn_rate(t, old = old.params, step = updates, G)
            new.params <- updater$new
            G <- updater$G

            # hard thresholding on alpha and beta parameters
            for(j in which(grepl("alpha_star_", param_names, fixed = TRUE))){
                if(new.params[1,j] < gamma_threshold[1]) new.params[1,j] <- gamma_threshold[1]
                if(new.params[1,j] > gamma_threshold[2]) new.params[1,j] <- gamma_threshold[2]
            }

            for(j in which(grepl("beta_star_", param_names, fixed = TRUE))){
                if(new.params[1,j] < gamma_threshold[1]) new.params[1,j] <- gamma_threshold[1]
                if(new.params[1,j] > gamma_threshold[2]) new.params[1,j] <- gamma_threshold[2]
            }

            alpha_star <- as.numeric(new.params[1, which(grepl("alpha_star_", param_names, fixed = TRUE))])
            beta_star <- as.numeric(new.params[1, which(grepl("alpha_star_", param_names, fixed = TRUE))])
            beta_1 <- as.numeric(new.params[1, which(grepl("beta_1", param_names, fixed = TRUE))])

            bayesLik <- bayes.logLik(categories, X, beta_1, alpha_star, beta_star, priors)
            pDIC <- 2 * (bayesLik - samps$simLikelihood)
            elpd <- bayesLik - pDIC
            DIC <- -2 * bayesLik + pDIC

            params[[t]] <- new.params
            elbo[[t-1]] <- data.frame(iter = t, elapsed = difftime(Sys.time(), start.time, units = "secs"), elbo = samps$ELBO,
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

        # check for convergence
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
        n = n,
        p = p,
        data = list(categories = categories, X = X),
        method = method
    ))
}

#' Fit a Bayesian Logistic Regression Model Using QMCVI
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
qmcvi_logitreg <- function(categories, X, method = "Naive", learn_rate = rate_constant(), mc_size = 1000, max_iter = 1000, min_iter = 1, priors = NULL, seed = 923,
                            converge = 1e-4, crit = "param", gamma_threshold = c(1e-1, 50), verbose = FALSE){
    p <- ncol(X)
    n <- nrow(X)

    priors <- unpack_priors(priors)

    seq_p <- 1:(p+1)

    # the total number of parameters given by the following:
    # - (p + 1) x 1 beta_1, alpha_star, and beta_star parameters
    coef_params <- paste0("beta_1_", seq_p)
    shape_params  <- paste0("alpha_star_", seq_p)
    rate_params <- paste0("beta_star_", seq_p)

    # initialize these variables by a random draw in the prior distributions4
    set.seed(923)

    params <- list(max_iter)
    elbo <- list(max_iter - 1)
    steps <- list(max_iter - 1)

    params[[1]] <- data.frame(t(c(
        rep(0, times = p + 1), # beta_1
        rep(2, times = p + 1), # alpha_star
        rep(3, times = p + 1)  # beta_star
    )))

    param_names <- c(coef_params, shape_params, rate_params)
    colnames(params[[1]]) <- param_names
    G <- 0
    converged <- FALSE
    start.time <- Sys.time()

    set.seed(seed)
    for(t in 2:max_iter){
        old.params <- params[[t - 1]]

        samps <- generate_qmc_samples(S = mc_size, categories, X, old.params, priors, method=method)
        updates <- generate_bbvi(samps)
        
        updater <- learn_rate(t, old = old.params, step = updates, G)
        new.params <- updater$new
        G <- updater$G

        # hard thresholding on ssq
        for(j in which(grepl("alpha_star_", param_names, fixed = TRUE))){
            if(new.params[1,j] < gamma_threshold[1]) new.params[1,j] <- gamma_threshold[1]
            if(new.params[1,j] > gamma_threshold[2]) new.params[1,j] <- gamma_threshold[2]
        }

        for(j in which(grepl("beta_star_", param_names, fixed = TRUE))){
            if(new.params[1,j] < gamma_threshold[1]) new.params[1,j] <- gamma_threshold[1]
            if(new.params[1,j] > gamma_threshold[2]) new.params[1,j] <- gamma_threshold[2]
        }

        alpha_star <- as.numeric(new.params[1, which(grepl("alpha_star_", param_names, fixed = TRUE))])
        beta_star <- as.numeric(new.params[1, which(grepl("alpha_star_", param_names, fixed = TRUE))])
        beta_1 <- as.numeric(new.params[1, which(grepl("beta_1", param_names, fixed = TRUE))])

        bayesLik <- bayes.logLik(categories, X, beta_1, alpha_star, beta_star, priors)
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
        n = n,
        p = p,
        data = list(categories = categories, X = X),
        method = method
    ))
}

#' Complete Joint Log-Likelihood
#'
#' @param categories A n x 1 array of (0,1) categories per observation
#' @param X A n x p matrix of observations for the p predictor variables
#' @param coeffs A (p + 1) x 1 array of beta coefficinents, including the intercept and slopes for p variables
#' @param alpha A list of alpha variables for the variances
#' @param priors The priors in list form, containing alpha0 and beta0
#' 
#' @return The complete joint likelihood p(y, zu, mu | mu0, tausq, sigsq, phi0)
logp.all <- function(categories, X, coeffs, alpha, priors){
    n <- nrow(X) # number of rows in the data
    p <- length(coeffs) - 1

    X <- cbind(X0 = rep(0, times = n), X)
    y  = 1/(1 + exp(-X %*% coeffs))

    alpha0 <- priors[["alpha0"]]
    beta0 <- priors[["beta0"]]

    y_star <- ifelse(y > 0.99999, 0.99999, y)
    y_star <- ifelse(y_star < 0.00001, 0.00001, y_star)
    log.data <- sum(categories * log(y_star) + (1 - categories) * log(1 - y_star))
    log.coeff <- dmvnorm(coeffs, mean = rep(0, times = p + 1), sigma = diag(1/alpha), log = TRUE)
    log.alpha <- sapply(alpha, function(x){
        dgamma(x, shape = alpha0, rate = beta0, log = TRUE)
    })

    logp <- sum(log.data) + sum(log.coeff) + sum(log.alpha)
    return(logp)
}

logp.beta <- function(categories, X, coeffs, alpha, priors){
    n <- nrow(X) # number of rows in the data
    p <- length(coeffs) - 1

    X <- cbind(X0 = rep(0, times = n), X)
    y  = 1/(1 + exp(-X %*% coeffs))

    alpha0 <- priors[["alpha0"]]
    beta0 <- priors[["beta0"]]

    y_star <- ifelse(y > 0.99999, 0.99999, y)
    y_star <- ifelse(y_star < 0.00001, 0.00001, y_star)
    log.data <- sum(categories * log(y_star) + (1 - categories) * log(1 - y_star))
    log.coeff <- dmvnorm(coeffs, mean = rep(0, times = p + 1), sigma = diag(1/alpha), log = TRUE)

    logp <- sum(log.data) + log.coeff
    return(logp)
}

logp.alpha <- function(categories, X, coeffs, alpha, priors){
    n <- nrow(X) # number of rows in the data
    p <- length(coeffs) - 1

    X <- cbind(X0 = rep(0, times = n), X)
    y  = 1/(1 + exp(-X %*% coeffs))

    alpha0 <- priors[["alpha0"]]
    beta0 <- priors[["beta0"]]

    log.coeff <- dmvnorm(coeffs, mean = rep(0, times = p + 1), sigma = diag(1/alpha), log = TRUE)
    log.alpha <- sapply(alpha, function(x){
        dgamma(x, shape = alpha0, rate = beta0, log = TRUE)
    })

    logp <- log.coeff + log.alpha
    return(logp)
}

#' Complete Mean-Field Variational Approximation to the Posterior
#'
#' @param coeffs A (p + 1) x 1 array of beta coefficinents, including the intercept and slopes for p variables
#' @param alpha A list of alpha variables for the variances
#' @param beta_1 The (p + 1) variational means for each beta coefficient
#' @param alpha_star The (p + 1) variational shape parameters for each alpha
#' @param beta_star The (p + 1) variational rate parameters for each alpha
#' 
#' @return The complete mean-field approximation q(z, mu | phi, m, s)
logq.all <- function(coeffs, alpha, beta_1, alpha_star, beta_star){
    p <- length(coeffs) - 1

    log.beta <- dmvnorm(coeffs, mean = beta_1, sigma = diag(1/alpha), log = TRUE)
    log.alpha <- sapply(seq_along(alpha), function(x){
        dgamma(alpha[x], shape = alpha_star[x], rate = beta_star[x], log = TRUE)
    })

    logq <- sum(log.beta) + sum(log.alpha)
    return(logq)
}

logq.alpha <- function(coeffs, alpha, beta_1, alpha_star, beta_star){
    p <- length(coeffs) - 1

    log.beta <- dmvnorm(coeffs, mean = beta_1, sigma = diag(1/alpha), log = TRUE)
    log.alpha <- sapply(seq_along(alpha), function(x){
        dgamma(alpha[x], shape = alpha_star[x], rate = beta_star[x], log = TRUE)
    })

    logq <- log.beta + log.alpha
    return(logq)
}

logq.beta <- function(coeffs, alpha, beta_1, alpha_star, beta_star){
    p <- length(coeffs) - 1

    log.beta <- dmvnorm(coeffs, mean = beta_1, sigma = diag(1/alpha), log = TRUE)

    logq <- log.beta
    return(logq)
}

#' Gradient Function for Coefficients, given beta1
#'
#' @param coeffs A (p + 1) x 1 array of beta coefficinents, including the intercept and slopes for p variables
#' @param alpha A list of alpha variables for the variances
#' @param beta_1 The (p + 1) variational means for each beta coefficient
#' 
#' @return The gradient function for mu, given m
grad.beta <- function(coeffs, alpha, beta_1){
    Sigma <- diag(1/alpha)
    grad <- as.numeric(solve(Sigma) %*% (coeffs - beta_1))
    return(grad)
}

#' Gradient Function for Coefficients, given alpha_star
#'
#' @param alpha A list of alpha variables for the variances
#' @param alpha_star The (p + 1) variational shape parameters for each alpha
#' @param beta_star The (p + 1) variational rate parameters for each alpha
#' 
#' @return The gradient function for mu, given m
grad.alpha_star <- function(alpha, alpha_star, beta_star){
    grad <- log(beta_star) - digamma(alpha_star) + alpha
    return(grad)
}

#' Gradient Function for Coefficients, given beta_star
#'
#' @param alpha A list of alpha variables for the variances
#' @param alpha_star The (p + 1) variational shape parameters for each alpha
#' @param beta_star The (p + 1) variational rate parameters for each alpha
#' 
#' @return The gradient function for mu, given m
grad.beta_star <- function(alpha, alpha_star, beta_star){
    grad <- alpha_star/beta_star - alpha
    return(grad)
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
unpack_parameters <- function(params){
    alpha_star  <- as.numeric(params[,grepl("alpha_star_",colnames(params))])
    beta_star <- as.numeric(params[,grepl("beta_star_",colnames(params))])
    beta_1 <- as.numeric(params[,grepl("beta_1_",colnames(params))])

    return(list(
        alpha_star = alpha_star,
        beta_star = beta_star,
        beta_1 = beta_1
    ))
}

#' Unpack List of Prior Distributions
#' 
#' For use in initializing the method.
#'
#' @param priors a list of prior specifications
#' 
#' @return An updated list of prior specifications, completed in case of missing priors
unpack_priors <- function(priors){
    
    if(is.null(priors[["alpha0"]])){
        priors[["alpha0"]] <- 1
    }

    if(is.null(priors[["beta0"]])){
        priors[["beta0"]] <- 1
    }

    return(priors)
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
generate_samples <- function(S = 1000, categories, X, params, priors, method="Naive"){
    n <- nrow(X)
    p <- ncol(X)

    this.params <- unpack_parameters(params)
    alpha_star <- this.params$alpha_star
    beta_star <- this.params$beta_star
    beta_1 <- this.params$beta_1

    elbograd <- list(S)
    ELBO <- list(S)
    likelihoods <- list(S)

    for(s in 1:S){
        # draw samples
        alpha <- sapply(1:(p+1), function(i){
            rgamma(1, shape = alpha_star[i], rate = beta_star[i])
        })

        coeffs <- as.numeric(rmvnorm(n = 1, mean = beta_1, sigma = diag(1/alpha)))
        

        # compute current ELBO
        ELBO[[s]] <- logp.all(categories, X, coeffs, alpha, priors) - logq.all(coeffs, alpha, beta_1, alpha_star, beta_star)

        # compute the simulated likelihood
        likelihoods[[s]] <- logp.all(categories, X, coeffs, alpha, priors)

        # compute ELBO gradient
        if(method %in% c("Naive","JS","JS+")){
            grad <- c(
                as.vector(grad.beta(coeffs, alpha, beta_1)),
                as.vector(grad.alpha_star(alpha, alpha_star, beta_star)),
                as.vector(grad.beta_star(alpha, alpha_star, beta_star))
            )

            diff <- rep(ELBO[[s]], times = length(grad))

            elbograd[[s]] <- grad * diff

        }else{
            # TODO: Create subroutine for RB, RB+ variants in case needed
            elbo.alpha <- as.vector(grad.beta(coeffs, alpha, beta_1)) * as.vector(logp.beta(categories, X, coeffs, alpha, priors) - logq.beta(coeffs, alpha, beta_1, alpha_star, beta_star))
            elbo.alpha_star <- as.vector(grad.alpha_star(alpha, alpha_star, beta_star)) * as.vector(logp.alpha(categories, X, coeffs, alpha, priors) - logq.alpha(coeffs, alpha, beta_1, alpha_star, beta_star))
            elbo.beta_star <- as.vector(grad.beta_star(alpha, alpha_star, beta_star)) * as.vector(logp.alpha(categories, X, coeffs, alpha, priors) - logq.alpha(coeffs, alpha, beta_1, alpha_star, beta_star))

            elbograd[[s]] <- c(elbo.alpha, elbo.alpha_star, elbo.beta_star)
        }

    }

    ELBO <- do.call('rbind', ELBO)
    elbograd <- do.call('rbind', elbograd)
    likelihoods <- do.call('rbind', likelihoods)
    
    list_data <- list(Z = elbograd, ELBO = mean(ELBO), simLikelihood = mean(likelihoods), S = S, method=method)
    return(list_data)
}

#' Perform Rejection Sampling With S Samples
#'
#' @param ELBO.Prev Previous ELBO Value
#' @param data The set of observations
#' @param K The number of clusters in the multivariate normal mixture
#' @param params The list of parameters drawn
#' @param priors The list of priors set for the mixture
#' @param method Which type of ELBO estimator to use, method = c("Naive","JS+","RB","RB+")
#' 
#' @return S samples from the corresponding mean-field variational approximations
generate_reject_samples <- function(S = 10, ELBO.Prev, tempered = 10, categories, X, params, priors, method="Naive"){
    n <- nrow(X)
    p <- ncol(X)

    this.params <- unpack_parameters(params)
    alpha_star <- this.params$alpha_star
    beta_star <- this.params$beta_star
    beta_1 <- this.params$beta_1

    elbograd <- list()
    ELBO <- list()
    likelihoods <- list()
    num_accept <- 0

    for(s in 1:S){
        # draw samples
        alpha <- sapply(1:(p+1), function(i){
            rgamma(1, shape = alpha_star[i], rate = beta_star[i])
        })

        coeffs <- as.numeric(rmvnorm(n = 1, mean = beta_1, sigma = diag(1/alpha)))

        # compute current ELBO
        ELBO.New <- logp.all(categories, X, coeffs, alpha, priors) - logq.all(coeffs, alpha, beta_1, alpha_star, beta_star)

        # compute the simulated likelihood
        Lik.New <- logp.all(categories, X, coeffs, alpha, priors)

        # compute ELBO gradient
        if(method %in% c("Naive","JS","JS+")){
            grad <- c(
                as.vector(grad.beta(coeffs, alpha, beta_1)),
                as.vector(grad.alpha_star(alpha, alpha_star, beta_star)),
                as.vector(grad.beta_star(alpha, alpha_star, beta_star))
            )

            diff <- rep(ELBO.New, times = length(grad))

            elbograd.new <- grad * diff

        }else{
            # TODO: Create subroutine for RB, RB+ variants in case needed

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
generate_qmc_samples <- function(S = 1000, categories, X, params, priors, method="Naive"){
    n <- nrow(X)
    p <- ncol(X)

    this.params <- unpack_parameters(params)
    alpha_star <- this.params$alpha_star
    beta_star <- this.params$beta_star
    beta_1 <- this.params$beta_1

    elbograd <- list(S)
    ELBO <- list(S)
    likelihoods <- list(S)
    
    sobol_gamma <- sobol(S, dim = length(alpha_star))
    sobol_gamma <- apply(sobol_gamma, 2, owen_scramble)

    sobol_norms <- sobol(S, dim = length(beta_1))
    sobol_norms <- apply(sobol_norms, 2, owen_scramble)

    for(s in 1:S){
        # draw samples
        alpha <- sapply(1:(p+1), function(i){
            as.numeric(qmc_rgamma(1, shape = alpha_star[i], rate = beta_star[i], sobol_samples = sobol_gamma[s,i]))
        })

        coeffs <- as.numeric(qmc_rmvnorm(n = 1, mean = beta_1, sigma = diag(1/alpha), sobol_samples = matrix(sobol_norms[s,], nrow = 1)))
        
        # compute current ELBO
        ELBO[[s]] <- logp.all(categories, X, coeffs, alpha, priors) - logq.all(coeffs, alpha, beta_1, alpha_star, beta_star)

        # compute the simulated likelihood
        likelihoods[[s]] <- logp.all(categories, X, coeffs, alpha, priors)

        # compute ELBO gradient
        if(method %in% c("Naive","JS","JS+")){
            grad <- c(
                as.vector(grad.beta(coeffs, alpha, beta_1)),
                as.vector(grad.alpha_star(alpha, alpha_star, beta_star)),
                as.vector(grad.beta_star(alpha, alpha_star, beta_star))
            )

            diff <- rep(ELBO[[s]], times = length(grad))

            elbograd[[s]] <- grad * diff

        }else{
            # TODO: Create subroutine for RB, RB+ variants in case needed

        }

    }

    ELBO <- do.call('rbind', ELBO)
    elbograd <- do.call('rbind', elbograd)
    likelihoods <- do.call('rbind', likelihoods)
    
    list_data <- list(Z = elbograd, ELBO = mean(ELBO), simLikelihood = mean(likelihoods), S = S, method=method)
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
bayes.logLik <- function(categories, X, beta_1, alpha_star, beta_star, priors){
    alpha <- alpha_star/beta_star

    logp <- logp.all(categories, X, beta_1, alpha, priors)

    return(logp)
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

#' Quasi-Monte Carlo Random Samples From the Gamma Distribution
#'
#' @param n The number of QMC samples to generate
#' @param shape The shape parameter alpha
#' @param rate The rate parameter beta
#' 
#' @return A vector of length n samples from the Gamma Distribution
qmc_rgamma <- function(n, shape, rate, sobol_samples = NULL){

    if(is.null(sobol_samples)){
        sobol_samples <- sobol(n, dim = 1)
    }
  
    samps <- qgamma(sobol_samples, shape = shape, rate = rate)
  
    return(samps)
}

logitreg_accuracy <- function(model, nsamp = 500, seed = 923, newdata = NULL){
    last_iter <- model$iterations
    last_parameters <- unpack_parameters(model$trace[last_iter,])
    p <- ncol(model$data$X)

    if(is.null(newdata)){
        t <- model$data$categories
        X <- cbind(X0 = 1, model$data$X)
    }else{
        t <- newdata$categories
        X <- cbind(X0 = 1, newdata$X)
    }

    # ROC-AUC Using Simulations
    set.seed(seed)

    roc_simul <- replicate(nsamp, {
        alpha <- sapply(1:(p+1), function(i){
            rgamma(1, shape = last_parameters$alpha_star[i], rate = last_parameters$beta_star[i])
        })

        coeffs <- as.numeric(rmvnorm(1, mean = as.numeric(last_parameters$beta_1), sigma = diag(1/alpha)))

        y <- 1/(1 + exp(-X %*% coeffs))
        ystar <- ifelse(y > 0.99999, 0.99999, y)
        ystar <- ifelse(ystar < 0.00001, 0.00001, ystar)

        roc_object <- suppressWarnings(suppressMessages(roc(t,ystar, quet=TRUE)))
        return(roc_object$auc)
    })

    summary_simul <- cbind(
        mean = mean(roc_simul, na.rm = TRUE),
        sd = sd(roc_simul, na.rm = TRUE)
    )

    # ROC-AUC Using Posterior Means
    coeffs <- as.numeric(last_parameters$beta_1)
    y <- 1/(1 + exp(-X %*% coeffs))
    ystar <- ifelse(y > 0.99999, 0.99999, y)
    ystar <- ifelse(ystar < 0.00001, 0.00001, ystar)
    roc_object <- suppressWarnings(suppressMessages(roc(t,ystar, quet=TRUE)))
    

    return(list(
        auc_simulations = summary_simul,
        auc_postmean = roc_object$auc,
        roc = roc_object
    ))
}

logitreg_predict <- function(model, nsamp = 500, seed = 923, newdata = NULL){
    last_iter <- model$iterations
    last_parameters <- unpack_parameters(model$trace[last_iter,])
    p <- ncol(model$data$X)

    if(is.null(newdata)){
        t <- model$data$categories
        X <- cbind(X0 = 1, model$data$X)
    }else{
        t <- newdata$categories
        X <- cbind(X0 = 1, newdata$X)
    }

    # ROC-AUC Using Posterior Means
    coeffs <- as.numeric(last_parameters$beta_1)
    y <- 1/(1 + exp(-X %*% coeffs))


    return(y)
}

logitreg_summary <- function(model){
    last_iter <- model$iterations
    last_elbo <- model$elbo[last_iter - 1,]
    last_parameters <- unpack_parameters(model$trace[last_iter,])

    coeftable <- do.call('cbind', last_parameters)

    return(last_elbo)
}