require(mvtnorm)

cartesian_index <- function(seq1, seq2){
    combs <- expand.grid(seq1,seq2)
    sapply(1:nrow(combs), function(x){
        return(paste0(combs[x,1], "_" ,combs[x,2]))
    })
}

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

bbvi_multimix <- function(data, clusters = 2, method = "Naive", learn_rate = rate_constant(), max_iter = 1000, mc_size = 1000, priors = NULL, seed = 923,
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
    set.seed(seed)
    init_locs <- as.matrix(data[sample(1:N, size = K), ])

    params <- list(max_iter)
    elbo <- list(max_iter - 1)
    steps <- list(max_iter - 1)

    params[[1]] <- data.frame(t(c(
        as.vector(init_locs),
        rep(10, times = length(var_params)),
        rep(1/K, times = length(prob_params))
    )))

    param_names <- c(mean_params, var_params, prob_params)
    colnames(params[[1]]) <- param_names
    G <- 0
    converged <- FALSE
    start.time <- Sys.time()

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
        elbo[[t-1]] <- data.frame(iter = t, elapsed = difftime(start.time, Sys.time(), units = "secs"), elbo = samps$ELBO,
                                    bayesLikelihood = bayesLik, simLikelihood = samps$simLikelihood, elpd = elpd, DIC = DIC)
        steps[[t-1]] <- updates

        # check for convergence
        if (t > 2){
            elbo.chg <- (samps$ELBO - elbo[[t-1]]$elbo[1])/abs(elbo[[t-1]]$elbo[1])
        }else{
            elbo.chg <- 1.0000
        }
        
        norm_old <- sum(old.params^2)
        norm_chg <- sum((new.params[1,] - old.params[1,])^2)
        lambda <- sqrt(norm_chg)/sqrt(norm_old)

        if((crit == "param" & lambda < converge) | crit == "elbo" & elbo.chg < converge){
            message(paste0("Algorithm converged after ", t, " iterations."))
            converged <- TRUE
            params[(t+1):max_iter] <- NULL
            elbo[t:(max_iter - 1)] <- NULL
            steps[t:(max_iter - 1)] <- NULL
            break
        }

        if(verbose) message(paste0("BBVI-",method,": Iteration ", t,
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
        elapsed.time = difftime(start.time, Sys.time(), units = "secs"),
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

#' Generate the BBVI Update Steps
#'
#' @param dat The generated samples drawn from `generate_samples`
#' 
#' @return A (N + 2) x (K) matrix of update steps, plus a column of parameter labels
generate_bbvi <- function(dat, method = NULL){
    if(is.null(method)) method = dat$method
    Z = dat$Z
    S = dat$S
    
    # Estimate the variance
    var_index <- sample(1:S, size=round(S/3), replace=TRUE)
    var_samples <- Z[var_index, ]
    variances <- apply(var_samples, 2, var)

    # Now estimate the mean
    means <- apply(Z, 2, mean)
    norm_est = sum(means^2)
    m = length(means) - 1

    if(method %in% c("JS+","RB+")){
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
#' @param eta The tuning parameter eta for performing AdaGrad
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

plot_multimix <- function(bbvi){
    N <- bbvi$N 
    K <- bbvi$K 
    p <- bbvi$p

    last_iter <- bbvi$iterations
    params <- bbvi$trace[last_iter, ]
    rownames(params) <- NULL
    this.params <- unpack_parameters(as.matrix(params), N, K, p)
    
    x <- seq(from = min(bbvi$data[,1]), to = max(bbvi$data[,1]), length = 100)
    y <- seq(from = min(bbvi$data[,2]), to = max(bbvi$data[,2]), length = 100)

    par(bg = "white")
    plot(bbvi$data)
    for(k in 1:K){
        mu <- as.numeric(this.params$m[k,])
        sigma <- this.params$ssq[[k]]
        f <- function(x, y) dmvnorm(cbind(x, y), mu, sigma)
        z <- outer(x, y, f)
        contour(x, y, z, add = TRUE, col = "red")
    }
}

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