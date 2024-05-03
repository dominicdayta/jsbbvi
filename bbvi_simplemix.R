#' Random Sample Generation
#'
#' @param n Size of mixture sample to be generated
#' @param props An array of K mixing proportions per component, can be unnormalized
#' @param means An array of K cluster means
#' @param vars An array of K cluster
#' @param seed Random seed (defaults to 923)
#' 
#' @return An array of n observations corresponding to the Gaussian mixture distribution
generate_mixture <- function(n, props, means, vars, seed=923){
    set.seed(seed)
    props <- props/sum(props)
    sizes <- round(n * props, 0)
    
    samp <- NULL
    for(i in seq_along(sizes)){
        samp <- c(samp, rnorm(sizes[i], mean=means[i], sd=sqrt(vars[i])))
    }

    return(samp)
}

#' Complete Joint Log-Likelihood For Naive and JS+ Estimator
#'
#' @param y An array of N observations containing the observed data
#' @param z An N x K matrix of one-hot vector encodings for the K components
#' @param mu An array of K cluster means
#' @param mu0 The prior mean for each mu \sim N(mu0, tausq), defaults to 0
#' @param tausq The prior variance for each mu \sim N(mu0, tausq), defaults to 1
#' @param sigsq The data variance for each y \sim N(zi^T mu, sigsq), defaults to 1
#' @param phi0 The prior probabilities for zi, zi \sim Categorical(phi0), defaults to 1/K
#' 
#' @return The complete joint likelihood p(y, zu, mu | mu0, tausq, sigsq, phi0)
logp.all <- function(y, z, mu, mu0 = 0, tausq = 1, sigsq = 1, phi0 = rep(1/length(mu), times = length(mu))){
    K <- length(mu)
    N <- length(y)

    term1 <- 0
    term2 <- 0
    for(i in 1:N){
        this.z <- as.numeric(z[i, ])
        term1 <- term1 + (y[i] - t(this.z) %*% mu)^2
        term2 <- term2 + sum(this.z * log(phi0))
    }

    logp <- N * log(1/sqrt(2 * pi * sigsq)) - (term1/(2 * sigsq)) + term2 + K * log(1/sqrt(2 * pi * tausq)) - (sum((mu - mu0)^2) / (2 * tausq))
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
logp.mu <- function(y, z, mu, mu0 = 0, tausq = 1, sigsq = 1, phi0 = rep(1/length(mu), times = length(mu))){
    N <- length(y)
    K <- length(mu)
    
    logp <- sapply(1:K, function(j){
        term1 <- 0
        for(i in 1:length(y)){
            if(z[i,j] == 0) next
            term1 <- term1 + log(1/sqrt(2 * pi * sigsq)) - (y[i] - z[i,j] * mu[j])^2 / (2 * sigsq)
        }

        .loglik <- term1
        .logprior <- log(1/sqrt(2 * pi * tausq)) - (mu[j] - m[j])^2 / (2 * tausq)

        return(.loglik + .logprior)
    })

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
logp.z <- function(y, z, mu, mu0 = 0, tausq = 1, sigsq = 1, phi0 = rep(1/length(mu), times = length(mu))){
    K <- length(mu)

    logp <- sapply(1:length(y), function(i){
        -1 * (y[i] - t(z[i,]) %*% mu)^2 / (2 * sigsq)
    })

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
    K <- length(mu)

    term1 <- 0
    for(i in 1:nrow(z)){
        term1 <- term1 + sum(z[i,] * log(phi[i, ]))
    }

    term2 <- sum(sapply(1:length(mu), function(j){
        log(1/sqrt(2 * pi * ssq[j])) - ((mu[j] - m[j])^2)/(2 * ssq[j])
    }))

    logq <- term1 + term2
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

    logq <- log(1/sqrt(2 * pi * ssq)) - (mu - m)^2 / (2 * ssq)

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
    grad <- (mu - m) / ssq
    return(grad)
}

#' Gradient Function for Mu, given ssq
#'
#' @param mu An array of K cluster means
#' @param m An array of K posterior means such that mu \sim q(m, ssq)
#' @param ssq An array of K posterior variances such that mu \sim q(m, ssq)
#' 
#' @return The gradient function for mu, given ssq
grad.mu_ssq <- function(mu, m, ssq){
    grad <- -1 / (2 * ssq) + (mu - m)^2/(2 * ssq^2)
    return(grad)
}

#' Gradient Function for z
#'
#' @param phi An N x K matrix of cluster probabilities such that z \sim q(phi)
#' @param phi An N x K matrix of cluster probabilities such that z \sim q(phi)
#' 
#' @return The gradient function for z
grad.z <- function(z, phi){
    grad <- sapply(1:nrow(z), function(i){
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
generate_samples <- function(S = 1000, y, phi, m, ssq, mu0 = 0, tausq = 1, sigsq = 1, 
                             phi0 = rep(1/length(m), times = length(m)), method="Naive"){

    Ancestor = rbind(
        phi,
        m,
        ssq
    )

    Z <- list(S)
    Samps <- list(S)
    ELBO <- numeric(S)
    for(s in 1:S){
        z.trans <- sapply(1:nrow(phi), function(i){
            x <- rmultinom(1, size=1, prob=phi[i,])
            z <- as.numeric(x[,1])
            return(z)
        })
        z <- t(z.trans)

        mu <- sapply(seq_along(m), function(j){
            rnorm(1, mean=m[j], sd=sqrt(ssq[j]))
        })

        Samps[[s]] <- data.frame(
            iter=s, 
            param=c(paste0("z", 1:nrow(z)), "mu"), 
            rbind(z, mu)
        )

        rownames(Samps[[s]]) <- NULL

        ELBO[s] <- as.numeric(logp.all(y=y, z=z, mu=mu, mu0=mu0, tausq=tausq, sigsq=sigsq, phi0=phi0) - logq.all(z=z, mu=mu, phi=phi, m=m, ssq=ssq))

        if(method %in% c("Naive","JS","JS+")){
            grads <- rbind(
                grad.z(z, phi),
                matrix(grad.mu_m(mu, m, ssq), nrow=1, byrow = TRUE),
                matrix(grad.mu_ssq(mu, m, ssq), nrow=1, byrow = TRUE)
            )

            diff <- as.numeric(logp.all(y=y, z=z, mu=mu, mu0=mu0, tausq=tausq, sigsq=sigsq, phi0=phi0) - logq.all(z=z, mu=mu, phi=phi, m=m, ssq=ssq))

            elbo <- grads * diff

        }else{

            elbo.z <- grad.z(z, phi)
            for(j in 1:ncol(elbo.z)){
                elbo.z[,j] <- as.numeric(elbo.z[,j]) * as.numeric(logp.z(y,z,mu, mu0, tausq, sigsq, phi0) - logq.z(z, mu, phi, m, ssq))
            }

            elbo.mu_m <- as.numeric(grad.mu_m(mu, m, ssq)) * as.numeric(logp.mu(y,z,mu, mu0, tausq, sigsq, phi0) - logq.mu(z,mu, phi, m, ssq))
            elbo.mu_ssq <- as.numeric(grad.mu_ssq(mu, m, ssq)) * as.numeric(logp.mu(y,z,mu, mu0, tausq, sigsq, phi0) - logq.mu(z,mu, phi, m, ssq))
            elbo <- rbind(
                elbo.z,
                matrix(elbo.mu_m, nrow=1),
                matrix(elbo.mu_ssq, nrow=1)
            )

        }

        for(j in 1:ncol(elbo)){
            elbo[, j] <- sapply(as.numeric(elbo[,j]), function(x){
                if(is.nan(x) || is.na(x) || is.infinite(x)) return(runif(1,-1,1))
                return(x)
            })
        }
        
        param <- c(paste0("z", 1:length(y)), "m", "ssq")
        Z[[s]] <- data.frame(iter=s, param=param, elbo)
    }

    list_data <- list(ELBO = mean(ELBO), Ancestor = Ancestor, MC = do.call('rbind', Samps), Z=do.call('rbind', Z), method=method)
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
    S = max(Z$iter)
    params = unique(Z$param)

    # Estimate the variance
    var_index <- sample(1:S, size=round(S/3), replace=TRUE)
    var_samples <- Z[Z$iter %in% var_index, ]
    var_list <- list(length(params))
    for(j in seq_along(params)){
        this.param <- var_samples[var_samples$param == params[j], ]
        this.param$iter <- NULL
        this.param$param <- NULL
        var.param <- matrix(apply(this.param, 2, var), nrow=1)
        var_list[[j]] <- data.frame(param=params[j], var.param)
    }
    var_est <- do.call('rbind', var_list)

    # Now estimate the mean
    mean_list <- list(length(params))
    for(j in seq_along(params)){
        this.param <- Z[Z$param == params[j], ]
        this.param$iter <- NULL
        this.param$param <- NULL
        mean.param <- matrix(apply(this.param, 2, mean), nrow=1)
        mean_list[[j]] <- data.frame(param=params[j], mean.param)
    }
    mean_est <- do.call('rbind', mean_list)

    norm_est = sum(mean_est[,-1]^2)
    m = nrow(mean_est) * (ncol(mean_est) - 1)

    if(method %in% c("JS+","RB+")){

        for(j in 2:ncol(mean_est)){
            vars <- as.numeric(var_est[,j])
            means <- as.numeric(mean_est[,j])
            mean_est[,j] <- positive(1 - ((m - 3) * vars)/norm_est) * means
        }
        
    }

    return(mean_est)
}

#' Constant Rate Function Constructor
#' 
#' Returns a function to use as a learning rate for the BBVI Update Steps
#'
#' @param rho The constant learning rate, at each step rho(t) = rho
#' 
#' @return Returns a function that the BBVI Update step uses for updating parameters
rate_constant <- function(rho = 1e-6){
    rate_function <- function(t, ...){
        return(rho)
    }

    return(rate_function)
}

#' Reciprocal Decay Rate Function Constructor
#' 
#' Returns a function to use as a learning rate for the BBVI Update Steps
#'
#' @param rho The constant learning rate, at each step rho(t) = rho/t
#' 
#' @return Returns a function that the BBVI Update step uses for updating parameters
rate_reciprocal <- function(rho = 1e-6){
    rate_function <- function(t, ...){
        return(rho/t)
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
rate_adagrad <- function(eta = 1){
    rate_function <- function(t, G){
        perturb <- runif(nrow(G), 0.001, 0.01)
        grads <- eta * 1/sqrt(diag(G) + perturb)
        return(grads)
    }

    return(rate_function)
}

simplemix.bbvi <- function(data, K, priors = list(), maxiter=100, init=NULL, learn_rate=rate_constant(1e-6), mc_size=1000, method="Naive",
                            verbose=FALSE, seed=923, thresh=1e-16, beta=0.9){
    
    tausq <- 1
    sigsq <- 1
    mu0 <- 0
    phi0 <- rep(1/K, times = K)

    if(! is.null(priors[["tausq"]])){
        tausq <- as.numeric(priors[["tausq"]])
    }

    if(! is.null(priors[["sigsq"]])){
        sigsq <- as.numeric(priors[["sigsq"]])
    }

    if(! is.null(priors[["mu0"]])){
        mu0 <- as.numeric(priors[["mu0"]])
    }

    if(! is.null(priors[["phi0"]])){
        phi0 <- as.numeric(priors[["phi0"]])
    }

    if(is.null(init)){
        phi <- t(sapply(data, function(x){
            y <- runif(K)
            return(y/sum(y))
        }))

        m <- rnorm(K, 0, 20)

        ssq <- runif(K, min=1, max=10)

        init <- data.frame(
            param = c(paste0("z", 1:nrow(phi)),"m","ssq"),
            rbind(
                phi,
                m,
                ssq
            )
        )

        rownames(init) <- NULL
    }

    params <- init
    paramlist <- list(maxiter)
    elbolist <- list(maxiter)

    if(verbose) message(paste0("Doing Black Box Variational Inference With Method ",method))
    start_time <- Sys.time()
    G <- 0

    orig.method <- method

    set.seed(seed)

    for(t in 1:maxiter){
        
        phi <- as.matrix(params[1:length(data),-1])
        m <- as.numeric(params[nrow(params) - 1,-1])
        ssq <- as.numeric(params[nrow(params), -1])

        samps <- generate_samples(mc_size, y = data, phi = phi, m = m, ssq = ssq, 
                                  mu0 = mu0, tausq = tausq, sigsq = sigsq, phi0 = phi0, method = method)
        updates <- generate_bbvi(samps)

        if(verbose && t %% 100 == 0) message(paste0("BBVI-",method,": Iteration ", t,
         " | Means at ", round(m[1],2), ", ", round(m[2],2), 
         " | Variances at ", round(ssq[1],2), ", ", round(ssq[2],2),
         " | ELBO: ", round(samps$ELBO,2)))

        old.params <- as.numeric(as.matrix(params[,-1]))
        elbo <- as.numeric(as.matrix(updates[,-1]))
        elbolist[[t]] <- data.frame(iter=t, elbo=samps$ELBO)

        #G <- G + elbo %*% t(elbo) # for AdaGrad
        G <- beta * G + (1 - beta) * (elbo %*% t(elbo)) # for AdaGrad

        new.params <- old.params + learn_rate(t, G) * elbo
        params[,-1] <- as.matrix(new.params, nrow = (length(data) + 2), ncol = (ncol(params) - 1))

        # hard thresholding on ssq
        params[nrow(params), -1] <- sapply(as.numeric(params[nrow(params), -1]), function(x){
            if(x < 0.01) return(0.01)
            if(x > 1e3) return(1e3)
            return(x)
        })

        # normalize each column in phi
        for(i in 1:length(data)){
            probs <- sapply(as.numeric(params[i,-1]), function(x){
                if(x < 0.01) return(0.01)
                if(is.na(x) || is.infinite(x) || is.nan(x)) return(0.01)
                return(x)
            })
            total_prob <- sum(probs, na.rm = TRUE)
            params[i,-1] <- probs/total_prob
        }

        paramlist[[t]] <- data.frame(iter=t, time=as.numeric(difftime(start_time, Sys.time(), units="mins")), params)

        # check for convergence
        if(t > 10){
            elbo.sub <- do.call('rbind', elbolist[t:(t-10)])
            elbo.mod <- lm(elbo ~ iter, data=elbo.sub)
            slope <- coef(elbo.mod)[2]
            if(abs(slope) < thresh){
                elbolist[(t+1):maxiter] <- NULL
                paramlist[(t+1):maxiter] <- NULL
                print(paste0("Convergence reach with ELBO change ", round(slope,4)))
                break
            }
        }
    }

    return(list(
        elbo=do.call('rbind', elbolist),
        paths=do.call('rbind', paramlist)
    ))
}

simplemix.logLik <- function(data, iterations){
    N <- length(data)
    K <- ncol(iterations) - 2
    niter <- max(iterations$iter)

    liks <- sapply(1:niter, function(k){
        params <- iterations %>%
            filter(iter == k)
        
        params$iter <- NULL
        params$param <- NULL
        params <- as.matrix(params)

        phi <- params[1:N, ]
        means <- as.numeric(params[N+1,])

        log_prob <- 0
        for(i in 1:N){
            probs <- 0
            for(j in 1:K){
                this_prop <- as.numeric(phi[i, j])
                probs <- probs + this_prop * dnorm(data[i], mean = means[j], sd = 1)
            }

            log_prob <- log_prob + log(probs)
        }

        return(log_prob)
    })

    lik_data <- data.frame(Iteration = 1:niter, logLik = liks)

    return(lik_data)    
}

simplemix.DIC <- function(data, iterations, S=10){
    N <- length(d)
    K <- ncol(iterations) - 2
    niter <- max(iterations$iter)

    liks <- sapply(1:niter, function(k){
        params <- iterations %>%
            filter(iter == k)
        
        params$iter <- NULL
        params$param <- NULL
        params <- as.matrix(params)

        phi <- params[1:N, ]
        means <- as.numeric(params[N+1,])

        pDIC_samps <- replicate(S, {
            phi2 <- t(sapply(1:N, function(i){
                rmultinom(1,1,prob = phi[i, ])
            }))

            means2 <- sapply(means, function(m){
                rnorm(1, mean = m, sd = 1)
            })

            log_prob <- 0
            for(i in 1:N){
                probs <- 0
                for(j in 1:K){
                    this_prop <- as.numeric(phi2[i, j])
                    probs <- probs + this_prop * dnorm(data[i], mean = means2[j], sd = 1)
                }

                log_prob <- log_prob + log(probs)
            }

            return(log_prob)
        })

        pDIC <- mean(pDIC_samps)
    })

    lik_data <- simplemix.logLik(data=d, trace.jsp)
    lik_data$pDIC <- 2 * (lik_data$logLik - liks)
    lik_data$DIC <- -2 * lik_data$logLik + 2 * lik_data$pDIC

    return(lik_data)
}


# logp.all(c(1,2), rbind(c(0,1),c(1,0)), c(-1,1))
# logp.all(c(1,2), rbind(c(0,1),c(1,0)), c(-1,1))
# logq.all(rbind(c(0,1),c(1,0)), c(-1,1), rbind(c(0.1,0.9),c(0.9,0.1)), c(0,0))