lagg <- function(x)lag(x, default = 1)

logLikelihood<- function(beta, Y, X, offset){

    pi <- plogis(crossprod(t(X), beta)  + offset) # P(Y|A,W)= expit(beta0 + beta1*X1+beta2*X2...)
    pi[pi==0] <- .Machine$double.neg.eps # to prevent taking the log of 0
    pi[pi==1] <- 1 - .Machine$double.neg.eps
    logLike <- sum(Y * log(pi) + (1 - Y) * log(1 - pi))

    return(-logLike)

}

grad <- function(beta, Y, X, offset){

    pi <- plogis(crossprod(t(X), beta) + offset) # P(Y|A,W)= expit(beta0 + beta1*X1+beta2*X2...)
    pi[pi==0] <- .Machine$double.neg.eps # for consistency with above
    pi[pi==1] <- 1 - .Machine$double.neg.eps
    gr <- crossprod(X, Y - pi) # gradient is -residual*covariates

    return(-gr)

}

glmc <- function(Y, X, offset, subset, pp = qp){
    Y <- Y[subset]
    X <- X[subset, , drop = FALSE]
    offset <- offset[subset]
    pp <- 0.99 * pp
    Y  <- (Y - pp) / (1 - 2 * pp)
    offset <- qlogis((plogis(offset) - pp) / (1 - 2 * pp))
    optim(par = rep(0, dim(X)[2]), fn = logLikelihood, gr = grad,
          Y = Y, X = X, offset = offset, method = "BFGS")$ par
}

glmclo <- function(Y, X, offset, subset, pp = qp){
    Y <- Y[subset]
    X <- X[subset, , drop = FALSE]
    offset <- offset[subset]
    pp <- 0.99 * pp
    Y      <- (Y - pp) / (1 - pp)
    offset <- qlogis((plogis(offset) - pp) / (1 - pp))
    optim(par = rep(0, dim(X)[2]), fn = logLikelihood, gr = grad,
          Y = Y, X = X, offset = offset, method = "BFGS")$par
}

glmcup <- function(Y, X, offset, subset, pp = qp){
    Y <- Y[subset]
    X <- X[subset, , drop = FALSE]
    offset <- offset[subset]
    pp <- 0.99 * pp
    Y      <- Y / (1 - pp)
    offset <- qlogis(plogis(offset) / (1 - pp))
    optim(par = rep(0, dim(X)[2]), fn = logLikelihood, gr = grad,
          Y = Y, X = X, offset = offset, method = "BFGS")$par
}

meanlist <- function(x){
    y <- sapply(x, I)
    if(is.null(ncol(y))) y <- matrix(y, ncol = length(y))
    rowMeans(y)
}

varlist  <- function(x)diag(var((t(sapply(x, I)))))


bound01 <- function(x, p = 1e-6){
    xx <- pmax(pmin(x, 1 - p), p)
    as.numeric(xx)
}

boundlo <- function(x, p = 0.01){
    xx <- pmax(x, p)
    as.numeric(bound01(xx))
}

boundup <- function(x, p = 0.01){
    xx <- pmin(x, 1 - p)
    as.numeric(bound01(xx))
}

npregbw2 <- function(formula, data, subset = rep(TRUE, nrow(data)), ...){

    xnames <- all.vars(formula[[3]])

    constantVars <- data %>% ungroup() %>% filter(subset) %>%
        select_(.dots = xnames) %>% var() %>%
        diag() == 0

    if(all(constantVars)) {
        formula <- update(formula, . ~ 1)
        environment(formula) <- environment()
        bw <- lm(formula, data = data, subset = subset)
    } else {
        form <- as.formula(paste0('. ~ ', paste(xnames[!constantVars],
                                                collapse = ' + ')))
        formula <- update(formula, form)
        environment(formula) <- environment()
        fit <- npregbw(formula, data = data, subset = subset, ...)
        bw <- npreg(fit)
    }

    return(bw)

}

