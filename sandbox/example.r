library(dplyr)
library(tidyr)
library(tmvtnorm)
library(survdr)

gAt <- function(U)plogis(- 2 * U)
gRt <- function(t,A,U)plogis(-4 + A + cos(t) * A - sqrt(t) * A * U)
ht  <- function(t,A,U)plogis(-3 + A - 2 * log(t) * U - 0.5 * A * U -
                             0.6 * sin(t) * (A + 1) * U)


genW <-  function(n) rtmvnorm(n, sigma = toeplitz((10:1) / 10),
                              lower = rep(-1.5, 10), upper = rep(1.5, 10))
genUh <-  function(W) {
    U <- abs(W[, 1])^(1/2) * abs(W[, 2])^(1/2) - abs(W[, 10])^(1/2) +
         cos((W[, 5])) - cos((W[, 6])) * cos(W[, 5])
    return(U)
}

genUg <-  function(W) {
    U <- abs(W[, 1])^(1/2) * abs(W[, 10])^(1/2) - abs(W[, 9])^(1/2) +
         cos((W[, 5])) - cos((W[, 7])) * cos((W[, 6]))
    return(U)
}

genData <- function(n, K, uni = TRUE){

    id <- rep(1:n, each =  K)
    m <- rep(1:K, n)

    W <- genW(n)
    Uh <- genUh(W)
    Ug <- genUg(W)

    gA1 <- gAt(Uh)
    A <- rbinom(n, 1, gA1)

    W <- data.frame(W)[id, , drop = FALSE]
    names(W) <- paste0('W', 1:ncol(W))
    A <- A[id]
    Ug <- Ug[id]
    Uh <- Uh[id]

    Rm <- rbinom(n * K, 1, gRt(m, A, Ug))
    Lm <- rbinom(n * K, 1, ht(m, A, Uh))

    df <- data.frame(id, m, A, W, Rm, Lm)

    df <- df %>% group_by(id) %>%
        mutate(tmpR = cumsum(Rm), tmpL = cumsum(Lm),
               tmpRR = cumsum(tmpR), tmpLL = cumsum(tmpL),
               Im = (tmpR == 0) * (tmpLL <= 1),
               Jm = (tmpRR <= 1) * (tmpLL <= 1),
               Lm = Im * Lm, Rm = Jm * Rm) %>%
        select(-tmpR, -tmpRR, -tmpL, -tmpLL)

    df <- ungroup(df)
    return(df)
}

df <- genData(100, 5)

Xnames <- paste('W', 1:10, sep = '')
Xsum <- paste(Xnames, collapse='+')
fitRw <- glm(as.formula(paste('Rm ~ m + A + ', Xsum)),
             data = select(df, -id, -Lm, -Im),
             subset = Jm == 1,
             family = binomial, maxit = 1000)
fitLw <- glm(as.formula(paste('Lm ~ m + A + ', Xsum)),
             data = select(df, -id, -Rm, -Jm),
             subset = Im == 1,
             family = binomial, maxit = 1000)
fitAw <- glm(as.formula(paste('A ~ ', Xsum)),
             data = select(df, -id, -Lm, -Im, -Rm, -Jm),
             subset = m == 1,
             family = binomial, maxit = 1000)
newdf <- df
newdf$A <- 1

df <- mutate(df,
             gAw   = predict(fitAw, newdata = newdf, type = 'response'),
             gRw   = predict(fitRw, newdata = newdf, type = 'response'),
             hw    = predict(fitLw, newdata = newdf, type = 'response'))

tmleest <- tmle(df, 5, 'gAw', 'gRw', 'hw', 0.01)
tmledrest <- tmledr(df, 5, 'gAw', 'gRw', 'hw', 0.01)
