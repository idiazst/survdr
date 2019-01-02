tmledr <- function(df, tau, gA, gR, h, qp){

    n <- length(unique(df$id))
    library(np)

    df <- df %>%
        mutate_(gR = gR, gA = gA, h = h) %>%
        group_by(id) %>%
        arrange(id, m) %>%
        mutate(Gm = boundlo(cumprod(1 - gR), qp),
               Sm = boundlo(cumprod(1 - h), qp),
               Cgm = boundlo(gA * Gm, qp),
               Stau = nth(Sm, tau),
               gA = boundlo(gA, qp),
               gR = boundup(gR, qp),
               h  = boundup(h, qp),
               Chm = Stau / Sm,
               M = sum(Sm * (m <= tau)))

    eRfit <- npregbw2((Rm - gR) / Cgm ~ m + Chm, data = df,
                    subset = df$A == 1 & df$Jm == 1)
    eLfit <- npregbw2((Lm - h) * Chm ~ m + Cgm, data = df,
                    subset = df$A == 1 & df$Im == 1)
    eAfit <- npregbw2((A - gA) / gA ~ M, data = df,
                    subset = df$m == 1)

    eR <- predict(eRfit, newdata = df)
    eA <- predict(eAfit, newdata = df)
    eL <- predict(eLfit, newdata = df)

    df <- df %>% bind_cols(eR = eR, eA = eA, eL = eL)

    qfit <- df %>% subset(m == 1) %>% npregbw2(A ~ Stau, data = .)

    df <- df %>% mutate(Dmm = 0, Bmm = 0, Umm = 0, Vmm = 0, HR = 0, HL2 = 0) %>%
        bind_cols(q = boundlo(predict(qfit, newdata = df), qp))

    for(k in 1:tau){

        dff <- df %>% mutate(Chk = Stau / ifelse(k == 1, 1, nth(Sm, k - 1)),
                             Cgk = gA * nth(Gm, k))
        dfit  <- npregbw2(Rm ~ m + Chk, data = dff, subset = dff$A == 1 & dff$Jm == 1)
        bfit  <- npregbw2(Lm ~ m + Chk, data = dff, subset = dff$A == 1 & dff$Im == 1)
        ufit  <- npregbw2(Rm ~ m + Cgk, data = dff, subset = dff$A == 1 & dff$Jm == 1)
        vfit  <- npregbw2(Lm ~ m + Cgk, data = dff, subset = dff$A == 1 & dff$Im == 1)

        Dmk <- paste0('Dm', k)
        Bmk <- paste0('Bm', k)
        Umk <- paste0('Um', k)
        Vmk <- paste0('Vm', k)
        Gmk <- paste0('Gm', k)

        df <- df %>% bind_cols(!!paste0('d', k) := predict(dfit, newdata = dff),
                               !!paste0('b', k) := predict(bfit, newdata = dff),
                               !!paste0('u', k) := predict(ufit, newdata = dff),
                               !!paste0('v', k) := predict(vfit, newdata = dff)) %>%
            mutate(!!Dmk := boundlo(cumprod(1 - !!as.name(paste0('d', k))), qp),
                   !!Bmk := boundlo(cumprod(1 - !!as.name(paste0('b', k))), qp),
                   !!Umk := boundlo(cumprod(1 - !!as.name(paste0('u', k))), qp),
                   !!Vmk := boundlo(cumprod(1 - !!as.name(paste0('v', k))), qp),
                   !!Gmk := (k == m) * Gm,
                   Dmm := Dmm + (k == m) * !!as.name(Dmk),
                   Bmm := Bmm + (k == m) * !!as.name(Bmk),
                   Umm := Umm + (k == m) * !!as.name(Umk),
                   Vmm := Vmm + (k == m) * !!as.name(Vmk),
                   HR := HR + (m > k) * nth(lagg(Vmm), k) / !!as.name(Vmk) *
                       nth(Umm, k) / !!as.name(Umk) * nth(eL, k) / nth(Gm, k),
                   HL2 := HL2 + (m <= k) * nth(lagg(Bmm), k) / lagg(!!as.name(Bmk)) *
                       nth(Dmm, k) / !!as.name(Dmk) * nth(eR, k) / nth(lagg(Sm), k))
    }

    ## Create auxiliary variables

    df <- df %>%
        mutate(HA = sum(Umm / boundlo(Gm * gA, qp) * lagg(Vmm) * eL),
               HR = HR * lagg(Gm) / boundlo(gA * Gm, qp),
               HL2 = lagg(Sm) * nth(Sm, tau) / Sm * (HL2  +
                                                     eA / boundlo(q * Dmm * lagg(Bmm), qp)),
               HL1 = nth(Sm, tau) / boundlo(Sm * Gm * gA, qp))

    iter <- 1
    crit <- TRUE

    while(crit) {

        fitA <- glmclo(df$A, matrix(df$HA), offset = qlogis(df$gA),
                                subset = df$m == 1, qp)
        fitR <- glmcup(df$Rm, matrix(df$HR), offset = qlogis(df$gR),
                                subset = df$A == 1 & df$Jm == 1, qp)
        fitL <- glmcup(df$Lm, cbind(df$HL1, df$HL2), offset = qlogis(df$h),
                                subset = df$A == 1 & df$Im == 1, qp)

        df <- df %>% select(-gR, -gA, -h, -Gm, -Sm) %>%
            bind_cols(gA = qp + (1 - qp) *
                          plogis(qlogis((df$gA - qp) / (1 - qp)) + fitA * df$HA),
                      gR = (1 - qp) * plogis(qlogis(df$gR / (1 - qp)) + fitR * df$HR),
                      h  = (1 - qp) * plogis(qlogis(df$h / (1 - qp)) +
                                                      fitL[1] * df$HL1 + fitL[2] * df$HL2)) %>%
            group_by(id) %>%
            arrange(id, m) %>%
            mutate(Gm = boundlo(cumprod(1 - gR), qp),
                   Sm = boundlo(cumprod(1 - h), qp),
                   HR = 0, HL2 = 0)

        for(k in 1:tau){

            Dmk <- paste0('Dm', k)
            Bmk <- paste0('Bm', k)
            Umk <- paste0('Um', k)
            Vmk <- paste0('Vm', k)
            Gmk <- paste0('Gm', k)

            df <- df %>%
                mutate(HR := HR + (m > k) * nth(lagg(Vmm), k) / !!as.name(Vmk) *
                           nth(Umm, k) / !!as.name(Umk) * nth(eL, k) / nth(Gm, k),
                       HL2 := HL2 + (m <= k) * nth(lagg(Bmm), k) / lagg(!!as.name(Bmk)) *
                           nth(Dmm, k) / !!as.name(Dmk) * nth(eR, k) / nth(lagg(Sm), k))

        }

        df <- df %>%
            mutate(HA = sum(Umm / boundlo(Gm * gA, qp) * lagg(Vmm) * eL),
               HR = HR * lagg(Gm) / boundlo(gA * Gm, qp),
               HL2 = lagg(Sm) * nth(Sm, tau) / Sm * (HL2  +
                                                     eA / boundlo(q * Dmm * lagg(Bmm), qp)),
               HL1 = nth(Sm, tau) / boundlo(Sm * Gm * gA, qp))


        tmle <- df %>% subset(m == tau) %>% pull(Sm) %>% mean()

        dfD <- summarise(df,
                         D  = -sum(A * Im * Stau * (Lm - h) / boundlo(gA * Gm * Sm,  qp) *
                                   (m <= tau)) + nth(Sm, tau) - tmle,
                         DR = -sum(A * Jm * HR * (Rm - gR) * (m <= tau)),
                         DL = -sum(A * Im * HL2 * (Lm - h) * (m <= tau)),
                         DA = -mean(HA * (A - gA)))


        crit <- dfD %>% select(D, DL, DR, DA) %>% summarise_all(mean) %>% abs() %>% max() >
            1e-4 * n^(-0.6) & iter < 20

        iter <- iter + 1

    }

    var <- with(dfD, var(D - DL - DR - DA)) / n

    return(list(estimate = c(tmle = tmle, sd = sqrt(var)),
                infun = with(dfD, D - DL - DR - DA)))

}

tmle <- function(df, tau, gA, gR, h, qp){

    n <- length(unique(df$id))

    df <- df %>%
        mutate_(gR = gR, gA = gA, h = h) %>%
        mutate(gA = boundlo(gA, qp),
               gR = boundup(gR, qp),
               h  = boundup(h, qp)) %>%
        group_by(id) %>%
        arrange(id, m)

    iter <- 1
    crit <- TRUE

    while(crit) {

        df <- df %>%
            mutate(Gm = boundlo(cumprod(1 - gR), qp),
                   Sm = boundlo(cumprod(1 - h), qp),
                   Stau = nth(Sm, tau),
                   HL = Stau / (Sm * Gm * gA))

        fitL <- glmcup(df$Lm, matrix(df$HL), offset = qlogis(df$h),
                                subset = df$A == 1 & df$Im == 1, qp)

        df <- df %>% select(-h, -Sm) %>%
            bind_cols(h  = (1 - qp) * plogis(qlogis(df$h / (1 - qp)) +
                                                      fitL * df$HL)) %>%
            group_by(id) %>%
            arrange(id, m) %>%
            mutate(Sm = boundlo(cumprod(1 - h), qp),
                   Stau = nth(Sm, tau))

        tmle <- df %>% subset(m == tau) %>% pull(Sm) %>% mean()

        dfD <- summarise(df,
                         D  = -sum(A * Im * Stau * (Lm - h) / (gA * Gm * Sm) *
                                   (m <= tau)) + nth(Sm, tau) - tmle)

        crit <- dfD %>% pull(D) %>% mean() %>% abs() > 1e-4 * n^(-0.6) & iter < 20

        iter <- iter + 1
    }

    var <- with(dfD, var(D)) / n

    return(list(estimate = c(tmle = tmle, sd = sqrt(var)),
                infun = with(dfD, D)))

}
