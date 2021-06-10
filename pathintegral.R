library(pracma)
library(parallel)

source("./UWerr.R")

simulateMCMCPathIntegral <- function(sim_idx, N_sweep, ep, h, L) {
    message(sprintf("Simulation - N_sweep: %s, ep: %s, h: %s, L: %s", N_sweep, ep, h, L))

    N_t  <- 250/ep 

    X2 <- replicate(N_sweep, 0) # To store <x^2 > from each path
    X4 <- replicate(N_sweep, 0) # To store <x^4 > from each path
    u <- replicate(N_t+1, 0) # Initialise a cold start
    U <- c()

    for (j in 1:N_sweep){
        # Acceptance rate
        accrate <- 0

        # Site visiting order
        order <- sample(2:N_t, N_t-1)

        # This loop constitutes 1 sweep
        for (i in 1:(N_t-1)){
            k <- order[i]
            rand <- rand()[,1]
            del <- 2*h*(rand-0.5) 
            uc <- u[k] + del 
            DS <- ep*del*((1+ep^2/2)*(uc+u[k])-(u[k-1] + u[k +1]))+L*ep^5*(uc^4-u[k]^4)
            rand <- rand()[,1]
            if (rand < exp(-DS)){
                u[k] <- uc 
                accrate = accrate + 1/(N_t-1)
            }
        }

        X2[j] <- mean(u^2) 
        X4[j] <- mean(u^4) 

        if (mod(j, 100) == 0){
            U <- c(U, ep*u)
        }
    }

    plaq.res2 <- uwerr(data=X2)
    plaq.res4 <- uwerr(data=X4)

    result <- list(
        N_sweep=N_sweep,
        ep=ep,
        h=h,
        L=L,
        N_t=N_t,
        accrate=accrate,
        X2=X2,
        X4=X4,
        U=U,
        UX2=plaq.res2,
        UX4=plaq.res4
    )
    saveRDS(result, sprintf("sim_results/sim_%s.rds", sim_idx))
}

sweeps <- c(seq(500, 10000, 500), 110000)
ep <- c(0.1, 0.5)
h <- c(3)
L <- c(0, 1, 10)
param_grid <- expand.grid(N_sweep=sweeps, ep=ep, h=h, L=L)

mclapply(1:nrow(param_grid), function(i) {
    with(param_grid[i, ], simulateMCMCPathIntegral(i, N_sweep, ep, h, L))
}, mc.cores=8)
