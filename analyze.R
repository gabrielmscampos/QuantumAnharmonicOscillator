analyze_simulation <- function(file) {
    UW <- readRDS(file)

    # Simulation parameters
    N_sweep <- UW$N_sweep
    ep <- UW$ep
    h <- UW$h
    L <- UW$L
    N_t <- UW$N_t

    # Acceptance rate
    accrate <- UW$accrate

    # <x2>
    plaq.res2 <- UW$UX2
    UWX2 <- plaq.res2$value
    UerrX2 <- plaq.res2$dvalue
    ererX2 <- plaq.res2$ddvalue
    tauX2 <- plaq.res2$tauint
    dtauX2 <- plaq.res2$dtauint

    # <x4>
    plaq.res4 <- UW$UX4
    UWX4 <- plaq.res4$value
    UerrX4 <- plaq.res4$dvalue
    ererX4 <- plaq.res4$ddvalue
    tauX4 <- plaq.res4$tauint
    dtauX4 <- plaq.res4$dtauint

    # Compute energy at ground level
    E_0 <- ep^2*UWX2 + 3*L*ep^4*UWX4
    ErrE_0 <- ep^2*UerrX2 + 3*L*ep^4*UerrX4

    return(list(
        N_sweep=N_sweep,
        ep=ep,
        h=h,
        L=L,
        N_t=N_t,
        accrate=accrate,
        E_0=E_0,
        ErrE_0=ErrE_0
    ))
}

# Analyze results
files <- list.files("sim_results", pattern = "\\.rds$", full.names = TRUE)
results <- lapply(files, analyze_simulation)
df <- do.call(rbind, lapply(results, as.data.frame))
df

# Plot the energy at the ground level
UW <- readRDS(files[6])
hist(UW$U, 50, prob=TRUE, border="black", col="peachpuff", xlab="x")
lines(density(UW$U), lwd=2, col="chocolate3")

# Plot the equilibrium at related to a N_sweep
plot(1:UW$N_sweep, UW$X2)
