###########################################################################
###########################################################################
###                                                                     ###
###                        POWER CALCULATION SNPS                       ###
###                                                                     ###
###########################################################################
###########################################################################

## packages needed --------------------------

## loaded directly via ::
## dplyr
## import
## purrr
## powerSurvEpi

library(ggplot2)
theme_set(theme_bw())

import::from(magrittr, "%>%")
import::from(magrittr, "%<>%")


##::::::::::::::::::::::::::::::::::
##  Variables need to be defined!!!  
##::::::::::::::::::::::::::::::::::

## number of SNPs --------------------------------------
p <- 60

## n: total number of subjects -------------------------
n <- 1103 # nrow(dat.single)  # total number of subjects

## alpha: type I error rate ----------------------------
alpha <- 0.05/p

## theta: postulated hazard ratio ----------------------
theta <- seq(0.1, 20, by = 0.05) 
## <<<<<<<<<<<<<<, set "by" very small if you want a smooth plot!

## p: proportion of subjects taking --------------------
## the value one for the covariate of interest.
p <- unique(c(0.33, 0.18, 0.25, 0.05, 0.13, 0.33, 0.30, 0.36, 0.44, 0.35, 0.14, 0.14, 0.41, 0.15, 0.27, 0.07, 0.05, 0.29, 0.02, 0.08, 0.29, 0.07))

## psi: proportion of subjects -------------------------
## had an event at of the disease of interest.
psi <- c(26, 45)/n

## rho2: square of the correlation between the ---------
## covariate of interest and the other covariate.
rho2 <- 0

## Assemble all parameter combinations ------------------
param_cmb <- expand.grid(n, theta, p, psi, rho2, alpha)
names(param_cmb) <- c("n", "theta", "p", "psi", "rho2", "alpha")

##////////////////////////////////////////////////////////////////
##                     Power calculation                        //
##////////////////////////////////////////////////////////////////

## Now we loop through all the possible parameter combinations. 

param_cmb %<>% dplyr::mutate(power = purrr::pmap(param_cmb, function(n, theta, p, psi, rho2, alpha) {
  powerSurvEpi::powerEpi.default(
    n = n,
    theta = theta,
    p = p,
    psi = psi,
    rho2 = rho2,
    alpha = alpha
  )
}) %>% unlist())


## /////////////////////////////////////////////
## visualise power
## /////////////////////////////////////////////


## turn psi (freq.event) into factor
param_cmb %<>% dplyr::mutate(psi = factor(psi))
levels(param_cmb$psi) <- format(as.numeric(levels(param_cmb$psi)), dig = 2)

## rename columns
param_cmb_plot <- param_cmb %>% 
  dplyr::rename(HR = theta, Power = power, Freq.risk = p, Freq.event = psi)

## plot
ggplot(data = param_cmb_plot, aes(HR, Power)) + 
  geom_line() + 
  facet_grid(Freq.event ~ Freq.risk, labeller=label_both) + 
  ggtitle(paste("alpha =", format(alpha, dig = 2), " and n =", n)) +
  geom_line(data = param_cmb_plot %>% dplyr::filter(Power > 0.8 & HR < 1), aes(HR, Power), color = "red", size = 1) +
  geom_line(data = param_cmb_plot %>% dplyr::filter(Power > 0.8 & HR > 1), aes(HR, Power), color = "red", size = 1) + 
  scale_y_continuous(breaks=seq(0, 1, 0.2)) 
