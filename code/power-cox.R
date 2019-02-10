###########################################################################
###########################################################################
###                                                                     ###
###                   POWER CALCULATION FOR STCS SNPS                   ###
###                                                                     ###
###########################################################################
###########################################################################

## Email Agnieszka --------------------------------------------------------
## Here are the SNPs we would like to calculate power for
## rs419598 rs16910526 rs5743611 rs4986790 rs4252125 rs5743618 rs5743708 rs5743810 
## rs1800450 rs1800972 rs4257674 rs1554013 rs2243250 rs5743604 rs3775291 rs16944 
## rs11713464 rs1800896 rs1800871 rs1143627 rs2287886 rs4804803 rs1800587 rs3921 
## rs17886395 rs1800629 rs10870077 rs72550870 rs7530511 rs2069705
## The SNPs are coded 0 "homozygous wt" 1"heterozygous" 2"homozygous mutant"
## recessive mode 0 1 vs 2
## dominant mode 0 vs 1 2

## I am sending you ready dataset for infection STCS_IMI.txt, where
## group is                          tx1vc==1
## censoring                        censm==36
## endpoint infection            ifmoldg2==2
## removed patients             ibmoldg2==1

## ------------------------------------------------------------------------


##////////////////////////////////////////////////////////////////
##                             Data                             //
##////////////////////////////////////////////////////////////////

# Read data ----------

dat <- read.table( "STCS_IMI.txt", sep = "\t", header = TRUE)
head(dat)



## Sort out date -----
## mix of french and english dates

f.convert.date <- function(x)
  {
    tmp.date <- as.character(x)
    id.point <- grep("\\.", tmp.date)

    date.point <- tmp.date[id.point]
    date.rest <- tmp.date[-id.point]


    ## switch french to engl
    month.engl <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")
    month.french <- c("janv", "feb", "mars", "apr", "may", "juin", "juil", "aug", "sept", "oct", "nov", "dec")

    for (i in 1:length(month.french))
      {
        id <- grep(month.french[i], date.point)
        date.point[id] <- sub(month.french[i], month.engl[i], date.point[id])
      }

    ## dates <- "15feb2007"
    date.rest2 <- as.Date(date.rest, "%d%b%Y")
    
    ## dates <- c("15.feb.07", "05.juin.85")
    date.point2 <- as.Date(date.point, "%d.%b.%y")
    ## control: sum(is.na(date.point2))


    x2 <- rep(NA, length(tmp.date))
    x2[id.point] <- as.character(date.point2)
    x2[-id.point] <- as.character(date.rest2)
    as.Date(x2)

  
  }

## Convert dates
dat$censdate3p5 <- f.convert.date(dat$censdate3p5)
dat$censor <- f.convert.date(dat$censor)
dat$date <- f.convert.date(dat$date)
dat$txdate <- f.convert.date(dat$txdate)

## Turn NAs into 0
dat$ifmoldg2[is.na(dat$ifmoldg2)] <- 0

## Create binary variables (instead of {1,2} becomes {FALSE, TRUE})
dat$ifmoldg2.2 <- dat$ifmoldg2 > 1
dat$ifmoldg2.12 <- dat$ifmoldg2 > 0
dat$immoldg2.2 <- dat$immoldg2 > 1
dat$immoldg2.12 <- dat$immoldg2 > 0

## Calculate time to event -----------
dat$time.2 <- ifelse(!dat$ifmoldg2.2, dat$censdate3p5 - dat$txdate, dat$date - dat$txdate)
dat$time.12 <- ifelse(!dat$ifmoldg2.12, dat$censdate3p5 - dat$txdate, dat$date - dat$txdate)

## Select SNPs -----------------------
snps <- c("rs419598", "rs16910526", "rs5743611", "rs4986790", "rs4252125", "rs5743618", "rs5743708",
          "rs5743810", "rs1800450", "rs1800972", "rs4257674", "rs1554013", "rs2243250", "rs5743604",
          "rs3775291", "rs16944", "rs11713464", "rs1800896", "rs1800871", "rs1143627", "rs2287886",
          "rs4804803", "rs1800587", "rs3921", "rs17886395", "rs1800629", "rs10870077", "rs72550870",
          "rs7530511", "rs2069705")
snps.r <- paste0(snps, "_R")
snps.d <- paste0(snps, "_D")
snps.all <- c(snps.r, snps.d)

## Remove duplicated lines ------------
dat.single <- dat[!duplicated(dat$id), ]


##////////////////////////////////////////////////////////////////
##                     Power calculation                        //
##////////////////////////////////////////////////////////////////

library(powerSurvEpi)
library(survival)


##::::::::::::::::::::::::::::::::::
##  Variables need to be defined!!!  
##::::::::::::::::::::::::::::::::::


## n: total number of subjects -------------------------
n <- 1103 # nrow(dat.single)  # total number of subjects

## alpha: type I error rate ----------------------------
alpha <- 0.05/length(snps.all)

## theta: postulated hazard ratio ----------------------
theta <- seq(0.1, 20, 0.05) ## <<<<<<<<<<<<<<, set the distance between the thetas very small if you want a smooth plot!

## p: proportion of subjects taking --------------------
## the value one for the covariate of interest.
p <- unique(c(0.33, 0.18, 0.25, 0.05, 0.13, 0.33, 0.30, 0.36, 0.44, 0.35, 0.14, 0.14, 0.41, 0.15, 0.27, 0.07, 0.05, 0.29, 0.02, 0.08, 0.29, 0.07))

## psi: proportion of subjects -------------------------
## died of the disease of interest.
psi <- c(26, 45)/n# apply(dat.single[, c("immoldg2.2", "immoldg2.12")], 2, function(x) sum(x, na.rm = TRUE)/nrow(dat.single)) #0.0009   # proportion of subjects with event

## rho2: square of the correlation between the ---------
## covariate of interest
##          and the other covariate.
rho2 <- 0

## Assemble the parameters -----------------------------
setting <- data.frame(matrix(NA, ncol = 7, nrow = length(alpha)*length(n)*length(theta) * length(p)
                             * length(psi) * length(rho2)))
names(setting) <- c("alpha", "n", "HR", "Freq.risk", "Freq.event", "rho2", "Power")

## Now we loop through all the possible parameter combinations. 
## (not the most elegant, but actually quite efficient
## could also be done with purrr package)
k <- 1
for (A in alpha)
{
    for (N in n)
    {
        for (T in theta)
        {
            for (P in p)
            {
                for (PSI in psi)
                {
                    for (RHO2 in rho2)
                    {
                        power <- powerEpi.default(n = N, theta = T, p = P, psi = PSI, rho2 = RHO2, alpha = A)
                        setting[k, ] <- c(A, N, T, P, PSI, RHO2, power)
                        k <- k + 1
                        cat(k, "\n")
                    } 
                }
            }
        }
    }
}




## /////////////////////////////////////////////
## table version
## /////////////////////////////////////////////
setting$Freq.risk <- as.numeric(as.character(setting$Freq.risk))

tabA <- setting[which(setting$Freq.event < 0.03), ]
tabA <- setting[order(setting$Freq.risk), ]
tabB <- setting[which(setting$Freq.event > 0.03), ]
tabB <- setting[order(setting$Freq.risk), ]

write.table(tabA, file = "results_power_freqA.txt", row.names = FALSE)
write.table(tabB, file = "results_power_freqB.txt", row.names = FALSE)

## /////////////////////////////////////////////
## graph version
## /////////////////////////////////////////////

library(reshape)
library(ggplot2)
theme_set(theme_bw())

setting$Freq.event <- factor(setting$Freq.event)
levels(setting$Freq.event) <- format(as.numeric(levels(setting$Freq.event)), dig = 2)

pdf("power.pdf", width = 14)
qplot(HR, Power,  data = setting, geom = "line") + facet_grid(Freq.event ~ Freq.risk,
                                         labeller=label_both) +  ggtitle(paste("alpha =",
                                         format(alpha, dig = 2), " and n =", n) ) +
    geom_hline(yintercept = 0.8, linetype = 2) + scale_y_continuous(breaks=seq(0, 1, 0.2))  
dev.off()

pdf("power_color.pdf", width = 14)
qplot(HR, Power,  data = setting, geom = "line") + facet_grid(Freq.event ~ Freq.risk,
                                         labeller=label_both) +  ggtitle(paste("alpha =",
                                         format(alpha, dig = 2), " and n =", n) ) +
    geom_line(data = setting[setting$Power > 0.8 & setting$HR < 1, ], aes(HR, Power), color = "red",
              size = 1) +
    geom_line(data = setting[setting$Power > 0.8 & setting$HR > 1, ], aes(HR, Power), color = "red",
              size = 1)+ scale_y_continuous(breaks=seq(0, 1, 0.2))  
dev.off()
