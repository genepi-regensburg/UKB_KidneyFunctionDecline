
source("/path/to/simulations_helpers.r")

# global parameters ----
set.seed(123)
nSims <- 10000
MAF <- 0.3
ncores <- 200
alpha <- 0.05
dir_out <- "/path/to/output/directory"

lengths <- c(4000, 3200, 2560, 2048, 8192)
mu1 <- 4
sd1 <- 1
t2 <- t3 <- t4 <- t5 <- 0.5

models <- list(
  "difference model" = c(paste0("y_diff ~ genotype") %>% as.formula(), "genotype"),
  "time model RI&RS" = c(paste0("y ~ time * genotype + (1+time|id)") %>% as.formula(), "time:genotype"),
  "age model RI&RS" = c(paste0("y ~ age * genotype + (1+age|id)") %>% as.formula(), "age:genotype"),
  "age model RI&RS uncorr." = c(paste0("y ~ age * genotype + (1+age||id)") %>% as.formula(), "age:genotype"),
  "age model RI-only" = c(paste0("y ~ age * genotype + (1|id)") %>% as.formula(), "age:genotype"),
  "BLUPs&LinReg" = c(paste0("y ~ age + (1+age|id)") %>% as.formula(), paste0("RS ~ genotype") %>% as.formula(), "genotype"),
  "age model RI&RS with singletons" = c(paste0("y ~ age * genotype + (1+age|id)") %>% as.formula(), "age:genotype")
)

# biobank UKB (cannot be shared due to data protection) ----
## age trajectories
# pheno <- read_file("/path/to/UKB/pheno/file.txt") # cannot be shared due to data protection

# empiricalTrajectories <- pheno %>%
#   arrange(id, age_c) %>% # Arrange by id and then by age to ensure the age trajectory is in order
#   split(.$id) %>% # Split the data frame into a list of data frames, each containing one individual's data
#   lapply(function(individual_df) individual_df$age_c) # For each individual's data frame, extract the age column as a vector

## parameters
# sigmaRE <- matrix(c(98.1, -32.9, -32.9, 52.1), 2, 2)
# intercept <- 100
# ageEffect <- -10
# trueInterceptEffect <- -0.25
# trueSlopeEffect <- -0.25 # or 0 for T1E simulations
# sigma <- 7

## run
# set.seed(1)

# res <- simulate(
#     nSims = nSims, ncores = ncores,
#     datatype = "empirical",
#     empiricalTrajectories = empiricalTrajectories, n = 20000,
#     sigmaRE = sigmaRE, MAF = MAF, intercept = intercept, ageEffect = ageEffect, trueSlopeEffect = trueSlopeEffect,
#     trueInterceptEffect = trueInterceptEffect, sigma = sigma, models = models, alpha = alpha, dir_out = dir_out, return_results=T)

# cohort KORA eGFR ----
## parameters
sigmaRE <- matrix(c(82.5, 2.1, 2.1, 10.4), 2, 2)
intercept <- 110
ageEffect <- -10
trueInterceptEffect <- -0.1
trueSlopeEffect <- -0.1 # or 0 for T1E simulations
sigma <- 3

## run
set.seed(3)

res <- simulate(
    nSims = nSims, ncores = ncores,
    datatype = "cohort",
    lengths = lengths, mu1 = mu1, sd1 = sd1, t2 = t2, t3 = t3, t4 = t4, t5 = t5,
    sigmaRE = sigmaRE, MAF = MAF, intercept = intercept, ageEffect = ageEffect, trueSlopeEffect = trueSlopeEffect,
    trueInterceptEffect = trueInterceptEffect, sigma = sigma, models = models, alpha = alpha, dir_out = dir_out, return_results = T)

# cohort KORA BMI 2.0 ----
## parameters
sigmaRE <- matrix(c(19.0, 1.2, 1.2, 1.4), 2, 2)
intercept <- 27
ageEffect <- 0.5
trueInterceptEffect <- 0.04
trueSlopeEffect <- 0.04 # or 0 for T1E simulations
sigma <- sqrt(1.6)

## run
set.seed(4)

res <- simulate(
    nSims = nSims, ncores = ncores,
    datatype = "cohort",
    lengths = lengths, mu1 = mu1, sd1 = sd1, t2 = t2, t3 = t3, t4 = t4, t5 = t5,
    sigmaRE = sigmaRE, MAF = MAF, intercept = intercept, ageEffect = ageEffect, trueSlopeEffect = trueSlopeEffect,
    trueInterceptEffect = trueInterceptEffect, sigma = sigma, models = models, alpha = alpha, dir_out = dir_out, return_results = T)
