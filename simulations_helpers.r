library(data.table)
library(dplyr)
library(lme4)
library(lmerTest)
library(mvtnorm)
library(doParallel)
library(foreach)

# functions ----
read_file <- function(a_infile, header="auto",sep="auto",fill=FALSE, ...){
	tblInRaw <- fread(a_infile, header=header, sep=sep)
	tblIn <- as.data.frame(tblInRaw,fill=fill)
	return(tblIn)
}

createAgeTrajectories <- function(lengths, mu1, sd1, t2, t3=0, t4=0, t5=0) {
  n <- sum(lengths)
  result_list <- vector("list", n)
  lengths_vec <- rep(1:5, lengths)
  num_timepoints <- sum(lengths != 0) # assuming no "gaps" in the lengths vector

  start_index <- 1
  for (i in 1:num_timepoints) {
    end_index <- start_index + lengths[i] - 1
    for (j in start_index:end_index) {
      result_list[[j]] <- numeric(i)
      result_list[[j]][1] <- rnorm(1, mu1, sd1)
      if (i > 1) result_list[[j]][2] <- result_list[[j]][1] + t2
      if (i > 2) result_list[[j]][3] <- result_list[[j]][2] + t3
      if (i > 3) result_list[[j]][4] <- result_list[[j]][3] + t4
      if (i > 4) result_list[[j]][5] <- result_list[[j]][4] + t5
    }
    start_index <- end_index + 1
  }

  result_list <- setNames(result_list, seq_along(result_list))
  return(result_list)
}

createData <- function(trajectories, sigmaRE, MAF, intercept, ageEffect, trueInterceptEffect, trueSlopeEffect, sigma=1, n=NULL){

    if(is.null(n)) n <- length(trajectories)

    # sample trajectories
    subjects <- sample(names(trajectories), n)
    sampleTrajectories <- trajectories[subjects]
    m <- length(unlist(sampleTrajectories))

    REs <- mvtnorm::rmvnorm(n, mean = c(0,0), sigma = sigmaRE)

    genotypes <- rbinom(n, 2, MAF)

    dataList <- lapply(seq_along(sampleTrajectories), function(i) {

        subjectTrajectory <- sampleTrajectories[[i]]
        len <- length(subjectTrajectory)
        first_age <- min(subjectTrajectory)
        time <- subjectTrajectory - first_age

        data.frame(
            id = rep(subjects[i], len),
            y = (intercept + REs[i, 1] + genotypes[i] * trueInterceptEffect) + subjectTrajectory * (ageEffect + REs[i, 2] + genotypes[i] * trueSlopeEffect) + rnorm(len, 0, sigma),
            age = subjectTrajectory,
            genotype = rep(genotypes[i], len),
            time = time
        )
    })

    data <- do.call(rbind, dataList)
    rownames(data) <- NULL

    return(data)

}

runModel <- function(df, model, f, trueSlopeEffect){

    if (!endsWith(model, "singletons")) {df <- df %>% group_by(id) %>% filter(n() >= 2) %>% ungroup()}

    if (!model %in% c("difference model", "BLUPs&LinReg")) {
        mod <- lmer(formula=f[[1]], data=df, control=lmerControl(optimizer = 'bobyqa'))
        sumStats <- summary(mod)$coefficients[f[[2]], c(1, 2, 5)]
    } else if (model == "difference model") {
        df_diff <- df %>% group_by(id) %>%
            summarise(y_diff = (last(y) - first(y)) / (last(age) - first(age)), genotype = first(genotype), age_first = first(age))
        mod <- lm(formula=f[[1]], data=df_diff)
        sumStats <- summary(mod)$coefficients[f[[2]], c(1, 2, 4)]
    } else {
        null_mod <- lmer(formula=f[[1]], data=df, control=lmerControl(optimizer = 'bobyqa'))
        blups <- coef(null_mod)$id %>% mutate(id = rownames(.) %>% as.character(), age = age - fixef(null_mod)["age"]) %>% rename(RI = `(Intercept)`, RS = age)
        df_BL <- df %>% group_by(id) %>% slice(1) %>% ungroup() %>% left_join(blups, by="id")
        mod <- lm(formula=f[[2]], data=df_BL)
        sumStats <- summary(mod)$coefficients[f[[3]], c(1, 2, 4)]
    }

    out <- data.frame(
        model = model,
        coef = sumStats[1],
        se = sumStats[2],
        p = sumStats[3]
    )

    out

}

simulate <- function(
    nSims, ncores,
    datatype = "cohort", # "cohort" or "empirical"
    lengths=NULL, mu1=NULL, sd1=NULL, t2=NULL, t3=NULL, t4=NULL, t5=NULL,
    empiricalTrajectories=NULL, n=NULL,
    sigmaRE, MAF, intercept, ageEffect, trueSlopeEffect, trueInterceptEffect, sigma=1,
    models,
    alpha=0.05,
    dir_out="", return_results=FALSE){

    # results loop
    registerDoParallel(cores = ncores)

    results <- foreach(i = 1:nSims, .combine = 'rbind', .packages = c("dplyr", "lme4", "lmerTest", "mvtnorm")) %dopar% {
        if(datatype == "cohort") {
            ageTrajectories <- createAgeTrajectories(lengths, mu1, sd1, t2, t3, t4, t5)
        } else if(datatype == "empirical") {
            ageTrajectories <- sample(empiricalTrajectories, n, replace = FALSE)
        }
        df <- createData(
            trajectories = ageTrajectories, sigmaRE = sigmaRE, MAF = MAF, intercept = intercept, ageEffect = ageEffect,
            trueSlopeEffect = trueSlopeEffect, trueInterceptEffect = trueInterceptEffect, sigma = sigma, n = NULL)
        res_list <- lapply(names(models), function(model) {
            runModel(df = df, model = model, f = models[[model]], trueSlopeEffect = trueSlopeEffect)})
        res <- do.call(rbind, res_list)
        return(res)}

    stopImplicitCluster()

    # summarize results

    summary_df <- results %>%
        group_by(model) %>%
        summarise(
            avg_coef = mean(coef, na.rm = TRUE),
            se_avg_coef = sd(coef, na.rm = TRUE) / sqrt(n()),
            ci_avg_coef = paste0(round(avg_coef - 1.96 * se_avg_coef, 3), "; ", round(avg_coef + 1.96 * se_avg_coef, 3)),
            power = sum(p < alpha, na.rm = TRUE) / n(),
            power_ci_lower = binom.test(sum(p < alpha), n(), 0.5)$conf.int[1],
            power_ci_upper = binom.test(sum(p < alpha), n(), 0.5)$conf.int[2],
            power_ci = paste0(power_ci_lower, "; ", power_ci_upper)
        ) %>%
        arrange(match(model, names(models)))

    path_out_summary <- file.path(
        dir_out,
        paste0(
            "powerSims_nSims", nSims,
            # "_mu1", mu1, "_sd1", sd1, "_t2", t2, "_t3", t3, "_t4", t4, "_t5", t5,
            datatype,
            "_sigmaRE", sigmaRE[1,1], "+", sigmaRE[2,2], "+", sigmaRE[2,1], "_MAF", MAF, "_intercept", intercept, "_ageEffect", ageEffect,
            "_trueSlopeEffect", trueSlopeEffect, "_trueInterceptEffect", trueInterceptEffect, "_sigma", sigma,
            "_alpha", alpha,
              ".txt"))

    path_out_raw <- file.path(
        dir_out,
        paste0(
            "powerSims_nSims", nSims,
            # "_mu1", mu1, "_sd1", sd1, "_t2", t2, "_t3", t3, "_t4", t4, "_t5", t5,
            datatype,
            "_sigmaRE", sigmaRE[1,1], "+", sigmaRE[2,2], "+", sigmaRE[2,1], "_MAF", MAF, "_intercept", intercept, "_ageEffect", ageEffect,
            "_trueSlopeEffect", trueSlopeEffect, "_trueInterceptEffect", trueInterceptEffect, "_sigma", sigma,
            "_alpha", alpha,
            "_rawData",
              ".txt"))

    fwrite(summary_df, path_out_summary, sep="\t")
    fwrite(results, path_out_raw, sep="\t")

    if(return_results) return(summary_df)

}
