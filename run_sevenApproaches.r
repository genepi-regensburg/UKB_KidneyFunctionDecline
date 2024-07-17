library(dplyr)
library(lmerTest)
library(ggplot2)
library(nephro)
if(packageVersion('nephro')<'1.4') stop("Please install R package nephro_v1.4.")

#########################################################
#########################################################
### 1. Raw data preparation

# We here assume input data: "raw_input.txt" with
##  - one row per person and crea measurement (long format; please also include persons with only one crea measurement)
##  - also including snp genotypes for the genetic variants of interest (long format)
## We assume the following input columns:
##   id: person id
##   crea_at_exam: serum crea measurement (in mg/dl) at the exam
##   sex: "m" for men, "w" for women
##   age_at_exam: age at exam
##   PC1, ..., PC20; 20 PCs (strechted to long format)
##   rs1, rs2, ... (rsids of the genetic variants): dosages of variants (long format)
## read data
tblin <- read.table("raw_input.txt",header=T)
arsid <- names(tblin)[grepl("^rs",names(tblin))] # grep array of rsids from tbl header

#########################################################
#########################################################
### 2. Validation R script
### 2.1 Additional phenotype steps

print("Additional phenotype steps ... ")

if(!all(tblin$sex %in% c("m","w"))) stop("Please code sex as m/w.")
tblin$sex01 <- ifelse(tblin$sex=="w",0,1)

# (1) obtain eGFR
tblin$eGFR <- CKDEpi_RF.creat(tblin$crea_at_exam, tblin$sex01, tblin$age_at_exam)

# (2) exclude value after egfr<15 observed
tblin <- tblin %>%
    arrange(id, age_at_exam) %>%
    group_by(id) %>%
    mutate(
        eGFRcrea = case_when(
            is.na(eGFR) ~ NA_real_,
            age == first(age) ~ eGFRcrea, # always keep first value (necessary until same-day-measurements are addressed)
            (eGFR < 15) & (!lag(cumany(eGFR < 15)) | is.na(lag(cumany(eGFR < 15)))) ~ eGFR,
            cumany(eGFR < 15) ~ NA_real_,
            TRUE ~ eGFR
        )
    ) %>%
    filter(!is.na(eGFR)) %>%
    ungroup()

# (3) center/scale age_at_exam
tblin$age_at_exam_c50_10 = (tblin$age_at_exam-50)/10

# (4) obtain residuals from age model
formula <- paste0("eGFR ~ sex + age_at_exam_c50_10 + ", paste(paste0("PC", 1:20), collapse = " + "), " + (1+age_at_exam_c50_10|id)")

mod <- lmer(
    formula = as.formula(formula),
    data = tblin,
    REML = TRUE,
    control = lmerControl(optimizer = 'bobyqa')
    )

res <- residuals(mod, type = "response")

# (5) remove outliers
tblin <- tblin %>%
    mutate(res = res) %>%
    filter(abs(res) <= 10*sd(res)) %>%
	select(-res)

# (6) winsorize
tblin$eGFR[tblin$eGFR>200] <- 200
tblin$eGFR[tblin$eGFR<15] <- 15

# (7) obtain age at baseline and time since baseline and rescale
tblin = tblin[order(tblin$id,tblin$age_at_exam),]
tblin_baseline = tblin[!duplicated(tblin$id),c("id","age_at_exam")]
names(tblin_baseline)[2] = "age_at_baseline"
tblin = merge(tblin,tblin_baseline,by="id")
tblin$time_since_baseline = tblin$age_at_exam - tblin$age_at_baseline

tblin$age_at_baseline_c50_10 = (tblin$age_at_baseline-50)/10
tblin$time_since_baseline_10 = (tblin$time_since_baseline)/10

# (8) save
fwrite(tblin,"pheno.txt",sep="\t",quote=F,row.names=F)

#########################################################
### 2.2 Statistical analyses

for(rsid in arsid) {

	print(paste0("Running decline models on ",rsid," ..."))

	tbl <- tblin
	# change column name to "G"
	idxg <- which(names(tbl)==rsid)
	names(tbl)[idxg] <- "G"

	isRemoved <- is.na(tbl$age_at_exam)|is.na(tbl$eGFR)|is.na(tbl$G)|is.na(tbl$sex)
	if(any(isRemoved)) tbl <- tbl[!isRemoved,]

	lsout <- list()

	lsout$n_all <- length(unique(tbl$id))
	lsout$eaf = sum(tbl$G)/2/lsout$n_all

	tbl_restricted <- tbl %>% group_by(id) %>% filter(n() > 1 & (last(time_since_baseline)-first(time_since_baseline) >= 1)) %>% ungroup()
	lsout$n <- length(unique(tbl_restriced$id))
	lsout$m <- nrow(tbl_restriced)

	## time model RI&RS (restricted)
	formula_time_model_rirs <- paste0("eGFR ~ sex + age_at_baseline_c50_10 + time_since_baseline_10*G + ", paste(paste0("PC", 1:20), collapse = " + "), " + (1 + time_since_baseline_10 | id)")

	lsout$time_model_rirs <- summary(lmer(
		formula = as.formula(formula_time_model_rirs),
		data = tbl_restricted,
		control = lmerControl(optimizer = 'bobyqa')))

	## age model RI&RS (restricted)
	formula_age_model_rirs <- paste0("eGFR ~ sex + age_at_exam_c50_10*G + ", paste(paste0("PC", 1:20), collapse = " + "), " + (1 + age_at_exam_c50_10 | id)")

	lsout$age_model_rirs <- summary(lmer(
		formula = as.formula(formula_age_model_rirs),
		data = tbl_restricted,
		control = lmerControl(optimizer = 'bobyqa')))

	## age model RI&RS uncorrelated (restricted)
	formula_age_model_rirs_uncorrelated <- paste0("eGFR ~ sex + age_at_exam_c50_10*G + ", paste(paste0("PC", 1:20), collapse = " + "), " + (1 + age_at_exam_c50_10 || id)")

	lsout$age_model_rirs_uncorrelated <- summary(lmer(
		formula = as.formula(formula_age_model_rirs_uncorrelated),
		data = tbl_restricted,
		control = lmerControl(optimizer = 'bobyqa')))

	## age model RI (restricted)
	formula_age_model_ri_only <- paste0("eGFR ~ sex + age_at_exam_c50_10*G + ", paste(paste0("PC", 1:20), collapse = " + "), " + (1 | id)")

	lsout$age_model_ri_only <- summary(lmer(
		formula = as.formula(formula_age_model_ri_only),
		data = tbl_restricted,
		control = lmerControl(optimizer = 'bobyqa')))

	## age model RS (unrestricted)
	lsout$age_model_rirs_unrestricted <- summary(lmer(
		formula = as.formula(formula_age_model_rirs),
		data = tbl,
		control = lmerControl(optimizer = 'bobyqa')))

	# (implicit) two-stage models (where eGFR decline is the phenotype of the second-stage model) ----

	## difference model
	tbl_diff <- tbl_restricted %>%
	  arrange(id, time_since_baseline) %>%
	  group_by(id) %>%
	  summarise(
		# Get the difference in eGFR between last and first observation
		eGFR_diff = eGFR[last(which(!is.na(eGFR)))] - eGFR[1],
		# Extract geno_id at first observation
		id = first(id),
		# Extract age at first observation
		age_at_baseline = first(age_at_baseline),
		# Extract sex at first observation
		sex = first(sex),
		# Compute difference in t between last and first observation
		t_diff = time_since_baseline[last(which(!is.na(time_since_baseline)))] - time_since_baseline[1],
		# Extract PC1 thru PC20 at first observation
		PC1 = first(PC1), PC2 = first(PC2), PC3 = first(PC3), PC4 = first(PC4), PC5 = first(PC5), PC6 = first(PC6), PC7 = first(PC7),
		PC8 = first(PC8), PC9 = first(PC9), PC10 = first(PC10), PC11 = first(PC11), PC12 = first(PC12), PC13 = first(PC13), PC14 = first(PC14),
		PC15 = first(PC15), PC16 = first(PC16), PC17 = first(PC17), PC18 = first(PC18), PC19 = first(PC19), PC20 = first(PC20),
		# Extract genotype at first observation
		G = first(G)
	  ) %>%
	  ungroup() %>%
	  mutate(
		eGFR_diff_yearly = eGFR_diff/t_diff,
		eGFR_diff_10_year = eGFR_diff_yearly*10)

	formula_difference_model <- paste0("eGFR_diff_10_year ~ G + ", paste(paste0("PC", 1:20), collapse = " + "))

	lsout$difference_model <- summary(lm(
		formula = as.formula(formula_difference_model),
		data = tbl_diff))

	## BLUPs&LinReg

	### first stage: estimating an age model RI&RS (150K) without genotype (null model)
	formula_age_model_rirs_nullModel <- paste0("eGFR ~ age_at_exam_c50_10 + sex + (1 + age_at_exam_c50_10 | id)")

	age_model_rirs_nullModel <- lmer(
		formula = formula_age_model_rirs_nullModel,
		data = tbl_restricted,
		control = lmerControl(optimizer = 'bobyqa'))

	blups <- coef(age_model_rirs_nullModel)$id %>%
		mutate(id = rownames(.) %>% as.integer(), age_at_exam_c50_10 = age_at_exam_c50_10 - fixef(age_model_rirs_nullModel)["age_at_exam_c50_10"]) %>%
		rename(random_intercept = `(Intercept)`, random_slope = age_at_exam_c50_10) %>%
		select(-sex)

	tbl_secondStage <- tbl_restricted %>%
	  arrange(id, time_since_baseline) %>%
	  group_by(id) %>%
	  summarise(
		id = first(id),
		# Extract age at first observation
		age_at_baseline = first(age_at_baseline),
		# Extract sex at first observation
		sex = first(sex),
		# Extract PC1 thru PC20 at first observation
		PC1 = first(PC1), PC2 = first(PC2), PC3 = first(PC3), PC4 = first(PC4), PC5 = first(PC5), PC6 = first(PC6), PC7 = first(PC7),
		PC8 = first(PC8), PC9 = first(PC9), PC10 = first(PC10), PC11 = first(PC11), PC12 = first(PC12), PC13 = first(PC13), PC14 = first(PC14),
		PC15 = first(PC15), PC16 = first(PC16), PC17 = first(PC17), PC18 = first(PC18), PC19 = first(PC19), PC20 = first(PC20),
		# Extract genotype at first observation
		G = first(G)
	  ) %>%
	  ungroup() %>%
	  left_join(blups, by="id")

	### second stage: estimating a linear regression model of BLUPs on genotype
	formula_linReg <- paste0("random_slope ~ G +", paste(paste0("PC", 1:20), collapse = " + "))

	lsout$linReg_model <- summary(lm(
		formula = as.formula(formula_linReg),
		data = tbl_secondStage))

	save(lsout, file = paste0(G,".models.Rdata"))

}
