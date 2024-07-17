library(data.table)
library(dplyr)
library(GMMAT)
library(MAGEE)
library(Matrix)
library(doParallel)
library(foreach)

file_pheno <- "path/to/pheno.txt" # pheno.txt as from run_sevenApproaches.r
names_geno <- c("/path/to/geno/for/chr", "/filename/continues.bgen") # must be all bgen files in one directory, e.g. path/to/chr1.bgen
bgen.samplefile <- "/path/to/samplefile.sample"
ncores <- 50
nperbatch <- 50
nclusters <- 4
dir_out <- "path/to/output/directory"

# helpers ----
read_file <- function(a_infile, header="auto",sep="auto",fill=FALSE, ...){
	tblInRaw <- fread(a_infile, header=header, sep=sep)
	tblIn <- as.data.frame(tblInRaw,fill=fill)
	return(tblIn)
}

# pheno ----
id <- "id"
pheno <- read_file(file_pheno)

# setup ----
formula <- paste0("eGFR ~ sex + age_at_exam_c50_10 + ", paste(paste0("PC", 1:20), collapse = " + "))
interaction <- "age_at_exam_c50_10"
random.slope <- "age_at_exam_c50_10"

# null model ----
mod <- glmmkin(
    fixed = formula,
    data = pheno,
    id = id,
    random.slope = random.slope,
    family = gaussian(link = "identity"))

# loop

## Set the number of parallel threads
cl <- makeCluster(nclusters)
registerDoParallel(cl)

## foreach loop
foreach(i = 22:1, .packages = c("data.table", "dplyr", "GMMAT", "MAGEE", "Matrix")) %dopar% {
  # geno
  geno.file <- paste0(names_geno[1], i, names_geno[2])

  # results
  outfile <- file.path(dir_out, paste0("results_chr", i, ".txt"))

  # SNP interactions
  glmm.gei(
    mod,
    interaction = interaction,
    geno.file = geno.file,
    outfile = outfile,
    bgen.samplefile = bgen.samplefile,
    covar.center = "none",
    geno.center = FALSE,
    meta.output = TRUE,
    missing.method = "omit",
    ncores = ncores,
    nperbatch = nperbatch)

  print(paste0("finished chromosome ", i))
}

## Stop the parallel threads
stopCluster(cl)