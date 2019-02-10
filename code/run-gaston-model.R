############################################################################
############################################################################
###                                                                      ###
###                  GWAS ANALYSIS USING GASTON                         ###
###                                                                      ###
############################################################################
############################################################################

##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## The input are several variables defined in
## src/setup.R + the variable args passed on from the bash file.
## The output are png's and txt's files in the DIR_SCRATCH folder.
## All intermediate files go into DIR_PROCESSED folder.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


##----------------------------------------------------------------
##                           1. Setup                           --
##----------------------------------------------------------------

## some default args that will be overwritten in the next line, but 
## are here for direct execution of the file when debugging
args <- c("some-path-to-datafile.dat", 
          "lmm", "t1", 4, "forreal")

## this reads the arguments from the bash-make-like file
## Rscript --vanilla run-gaston-model.R  ${PATHOGEN[$index]} $MODEL $TAXID $CORES $DEBUGGING ${THRESH[$index]}
args <- commandArgs(trailingOnly = TRUE)

## assign args to variables
PATHOGEN <- args[1]
model_method <- args[2]
taxid <- args[3]
n_cores <- as.numeric(args[4]) ## number of threads used
debugging <- ifelse(args[5]=="debugging", TRUE, FALSE) ## if "debugging" make datasets small for testing
outcome_thresh <- as.numeric(args[6]) ## threshold for outcome frequency (should be 0.3 for AA and 0.1 for genes)
  
## setting the R code directory
DIR_SRC <- here::here("src")

## This R script creates all the variables
## takes counter_for_pathogen_outcome as input!
source(glue::glue("{DIR_SRC}/setup.R"))

## Here is a list of the variables it creates that we are recycling later
## - HOST: host file string, with path
## - NAM: thats the PATHOGEN file string, but without path
## - PATHOGEN: same as NAM but no path
## - DIR_SRC: where this file is stored
## - DIR_SCRATCH: where all results are stored
## - DIR_OUTPUT: where all processed raw data are stored

## add string to debugging (so files can be easily deleted)
if(debugging) NAM <- glue::glue("{NAM}_debugging")

## read in the plan generated in DIR_SRC/plan.R
plan <- read_delim(glue::glue("{DIR_SCRATCH}/plan_{NAM}.txt"), delim = " ")

##----------------------------------------------------------------
##                         2. Host Data                         --
##----------------------------------------------------------------
## see also https://onlinelibrary.wiley.com/doi/epdf/10.1002/mpr.1608
## Data is in binary plink format

source(glue::glue("{DIR_SRC}/prep-host.R"))

## QC:
HOST_QC <-
  glue::glue("{HOST}_QC") %>% str_replace(DIR_HOST, DIR_PROCESSED)

if (debugging) {
  HOST_QC <- glue::glue("{HOST_QC}_debugging")
  ## apply callrate snps
  system(
    glue::glue(
      "{PLINK} --bfile {HOST} --geno {callrate.snps.thresh} --from-bp {from_to_debugging[1]} --to-bp {from_to_debugging[2]} --chr {chr_debugging} --maf {maf.thresh} --hwe {hwe.thresh} --make-bed --out {HOST_QC}_intermediate"
    )
  )
  ## apply callrate individuals
  system(
    glue::glue(
      "{PLINK} --bfile {HOST_QC}_intermediate --from-bp {from_to_debugging[1]} --to-bp {from_to_debugging[2]} --chr {chr_debugging} --mind {callrate.inds.thresh} --maf {maf.thresh} --hwe {hwe.thresh} --make-bed --out {HOST_QC}"
    )
  )
  
} else {
  ## apply callrate snps
  system(
    glue::glue(
      "{PLINK} --bfile {HOST} --geno {callrate.snps.thresh} --maf {maf.thresh} --hwe {hwe.thresh} --make-bed --out {HOST_QC}_intermediate"
    )
  ) 
  ## apply callrate individuals
  system(
    glue::glue(
      "{PLINK} --bfile {HOST_QC}_intermediate --mind {callrate.inds.thresh} --maf {maf.thresh} --hwe {hwe.thresh} --make-bed --out {HOST_QC}"
    )
  ) 
}


## load post-qc files ---
data_host_raw <-
  load_prepare_host(basename = HOST)

## load post-qc files ---
data_host <-
  load_prepare_host(basename = HOST_QC)

## plotting ---
if (!file.exists(glue::glue("{DIR_SCRATCH}/desc_postqc_host.png")) &
    !file.exists(glue::glue("{DIR_SCRATCH}/desc_preqc_host.png"))) {
  summarise_host(
    data = data_host_raw,
    path_out = glue::glue("{DIR_SCRATCH}/desc_preqc_host.png")
  )
  summarise_host(
    data = data_host,
    path_out = glue::glue("{DIR_SCRATCH}/desc_postqc_host.png")
  )
}

## data_host and data_host_raw are of class bed.matrix


##----------------------------------------------------------------
##                       3. Pathogen Data                       --
##----------------------------------------------------------------

source(glue::glue("{DIR_SRC}/prep-pathogen.R"))

## define coverage filenames ---
files_coverage <- case_when(
  str_detect(PATHOGEN, "t1") ~ as.character(glue::glue("{DIR_PATHOGEN}/NC_009334_with_t1_alts.covstats.dat")),
  str_detect(PATHOGEN, "t2") ~ as.character(glue::glue("{DIR_PATHOGEN}/NC_007605_with_t2_alts.covstats.dat"))
)

## files_coverage provides us with a quality measure of the pathogen files. 

## load raw files ---
data_pathogen_raw <-
  load_prepare_pathogen(files_coverage = files_coverage, files_data = PATHOGEN, debugging = debugging)


## apply QC to pathogen data ---
data_pathogen <- apply_qc_pathogen(data_pathogen_raw)
## data_pathogen is a list!
## data_pathogen$data 
## data_pathogen$outcome


## summarise data as a visualisation ---

## pre qc summarising
summarise_pathogen(
  data = data_pathogen_raw,
  path_out = glue::glue("{DIR_SCRATCH}/desc_pathogen_preqc_{NAM}"),
  filename = NAM,
  format = "png"
)

## post qc summarising
summarise_pathogen(
  data = data_pathogen,
  path_out = glue::glue("{DIR_SCRATCH}/desc_pathogen_postqc_{NAM}"),
  filename = NAM,
  format = "png"
)


## create extra column in plan that indicates whether pathogen will be tested or not
## filter for outcome_thresh first
plan <- plan %>% mutate(include = (outcome.freq >= outcome_thresh & outcome.freq <= (1-outcome_thresh)))

## we have to take the intersect between the QC done before for outcome, and the QC based on the thresholding for frequencies)
outcome <- intersect(data_pathogen$outcome, plan$outcome[plan$include])

## special case for when debugging
if(!debugging) {
  outcome <- intersect(outcome, plan$outcome[plan$include])
}

## load the genetic variants of the pathogen 
## to calculate the GRM matrix
data_pathogen_grm <-
  apply_qc_pathogen(
    load_prepare_pathogen(files_coverage = files_coverage, files_data = PATHOGEN_FORGRM, debugging = debugging)
  )

## data_pathogen_grm is a list!
## data_pathogen_grm$data 
## data_pathogen_grm$outcome


##---------------------------------------------------------------
##                        4. Covariates                        --
##---------------------------------------------------------------

source(glue::glue("{DIR_SRC}/prep-covar.R"))

file_clinical <-
  glue::glue("{DIR_COVAR}/data-covar.covar")
## this file still has all the PCs in from older data. do not use!!!
## our PCs will be generated later.

data_covar <- load_prepare_covar(file_clinical = file_clinical,
                                 file_pca = NULL)
data_covar %<>% mutate(age.log10 = log10(age), rna.log10 = log10(rna), cd4.log10 = log10(cd4))

## summarise covariables in a visualisation
summarise_covar(
  data_covar,
  vars = c("rna", "sex", "age", "cd4", "rna.log10", "age.log10", "cd4.log10"),
  path_out = glue::glue("{DIR_SCRATCH}/desc_covar"),
  format = "png"
)



##---------------------------------------------------------------
##                           5. GRMs                           --
##---------------------------------------------------------------

## Right now, we calculate the GRMs for the pathogen and the host on a variant basis
## TODO: Apply LOCO for host and leave AA out for pathogen

## from: https://stackoverflow.com/questions/32675820/how-to-handle-missing-values-in-crossprod-in-r
tcrossprod.replacena <- function(x, val = 0) {
  tcrossprod(replace(x, is.na(x), val))
}

## this only works because the pathogen is relatively small (12400 columns)
grm_pathogen <-
  tcrossprod.replacena(scale(data_pathogen_grm$data[, data_pathogen_grm$outcome])) /
  nrow(data_pathogen_grm$data) ## same as in gaston vignette
colnames(grm_pathogen) <-
  rownames(grm_pathogen) <- data_pathogen_grm$data$id

## if we wanted the same method as in gaston
#scale(data_pathogen[,outcome]) or
#.Call("gg_Kinship_w", PACKAGE = "gaston", data_pathogen[,outcome])

grm_host <- gaston::GRM(data_host)

## what variants are in the gene that we are looking at?
##variants_notinoutcome <- gene2variant_translation %>% filter(gene != "BNRF1") %>% select(variant)

## add PCs from pathogen
replacena <- function(x, val = 0) {
  (replace(x, is.na(x), val))
}
X <- replacena(data_pathogen_grm$data[, data_pathogen_grm$outcome])
PC <- prcomp(X, scale = TRUE)$x %>% as.data.frame() %>% select(PC1, PC2) %>% mutate(id = data_pathogen_grm$data$id)

## add PCs to data_covar
data_covar <- data_covar %>% select(id, rna, sex, age, cd4, rna.log10, cd4.log10, age.log10) %>% right_join(PC, by = c("id"="id"))

## special case for debugging
if(debugging) {
  str_covar <- c("rna.log10", "sex", "age.log10", "cd4.log10")
} else {
  str_covar <- c("rna.log10", "sex", "age.log10", "cd4.log10", "PC1", "PC2")
}

##----------------------------------------------------------------
##                       6. LM-Covariates                       --
##----------------------------------------------------------------
## Model with only covariates (no SNPs), called "null model"

type_outcome <- plan$model_outcome[1]

tmp <- data_covar %>% right_join(data_pathogen$data[, c("id", data_pathogen$outcome)])

## define GLM family
if (type_outcome == "binary") FAMILY <- "binomial"
if (type_outcome == "quantitative") FAMILY <- "gaussian"

## apply glm to only covariates
fit.null <- purrr::map_df(outcome, function(x) 
  glm(as.formula(paste(x, "~", paste(str_covar, collapse ="+"))), data = tmp, family = FAMILY) %>% broom::tidy() %>% mutate(outcome = x))
write_delim(fit.null, path = glue::glue("{DIR_SCRATCH}/lm_null_{NAM}.txt"), delim = " ")

## update plan
plan[, "str.out.null"] <- glue::glue("{DIR_SCRATCH}/lm_null_{NAM}.txt")
plan[, "covars.null"] <- paste(str_covar, collapse ="+")

## write out results
write_delim(plan, path = glue::glue("{DIR_SCRATCH}/plan_{NAM}.txt"), delim = " ")

## remove unnecessary objects
rm(list = c("data_pathogen_raw", "data_pathogen_grm"))



##----------------------------------------------------------------
##                            7. LMM                            --
##----------------------------------------------------------------

## here we run the actual analysis
## important: make sure that all datasets have synchronised ids
## 
## the function run_analysis has several arguments that define 
## what kind of data are used (data_rhs, data_outcome, data_covar), 
## what kind of model should be used (model_method)
## what kind of GRM should be used (grm_rhs, grm_outcome)
## what kind of strings shoudl be inherited(str_covar, str_outcome, str_id)
## what the results should look like (str_write, return, dir_out)

source(glue::glue("{DIR_SRC}/run-model.R"))

## using multicores:
## - furrr did not work
## - future did not work
## - gaston.utils::association.test.parallel works

## here we are using a for loop and for each index k, the GWAS is parallelised
## through the function gaston.utils::association.test.parallel.

for (k in 1:length(outcome)) {
    run_analysis(
    data_rhs = data_host,
    data_outcome = data_pathogen$data[, c("id", outcome[k])],
    data_covar = data_covar,
    str_covar = str_covar,
    str_outcome = outcome[k],
    str_id = "id",
    dir_out = DIR_SCRATCH,
    dir_processed = DIR_PROCESSED,
    return = FALSE,
    str_write = glue::glue("{NAM}"),
    grm_rhs = grm_host,
    grm_outcome = NULL,
    method = "lmm",
    type_outcome = type_outcome,
    n_pcs = 0,
    dir_plan = glue::glue("{DIR_SCRATCH}/plan_{NAM}.txt"),
    n_cores = n_cores
  )
}


##----------------------------------------------------------------
##                          8. Report                           --
##----------------------------------------------------------------
## analyse_descriptive(data_pathogen, data_host = NULL)
## run this externally in reports/reports.Rmd, till drake is set up
## rmarkdown::render(here::here("reports", "report.Rmd"))
