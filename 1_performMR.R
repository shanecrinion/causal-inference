#!/usr/bin/env Rscript 

# Shane Crinion / shanecrinion@gmail.com / s.crinion1@nuigalway.ie
# 16-3-2021
# Command line 2-sample Mendelian Randomisation, suitable for local or MRBase data. 

# usage: ./1_performMR.R -m mrbase -r 0.001 -e "ADHD" -o "Chronotype"

# preprocessing: 
## column names must be SNP, pos, effect_allele, other_allele, eaf, beta
### sed -a -e '1s/BP/pos/' -e '1s/A1/effect_allele/' -e '1s/A2/other_allele/' -e '1s/FRQ_A_67390/eaf/' -e '1s/OR/beta/' -e '1s/SE/se/' -e '1s/\tP\t/\tpval\t/' PGC3-cp.tsv

##### 1. Analysis and script flag set up ####

# import required libraries
suppressPackageStartupMessages({
  library(optparse) # import via flags
  library(TwoSampleMR) # perform MR
  })

# flag list
option_list <- list(
  make_option(c("-i", "--instruments"), type="character", action="store",
              default = NULL, help="Method for instrumental data - enter 'toploci' 'sumstats' or 'mrbase'."),
  make_option(c("-s", "--snps"), type="character", action="store",
              default = NULL, help="File location for toploci or sumstats."),
  make_option(c("-o", "--outcome"), type="character", action="store",
              default = NULL, help="Method for outome data - enter 'sumstats' or 'mrbase'."),
  make_option(c("-f", "--file"), type="character", action="store",
              default=NULL, help="File location for sumstats of outcome."),
  make_option(c("-r", "--r2"), type="double",
              default = 0.001, help= "R2 threshold for clumping."),
  make_option(c("-e", "--exp"), type="character",action = "store",
              default="Exposure", help="Name of exposure."),
  make_option(c("-d","--disorder"),type ="character", action="store",
              default = "Outcome", help="Name of outcome."),
  make_option(c("-b", "--beta"), type = "character", 
              default = NULL, help="Required for sumstats and toploci only. Enter 'exposure', 'outcome' or 'both'."))

# flag parser 
opt_parser <- OptionParser(option_list = option_list,add_help_option = F)
opt <- parse_args(opt_parser)

# set up output dir
if ("results" %in% list.files()){
  print("results dir exists")
} else {
  print("creating results directory")
  dir.create("results")
}

output.dir <- tolower(paste0("results/exp.", gsub(" ", "", opt$e), ".out.", gsub(" ", "", opt$d)))

if (output.dir %in% list.files()){
  message("writing to output directory ", 
          output.dir)
} else {
  message("creating and writing to output directory: ", 
          output.dir)
  dir.create(output.dir)
}

# set up log file
con <- file(paste0(output.dir,"/out.log"))
sink(con, append=FALSE)
sink(con, append=FALSE, type="message")

# Message for user:
message("....")
message("Instrument method: ", opt$i)
message("Outcome method: ", opt$o)
message("Exposure: " , opt$e)
message("Outcome: ", opt$d)
message("...")

print(opt)

# set up warning messages 
if (is.null(opt$i) | is.null(opt$o) | is.null(opt$e) | is.null(opt$d)){
  print_help(opt_parser)
  stop("At least one necessary flag (-i, -e, -o) is missing", call.=FALSE)}

##### 2. Exposure set-up ####

# load mrbase available outcomes 
if (opt$o  == "mrbase" | opt$i == "mrbase"){
  suppressMessages(ao <- available_outcomes())
}

# input data method selection
  if (opt$i ==  "toploci") { # Method 1. Read file containing top independent loci reported by GWAS
    message("Reading top loci file:")
    exposure_data <- read.csv(opt$s, header = T) 
    exposure_data <- format_data(exposure_data, type = "exposure")
    } else if (opt$i == "sumstats") { # Method 2. Extract top SNPs from sumstats
    message("Reading and extracting instrumental SNP P < 5e-8 from exposure summary data...")
    data <- read.table(opt$s, sep="\t", header=T)
    data <- data[data$pval < 5e-8,]
    message("Formatting SNPS...")
    exposure_data <- format_data(data)
    message("Clumping instrumental SNPs..can take a while...")
    exposure_data <- clump_data(exposure_data, clump_r2 = opt$r)
    } else if (opt$i == "mrbase") { # Method 3. Extract SNPs from MRBase
    message("Importing data using MRBase")
    data <- ao[ao$trait == opt$e,]
    data <- data[data$population=="European",]
    id <- data$id
    exposure_data <- extract_instruments(outcomes = id, r2 = opt$r, force_server = T)
    }


# preview data for user
message("Exposure data imported.... ")
print(head(exposure_data))
message("...")
message("Number of instruments imported: ")
message(nrow(exposure_data))
message("...")

##### 3. Outcome set-up ####
if (opt$o ==  "sumstats") { # Method 1. Read sumstats file 
  message("Reading outcome summary stats and extracting instrumental SNPs:") 
  outcome_data <- read_outcome_data(
    snps = exposure_data$SNP, 
    filename = opt$f, 
    sep = "\t")
  outcome_data$outcome <- opt$d
} else if (opt$o == "mrbase"){ # Method 2. Read mrbase data
  data <- ao[ao$trait == opt$d,]
  message("Preview of data available for disorder:")
  print(data)
  data <- data[data$population=="European",]
  id.out <- data$id[is.na(data$id)==FALSE]
  outcome_data <- extract_outcome_data(exposure_data$SNP, id.out, proxies = T)
}

# preview data for user
message("Outcome data imported...")
print(outcome_data[1:6, 1:6])
message("...")
message("Number of instruments available in outcome data:")
print(nrow(outcome_data))
if (nrow(outcome_data) < 3){
  message("Not enough instruments. If you continue, MR will not be robust. Edit the outcome name eg. can search 'Autism|ASD|Autism Spectrum Disorder|autism' instead of 'Autism'")
} else if (length(id.out) >= 3) {
  message("Enough instruments to perform MR.")
}

#### 4. Clean & Harmonize Datasets ####
message("...")
message("Converting OR -> Beta")
# convert OR to beta if necessary
if (is.null(opt$b)){
  exposure_data$beta <- exposure_data$beta
  message("OR -> beta conversion not required")
} else if (opt$b == "both") {
  exposure_data$beta <- log(exposure_data$beta)
  outcome_data$beta.outcome <- log(outcome_data$beta.outcome)
} else if (opt$b == "exposure") {
  exposure_data$beta <- log(exposure_data$beta)
} else if (opt$b == "outcome") {
  outcome_data$beta.outcome <- log(outcome_data$beta.outcome)
} 


# harmonize data
message("Harmonising data to flip the effect to co-ordinate between exposure and outcome...")

dat <- harmonise_data(exposure_data,
                      outcome_data, 
                      action = 2)

#### 5. Estimate causal effect of exposure on outcome ####
message("Performing MR (IVW, Egger regression, penalised weighted median, weighted median and IVW radial)") 
mr_results <- mr(dat, method_list=c("mr_ivw", 
                                  "mr_egger_regression",
                                  'mr_penalised_weighted_median',
                                    "mr_weighted_median", 
                                  "mr_ivw_radial")) # heterogeneity tests = all but median tests

# add OR
mr_results$OR <- exp(mr_results$b)
# add CI 
mr_results$CI_lower <- 
  exp(mr_results$b-(1.96*mr_results$se))
mr_results$CI_upper <- 
  exp(mr_results$b+(1.96*mr_results$se))

# sensitivity tests
het <- mr_heterogeneity(dat) # Get heterogeneity statistics - maybe unnecessary?
pleiotropy <- mr_pleiotropy_test(dat) # erforms MR Egger and returns intercept values.
res_single <- mr_singlesnp(dat) # Perform mr_wald_ration on each SNP individually. Compare to all mr_ivw and mr_egger_regression 


message("Writing files for MR results and sensitivity tests")
write.csv(x = mr_results, file = paste0(output.dir, "/mrresults.csv"))
write.csv(x = het, file=paste0(output.dir, "/heterogeneity.csv"))
write.csv(x = pleiotropy, file=paste0(output.dir, "/pleiotropy.csv"))
write.csv(x = res_single, file=paste0(output.dir, "/singlesnpMR.csv"))
message("...")



#### Visualize causal effect of exposure on outcome ####
message("Performing visualisations...")
library(ggplot2)
# Generate a scatter plot comparing the different methods
png(paste0(output.dir, "/scatterplot.png")) 
mr_scatter_plot(mr_results, dat)[[1]] + 
  theme_classic() + 
  theme(axis.line.y = element_blank(), 
        axis.line.x = element_blank(), 
        legend.position = "top")
dev.off()

# Generate a funnel plot to check asymmetry
png(paste0(output.dir, "/funnelplot.png"))
mr_funnel_plot(res_single)[[1]] +
theme_classic() + 
  theme( legend.position = "top")
dev.off()

# Run a leave-one-out analysis and generate a plot to test whether any one SNP is driving any pleiotropy or asymmetry in the estimates
res_loo <- mr_leaveoneout(dat)
png(paste0(output.dir, "/loo.png"))
mr_leaveoneout_plot(res_loo)[[1]] +
  theme_classic() + 
  theme( legend.position = "top")
dev.off()


message("Analysis complete...")

# Restore output to console
sink() 
sink(type="message")

## holes - not necessarily issues but won't printing errors
# beta conversion - message for "sumstats or toploci not being used - beta entry being ignored"

