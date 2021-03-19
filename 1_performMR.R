#!/usr/bin/env Rscript 

# Shane Crinion / shanecrinion@gmail.com
# 16-3-2021
# Compact and versatile 2SMR suitable for local or MRBase data. 

# usage: ./1_performMR.R -m toploci d toploci_chronotype.csv

# preprocessing: 
## column names must be SNP, pos, effect_allele, other_allele, eaf, beta
### sed -i -e '1s/BP/pos/' -e '1s/A1/effect_allele/' -e '1s/A2/other_allele/' -e '1s/FRQ_A_67390/eaf/' -e '1s/OR/beta/' -e '1s/SE/se/' -e '1s/\tP\t/\tpval\t/' PGC3-cp.tsv

# import required libraries
suppressPackageStartupMessages({
  library(optparse) # import via flags
  library(TwoSampleMR) # perform MR
  })

# create flag list
option_list <- list(
  make_option(c("-m", "--methodinst"), type="character", action="store",
              default = NULL, help="Enter 'toploci' 'sumstats' or 'mrbase'"),
  make_option(c("-d", "--datainst"), type="character", action="store",
              default = NULL, help="File location for toploci or sumstats"),
  make_option(c("-z", "--dataout"), type="character", action="store",
              default=NULL, help="File location for sumstats of outcome"),
  make_option(c("-r", "--r2"), type="double",
              default = 0.001, help= "R2 threshold for clumping"),
  make_option(c("-e", "--exp"), type="character",action = "store",
              default="Exposure", help="Name of exposure"),
  make_option(c("-o","--out"),type ="character", action="store",
              default = "Outcome", help="Name of outcome"),
  make_option(c("-b", "--beta"), type = "logical", 
              default = TRUE, help="Is SNP effect measured as beta? Enter TRUE for beta or FALSE for odds ratio"))

# flag parser 
opt_parser <- OptionParser(option_list = option_list,add_help_option = F)
opt <- parse_args(opt_parser)

# Message for user:
message("....")
message("Input Method: ", opt$m)
message("Exposure: " , opt$exp)
message("Outcome: ", opt$out)
message("...")


# set up warning messages 
if (is.null(opt$m) | is.null(opt$e) | is.null(opt$o)){
  print_help(opt_parser)
  stop("At least one necessary flag (-m, -e, -o) is missing", call.=FALSE)}

##### 1. Import instruments ########
  if (opt$m ==  "toploci") { # Method 1. Read file containing top independent loci reported by GWAS
    message("Reading top exposure loci file:")
    exposure_data <- read.csv(opt$d) 
    } else if (opt$m == "sumstats") { # Method 2. Extract top SNPs from sumstats
    message("Reading and extracting SNP P < 5e-8 from exposure summary data...")
    data <- read.table(opt$d, sep="\t", header=T)
    data <- data[data$pval < 5e-8,]
    message("Formatting SNPS...")
    exposure_data <- format_data(data)
    message("Clumping SNPs..can take a while...")
    exposure_data <- clump_data(exposure_data, clump_r2 = opt$r)
    } else if (opt$m == "mrbase") { # Method 2. Extract SNPs from MRBase
    message("Importing data using MRBase")
    ao <- available_outcomes()
    data <- ao[ao$trait == opt$e,]
    PMID <- readline(prompt="Enter PMID (or press Enter is NA): ")
    Population <- readline(prompt = "Enter Population eg. Mixed, European or Asian")
    if (PMID == "") {
      data <- data[data$population==Population,]
    } else {
       data <- data[data$population==Population & data$PMID==PMID,] }
    id <- data$id
    exposure_data <- extract_instruments(outcomes = id, r2 = opt$r)
    }


message("Exposure data imported...")
print(head(exposure_data))


