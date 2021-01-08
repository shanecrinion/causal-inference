library(tidyverse)
library(dplyr)

# exposure = sz
# outcome = chronotype

## Input outcome data
# read in scz data
data_scz <- read.delim('~/Desktop/files/data/PGC3_SCZ_wave3_public.v2.tsv.gz')
# rename the headings


# find overlapping snps
#overlap <- intersect(data_scz$SNP, data_chronotype$SNP)
#extract the subsets
#data_scz.subset <- data_scz[data_scz$SNP %in% overlap,]
# data_chronotype.subset <- data_chronotype[data_chronotype$SNP %in% overlap,]

# Test effect of Chronotype (X) on Schizophrenia
# limit to required columns for simplicity



## Data Harmonisation
# a. How many SNPs do authors describe as being independently associated with Schizophrenia?
# get snps associated with schizophrenia from newest paper
# make sure they're not associated with chronotype
# extract the snps with highest p-value
length(which(data_scz$P <= 5E-8)) #22345
head(data_scz[which(data_scz$P <= 5E-8),])
# later - get the list of snps from PGC3 paper

# b.  What information does this file contain that are needed for Mendelian randomization analyses?
colnames(data_scz)

# c. How many are associated within only Europeans?
# 80% European and 20% East Asian

# d. Why might it be best to use the SNPs that have been identified as being associated with BMI in Europeans only?
# SNPs are inherited in the same LD blocks and does not cause false positives

### 2. Read in the second sheet of this file to get the estimates of the SNPs associated with BMI in Europeans.
# we don't have info for the ethnicity

# a. Check that these are all associated with BMI at a conventional level of genome-wide significance
# i will limit to only SNPs that are GWS 
data_scz <- data_scz[which(data_scz$P <= 5E-8),]
dim(data_scz)

# b. Are all of these SNPs "good instruments"? What else might we want to check to see if they are strongly and independently associated with BMI? e.g., what happens if you look at the SNPs per chromosome?
# this explains why there are so many GWS SNPs (i think)
# we need to check whether the snps are in high LD - if they are not then they are good insturments
# we can do this using LDlink or SNPsnap
# As there are so many, I will get the list of associated SNPs from the paper
# import the list of top SNPs from Supp table 3 (after converting to CSV)
top_snps <- read.csv("~/Desktop/files/data/PGC3-supp/Supplementary Table 3 - Combined discovery-replication loci.csv")
# remove commas in values and make decimal
snp_vals <- c("top.freq", "top.info","top.P","top.OR", "top.SE")
for(i in snp_vals){
 top_snps[,i] <-   as.numeric(gsub(",", ".", gsub("\\.", "", top_snps[,i])))
}

### 3.	We're going to make sure the effect allele is the allele that increases schizophrenia risk using the effect allele and beta column 
# Browse the data, are all SNPs coded so that the effect allele increases BMI? 
# a.	Are all SNP effects in the same direction? 
top_snps$top.OR <- log(top_snps$top.OR) # got log odds ratio as required later on

#top_snps[,c("top.index","top.alleles..A1.A2.","top.OR","top.SE","top.P"),] 


# b. What would we need to do in order to make sure that all effect estimates were in the same direction?
library(tidyr)
top_snps <- top_snps %>% 
  separate(top.alleles..A1.A2., c("A1", "A2"), sep = "/")
head(top_snps)

# how many snps are not in positive direction?
length(which(top_snps$top.OR < 0))

top_snps_test <- top_snps
names(top_snps)

swapdir <- function(data, OR, A1, A2, A1_freq){
  for (i in 1:nrow(data)){
    # first check if beta direction matches
    if (data[i,OR] < 0){
      # change direction 
      data[i,OR] <- -1 * data[i,OR]
      # extract allele info for swapping
      new_ref <- data[i,A2] ; new_alt <- data[i,A1]
      new_freq <- 1 - data[i,A1_freq] 
      # swap values 
      data[i,A1] <- new_ref; data[i,A2] <- new_alt
      data[i,A1_freq] <- new_freq
    }}}

swapdir(top_snps_test, "top.OR", "A1", "A0", "top.freq")

top_snps[1:10,1:10]
top_snps_test[1:10,1:10]

for (i in 1:nrow(top_snps_test)){
  # first check if beta direction matches
  if (top_snps_test[i,"top.OR"] < 0){
    # change direction 
    top_snps_test[i,"top.OR"] <- -1 * top_snps_test[i,"top.OR"]
    # extract allele info for swapping
    new_ref_exp <- top_snps_test[i,"A0"] ; new_alt_exp <- top_snps_test[i,"A1"]
    new_freq_exp <- 1 - top_snps_test[i,"top.freq"] }}


###################
# supp info says SE = "Standard error of the meta-analyzed coefficients(i.e. SE of log(OR))from the logistic regressionreported by METAL".
 ## supp info says OR = "Odds ratiofrom the meta-analysisusing A1 allele as the baseline. This is the exp(beta), where beta is the effect estimate as reported by METAL". 
  ## i might need to convert OR to log to match with chronotype data but I don't know how
# ** not sure if it's right to use log here **


################################################################# PART 2 #################################################################

############################################################					
#SNP LOOKUP IN A GWAS OF CHRONOTYPE
############################################################


# a. How many SNPs are in this truncated CARDIOGRAM dataset?
# Input exposure data
data_chronotype <- read.delim('~/Desktop/files/data/chronotype_raw_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt.gz')

# how many are HWE
dim(data_chronotype[which(data_chronotype$HWE_P < 5e-8),])
dim(data_chronotype)

# b. Does this file contain everything that is required to perform a two-sample Mendelian randomization analysis?
colnames(data_chronotype) # yes: SNP name, ref & other allele, p-val, Log odds ratio and standard error

# convert beta to log odds ratio (OR = exp(beta))
data_chronotype$BETA <- log(exp(data_chronotype$BETA))
# convert standard error (of beta) to log standard error of OR 
data_chronotype$SE <- log(exp(data_chronotype$SE))


### 2. How many of the SCZ SNPs are included in the chronotype dataset?
SCZ_SNPs <- top_snps$top.index
SCZ_SNPs <- as.vector(SCZ_SNPs)
#matches <- unique(grep(paste(SCZ_SNPs, collapse="|"), data_chronotype$SNP, value=TRUE))
matches <- Reduce(intersect, list(SCZ_SNPs, data_chronotype$SNP)) # 262

# View the data from Chronotype for our SCZ SNPs
chronotype_SCZ <- data_chronotype[which(data_chronotype$SNP %in% SCZ_SNPs),]
dim(chronotype_SCZ)
head(chronotype_SCZ)

# 3. Merge the GIANT and CARDIOGRAM SNP summary associations
# First, make sure the column headings are easy to understand (i.e., add "BMI" and "CHD" onto the respective datasets)
colnames(top_snps) <- paste("SCZ", colnames(top_snps), sep = "_")
head(top_snps)
colnames(data_chronotype) <- paste("CRTPE", colnames(data_chronotype), sep = "_")
head(data_chronotype)

merged <- merge(top_snps,data_chronotype, by.x="SCZ_top.index", by.y="CRTPE_SNP")
dim(merged) 
head(merged)

# remove unnecessary columns
merged <- within(merged, 
                 rm("SCZ_top.info","SCZ_merge.LEFT", 
                    "SCZ_merge.RIGHT","SCZ_ENSEMBL.top5.genes",
                 "SCZ_ENSEMBL.genes..all..ENSEMBL.IDs.",
                 "SCZ_ENSEMBL.genes..all..clear.names."))



################################################################# PART 3 #################################################################

#####################################################################	
#HARMONIZING THE EFFECT ALLELES IN THE GIANT AND CARDIOGRAM DATASETS#
#####################################################################

### 1. Make sure that the effect alleles in the Chronotype and SCZ datasets are the same. We want the cardiogram effect allele to be the allele that increases SCZ risk
# But be careful of palindromic SNPs or SNPs on different strands
# First we need to see whether the effect alleles are the same
# Browse the data
merged[,c("SCZ_top.index", "SCZ_A1", "CRTPE_ALLELE1","SCZ_A2","CRTPE_ALLELE0", "SCZ_top.freq","CRTPE_A1FREQ")]

# a.	How can we tell if the Chronotype and SCZ SNPs are coded using the same reference strand?
# they have the same reference allele 

# check that all the alleles are the same - they are
x <- 0
for (i in 1:nrow(merged)){
  if(merged[i,"SCZ_A1"] == merged[i,"CRTPE_ALLELE1"]){
    x <- x + 1
  }
print(x)}


# b.	Are Chronotype and the SCZ SNPs coded using the same reference strand? 
# yes

# c.	Are there any palindromic SNPs? 
palindromic_at<-subset(merged,SCZ_A1 %in% "A" & SCZ_A2 %in% "T")
palindromic_ta<-subset(merged,SCZ_A1 %in% "T" & SCZ_A2 %in% "A")
palindromic_gc<-subset(merged,SCZ_A1 %in% "G" & SCZ_A2 %in% "C")
palindromic_cg<-subset(merged,SCZ_A1 %in% "C" & SCZ_A2 %in% "G")
dim(palindromic_at)
dim(palindromic_ta)
dim(palindromic_gc)
dim(palindromic_cg)

# d.	How can we tell whether the effect alleles are the same in both datasets for palindromic SNPs
## (i.e., the allele that increases SCZ is the same as the reference allele in Chronotype)?

### 2.	Make sure the Chronotype log odds ratio reflects the allele that increases SCZ in the PGC data
# First, find the positions of SNPs with different effect alleles
# the log odds are not correct so I will flip the reference allele if beta is negative

merged_test <- merged



# ---------------------------
## Limitations: PGC data is European + East Asian ancestry. 


