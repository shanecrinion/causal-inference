# Mendelian Randomisation
Mendelian Randomisation analysis to explore the causal relationship between a potentially modifiable risk factor and an outcome. This script compacts the input and analysis steps in TwoSampleMR. It allows command line input of the instrumental data for 3 approaches - 1) top independent loci, 2) full summary statistics or 3) Using MRBase data. 

<!-- GETTING STARTED -->
## Getting Started
### Packages:
- optparse
- TwoSampleMR

## Preprocessing:
If using approaches 1 or 2, you need to match format column names correctly. The column names must include pos, effect_allele, other_allele, eaf, beta and pval.

Column names can be edited with sed:
sed -i -e '1s/BP/pos/' -e '1s/A1/effect_allele/' -e '1s/A2/other_allele/' -e '1s/FRQ_A_67390/eaf/' -e '1s/OR/beta/' -e '1s/SE/se/' -e '1s/\tP\t/\tpval\t/' PGC3-cp.tsv

That says, on line 1 only, substitute BP with pos, A1 with effect_allele, etc. Note that '\t' is included for the 'P' column to avoid replacing 'SNP' with 'SNPval'. 

<!-- USAGE EXAMPLES -->
## Flags:
- -m/--methodinst\*:  The approach being used for MR analysis. Depending on your desired approach, enter toploci, sumstats or mrbase.  
- -d/--datainst: The file location for exposure data.
- -r/--r2: Cutoff for R2 to detect linkage disequilibrium and clump SNPs. Default value is 0.1 but R2 < 0.01, 0.001 have also been used in the literature.
- -e/--exp\*: Name of the exposure.
- -o/--out\*: Name of the outcome.
- -b/--beta\*: Is SNP effect measured as odds ratio or beta? TRUE=beta.
\*: compulsory flags.


## Usage:
./1_performMR.R -m sumstats -d /home/shane/Desktop/files/data/PGC3-supp/PGC3-cp.tsv -r 0.001 -e "Schizophrenia" -o "Chronotype" -b TRUE


<!-- CONTACT -->
## Contact

<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements

