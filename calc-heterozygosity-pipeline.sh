## this is the master pipeline with which to calculate heterozygosity per individual 
# (p's and q's) for diversity


# last edit: October 5th 2020 

> log onto porygon ssh dazarate@porygon-test.ucsd.edu
> cd ~/nucleotide.diversity
> cd heterozygosity  # output directory where all data will be sent 


#######################################################################################

# Step 1: Calculate allele frequencies across all populations and reference groups
# per scaffold
# PROGRAM: ANGSD
# SCRIPT:  genoLike_omniSNPs.sh ## adopted from Erin Calfee (UC Davis)
# INPUT:
	(1) BAMs of all populations and reference groups
	(2) list of all scaffolds
	
# OUTPUT
	(1) genotype likelihood file (beagle.gz)
	(2) allele freq file (mafs.gz)
	(3) log files (.arg)
	(4) Output to dir SNPs


# SCRIPT = genotype_likelihood_all_SNPs.sh
# Run on a scaffold by scaffold basis, in parallel, by providing a list of scaffolds 
# TO RUN: 
parallel \
--noswap
 --joblog logs/geno_snps_log \
 --jobs 7 \
 ./script.sh $1 {1} $3 $4 $5 :::: scaffolds.txt

# to run 
# use screen 
parallel --noswap --joblog /logs/.log --jobs 7 ./script.sh $BAM_FILE {1} $MIN_IND $MAX_DEPTH $DIR_OUT :::: scaffolds.txt
 
# This step took 19 hours. 
# log file is in logs dir : genoLike1.log 

#!/bin/bash

# this script uses ANGSD to calculate genotype likelihoods (GL)
# at all putative SNPs for a region of the genome

# command line arguments:
# note: all paths relative to bees/geno_lik_and_SNPs/
BAMS_LIST=$1 # list of paths to bee bams
SCAFFOLD=$2
MIN_IND=$3 # minimum number of individuals with data for total sample to keep a site
MAX_DEPTH=$4 # maximum depth for total sample to keep a site
DIR_OUT=$5 # output file goes in this directory, labelled by  scaffold.

# also requires honeybee reference genome (indexed by samtools faidx)

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
set –o pipefail
set –o errexit
set –o nounset

# make directory to store output (if doesn't yet exist)
# mkdir -p $DIR_OUT

echo "calling variants and GL using ANGSD on BAMS for hilo genomic region: "${SCAFFOLD}

angsd.v930 -out ${DIR_OUT}"/"${SCAFFOLD} \
-r ${SCAFFOLD}: \
-ref "/media/data/dazarate/European_genome_scaffolds/Amel_4.5_scaffolds.fa" \
-bam ${BAMS_LIST} \
-remove_bads 1 \
-minMapQ 30 -minQ 20 \
-doMajorMinor 2 \
-doCounts 1 -minMaf 0.05 -doMaf 8 \
-GL 1 -doGlf 2 \
-P 1 \
-setMaxDepth ${MAX_DEPTH} \
-minInd ${MIN_IND}

echo "all done!"

# settings:
# -r specifies a genomic region to work on, e.g. scaffoldName: specifies the whole scaffold
# -remove_bads removes reads with flags like duplicates
# -doMajorMinor 2: infer major and minor from allele counts
# -bam list of bams to include
# -GL 1: use samtools genotype likelihood method
# -doGlf 2: output beagle likelihood file
# -minMapQ 30 -minQ 20: filter out sites with low mapping quality or base/BAQ quality
# (I pre-computed BAQ scores and replaced quality with minimum of BAQ/base quality,
# so this is equivalend to -baq 2 option here)
# we use AlleleCounts method for MAF (-doMaf 8) with -doCounts
# which doesn't consider base quality and ignores all reads with non-target alleles
# but also doesn't rely on HW equlibrium assumptions to get ML allele freq from genotype freqs
# -minMaf x: and then do a cutoff to only include variant sites with >x minor allele freq.
# -minInd N: only keep sites with information (at least one read) from N individuals (skipped for now)
# -P n means use n threads/nodes for each angsd task (here task=chromosome; then merges threads within-chrom)
# -setMaxDepth -setMinDepth filters out sites where total depth is below or exceeds some threshold

#######################################################################################



# Step 2: Extract human-readable allele frequencies of SNPS from binary file and index

> Use for loop to extract allele freq columns from file and index in angsd 
> create var.sites file for use in allele frequency per population calculations 

# run in one line 
mkdir ../var.sites # make var.sites dir to output results
cp ../scaffolds.txt . # make sure scaffolds.txt is in the working dir 
for i in $(cat scaffolds.txt); do zcat ./${i}.mafs.gz | tail -n +2 | awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' > ../var.sites/${i}.var.sites; sleep 2s; angsd.v930 sites index ../var.sites/${i}.var.sites; done

# check a few files and see if the head is alright 
head <file.var.sites> 
NC_037638.1	4339	G	A
NC_037638.1	4356	A	G
NC_037638.1	4371	G	A
NC_037638.1	4399	T	G
NC_037638.1	4428	G	A
NC_037638.1	4483	G	A 


# takes about 30  min to run across all 

#######################################################################################


# Step 3 : Calculate allele frequencies per population at SNPs across all regions
# SCRIPT: pop_allele_freq.sh (housed in heterozygosity dir)
# INPUT:
	(1) bam list of each population (.pop files)
	(2) variable sites (SNP) list calculated from step 1 across all pops and ref groups
	(3) list of scaffolds
	(4) list of regions (same as scaffolds)

# OUTPUT: 
 (1) .mafs.gz file (labeled as $POP.$SCAFFOLD.AF.mafs.gz)
 (2) arg file 
 (3) into $pop directory 
-------------------------------------- begin script : pop_allele_freq.sh

#!/bin/bash -l

# this script calculates allele frequencies for a population and list of bams

# to run: 
# use parallel & screen
# basic command
./allele_freq.sh $POP $SCAFFOLD 

# to run a small piece to test:
./allele_freq.sh san.diego Group1.1

# parallel command
# the link option will allow the scaffolds to pair the first item from scaffolds.txt to first item of scaffolds.txt, instead of running all the combinations between the two. 
parallel --link --noswap --joblog /logs/.log --jobs 7 ./script.sh $POP {1} {2} :::: scaffolds.txt :::: scaffolds.txt
# example
parallel --link --noswap --joblog logs/san.diego.log --jobs 7 ./pop_allele_freq.sh san.diego {1} {2} :::: scaffolds.txt :::: scaffolds.txt

POP="${1}"
SCAFFOLD="${2}"
REGION="{3}"
#PREFIX="${3}"
#BEE_ID_FILE="../bee_samples_listed/byPop/${POP}.list"
SNP_FILE="./var.sites/${SCAFFOLD}.var.sites"
#REGIONS_FILE="../geno_lik_and_SNPs/results/${PREFIX}/variant_sites/${SCAFFOLD}.regions"
DIR_OUT="./${POP}"
BAM_LIST="./${POP}.pop"

# general bash script settings to make sure if any errors in the pipeline fail
# then it’s a ‘fail’ and it passes all errors to exit and allows no unset variables
# i dont think these work in a script - ask Jude about them when / if they fail again
set –o pipefail
set –o errexit
set –o nounset

REF="/media/data/dazarate/European_genome_scaffolds/Amel_4.5_scaffolds.fa"

echo "Calculating allele frequencies for POP: "$POP" across all ancestries"

# make directory to store output (if doesn't yet exist)

mkdir -p "$DIR_OUT"

echo "bam list:" $BAM_LIST

echo "finding site allele frequencies"
angsd.v930 -out "${DIR_OUT}/${POP}.${SCAFFOLD}.AF" \
-bam "$BAM_LIST" \
-ref "$REF" \
-r "$REGION": \
-underFlowProtect 1 \
-remove_bads 1 \
-minMapQ 30 -minQ 20 \
-doMajorMinor 3 \
-sites "${SNP_FILE}" \
-doCounts 1 \
-doMaf 8 \
-P 1

echo "done getting allele frequencies!"


# options
# basic quality filtering for reads
# -doMajorMinor 3: takes major & minor allele from sites file
# ANGSD calculates freq of minor allele (ignoring all other non maj-min alleles seen)
# -doMaf 8 is an unbiased estimator of pop freq using straight counts and no HWE assumption; -doMaf 8 with -doCounts 1
# -minMapQ 30 -minQ 20: filter out sites with low mapping quality or base/BAQ quality
# (I pre-computed BAQ scores and replaced quality with minimum of BAQ/base quality,
# so this is equivalend to -baq 2 option here)

# underFlowProtect is necessary for large #s of bams

-------------------------------------- end script : pop_alele_freq.sh

#######################################################################################
                                                                                                                    

# Step 4: Use an R script to combine the sites file produced and the mafs.gz file previously produced, merge these two files 
# with also placing NAs for any missing data between the two files 

> Use the R_combine.R script from Erin Calfee 

> Need 
  (1) scaffolds.txt 
  (2) parallel / screen 

Rscript merge_alleleFreq_sites.R $POP $SCAFFOLD 
  # example
  parallel Rscript ./merge_alleleFreq_sites.R san.diego {1} :::: scaffolds.txt 

-------------------------------------- start script : merge_alleleFreq_sites.R

#!/usr/bin/env Rscript
# this script loads population allele frequencies 
# for a pop, and combines them with the sites file,
# producing a new file with a frequency for every site (or NA if no data)
# and number of individuals with data for each site

library(dplyr)

# to run:
Rscript <SCRIPT.R> $POP $SCAFFOLD

# arguments
args = commandArgs(trailingOnly=TRUE)

POP = args[1]
#PREFIX = args[2]
#ANCESTRY = args[3]
SCAFFOLD = args[2]


# paths to input and output directories
sites_file = paste0("./var.sites/", SCAFFOLD, ".var.sites") # points to var.sites files prodcued from all pops/refs
path_allele_freqs = paste0("./POP") # points to where mafs.gz files for the pop are stored 


# read in sites file and set column names 
sites0 <- read.table(sites_file, 
                    stringsAsFactors = F, header = F) %>%
  data.table::setnames(c("scaffold", "pos", "major", "minor"))

# read in the allele freqs table, designate character of column classes, and then merge both sites and mafs files, putting NAs where missing data exists 
allele_freq <- read.table(paste0(path_allele_freqs, "/", POP, ".", SCAFFOLD, ".", "AF", ".mafs.gz"),
                          stringsAsFactors = F, header = T,
                          colClasses = c("character", "numeric", 
                                         "character", "character",
                                         "character", "numeric",
                                         "integer")) %>%
  left_join(sites0, ., by = c("scaffold"="chromo", "pos"="position",
                              "major", "minor"))

# separate out allele freqs, also if MAF > 1 then it will change it to 1 and if MAF < 1 then it will just print MAF 
# MAF is "p" and then this will have column name "pop" and then we will export this 
f <- allele_freq %>%
  mutate(p = ifelse(phat > 1, 1, phat)) %>% # gets rid of weird angsd rounding error
  dplyr::select(p) %>%
  data.table::setnames(POP)
write.table(f, paste0(path_allele_freqs, "/", POP, ".", SCAFFOLD, ".freqs.txt"),
            col.names = T, row.names = F, quote = F)

# separate out number of individuals with data and write to file
# ifelse statement checks to see is nInd has NA and if it does it returns a 0 and if it doesnt it returns nIND value 
n <- allele_freq %>%
  mutate(n = ifelse(is.na(nInd), 0, nInd)) %>% # NA means no individuals with data
  dplyr::select(n) %>%
  data.table::setnames(POP)
write.table(n, paste0(path_allele_freqs, "/", POP, ".", SCAFFOLD,".nInd"),
            col.names = T, row.names = F, quote = F)



-------------------------------------- end script : merge_alleleFreq_sites.R

#######################################################################################
Step 5: Concatenate all files from one population together. 

## let me try to do this first with just one population at a time and then 
# i can try doing this with all the pops 

# for one pop

cat *freqs.txt > pop.all.freqs.txt 


# for multiple pops

so I just used shell/bash to 'paste' the different population files together after 
running the R script to combine them each with the sites file

e.g. 

for j in {1..16}; do paste $(for i in $(cat my_pops.list); do echo 
results/Group$j/$i.freqs.txt; done) > results/pops_included_every_SNP.freqs.txt; done


so the inner loop goes population by population to pastes together every population's 
results for one chromosome ... then the outer loop goes chr by chr ('Group1', 'Group2' etc) and
concats the results together to output just 1 file with all pops and all chromosomes

# download the files to local
scp dazarate@porygon-test.ucsd.edu:~/nucleotide.diversity/heterozygosity/allFreqs/*txt ~/Documents

#######################################################################################

# Step 6 : Calculate 2pq with the small sample correction 

# concatenate nInd files into one file, make sure that it correponds to the Freq file 
cat *.nInd > c.all.nInd.txt 

# download the files to local
scp dazarate@porygon-test.ucsd.edu:~/nucleotide.diversity/heterozygosity/all_nInd/* ~/Documents

# load the nInd files into R using 
calculating.pi.R script 
small_sample_correction_pi.R 

# small sample correction pi function 
het_small_sample_correction <- function(p, n){ # here n is the number of haplotypes observed
  2*(p - p^2*(n/(n-1)) + p*(1/(n-1)))}



#######################################################################################

# STEP 6.5 : Calculate qs and 2pqs from ps (if not using small sample size correction)

t <- f %>%
  mutate(q = 1 - p) %>% 
  mutate(x2pq = 2 * p * q)
# Step 6: Take Average of 2pq 

pi = average across sites of 2*p*(1-p)

#######################################################################################
# Step 7: Weight 2pq by multiplying by 
# use calculating_pi.R script in R for this 
(total # snps in the genome)/(total # of positions in the genome)

# Finding total number of positions in the genome 

from Amel_4.5_scaffolds.fa.fai index file of reference genome
load into R 

# Finding total # snps used 
sum of $v2 column in Amel...fai file 

#######################################################################################

# STEP 8: Jackkniving the pi estimates via chromosome 

# this step needs to start in the ~/nucleotide/heterozygosity/$POP dir on porygon 
# I need to remove one chromosome at a time and cat everything 
 
 mkdir freq 
 mkdir nInd

 # copy all freqs to /freq dir 

# use jude.jackknife.sh script to create jackknife files for both freq and nInd in both jacks directories 
# jude.jackknife.sh should be located in $POP root directory 
# to run 

./jude.jackknife.sh $POP $TYPE

# where type is either freq or nInd

################ use jude.jackknife.sh 

#!/bin/bash
# to run ./script {1} {2} 
# where {1} is  pop e.g. san.diego
# where {2} is type e.g. nInd or freq


POP="${1}"
TYPE="${2}"


cd ${TYPE}
mkdir -p ./jacks

shopt -s extglob

# just start with a string that doesn't match anything
globString="jacks"

for g in {1..15}
do
    globString="${globString}|*Group${g}.*"
    filename="./jacks/"${POP}".jackknife${g}."${TYPE}".txt"
    echo "working on ${filename}"
    cat !(${globString}) > ${filename}
done

################ end script 


# Step 9: Create all the Jackknives for the nInd files as well using the same script to create the jacks for the freqs files 
# within each POP subdirectory, there is a freqs and a nInd directory, within each is housed jacks directories that contain the final 15 jack files 

# also add the full complete original file (used to calculate pi initially with all 16 chromos) to the jacks directories for both freq & nInd and call it the 16th jack 
# all the files should be in the heterozygosity/allFreqs and /allnInd dirs 
 cp ./mexico.all.freqs.txt ../mexico/freq/jacks/mexico.jackknife16.freq.txt

# then use the "jackknife.pi" non-interactive script to read in the nInd and freq jackknife files. Use parallel -link option to read them in in appropriate order

# make 
mkdir corrected_pi 
mkdir logs 

# run it from the $POP root subdir 
parallel --link --joblog ./logs/a2.jackknife.log Rscript ./pi_jackknife.R a {1} {2} :::: chromo.txt :::: jack_positions.txt 


################################# BEGIN pi_jackknife.R script 


# Rscript ./script $POP $CHROMOSOME_NUMBER
# Rscript ./script san.diego 1 

# to run in parallel 
# link option requiered to make sure the items from the lists are matched line by line 
# parallel --link --noswap --joblog /logs/.log --jobs 7 ./script $POP $CHROMSOME_NUMBER $JACK_POS :::: chromo.txt :::: jack_positions.txt 

# the paths specified in this script assume that the script is being launched from the POP subfolder
# path: /media/data/dazarate/nucleotide.diversity/heterozygosity/san.diego
# or ~/nucleotide.diversity/heterozygosity/san.diego

# arguments
args = commandArgs(trailingOnly=TRUE)

POP = args[1]
CHROMO = args[2]
JACK_POS = args[3]

# read the nInd and Freq file in, these variables should stand in for the full path of these files 
POP.freq =  paste0("./freq/jacks/", POP, ".jackknife", CHROMO, ".freq.txt")
POP.nInd =  paste0("./nInd/jacks/", POP, ".jackknife", CHROMO, ".nInd.txt")

# read in the table given the path provided above 
POP.freq.table <- read.table(POP.freq,  stringsAsFactors = F, header = T)
POP.nInd.table <- read.table(POP.nInd,  stringsAsFactors = F, header = T)

# remove all the "POP" headers that creeped into it during the concatenate
# write a function to perform this step

clean_names <- function(x, ...) {
  x[!x[1] == POP,]
}

# convert the list to a data frame and then remove headers

# for nInd
POP.nInd <- as.data.frame(POP.nInd.table)
POP.nInd2 <- as.data.frame(clean_names(POP.nInd))
# for freq
POP.freq <- as.data.frame(POP.freq.table)
POP.freq2 <- as.data.frame(clean_names(POP.freq))


# merge the nInd and freq data frame per population 
# first rename the columns to match 
X = setNames(POP.nInd2,c('head'))
Y = setNames(POP.freq2,c('head'))
# then bind
POP_frame <- cbind(X, Y)
# create a third column with 2*nInd
# change columns back to actual names 
POP_frame = setNames(POP_frame,c('nInd', 'freq'))
# convert class from character to numeric 
POP_frame[] <- lapply(POP_frame, function(POP_frame) as.numeric(POP_frame))
# multiply nInd by 2 to get haplotype number 
POP_frame$haplotype <- POP_frame$nInd * 2

# use the small sample size correction function to output a new column with the corrected 2pq (pi) value 

# small sample size correction function 
het_small_sample_correction <- function(df, col_name, p, n){ # here n is the number of haplotypes observed
  df[[col_name]] <- 2*(df[[p]] - df[[p]]^2*(df[[n]]/(df[[n]]-1)) + df[[p]]*(1/(df[[n]]-1)))
  df
}

# apply small sample correction to POP frame and produce a correction column 
POP_corr <- het_small_sample_correction(POP_frame,"correction","freq","haplotype")

# take average of the corrective measure 
POP_corr_mean <- mean(POP_corr$correction,na.rm=TRUE)

# Weight by total SNPS / total Positions
# to find total positions 
# establish the total_position number 
total_pos <- as.numeric(JACK_POS)

# count the number of rows in the data frame 
count <- nrow(POP_corr)

# produce the weight 
weight <- count / total_pos

# multiply mean by weight to achieve corrected pi 
POP_corr_pi <- POP_corr_mean * weight

# write pop_corr_pi to a file
# again, the POP within the string is going to be read like POP unless I do something to escape it 
# let's see if this works : 
write.table(POP_corr_pi, paste0("./corrected_pi", POP, CHROMO, ".jackknife.corrected.pi"),
            col.names = F, row.names = F, quote = F)


############################ end pi_jackknife.R script 

cat all the corrected pi values into a list 
cat * > m.cat.jacks.txt

# write an Rscript that will calculate standard error from the cat.jacks.txt file 
# run from $POP directory 
# will output two estimates of standar error using two different techniques

# Rscript ./standard.error.R $POP

############################### start standard.error.R script 
#!/usr/bin/env Rscript 
# this script takes a list of numbers and calculates standard error from it 
# run from root $POP directory 

args = commandArgs(trailingOnly=TRUE)

POP = args[1]


#read in table
POP_pi_jacks = paste0("./corrected_pi/", POP, ".cat.jacks.txt")

# read in the table given the path provided above 
POP_pi_jacks_table <- read.table(POP_pi_jacks,  stringsAsFactors = F, header = F)


# calculate standard error via package 
install.packages("plotrix")
library(plotrix)
POP_standard_error1 <- std.error(POP_pi_jacks_table,na.rm)

write.table(POP_standard_error1, paste0("./corrected_pi/", POP, ".standard.error1.txt"),
            col.names = F, row.names = F, quote = F)


# calculate standard error via function 
std <- function(x) {sd(x)/sqrt(length(x))}

# convert table column to vector 
vec1 <- POP_pi_jacks_table$V1

# run function 
POP_standard_error2 <- std(vec1)

write.table(POP_standard_error2, paste0("./corrected_pi/", POP, ".standard.error2.txt"),
            col.names = F, row.names = F, quote = F)

############################### end standard.error.R script 

# should output two values to corrected_pi directory under standard error labels 
