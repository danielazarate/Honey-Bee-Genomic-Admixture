# WORKFLOW DOCUMENT
## GENOMIC ADMIXTURE IN 60 HONEY BEES SAMPLED FROM 4 DISTINCT LOCALITIES (SAN DIEGO, MEXICO, COSTA RICA, PANAMA)
### LAST MODIFIED NOV 16, 2020
#### AUTHOR: DANIELA ZARATE 

#######################################################################################################################

####### KEY COMMANDS MAGIC BOX ####### 


ssh dazarate@tscc-login.sdsc.edu # login Triton Supercomputer Cluster (TSCC)
qsub -q hotel -l nodes=1:ppn=16 -l walltime=25:00:00 script.sh 
qstat -u dazarate # check on status of the job 
scp -C *fastq.gz dazarate@tscc-login.sdsc.edu:/oasis/tscc/scratch/dazarate/harpur-honeypot/harpur-bams
yqd --user=id --queue
ssh dazarate@porygon-test.ucsd.edu
gstatement -u <user_name> -s 2019-01-01 -e 2019-02-28
#PBS -l nodes=tscc-1-1:ppn=2 # to specify a node 


####### KEY COMMANDS MAGIC BOX ####### 


#######################################################################################################################

## trim for quality and length with popoolation_1.2.2 basic perl trimming script

## maybe I should change the fastq type from sanger to illumina 
# quality 24 length 40 

## example perl trimmer
## note we should change sanger to illumina
perl /projects/ps-dazarate/popoolation_1.2.2/basic-pipeline/trim-fastq.pl \
--input1 /projects/ps-dazarate/Wallberg.Ref.HoneyBees/scutellata.B4.raw.reads.fastq \
--fastq-type sanger \
--output /oasis/tscc/scratch/dazarate/scutellata.B4.trim \
--quality-threshold 25 --min-length 40

# Index reference genome in BWA and map the individual reads to reference 
module load bwa      
bwa index -a is /projects/ps-dazarate/European_genome_scaffolds/Amel_4.5_scaffolds.fa
bwa mem ${REF} ${TRIM_1}.trim > ${TRIM_1}.map 

# Write as a BAM file (binary)
module load samtools 
samtools view -bS ${FILE}.map > ${FILE}.bam

# Remove low quality reads
module load samtools
samtools view -q 20 -bS ${FILE}.bam > ${FILE}.filtered.bam

# Sort the BAM files individually (for pure bees only 1 file)
samtools sort ${FILE}.filtered.bam > ${FILE}.sort.bam

# index 
samtools index ${FILE}.bam 

# create BAM files for all admix pops and source pops (Wallberg et al. 2014) 70 honey bees: 10 A, 10 M, 10 C, 10 O 
# once all the BAM files are created 

#######################################################################################################################


#!/bin/bash
# use ANGSD to generage BEAGLE genotype likelihood file from BAM files of all 130 honey bees 

NGSADMIX_BOX=/media/data/dazarate/ngsadmix_box

angsd.v930 \
-GL 1 \
-out ${NGSADMIX_BOX}/all.bees \
-nThreads 20 \
-doGlf 2 \
-doMajorMinor 1 \
-SNP_pval 1e-6 \
-doMaf 1 \
-minMapQ 40 \
-minQ 20 \
-rmTriallelic 


-bam ${NGSADMIX_BOX}/list.all.bees.txt

# using only invariant sites 
# -doGlf 2 : Beagle haplotype imputation can be performed directly on genotype likelihoods. To generate beagle input file use this 
# -doMajorMinor 1 : infer major and minor allele 
# -GL 1 : estimate genotype likelihood, SAMTOOLS model (4 available: GATK, SOAPsnp, SYK)
# -SNP_pval 1e-6 : We test for polymorphic sites and only output the ones with are likelhood ratio test p-value<1e-6
# -minMapQ : Discard reads with mapping quality below
# -minQ : Discard bases with base quality below
# -trim : Number of based to discard at both ends of the reads
# -remove_bads : Discard 'bad' reads, (flag >=256)
# -setMinDepthInd : Only use site if atleast minInd of samples has this minimum depth
# -setMaxDepth : If total depth is larger then site is removed from analysis
# -setMaxDepthInd : If depth persample is larger then individual is removed from analysis (from site).
# -setMinDepth : If total depth is smaller then site is removed from analysis
# -rmTriallelic 0.000000 : (Remove sites with a pvalue lower) Erin suggests: .05 and .001


## this seems to take 15 hours for all the bees together 

#######################################################################################################################

# use NGSADMIX to estimate ancestry 

#!/bin/bash

ngsadmix \
 -likes /media/data/dazarate/all.bees.beagle.gz \
 -K 4 -P 8 \
 -o /media/data/dazarate/all.bees.K4 \
 -minMaf 0.05 \
 -minInd 122


#-likes beagle file of genotype likelihoods
# -K number of clusters
# -o prefix of output file names
# -P Number of threads used
# -misTol Tolerance for considering a site as missing. Default = 0.05.
# To include high quality genotypes only, increase this value (for example, 0.9)
# -minInd Minumum number of informative individuals. Default = 0
# It only keeps sites where there is at least x # of individuals with NGS data
# -minMaf Minimum minor allele frequency. Default = 5%
# -minLrt Minimum likelihood ratio value for maf>0. Default = 0

## Input file has dim (AFTER filtering): nsites=2,726,805  

