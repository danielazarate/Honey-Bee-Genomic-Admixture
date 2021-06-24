# The Master Calculating Pi per Population Script 

### Last modified: 14 September 2020
# to run 

#######################################################################################
# Step 1: Calculate SAF Global Estimate using ANGSD per scaffold basis 
# -doSaf 1 : Calculate the Site allele frequency likelihood based on individual genotype likelihoods assuming HWE
#  for now, we can use the REF as ANC...ancestral state needs to be supplied for the full SFS
## but you can use the -fold 1 to estimate the folded SFS and then use the reference as ancestral.
# use the entire population for this

# To Run: script needs to be in the same directory as pop and scaffold files 
parallel --noswap --joblog logs/m_pop_saf_globalEstimate.log --jobs 4 ./saf.global.estimate.sh $bam_list {1} $minINd :::: scaffolds.txt 
parallel --noswap --joblog logs/o_pop_saf_globalEstimate.log --jobs 4 ./scripts/saf.global.estimate.sh san.diego.pop {1} 8 :::: scaffolds.txt
parallel --noswap --joblog logs/san_pop_saf_globalEstimate.log --jobs 4 ./saf.global.estimate.sh san.diego.pop {1} 8 :::: scaffolds.txt

# Needed:
INPUT=$BAM_FILE # list of paths of bam files
SCAFFOLD # list of all scaffolds

> use parallel to run a list of scaffolds per population 
> use saf.global.estimate.sh Script

> Actual Script Below:

#!/bin/bash
# this script calculates global SFS for a population of individuals
# for eventual nucleotide diversity calculation


# define your variables
POP=$1 # first argument is a list of BAMS
SCAFFOLD=$2
minInd=$3 # change according to number of individuals
# reference genome
REF=/media/data/dazarate/European_genome_scaffolds/Amel_4.5_scaffolds.fa



echo "calculating site allele frequency likelihood global estimate based on individual genotype likelihoods"

OUTDIR=/media/data/dazarate/nucleotide.diversity/pi/results
OUT=$POP.$SCAFFOLD.out

echo "out is $OUT"
echo "pop is $POP"
echo "minInd is $minInd"
echo "scaffold is ${SCAFFOLD}"

angsd.v928 \
-bam $POP \
-r ${SCAFFOLD}: \
-doSaf 1 \
-anc $REF \
-GL 2 \
-P 10 \
-minQ 20 \
-out $OUTDIR/$OUT \
-minMapQ 30 \
-minInd $minInd

OUTPUTS= .out.arg, .out.saf.pos.gz, .out.saf.idx

########################################################################################

# Step 2 : Merge SAF files for the entire population for all scaffolds
# fast, doesn't need screen 

### changed angsd version here to see if this fixes merge error 
singularity exec /media/data/dazarate/singularity.images/angsd.v928.sif /opt/angsd/angsd/misc/realSFS cat -outnames san.pop.merged san.diego*saf.idx
# the outputs will automatically have SAF extensions 

./realSFS cat a.pop*saf.gz
	-> This will cat together .saf files from angsd, both a .saf.gz and a .saf.idx file 

> cp m.pop.merged.saf* ~/nucleotide.diversity/pi/mergedSAF  



########################################################################################

# Step 3 : Calculate SFS from SAF estimate, fold if not providing ancestral state

> cd ~/nucleotide.diversity/pi/mergedSAF 
> use saf2sfs.sh script 
> input is .saf.idx file
> to run: ./saf2sfs.sh $MERGED_SAF

#!/bin/bash
# script to fold the SAF to create the SFS using realSFS from angsd


SAF_INPUT=$1 # first argument is the SAF

OUTDIR=/media/data/dazarate/nucleotide.diversity/pi/SFS
OUT=$OUTDIR/$SAF_INPUT

echo "SAF_INPUT is $SAF_INPUT"
echo "SHORT_SAF is $SHORT_SAF"
echo "out is $OUT"

echo "Obtaining maximum likelihood estimate of the SFS using realSFS"
singularity exec /media/data/dazarate/singularity.images/angsd.v928.sif /opt/angsd/angsd/misc/realSFS $SAF_INPUT -fold 1 -P 10 > $OUT.sfs


########################################################################################
# Step 4: Calculate thetas from SFS using ANGSD. 

$ angsd -bam $bam -out $out -doThetas 1 -doSaf 1 -pest $out.sfs -anc $ANC -GL 2

# Need
> list of BAMS
> SFS file
> Reference Genome
> script = theta.caller.sh

> to run :

> cd ~/nucleotide.diversity/pi/thetas
> ./theta.caller.sh $BAM $SFS
> ./theta.caller.sh ../o.pop ../SFS/o.pop.merged.saf.idx.sfs 

# actual script:

#!/bin/bash


BAM=$1
SFS=$2
#SCAFFOLD=$3
OUT_DIR=/media/data/dazarate/nucleotide.diversity//pi/thetas
OUT=$OUT_DIR"/"${BAM}.out

REF=/media/data/dazarate/European_genome_scaffolds/Amel_4.5_scaffolds.fa


angsd.v928 \
-bam $BAM \
-out $OUT \
-doThetas 1 \
-doSaf 1 \
-pest $SFS \
-anc $REF \
-GL 2

> outputs .thetas.gz and a .thetas.idx file 
> san diego pop took 4/5 hours 

########################################################################################

# Step 5:  Use thetaStat to print summary statistics 
./misc/thetaStat do_stat out.thetas.idx

> thetaStat do_stat <thetas.idx> -outnames <filename>
> can decide to use window sizes 
> thetaStat do_stat <thetas.idx> x -win 50000 -step 10000  -outnames <filename>

> singularity exec /media/data/dazarate/singularity.images/angsd.v928.sif /opt/angsd/angsd/misc/thetaStat do_stat san.diego.thetas.idx -outnames san.diego
> produces pestPG file


########################################################################################
# Step 6: Calculate genome-wide estimate of pi per population 

> download pestPG files 
> scp dazarate@porygon-test.ucsd.edu:~/nucleotide.diversity/pi/pestPG/san.diego.pestPG ~/Documents/SCRIPTS/pi_calculations/
> use R to calculate average of tP (and tW if you so choose) and average of nSites 
> use pi.calculation.R script
> pi = (average of tP / average of nSites) # is this right? I think I have to divide by total sites in the genome ? 



######################################################################################

# Step 7: I Believe I have to account for small sample size here ??? 

########################################################################################
# Step 7: Jackknife the pi estimates 

> use jackknife.R script in R 
