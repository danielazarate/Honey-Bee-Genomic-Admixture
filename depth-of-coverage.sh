#!/bin/bash
# this script uses samtools, parallel, awk, bash script to calculate depth of coverage from honey bee BAM files 

# to run in parallel
# use screen
# beebams.txt = list to paths of all genomes
# size.metrics.txt = list of sizes of all genomes 
# bash script is this one. 

# to correct using total size of the genome (including all uncovered bases): 
# calculate total size of genome for all individuals, using samtools, awk, parallel and outputting to a file so that it appends values to one column. 

parallel --noswap --joblog logs/depth.log --jobs 1 -k ./size-genome-calc.sh {1} :::: list.all.bees.txt

############### mini script: size-genome-calc.sh ############# 
#!/bin/bash

BEE="${1}"
DIR_OUT="/media/data/dazarate/depth-of-coverage"

mkdir -p "$DIR_OUT"

echo "samtools is calculating genome size of ${BEE}" 
samtools view -H ${BEE} | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}' >> ${DIR_OUT}/size.metrics.txt
echo "ending script, all genome size values outputted to ${DIR_OUT}"

############### mini script ############# 

# then, calculate depth of coverage, correcting for total genome size (including both covered and uncovered bases)

# parallel --noswap --link --joblog /logs/.log --jobs 1 -k ./depth-cov.sh {1} {2} :::: beebams.txt :::: size.metrics.txt

############### mini script: depth-cov.sh ############# 
#!/bin/bash

BEE="${1}"
GENOME_SIZE="${2}"

DIR_OUT="/media/data/dazarate/depth-of-coverage"

echo "samtools is calculating depth of coverage for ${BEE}, correcting for genome size" 
samtools depth $BEE |  awk '{sum+=$3} END { print sum/${GENOME_SIZE}}' >> ${DIR_OUT}/depth.txt

echo "ending script"

############### mini script ############# 


# If you want to include regions that were not covered in this calculation you need ot use: samtools depth -a


# to just do a quick back of the hand depth of coverage calculation, divide by NR (positions, as counted by line count) 

parallel --noswap --joblog logs/depth.log --jobs 1 -k ./rough-depth.sh {1} :::: list.all.bees.txt

############### mini script: rough-depth.sh ############# 
#!/bin/bash

BEE="${1}"
ß
DIR_OUT="/media/data/dazarate/depth-of-coverage"

echo "samtools is calculating depth of coverage for ${BEE}, a rough backhand estimate" 
samtools depth $BEE |  awk '{sum+=$3} END { print sum/NR}' >> ${DIR_OUT}/rough-depth.txt

echo "ending script"

############### mini script ############# 

# can then export to local excel and do quick averages/stats 
