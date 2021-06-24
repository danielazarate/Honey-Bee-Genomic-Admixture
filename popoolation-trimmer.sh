#!/usr/bin/env bash
# script to trim raw fastq reads for quality and length using popoolation_1.2.2 perl script
# Created: January 19, 2021

# to run: ./script.sh anatoliaca (one argument, giving the bee subspecies) 

set -Eeuo pipefail

BEE_SUBSPECIES=${1}

PROJECTS="/projects/ps-dazarate"
OASIS="/oasis/tscc/scratch/dazarate"

WALLBERG="/projects/ps-dazarate/Wallberg.HoneyBees"
HARPUR="/projects/ps-dazarate/Harpur.HoneyBees"

POPOOL="/projects/ps-dazarate/popoolation_1.2.2/basic-pipeline"

QUALITY=25
MIN_LENGTH=40

DIR_OUT="/oasis/tscc/scratch/dazarate/trimmed-reads"
echo "making directory ${DIR_OUT}"
mkdir -p "$DIR_OUT"

echo "Popoolation_1.2.2 is trimming single-pair Wallberg honeybee ${BEE_SUBSPECIES} for quality ${QUALITY} and length ${MIN_LENGTH} ${crossbones}"

for g in {1..10}
do
perl ${POPOO}/trim-fastq.pl \
--input1 ${WALLBERG}/${BEE_SUBSPECIES}.0$g.fastq \
--fastq-type sanger \
--output ${DIR_OUT}/${BEE_SUBSPECIES}.$g.trim \
--quality-threshold ${QUALITY} --min-length ${MIN_LENGTH}
done

echo "trimming is completed, all trimmed FASTQ sequences outputted to ${DIR_OUT}"
echo "great job champion!"
