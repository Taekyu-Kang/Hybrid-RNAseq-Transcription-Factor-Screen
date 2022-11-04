#!/bin/bash

DIR=$1 #directory with paired end files
OUTDIR=$2


for i in $(ls $DIR*.fq | rev | cut -d"/" -f1 | rev | cut -d"_" -f1 | uniq)  #grabs all unique file identifiers after last '/' and before last '_'
do
	echo ${i}
	samtools index ${OUTDIR}${i}/sorted.bam ${OUTDIR}${i}/sorted.bam.bai &
done
