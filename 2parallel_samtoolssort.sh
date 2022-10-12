#!/bin/bash

DIR=$1 #directory with paired end files
OUTDIR=$2


for i in $(ls $DIR*.fa | rev | cut -d"/" -f1 | rev | cut -d"_" -f1 | uniq)  #grabs all unique file identifiers after last '/' and before last '_'
do
	echo ${i}
	samtools sort -O BAM -o ${OUTDIR}${i}/sorted.bam -T ${OUTDIR}${i}/temp ${OUTDIR}${i}/${i}.sam &
done
