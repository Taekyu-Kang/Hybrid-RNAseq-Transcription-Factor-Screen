#!/bin/bash

DIR=$1 #directory with paired end files
OUTDIR=$2
QUAL=$3 #minimum MAPQ to keep


for i in $(ls $DIR*.fa | rev | cut -d"/" -f1 | rev | cut -d"_" -f1 | uniq)  #grabs all unique file identifiers after last '/' and before last '_'
do
	echo ${i}
	samtools view -q ${QUAL} -O SAM ${OUTDIR}${i}/${i}.sam -o ${OUTDIR}${i}/unique_q${QUAL}.sam &
done
