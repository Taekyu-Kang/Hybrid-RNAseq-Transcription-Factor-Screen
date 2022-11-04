#!/bin/bash

DIR=$1 #directory with paired end files
OUTDIR=$2
LIBTYPE=$3 #library type, fr-unstranded, fr-firststrand, fr-secondstrand
INDEX=$4 #location of bowtie index, ending with index name


for i in $(ls $DIR*.fq | rev | cut -d"/" -f1 | rev | cut -d"_" -f1 | uniq)  #grabs all unique file identifiers after last '/' and before last '_'
do
	echo ${i} #${OUTDIR}${i}/ ${LIBTYPE} ${INDEX} ${DIR}${i}_1.fq ${DIR}${i}_2.fq
	tophat -o ${OUTDIR}${i}/ -N 0 --library-type ${LIBTYPE} ${INDEX} ${DIR}${i}_1.fq ${DIR}${i}_2.fq &
done
