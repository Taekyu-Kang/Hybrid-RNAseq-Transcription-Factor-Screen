#!/bin/bash

DIR=$1 #directory with paired end files
OUTDIR=$2
LIBTYPE=$3 #library type, fr-unstranded, fr-firststrand, fr-secondstrand
INDEX=$4 #location of bowtie index, ending with index name


for i in $(ls $DIR*.fq.gz | rev | cut -d"/" -f1 | rev | cut -d"_" -f1 | uniq)  #grabs all unique file identifiers after last '/' and before last '_'
do
	echo ${i}
	tophat -o ${OUTDIR}${i}/ -N 0 --library-type ${LIBTYPE} ${INDEX} ${DIR}${i}_1.fq.gz ${DiR}${i}_2.fq.gz &
done
