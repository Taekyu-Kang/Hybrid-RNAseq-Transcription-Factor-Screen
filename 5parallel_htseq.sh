#!/bin/bash

DIR=$1 #directory with paired end files
OUTDIR=$2
GTF=$3 # GTF file for annoatation
STRAND=$4 #strandedness of library: no, yes, reverse (fr-unstranded,fr-secondstrand,fr-firststrand on tophat)
SORT=$5 #sorting of reads: pos, name

for i in $(ls $DIR*.fa | rev | cut -d"/" -f1 | rev | cut -d"_" -f1 | uniq)  #grabs all unique file identifiers after last '/' and before last '_'
do
	echo ${i}
	python /bigrock_home/tkang/HTSeq/scripts/count.py -s ${STRAND} -r ${SORT} ${OUTDIR}${i}/unique.sam ${GTF} > {OUTDIR}${i}/${i}_htseq.out &
done
