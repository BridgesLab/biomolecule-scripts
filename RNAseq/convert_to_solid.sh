#!/bin/bash
#This script converts all of the fastqcsanger files in the directory into qual and csfasta files for bowtie

for file in `ls *.fastqcssanger` ; do
    echo "python csfastq2solid.py $file" > ${file%%.*}.sh
    qsub -cwd ${file%%.*}.sh
done
