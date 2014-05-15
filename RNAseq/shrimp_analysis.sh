#!/bin/bash

#This script cycles through all fastqcsanger files in a directory first generating a sam file then converting that file into a bam file
#The name of the file is taken in and used for the output files
REFERENCE="reference-annotations/Homo_sapiens.hg19.fa"
echo "Using reference genome file in $REFERENCE"

rm -r shrimp_out
mkdir shrimp_out
stat shrimp_out

for file in `ls *.fastqcssanger` ; do
    #this runs the colorspaced gmapper then pipes the result through samtools to generate a bam file
    echo "gmapper-cs -Q $file $REFERENCE | samtools view -bS > shrimp_out/${file%%.*}.bam" > ${file%%.*}.sh
    echo "Processing ${file%%.*}"
    #submits the job to the queue 
    qsub -cwd ${file%%.*}.sh
    rm ${file%%.*}.sh
done
