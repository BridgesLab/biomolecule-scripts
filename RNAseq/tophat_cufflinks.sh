#Before running this script, the following links need to be established:
#reference-annotations folder is symlinked
#presumes the fasta files are in the directory called sequence-files

#this is the GTF file for the assembly
GTF="reference-annotations/Homo_sapiens.GRCh37.69.gtf"
#this is the name (without the .fa) of the index built with bowtie2-build
REFERENCE="reference-annotations/Homo_sapiens.GRCh37.69"

echo "This alignment is using bowtie with the $REFERENCE genome"
echo "The alignments and assemblies are made using $GTF"

#remove existing alignments and assemblies
rm -r tophat_out
rm -r cufflinks_out

mkdir tophat_out
mkdir cufflinks_out

for file in `ls sequence-files/*.fasta | xargs -n1 basename` ; do
  #run tophat alignment
  #this initializes the output directories to avoid a filesystem detection problem in glusterfs
  echo "stat tophat_out cufflinks_out" > ${file%%.*}.sh
  #this uses multiple processors (-p 11) and to first map to known sequences (-G) before matching other sequences
  echo "tophat2 -p 11 -G $GTF -o tophat_out/${file%%.*} $REFERENCE $sample.fa" >> ${file%%.*}.sh
  #The options are to use multiple cores (-p 11), to do a soft RABT assembly (-g) and to specify the output directory (-o)
  echo "cufflinks -p 11 -u -b -g $GTF --max-bundle-frags 9999999999 -o cufflinks_out/${file%%.*} tophat_out/${file%%.*}/accepted_hits.bam" >> $sample.sh
  echo "cufflinks_out/${file%%.*}/transcripts.gtf" >> assemblies.txt
  qsub -cwd $sample.sh
  rm $sample.sh
done

