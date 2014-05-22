#Invoke this in a directory containing tophat_out as:
#Rscript counts_tables.R
#This should generate csv files for both exon and transcript counts.
#THis presumes you used the current ensembl mouse genome release
require(GenomicFeatures)
require(biomaRt)
#make a database of transcripts from the ensembl assembly
#first get the current release from ensembl
txdb <- makeTranscriptDbFromBiomart(biomart="ensembl",dataset = 'mmusculus_gene_ensembl')
#make exon and transcript annotation objects
exons <- exons(txdb, columns=c('gene_id', 'tx_id', 'tx_name', 'tx_chrom', 'tx_strand', 'exon_id', 'exon_name', 'cds_id', 'cds_name', 'cds_chrom', 'cds_strand', 'exon_rank'))
transcripts <- transcripts(txdb, columns=c('gene_id', 'tx_id', 'tx_name', 'exon_id', 'exon_name', 'exon_chrom', 'exon_strand', 'cds_id', 'cds_name', 'cds_chrom', 'cds_strand', 'exon_rank'))
require(GenomicRanges)
require(Rsamtools)

#set list of sample ids as a vector
sample_ids = dir('tophat_out')

transcript.countsTable <- data.frame(
  row.names = as.vector(unlist(elementMetadata(transcripts)['tx_name'])))
exon.countsTable <- data.frame(
  row.names = as.vector(unlist(elementMetadata(exons)['exon_name'])))

#this forloop iterates over the sample_ids and generates exon and transcript counts for each sample_id
for(sample_id in sample_ids) {
  #read alignment
  align <- readGAlignmentsFromBam(sprintf("tophat_out/%s/accepted_hits.bam", sample_id))
  #count the overlapping reads for the transcripts
  transcript.counts <- countOverlaps(transcripts, align)
  #reassign to a specific transcript.counts object.
  assign(sprintf("transcript.counts.%s", sample_id), transcript.counts)
  #add this column to the countsTable
  transcript.countsTable <- cbind(transcript.countsTable,transcript.counts)
  remove(transcript.counts)
  #count the overlapping reads for the exons
  exon.counts <- countOverlaps(exons, align)
  #reassign to a specific transcript.counts object.
  assign(sprintf("exon.counts.%s", sample_id), exon.counts)
  #add this column to the countsTable
  exon.countsTable <- cbind(exon.countsTable, exon.counts)  
  #remove the align, transcript.counts and exon.counts objects for the next loop
  remove(align)
  remove(exon.counts)
  }

summary(transcript.countsTable)
summary(exon.countsTable)

#write these two counts tables to csv files.
write.csv(transcript.countsTable, "transcript_counts_table.csv")
write.csv(exon.countsTable, "exon_counts_table.csv")
