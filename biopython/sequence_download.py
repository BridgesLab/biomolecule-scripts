#!/usr/bin/python
'''This script takes a results file, generated by blast_search.py and downloads all the FASTA sequences for these results.

It is invoked simply as:

    python sequence_download.py
    
and it generates a file named blast_results.fasta containing all the sequences in FASTA format.
'''

import csv
from collections import OrderedDict

from Bio import Entrez, SeqIO

Entrez.email = "dave.bridges@gmail.com"

datafile = open('results.csv', 'r')
datareader = csv.reader(datafile)
data = []
for row in datareader:
    data.append(row)

gi_numbers = []    
for x in range(0,len(data)):
    gi_numbers.append(data[x][3])
unique_gi_numbers = OrderedDict.fromkeys(gi_numbers).keys()

records = []        
for gi in unique_gi_numbers:
    handle = Entrez.efetch(db="protein", id=gi, rettype="fasta", retmode="text")
    records.append(SeqIO.read(handle, "fasta"))
    
with open('blast_results.fasta', 'a') as outfile:
    for record in records:
        outfile.write(record.format("fasta"))

