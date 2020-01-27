#!/usr/bin/env python

'''
Purpose:
    This simple script identifies all genes in a GFF file with multiple exons,
    from a list of gene sequence IDs (column 1 in the GFF file).
    Often, intronic sequences are not listed in a GFF file, thus multiple
    exons are used as a way to determine if a gene has an intron. This method
    may be inherently flawed, thus careful annotation of these genes is 
    required.
Note: 
    The SeqId file must contain the SeqID in the first column, one SeqId per
    line.
    Example:
        SeqId1
        SeqId2
        SeqId3
Usage: 
    python SeqID.txt GFFfile.gff
Output: 
    Detailed txt and GFF files of multi-exonic genes and single exon genes,
    seperately.
'''

import gffutils
import sys
import os

SeqIdfile = sys.argv[1]
GffFile = sys.argv[2]

# Gather all SeqIds into unique list
SeqIdList = []
with open(SeqIdfile,'r') as SiFile:
    for line in SiFile:
        cols = line.split('\t')
        SeqIdList.append(str(cols[0]))
    SiFile.close()

# Insure that all seqIds are unique
SeqIdSet = set(SeqIdList)
SeqIdList = list(SeqIdSet) 
   
# Create GFF DB
DbName = str(GffFile) + ".db"
if not os.path.isfile(DbName): 
    GffDb = gffutils.create_db(GffFile, dbfn=DbName, force=True, 
                                keep_order=True, merge_strategy='merge', 
                                sort_attribute_values=True, 
                                force_gff = True)
elif os.path.isfile(DbName):
    GffDb = gffutils.FeatureDB(DbName)



# Query GffDb with unique SeqIDs
IntronicList = []
IntronGeneList = []
IntronDetailedList = []
ExonDetailedList = []
OneExonList = []
for SeqId in SeqIdList:
    gene = GffDb[SeqId]
    Count = len(
            list(GffDb.children(gene, featuretype='exon', order_by='start')))
    if Count >= 1:
        IntronGeneList.append(gene)
        IntronicList.append(SeqId)
        IntronDetailedList.append(SeqId)
        IntronDetailedList.append(Count)
        for i in GffDb.children(gene, featuretype='exon', order_by='start'):
            IntronDetailedList.append(i)
    elif Count == 0:
        OneExonList.append(SeqId)
        ExonDetailedList.append(gene)

# Print out to file
def printout(FileName, Data):
    print "printing out... " + str(FileName)
    with open(FileName, "w") as fileout:
        for Item in Data:
            fileout.write(str(Item) + "\n")
        fileout.close()

# Print out intronic lists
printout("IntronSeqIds.txt", IntronicList)
printout("IntronicGenes.gff", IntronGeneList)
printout("DetailedIntrons.txt", IntronDetailedList)

# Print out one-exon list
printout("OneExonSeqIds.txt", OneExonList)
printout("OneExonGenes.gff", ExonDetailedList)
