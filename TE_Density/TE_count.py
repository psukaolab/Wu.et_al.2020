#!/usr/bin/env python
'''
Purpose:
This script, using RepeatMasker to generate GFF files containing repeat 
information, plots these transposable elements in relation to the position in a
scaffold or contig. Both window size (bp), and step size (bp) are user defined.

Note: There are two main options: use RepeatMasker to generate a GFF file, OR
provide a GFF file in the exact format that RepeatMasker does.

Note: If RepeatMasker option is used, it must be available to your environment.

Usage:
    python TE_count.py FastaFile -r True
    **OR**
    python TE_count.py FastaFile -gff GffFile.gff

Default optional parameters:
    -s, Step Size, default = 250
    -w, Window Size, default = 1000
    -gff, User Provided Repeatmasker GFF file, default="empty.gff"
    -r, RepeatMasker, default = False
'''

import gffutils
import argparse

import sys
import os

from Bio import SeqIO
import TEcount_modules

# Capture command line args, with or without defaults
if __name__ == '__main__':
    # Parse the arguments
    LineArgs = TEcount_modules.parseArguments()

# Populate vars with args
FastaFile = LineArgs.FastaFile
GffFile = LineArgs.gff
Step = LineArgs.s
WinSize = LineArgs.w
RepMask = LineArgs.r

# Run RepeatMasker
if RepMask == True:
    print "How many threads for RepeatMasker?"
    threads = raw_input('Threads? (numeric)... ')
    print "Which species for RepeatMasker?"
    species = raw_input('species? (exact string)... ')
    try:   
        command = ('RepeatMasker -s -nolow -gff -pa ' 
                    + str(threads) + ' -species "' + str(species)
                    + '" ' + str(FastaFile))
        os.system('%s' % command)
        GffFile = str(FastaFile) + ".out.gff"
    except IndexError as e:
        print 'Error: '+ str(e)
        sys.exit('Exiting script...')
elif RepMask == False and GffFile == "empty.gff":
    sys.exit('Enter correct gff file of run RepeatMasker.\nExiting script...')

# Edit RepeatMasker GFF (not compatible with GFFUtils)
Check = TEcount_modules.gff_check(GffFile)
if Check == True:
    tempfile = "temp.gff"
    GffFile = TEcount_modules.edit_gff(GffFile, tempfile)
elif Check == False:
    sys.exit('Exiting script...')

# Create temporary GFF DB
GffDb = gffutils.create_db(tempfile, dbfn='tempGFF.db', force=True, 
                            keep_order=True, merge_strategy='merge', 
                            sort_attribute_values=True, 
                            force_gff = True)

# Step through each sequence, initate classes, and
# query GFF database and generate plots
for SeqRecord in SeqIO.parse(FastaFile, "fasta"):
    print SeqRecord.id

    # Determine if sequences and args are acceptable
    TEcount_modules.arg_seqcheck(SeqRecord, WinSize, Step)

    # Initialize TESeqRec class
    TESeqRecord = TEcount_modules.TESeqRec(SeqRecord.id, 
                                                len(SeqRecord), WinSize, Step)

    # Compute window ranges, get TE counts by bin
    ForGGplot = TESeqRecord.get_tes(GffDb)

    # Generate figure
    TotalTE = sum(ForGGplot[2])
    if TotalTE <= 0:
        print("Warning, " + str(TotalTE) 
                + " TEs found for " + str(SeqRecord.id) 
                + ", consider changing parameters")
    else: 
        TEcount_modules.gen_fig(SeqRecord.id, ForGGplot)

# Clean up files
print "Cleaning Files..."
os.remove(tempfile)
os.remove("tempGFF.db")
