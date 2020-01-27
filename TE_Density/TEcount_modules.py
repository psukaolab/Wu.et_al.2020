#!/usr/bin/env python
'''
module for TE_count script
'''
import Bio.SeqUtils
import gffutils
import argparse
import os
import sys
import re

import pandas as pd
import numpy as np
from ggplot import *

# Retrieve args from command line, usage help
def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("FastaFile", help="Fasta File.", type=str)

    # Optional arguments
    parser.add_argument("-s", help="Step Size.", type=int, default=250)
    parser.add_argument("-w", help="Window Size.", type=int, default=1000)
    parser.add_argument("-r", help="Run RepeatMasker.", 
                        type=bool, default=False)
    parser.add_argument("-gff", help="Matching Repeatmasker GFF File.",
                        type=str, default="empty.gff")

    # Parse arguments
    args = parser.parse_args()

    return args

# Class contains TE by window, sequence length
class TESeqRec:
    def __init__(self, SeqId, SeqLength, WinSize, Step):
        self.SeqId = SeqId
        self.SeqLength = SeqLength
        self.WinSize = WinSize
        self.Step = Step

        # Pre-Determine the number of TE ranges
        NumOfChunks = ((self.SeqLength-self.WinSize)/self.Step)+1
        self.NumOfChunks = NumOfChunks

    # Create Ranges formatted for GFFutils
    def get_tes(self, GffDb):
        Region = ""
        BpBin = []
        SeqIdList = []
        TEcounts = []
        ForGGplot = []
        for i in range(0, self.NumOfChunks*self.Step, self.Step):
            WinStart = i
            WinStop = i + self.WinSize
            Region = str(self.SeqId) + ":" + str(WinStart) + "-" + str(WinStop)
            TEhits = GffDb.region(strand=None, region=Region, 
                                    completely_within=False)
            hits = len(list(TEhits))
            SeqIdList.append(self.SeqId)
            BpBin.append(WinStop)
            TEcounts.append(hits)
        ForGGplot.append(SeqIdList)
        ForGGplot.append(BpBin)
        ForGGplot.append(TEcounts)
        return ForGGplot 

# GFF check (insure Repeatmasker format)
def gff_check(GffFile):
    cond1 = False
    cond2 = False
    cond3 = False
    with open(GffFile, "r") as GFF:
        for line in GFF:
            SplitLine = line.split("\t") 
            if len(SplitLine) == 9:
                cond1 = True
                if SplitLine[2] == "similarity":
                    cond2 = True
                if "Motif:" in SplitLine[8]:
                    cond3 = True
                    break
    if cond1 and cond2 and cond3:
        return True
    else:
        print "This RepeatMasker GFF file may not be compatible\n"
        print ("Consider replacing the last (9th) column, or renaming "  
                + "the 3rd column to 'similarity'")
        while True:
            Reply = raw_input('Continue Anyway? (Y/N)... ')
            if Reply == 'Y':
                return True
            elif Reply == 'N':
                sys.exit('Exiting script...')
            else:
                print 'Incorrect input. Use Y/N'

# GFF edit repeatmasker
def edit_gff(GffFile, tempfile):
    try:
        pattern = re.compile(r'\"Motif\:.+\"')
        with open(GffFile, "r") as GFF:
            with open(tempfile, "w") as tempGFF:
                for line in GFF:
                    SplitLine = line.split("\t")
                    if len(SplitLine) < 9:
                        tempGFF.write(str(line))
                    if len(SplitLine) == 9:
                        result = pattern.search(str(SplitLine[8]))
                        if result:
                            string = result.group(0)
                            newstr = string.replace("\"", "")
                            newstr = re.sub('Motif:', 'Motif=', newstr)
                            del SplitLine[-1]
                            SplitLine.append(newstr)
                            tempGFF.write("\t".join(SplitLine) + "\n")
    except IndexError as e:
        print 'Error: '+ str(e)
        sys.exit('error formatting GFF file...')

# Insure WinSize, SeqRecord, and Step are correct lengths
def arg_seqcheck(SeqRecord, WinSize, Step):
    if not ((type(WinSize) == type(0)) and (type(Step) == type(0))):
        raise Exception("**NOTE type(WinSize) and type(Step) must be int.")
    if Step > WinSize:
        raise Exception("**NOTE Step must not be larger than WinSize.")
    if WinSize > len(SeqRecord.seq):
        raise Exception('''**NOTE WinSize must not 
                        be larger than Sequence length.''')

# Generate GGplot of TE density
def gen_fig(SeqId, ForGGplot):
    OutFile = str(SeqId) + ".png"
    try:
        TEdf = pd.DataFrame({
                        'SeqId': ForGGplot[0],
                        'xval': ForGGplot[1],
                        'yval': ForGGplot[2]})

        p = ggplot(aes(x='xval', y='yval'), data=TEdf) \
            + geom_line() \
            + ylab("RepeatMasker TE Density") \
            + xlab("Position (bp) in Scaffold " + str(SeqId)) \
            + ggtitle(str(SeqId) 
                + " Transposable Element Density as predicted by RepeatMasker")
        p.save(OutFile)
    except IndexError as e:
        print 'Error: '+ str(e)
        sys.exit('Exiting script...')
