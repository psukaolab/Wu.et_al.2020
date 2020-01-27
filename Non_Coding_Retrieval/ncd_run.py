#!/usr/bin/env python
'''
Purpose:
non-coding sequence retrieval surrounding a gene or sequence of interest
from large genomic scaffolds. 

This script uses BLAST, so the  exact location of the gene is not required. 
Because there may be multiple matching HSPs (high scoring segment pairs) 
produced by BLAST for each gene of interest, this script outputs all HSPs of a 
minimum percent similarity (provided by the user).

Prerequisites:
ncd_run.py requires ncdmodules.py in the same directory
this script also requires Biopython and ncbi-blast+

Example:
python ncd_run.py protocadherinBs.fasta human_scaffolds.fasta 90

Output:
Three fasta files per genomic region of interest: *_upstream.fasta, 
*_downstream.fasta, and *_fullseq.fasta
*=These files are named after the genomic scaffold queried, so a brief name in
the fastafile will be desired, i.e. >CHR1, >CHR2, >scaffold34234, etc.

The sequences in fasta format are named after the corresponding matching gene, 
the percent similarity to that gene, and the location of the gene (based on the
scaffold queried)

example:
>protocadherinB-1|99%|1600000
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG


Usage:
python ncd_run.py Genes.fasta genome.fasta perID

where
	genes.fasta = list of genes in 5->3 orientation
	genome.fasta = genome of interest containing gene of interest
    perID = minimum percent similarity to the gene

NOTE: percent similarity is written without "%"
percent similarity represents the similarity of the HSP produced by BLAST
'''

import sys
import os
import imp
import subprocess
import shlex
import fnmatch
import tempfile

import Bio
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

#Alter this variable to change the bps retrieved upstream & downstream
SeqLen = int(10000)

#find requisite libraries
def findlibs():
    try:
        imp.find_module('Bio')
    except ImportError:
        print 'install Bio libraries (Biopython)'
        sys.exit('Exiting script...')

#making temporary blast db
def makeblastdb(FastaFile):
    DBname = FastaFile.replace('.fasta','.db')
    try:
        Command = ('makeblastdb -in '
                    + str(FastaFile)
                    + ' -input_type fasta -dbtype nucl -parse_seqids -out ' 
                    + str(DBname))
        Args = shlex.split(Command)
        subprocess.Popen(Args, stdout=subprocess.PIPE)
    except OSError as e:
        print "Error({0}): {1}".format(e.errno, e.strerror)
        print "Insure ncbi-blast+ is installed and available to environment"

#delete temp DB
def delete_db(DBname):
    File_Del = str(DBname) + '.*'
    for file in os.listdir('.'):
        if fnmatch.fnmatch(file, str(File_Del)):
            print file
            os.remove(str(file))

#Checking temporary blastdb
def check_blastdb(FastaFile):
    DBname = FastaFile.replace('.fasta','.db')
    while True:
        if os.path.isfile(str(DBname)+'.nhr'):
            print ('blast DB named... ' 
                    + str(DBname) 
                    + ' already exists. Overwrite?')
            Reply = raw_input('(Y/N)... ')
            if Reply == 'Y':
                makeblastdb(FastaFile)
                return DBname
            elif Reply == 'N':
                sys.exit('Exiting script...')
            else:
                print 'Incorrect input. Use Y/N'
        else:
            makeblastdb(FastaFile)
            return DBname

#Parse blast hits and Records for sequences
def parse_seq(Blast_Record, Hsp, Record, SeqLen):
    # library: from Bio import Seq for reverse complimentation    
    QueryStart = int(Hsp.query_start)
    QueryStop = int(Hsp.query_start) + int(Hsp.align_length)

    #Get Coordinates
    if int(Hsp.query_start) - SeqLen > 0:
        MaxStart = int(Hsp.query_start) - SeqLen
    else:
        MaxStart = 1
    if int(Blast_Record.query_length) - int(QueryStop + SeqLen) > 0:
        MaxStop = int(QueryStop + SeqLen)
    else:
        MaxStop = int(Blast_Record.query_length)

    #Extract (with gene/or seq of interest) 
    #upstream and downstream regardless of orientation
    Seq = Record.seq
    FullSeq = Seq[MaxStart:MaxStop]
    UpSeq = Seq[MaxStart:QueryStop]
    DownSeq = Seq[QueryStart:MaxStop]

    #Output 5'->3' oriented sequences, 
    #NOTE:We are assuming gene/seq of interest is in 5'->3' orientation
    Cond1 = int(Hsp.sbjct_start) >= int(Hsp.align_length)
    Cond2 = int(Hsp.sbjct_start) >= 1
    Cond3 = int(Hsp.sbjct_start) < int(Hsp.align_length)
    if Cond1: #TRUE
        #revcomp sequences
        UpOut = DownSeq.reverse_complement()
        DownOut = UpSeq.reverse_complement()
        FullOut = FullSeq.reverse_complement()
        return UpOut, DownOut, FullOut
    elif Cond2 and Cond3: # FALSE
        print("Partial match, may need to check strand orientation: " 
                + str(Record.id))
        #nothing happens, no revcomp
        return UpSeq, DownSeq, FullSeq

#remove characters 
def remove(filename, deletechars='\/:*?"<>|'):
    for c in deletechars:
        filename = filename.replace(c,'')
    return filename;

#Print out results by blast Record
def print_out(SeqList, Headers, FlPt1, FlPt2):
    FlPt1 = remove(FlPt1)
    OutFile = str(FlPt1) + str(FlPt2)
    OutHandle = open(OutFile, 'ab+')
    for x in range(0, len(SeqList)):
        OutHandle.write('>%s\n'%(Headers[x]))
        OutHandle.write('%s\n'%(SeqList[x]))
    OutHandle.close()

#BLAST genomic seq against temp gene database
def blast(FastaFile, BlastDB, perID, SeqLen):
    '''
    libraries:
        from Bio import SeqIO
        from Bio.Blast.Applications import NcbiblastnCommandline
        from Bio.Blast import NCBIXML
    '''    
    Fasta_Handle = open(FastaFile, "r")

    for Record in SeqIO.parse(Fasta_Handle, "fasta"):
        #generate temporary fasta file input, and BLASTxml output
        TempFasta = tempfile.NamedTemporaryFile()
        TempFasta.write(">%s\n%s\n" % (Record.id, Record.seq))
        TempBlastXML = tempfile.NamedTemporaryFile()

        #BLAST Record
        Blast_Command = NcbiblastnCommandline(
            query=TempFasta.name, 
            db=BlastDB, evalue=1e-10, 
            out=TempBlastXML.name, outfmt=5)
        std_output, err_output = Blast_Command()
        TempFasta.close()
        Result_Handle = open(TempBlastXML.name)
        Blast_Records = NCBIXML.parse(Result_Handle)

        #lists
        UpList = []
        DownList = []
        FullSeqList = []
        Headers = []

        #loop over Records, check perID
        for Blast_Record in Blast_Records:
            for Alignment in Blast_Record.alignments:
                for Hsp in Alignment.hsps:
                    Hsp_perID = ((float(Hsp.positives)
                                /float(Hsp.align_length)) 
                                * 100)
                    if Hsp_perID >= int(perID):
                        #call seq function
                        UpStream, DownStream, FullSeq = parse_seq(
                            Blast_Record, Hsp, Record, SeqLen)
                        #create list of seqs and Headers
                        UpList.append(UpStream)
                        DownList.append(DownStream)
                        FullSeqList.append(FullSeq)
                        #create header
                        Sbjct_Name = Alignment.title
                        Sbjct_Edit = Sbjct_Name.replace(
                            ' No definition line', '')
                        Header_String =  (str(Sbjct_Edit)
                                            + '|' + str(round(Hsp_perID,1)) 
                                            + '%' + '|' + str(Hsp.query_start))
                        Headers.append(Header_String)

        #print out and close
        print_out(UpList, Headers, Record.id, '_upstream.fasta')
        print_out(DownList, Headers, Record.id, '_downstream.fasta')
        print_out(FullSeqList, Headers, Record.id, '_fullseq.fasta')
        Result_Handle.close()

#parameter input
try:
    GeneFile = sys.argv[1]
    GenomicFile = sys.argv[2]
    PerID = sys.argv[3]
except IndexError as e:
    print 'Error: '+ str(e)
    print "Incorrect parameter usage, see instructions"
    sys.exit('Exiting script...')

#Check for biopython libs
findlibs()

#Check if DB already exists. Also calls makeblastdb
DBname = check_blastdb(GeneFile)
print "Making temp DB file..."
print DBname

#performs BLAST, calls parse_seq and print_out
print 'printing out files... '
blast(GenomicFile, DBname, PerID, SeqLen)

#delete temp DB
print 'deleting temp DB files... '
delete_db(DBname)
