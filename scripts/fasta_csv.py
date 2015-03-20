#!/usr/bin/python
import sys

def help_screen():
    print """
----------------------------------------------------------------
FASTA TO CSV: Turn a FASTA into a CSV file!
by D. Gonzales
               
This program takes a multi-fasta file and converts it to
a csv file.

Input:
>seq1
ATGC
>seq2
ATTC

Output:
seq1,ATGC
seq2,ATTC

Usage: fasta_csv.py [-h] [-i in] [-o out]
[-h] = This helpful help screen.
[-i] = FASTA input file. 
[-o] = CSV output file.

*Labels should not contain white spaces.
*Requirements: python 2.7, sys, getopt, biopython. 
----------------------------------------------------------------"""
    sys.exit()

	
def fasta_csv(infile,outfile):	
	from Bio import SeqIO

	infile=SeqIO.index(infile,'fasta')
	keys=infile.keys()
	outfile=open(outfile,'w')

	outfile.write('key,seq\n')
	for key in keys:
		outfile.write(key+','+str(infile[key].seq)+'\n')
    
	outfile.close()


def main():
    import getopt
    try:
        opts,args=getopt.getopt(sys.argv[1:],'i:o:h')
    except getopt.GetoptError as err:
        print str(err)
        print 'Usage: fasta_csv.py [-h] [-i in] [-o out]'
        sys.exit()
    
    for o,a in opts:
        if o=='-h':
            help_screen()
        if o=='-i':
            infile=str(a)
        if o=='-o':
            outfile=str(a)         
    
    options=[opt[0] for opt in opts]
    for option in ['-i','-o']:
        if option not in options:
            print 'missing arguement'
            print 'Usage: fasta_csv.py [-h] [-i in] [-o out]'
            sys.exit()      
    
    fasta_csv(infile,outfile)

if __name__=='__main__':
    main()
    sys.exit()
