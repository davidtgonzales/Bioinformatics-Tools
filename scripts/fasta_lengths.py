#!/usr/bin/python
import sys

def help_screen():
    print """
----------------------------------------------------------------
FASTA LENGTHS: Extract lengths of sequences from a FASTA file!
by D. Gonzales
               
This program takes a multi-fasta file and extracts the lengths
of each sequence and saves it as a csv file.
       
Usage: fasta_lengths.py [-h] [-i in] [-o out]
[-h] = This helpful help screen.
[-i] = FASTA input file. 
[-o] = CSV output file.

*Labels should not contain white spaces.
*Requirements: python 2.7, sys, getopt, biopython. 
----------------------------------------------------------------"""
    sys.exit()

	
def fasta_lengths(infile,outfile):	
	from Bio import SeqIO

	infile=SeqIO.index(infile,'fasta')
	keys=infile.keys()
	outfile=open(outfile,'w')

	outfile.write('key,length\n')
	for key in keys:
		outfile.write(key+','+str(len(infile[key].seq))+'\n')
    
	outfile.close()


def main():
    import getopt
    try:
        opts,args=getopt.getopt(sys.argv[1:],'i:o:h')
    except getopt.GetoptError as err:
        print str(err)
        print 'Usage: fasta_lengths.py [-h] [-i in] [-o out]'
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
            print 'Usage: fasta_lengths.py [-h] [-i in] [-o out]'
            sys.exit()      
    
    fasta_lengths(infile,outfile)

if __name__=='__main__':
    main()
    sys.exit()
