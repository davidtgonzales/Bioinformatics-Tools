#!/usr/bin/python
import sys

def help_screen():
    print """
----------------------------------------------------------------
FASTA EXTRACT: Extract sequences from a multi-fasta file!
by D. Gonzales
               
This program takes a multi-fasta file and a list of read IDs and
extracts the sequences of the IDs into a new file.
       
Usage: fasta_extract.py [-h] [-i in] [-l list] [-o out]
[-h] = This helpful help screen.
[-i] = FASTA input file.
[-l] = List of sequence IDs. 
[-o] = FASTA output file.

*Labels should not contain white spaces.
*Requirements: python 2.7, sys, getopt, biopython. 
----------------------------------------------------------------"""
    sys.exit()

def fasta_extract(infile,inlist,outfile):	
	from Bio import SeqIO

	infile=SeqIO.index(infile,'fasta')
	inlist=open(inlist,'r').readlines()
	outfile=open(outfile,'w')
	for item in inlist:
		item=item.strip().split(',')[0]
		outfile.write('>'+item+'\n'+str(infile[item].seq))
		if item!=inlist[len(inlist)-1]:
			outfile.write('\n')

	outfile.close()


def main():
    import getopt
    try:
        opts,args=getopt.getopt(sys.argv[1:],'i:l:o:h')
    except getopt.GetoptError as err:
        print str(err)
        print 'Usage: fasta_extract.py [-h] [-i in] [-l list] [-o out]'
        sys.exit()
    
    for o,a in opts:
        if o=='-h':
            help_screen()
        if o=='-i':
            infile=str(a)
        if o=='-l':
            inlist=str(a)
        if o=='-o':
            outfile=str(a)         
    
    options=[opt[0] for opt in opts]
    for option in ['-i','-l','-o']:
        if option not in options:
            print 'missing arguement'
            print 'Usage: fasta_extract.py [-h] [-i in] [-l list] [-o out]'
            sys.exit()      
    
    fasta_extract(infile,inlist,outfile)

if __name__=='__main__':
    main()
    sys.exit()
