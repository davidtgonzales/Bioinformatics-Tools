#!/usr/bin/python
import sys

def help_screen():
    print """
----------------------------------------------------------------
FASTA GC: Calculate the nucleotide percentages of your reads!
by D. Gonzales
               
This program takes a fasta file and calculates the ATGC percentages
for each read and stores these values in a CSV file format:
                   
ID,A,T,G,C
read_1,20,20,30,30
read_2,10,30,20,40

Usage: fasta_gc.py [-h] [-i in] [-o out]
[-h] = This helpful help screen.
[-i] = FASTA input file. 
[-o] = CSV output file.

*Labels should not contain white spaces.
*Requirements: python 2.7, sys, getopt, biopython. 
----------------------------------------------------------------"""
    sys.exit()

def fasta_gc(infile,outfile):
    from Bio import SeqIO
    
    outfile=open(outfile,'w')
    outfile.write('id'+','+'A'+','+'T'+','+'G'+','+'C')
    infile=SeqIO.index(infile,'fasta')
    keys=infile.keys()
    for key in keys:
        seq=str(infile[key].seq).upper()
        A=str(round(float(seq.count('A'))/float(len(seq))*100,2))
        T=str(round(float(seq.count('T'))/float(len(seq))*100,2))
        G=str(round(float(seq.count('G'))/float(len(seq))*100,2))
        C=str(round(float(seq.count('C'))/float(len(seq))*100,2))
        outfile.write('\n'+key+','+A+','+T+','+G+','+C)

    outfile.close()

def main():
    import getopt
    try:
        opts,args=getopt.getopt(sys.argv[1:],'i:o:p:h')
    except getopt.GetoptError as err:
        print str(err)
        print 'Usage: fasta_gc.py [-h] [-i in] [-o out]'
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
            print 'Usage: fasta_gc.py [-h] [-i in] [-o out]'
            sys.exit()      
    
    fasta_gc(infile,outfile)
    
if __name__=='__main__':
    main()
    sys.exit()