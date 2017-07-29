#!/usr/bin/python
import sys

def help_screen():
    print """
----------------------------------------------------------------
FASTA PREFIX: Add a prefix to the labels of a fasta file!
by D. Gonzales
               
This program takes a multi-fasta file and rewrites the labels 
on each sequence with an input prefix:

>gi|30018278|ref|NC_004722.1| Bacillus cereus ATCC 14579, complete genome
TAGCCACTTTTTTTTGATATTATAGTTGTGTTTTCACTTTGAATAAGTTTTCCACATCTTTATCTTATCC
ACAATTTGTGTATAACATGTGGACAGTTTTAATCACATGTGGGTAAATAGTTGTCCACATTTGCTTTTTT

>PREFIX|gi|30018278|ref|NC_004722.1| Bacillus cereus ATCC 14579, complete genome
TAGCCACTTTTTTTTGATATTATAGTTGTGTTTTCACTTTGAATAAGTTTTCCACATCTTTATCTTATCC
ACAATTTGTGTATAACATGTGGACAGTTTTAATCACATGTGGGTAAATAGTTGTCCACATTTGCTTTTTT

       
Usage: fasta_prefix.py [-h] [-i in] [-p prefix] [-o out]
[-h] = This helpful help screen.
[-i] = FASTA input file.
[-p] = Prefix input. 
[-o] = FASTA output file.

*Labels should not contain white spaces.
*Requirements: python 2.7, sys, getopt. 
----------------------------------------------------------------"""
    sys.exit()

def fasta_prefix(prefix,infile,outfile):
    infile=open(infile,'r')
    outfile=open(outfile,'w')
    start=0
    for line in infile:
        if line[0]=='>':
            outfile.write('>'+prefix+'|'+line[1:])
        else:
            outfile.write(line)
            
    outfile.close()

def main():
    import getopt
    try:
        opts,args=getopt.getopt(sys.argv[1:],'p:i:o:h')
    except getopt.GetoptError as err:
        print str(err)
        print 'Usage: fasta_prefix.py [-h] [-i in] [-p prefix] [-o out]'
        sys.exit()
    
    for o,a in opts:
        if o=='-h':
            help_screen()
        if o=='-p':
            prefix=str(a)
        if o=='-i':
            infile=str(a)
        if o=='-o':
            outfile=str(a)         
    
    options=[opt[0] for opt in opts]
    for option in ['-p','-i','-o']:
        if option not in options:
            print 'missing arguement'
            print 'Usage: fasta_prefix.py [-h] [-i in] [-p prefix] [-o out]'
            sys.exit()      
    
    fasta_prefix(prefix,infile,outfile)

if __name__=='__main__':
    main()
    sys.exit()
