#!/usr/bin/python
import sys

def help_screen():
    print """
----------------------------------------------------------------
FASTA LABEL SPACE REMOVER: Remove the spaces and commas in the 
labels of a fasta file!
by D. Gonzales
               
This program removes all the spaces in the labels of a fasta file:

>gi|30018278|ref|NC_004722.1| Bacillus cereus ATCC 14579, complete genome
TAGCCACTTTTTTTTGATATTATAGTTGTGTTTTCACTTTGAATAAGTTTTCCACATCTTTATCTTATCC
ACAATTTGTGTATAACATGTGGACAGTTTTAATCACATGTGGGTAAATAGTTGTCCACATTTGCTTTTTT

>PREFIX|gi|30018278|ref|NC_004722.1|BacilluscereusATCC14579completegenome
TAGCCACTTTTTTTTGATATTATAGTTGTGTTTTCACTTTGAATAAGTTTTCCACATCTTTATCTTATCC
ACAATTTGTGTATAACATGTGGACAGTTTTAATCACATGTGGGTAAATAGTTGTCCACATTTGCTTTTTT

       
Usage: fasta_label_space_remover.py [-h] [-i in] [-o out]
[-h] = This helpful help screen.
[-i] = FASTA DNA input file.
[-o] = FASTA protein output file.

*Labels should not contain white spaces.
*Requirements: python 2.7, sys, getopt. 
----------------------------------------------------------------"""
    sys.exit()

def fasta_label_space_remover(infile,outfile):
    infile=open(infile,'r')
    outfile=open(outfile,'w')
    start=0
    for line in infile:
        if line[0]=='>':
            outfile.write(''.join(line.split()).replace(',','')+'\n')
        else:
            outfile.write(line.upper())
            
    outfile.close()

def main():
    import getopt
    try:
        opts,args=getopt.getopt(sys.argv[1:],'i:o:h')
    except getopt.GetoptError as err:
        print str(err)
        print 'Usage: fasta_label_space_remover.py [-h] [-i in] [-o out]'
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
            print 'Usage: fasta_label_space_remover.py [-h] [-i in] [-o out]'
            sys.exit()      
    
    fasta_label_space_remover(infile,outfile)

if __name__=='__main__':
    main()
    sys.exit()
