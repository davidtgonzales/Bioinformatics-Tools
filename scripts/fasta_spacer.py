#!/usr/bin/python
import sys

def help_screen():
    print """
----------------------------------------------------------------
FASTA SPACER: Add or remove white spaces in your fasta files!
by D. Gonzales
               
This program takes a multi-fasta file and adds or removes a white 
space in between each seqeunce.
       
Usage: fasta_spacer.py [-h] [-s str] [-i in] [-o out]
[-h] = This helpful help screen.
[-s] = Option to add or remove white spaces (use 'add' or 'rem').
[-i] = FASTA input file. 
[-o] = FASTA output file.

*Labels should not contain white spaces.
*Requirements: python 2.7, sys, getopt. 
----------------------------------------------------------------"""
    sys.exit()

def fasta_spacer(space,infile,outfile):
    infile=open(infile,'r')
    outfile=open(outfile,'w')
    start=0
    if space=='add':
        for line in infile:
            if line[0]=='>' and start!=0:
                outfile.write('\n'+line)
            else:
                outfile.write(line)
                start=1
    
    if space=='rem':
        for line in infile:
            if line.strip()!='' and start!=0:
                outfile.write('\n'+line.strip())
            elif line.strip()!='' and start==0:
                outfile.write(line.strip())
                start=1
    
    outfile.close()            
    
def main():
    import getopt
    try:
        opts,args=getopt.getopt(sys.argv[1:],'s:i:o:h')
    except getopt.GetoptError as err:
        print str(err)
        print 'Usage: fasta_spacer.py [-h] [-s str] [-i in] [-o out]'
        sys.exit()
    
    for o,a in opts:
        if o=='-h':
            help_screen()
        if o=='-s':
            space=str(a)
            if space!='add' and space!='rem':
                print 'option error'
                print 'Use \'add\' or \'rem\' only for -s option'
                sys.exit()
        if o=='-i':
            infile=str(a)
        if o=='-o':
            outfile=str(a)         
    
    options=[opt[0] for opt in opts]
    for option in ['-s','-i','-o']:
        if option not in options:
            print 'missing arguement'
            print 'Usage: fasta_spacer.py [-h] [-s str] [-i in] [-o out]'
            sys.exit()      
    
    fasta_spacer(space,infile,outfile)

if __name__=='__main__':
    main()
    sys.exit()