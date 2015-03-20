#!/usr/bin/python
import sys

def help_screen():
    print """
----------------------------------------------------------------
FASTA DEREPLICATOR: Remove sequence replicates in a fasta file!
by D. Gonzales
               
This program takes a fasta file and removes identical sequences 
that appear more than once, even as a substring.
       
Usage: fasta_dereplicator.py [-h] [-i in] [-o out]
[-h] = This helpful help screen.
[-i] = Input FASTA file.
[-o] = Output FASTA file.

*Requirements: python 2.7, sys, getopt. 
----------------------------------------------------------------"""

def fasta_dereplicator(infile,outfile):
    infile=open(infile,'r')
    contigs=[]
    for line in infile:
        if line[0]!='>':
            contigs.append(line.strip().upper())

    derep1=list(set(contigs))
    joined='#'.join(derep1)
    derep2=[]
    for item in derep1:
        if joined.count(item)==1:
            derep2.append(item)

    outfile=open(outfile,'w')
    index=0
    for item in derep2:
        if index==0:
            outfile.write('>contig_'+str(index)+'\n'+item)            
        else:
            outfile.write('\n'+'>contig_'+str(index)+'\n'+item)
        index=index+1

    outfile.close()


def main():
    import getopt
    try:
        opts,args=getopt.getopt(sys.argv[1:],'i:o:h')
    except getopt.GetoptError as err:
        print str(err)
        print 'Usage: fasta_dereplicator.py [-h] [-i in] [-o out]'
        sys.exit()

    for o,a in opts:
        if o in ('-h'):
            help_screen()
        if o=='-i':
            infile=str(a)
        if o=='-o':
            outfile=str(a)

    for option in ['-i','-o']:
        if option not in [opt[0] for opt in opts]:
            print 'missing arguement'
            print 'Usage: fasta_dereplicator.py [-h] [-i in] [-o out]'
            sys.exit() 
    
    fasta_dereplicator(infile,outfile)
            
if __name__=='__main__':
    main()
    sys.exit()
