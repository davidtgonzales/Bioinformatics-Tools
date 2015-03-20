#!/usr/bin/python
import sys

def help_screen():
    print """
----------------------------------------------------------------
FASTA TRANSLATOR: Translate a fasta file according to a given frame!
by D. Gonzales
               
This program takes a multi-fasta file and a list of IDs and frames
to translate DNA to protein sequences.
       
Usage: fasta_translator.py [-h] [-i in] [-l list] [-o out]
[-h] = This helpful help screen.
[-i] = FASTA DNA input file.
[-l] = List of IDs and frames in csv format (ID,frame). 
[-o] = FASTA protein output file.

*Labels should not contain white spaces.
*Requirements: python 2.7, sys, getopt. 
----------------------------------------------------------------"""
    sys.exit()

def fasta_translator(infile,inlist,outfile):	
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.Alphabet import IUPAC

	infile=SeqIO.index(infile,'fasta')
	inlist=open(inlist,'r').readlines()
	outfile=open(outfile,'w')
	for item in inlist:
		seqid=item.strip().split(',')[0]
		frame=int(item.strip().split(',')[1])
		if frame>0:
			protein=str(infile[seqid].seq[frame-1:].translate())
			outfile.write('>'+item.strip()+'\n'+protein)

		elif frame<0:
			protein=str(infile[seqid].seq.reverse_complement()[(-1*frame)-1:].translate())
			outfile.write('>'+item.strip()+'\n'+protein)


		if item!=inlist[len(inlist)-1]:
			outfile.write('\n')

	outfile.close()

def main():
    import getopt
    try:
        opts,args=getopt.getopt(sys.argv[1:],'i:l:o:h')
    except getopt.GetoptError as err:
        print str(err)
        print 'Usage: fasta_translator.py [-h] [-i in] [-l list] [-o out]'
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
            print 'Usage: fasta_translator.py [-h] [-i in] [-l list] [-o out]'
            sys.exit()      
    
    fasta_translator(infile,inlist,outfile)

if __name__=='__main__':
    main()
    sys.exit()
