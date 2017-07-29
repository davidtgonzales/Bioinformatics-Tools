#!/usr/bin/python
import sys

def help_screen():
    print """
----------------------------------------------------------------
FASTA CUTTER: Cut your fasta files into a specified length!
by D. Gonzales
               
This program takes a fasta file and cuts the sequences into a
given length, l. Sequences shorter than the specified length
will be taken whole. The last cut of each sequence will be the
last l bp window. For example, with a specified length of 5 bp:
                   
INPUT          OUTPUT
>seqA          >seqA_0
ATGCA          ATGCA
>seqB          >seqB_0
ATGCATgca      ATGCA
               >seqB_4
               ATgca
       
Usage: fasta_cutter.py [-h] [-l int] [-p str] [-d str] [-i in] [-o out]
[-h] = This helpful help screen.
[-l] = Length to cut the sequences.
[-p] = Prefix to the label (optional).
[-d] = Delimiter label for the prefix and start index (0).
[-i] = FASTA input file. 
[-o] = FASTA output file.

*Labels should not contain white spaces.
*Requirements: python 2.7, sys, getopt, biopython. 
----------------------------------------------------------------"""
    sys.exit()


def fasta_cutter(length,prefix,delimiter,infile,outfile):
    from Bio import SeqIO

    outfile=open(outfile,'w')
    infile=SeqIO.index(infile,'fasta')
    keys=infile.keys()
    
    start=0
    for key in keys:
        index=0
        window=0
        while window<len(infile[key]):             
            if window+length<len(infile[key]):
                if start==0:
                    outfile.write('>'+prefix+key+delimiter+str(window))
                    start=1
                else:
                    outfile.write('\n>'+prefix+key+delimiter+str(window))                  
                outfile.write('\n'+str(infile[key].seq[window:window+length]))
            else:
                if len(infile[key])-length<0:
                    if start==0:
                        outfile.write('>'+prefix+key+delimiter+str(window))
                        start=1
                    else:
                        outfile.write('\n>'+prefix+key+delimiter+str(window))                      
                    outfile.write('\n'+str(infile[key].seq))
                else:
                    outfile.write('\n>'+prefix+key+delimiter+str(len(infile[key])-length))
                    outfile.write('\n'+str(infile[key].seq[len(infile[key])-length:]))
            index=index+1	
            window=window+length
    outfile.close()
    
    
def main():
    import getopt
    try:
        opts,args=getopt.getopt(sys.argv[1:],'l:d:i:o:p:h')
    except getopt.GetoptError as err:
        print str(err)
        print 'Usage: fasta_cutter.py [-h] [-l int] [-p str] [-d str] [-i in] [-o out]'
        sys.exit()
    
    prefix=''
    for o,a in opts:
        if o=='-h':
            help_screen()
        if o=='-l':
            length=int(a)
        if o=='-p':
            prefix=str(a)
        if o=='-d':
            delimiter=str(a)
        if o=='-i':
            infile=str(a)
        if o=='-o':
            outfile=str(a)         
    
    if prefix!='':
        prefix=prefix+delimiter
    options=[opt[0] for opt in opts]
    if '-p' in options:
        options.remove('-p')
    for option in ['-l','-d','-i','-o']:
        if option not in options:
            print 'missing arguement'
            print 'Usage: fasta_cutter.py [-h] [-l int] [-p str] [-d str] [-i in] [-o out]'
            sys.exit()      
    
    fasta_cutter(length,prefix,delimiter,infile,outfile)
    
    
if __name__=='__main__':
    main()
    sys.exit()