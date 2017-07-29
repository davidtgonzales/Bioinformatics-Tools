#!/usr/bin/python
import sys

def help_screen():
    print """
----------------------------------------------------------------
HMM EXTRACT: Extract data from hmmsearch results!
by D. Gonzales
               
This program takes the result of HMMSEARCH runs and summarizes 
it into a new csv file.
       
Usage: hmm_extract.py [-h] [-i in] [-o out]
[-h] = This helpful help screen.
[-i] = HMMSEARCH output file.
[-o] = CSV output file.

*Labels should not contain white spaces.
*Requirements: python 2.7, sys, getopt, biopython. 
----------------------------------------------------------------"""
    sys.exit()

def hmm_extract(infile,outfile):
    infile=open(infile,'r').readlines()
    outfile=open(outfile,'w')
    outfile.write('domain'+','+'query'+','+'score'+','+'cev'+','+'domain_seq'+','+'alignment'+','+'query_seq'+','+'start'+','+'end'+'\n')
    count=0
    for line in infile:
        if line[0:6]=='Query:':
            domain=line.split()[1]
        if line[0:2]=='>>':
            query=line.strip().strip('>>')
        if line[0:4]=='  ==':
            score=float(line.split()[4])
            cev=float(line.split()[8])
            domain_seq=infile[count+1].split()[2]
            alignment=infile[count+2].strip()
            query_seq=infile[count+3].split()[2]
            start=int(infile[count+3].split()[1])
            end=int(infile[count+3].split()[3])        
            outfile.write(domain+','+query+','+str(score)+','+str(cev)+','+domain_seq+','+alignment+','+query_seq+','+str(start)+','+str(end)+'\n')
        count=count+1

    outfile.close()


def main():
    import getopt
    try:
        opts,args=getopt.getopt(sys.argv[1:],'i:o:h')
    except getopt.GetoptError as err:
        print str(err)
        print 'Usage: hmm_extract.py [-h] [-i in] [-o out]'
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
            print 'Usage: hmm_extract.py [-h] [-i in] [-o out]'
            sys.exit()      
    
    hmm_extract(infile,outfile)

if __name__=='__main__':
    main()
    sys.exit()
