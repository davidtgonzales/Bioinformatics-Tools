#!/usr/bin/python
import sys

def help_screen():
    print """
----------------------------------------------------------------
FASTA COVERAGE SORTER: Extract and divide sequences according
to coverage values!
by D. Gonzales
               
This program takes a multi-fasta file, a list of read IDs with
their corresponding coverage values, and a list of bin coverage
ranges and then sorts the sequences into separate fasta files 
according to a set coverage bins.

coverage csv format:
contig,size,cov
seq1,100,30,300
seq2,50,50,2500

bin sizes csv format:
bin,min,max
bin1,0,100
bin2,100,200
bin3,200,300
       
Usage: fasta_coverage_sorter.py [-h] [-f fasta] [-v cov] [-b bin] [-o out]
[-h] = This helpful help screen.
[-f] = FASTA input file.
[-v] = csv of coverage values.
[-b] = csv of bin sizes of coverage values.
[-o] = FASTA output file base name.

*Labels should not contain white spaces.
*Requirements: python 2.7, sys, getopt, biopython. 
----------------------------------------------------------------"""
    sys.exit()

def fasta_coverage_sorter(infasta,incov,inbin,outfile):
    from Bio import SeqIO    

    infasta=SeqIO.index(infasta,'fasta')   
    incov=open(incov,'r').readlines()[1:]
    inbin=open(inbin,'r').readlines()[1:]

    coverage=[]    
    for item in incov:
        coverage.append(item.strip().split(','))

    binning=[]
    for item in inbin:
        temp=item.strip().split(',')
        temp.append([])
        binning.append(temp)
    
    for item in coverage:
        for ranges in binning:
            if float(item[2]) >= float(ranges[1]) and float(item[2]) < float(ranges[2]):
                ranges[3].append(item[0])
                break             

    for item in binning:
        outbin=outfile+'_'+item[0]
        records=(infasta[id] for id in item[3])
        SeqIO.write(records,outbin,'fasta')

def main():
    import getopt
    try:
        opts,args=getopt.getopt(sys.argv[1:],'f:v:b:o:h')
    except getopt.GetoptError as err:
        print str(err)
        print 'Usage: fasta_coverage_sorter.py [-h] [-f fasta] [-v cov] [-b bin] [-o out]'
        sys.exit()
    
    for o,a in opts:
        if o=='-h':
            help_screen()
        if o=='-f':
            infasta=str(a)
        if o=='-v':
            incov=str(a)
        if o=='-b':
            inbin=str(a)
        if o=='-o':
            outfile=str(a)         
    
    options=[opt[0] for opt in opts]
    for option in ['-f','-v','-b','-o']:
        if option not in options:
            print 'missing arguement'
            print 'Usage: fasta_coverage_sorter.py [-h] [-f fasta] [-v cov] [-b bin] [-o out]'
            sys.exit()      
    
    fasta_coverage_sorter(infasta,incov,inbin,outfile)

if __name__=='__main__':
    main()
    sys.exit()

