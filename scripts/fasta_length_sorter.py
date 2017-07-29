#!/usr/bin/python
import sys

def help_screen():
    print """
----------------------------------------------------------------
FASTA LENGTH SORTER: Sort your sequences by length!
by D. Gonzales
               
This program takes a multi-fasta file and sorts them into separate
fasta files according to length.
       
Usage: fasta_length_sorter.py [-h] [-i in] [-l lower limit] [-u upper limit] [-o out]
[-h] = This helpful help screen.
[-i] = FASTA input file. 
[-l] = lower limit.
[-u] = upper limit.
[-o] = FASTA output prefix.

*Labels should not contain white spaces.
*Requirements: python 2.7, sys, getopt, biopython. 
----------------------------------------------------------------"""
    sys.exit()

	
def fasta_length_sorter(infile,lower,upper,outprefix):	
    from Bio import SeqIO
    infile=SeqIO.index(infile,'fasta')
    keys=infile.keys()
    outfile_lower=open(outprefix+'_'+str(lower)+'_lower.fasta','w')
    outfile_mid=open(outprefix+'_'+str(lower)+'-'+str(upper)+'_mid.fasta','w')
    outfile_upper=open(outprefix+'_'+str(upper)+'_upper.fasta','w')
    
    count_lower=0
    count_mid=0
    count_upper=0
    for key in keys:
        if count_lower!=0:
            outfile_lower.write('\n')
            count_lower=count_lower+1
        if count_mid!=0:
            outfile_mid.write('\n')
            count_mid=count_mid+1
        if count_upper!=0:
            outfile_upper.write('\n')
            count_upper=count_upper+1
        if len(infile[key].seq)<=int(lower):
            outfile_lower.write('>'+key+'\n'+str(infile[key].seq)+'\n')
    	elif len(infile[key].seq)>int(lower) and len(infile[key].seq)<int(upper):
            outfile_mid.write('>'+key+'\n'+str(infile[key].seq)+'\n')
        elif len(infile[key].seq)>=int(upper):
            outfile_upper.write('>'+key+'\n'+str(infile[key].seq)+'\n')

    outfile_lower.close()
    outfile_mid.close()
    outfile_upper.close()

def main():
    import getopt
    try:
        opts,args=getopt.getopt(sys.argv[1:],'i:l:u:o:h')
    except getopt.GetoptError as err:
        print str(err)
        print 'Usage: fasta_length_sorter.py [-h] [-i in] [-l lower limit] [-u upper limit] [-o out]'
        sys.exit()
    
    for o,a in opts:
        if o=='-h':
            help_screen()
        if o=='-i':
            infile=str(a)
        if o=='-l':
            lower=str(a)
        if o=='-u':
            upper=str(a)
        if o=='-o':
            outprefix=str(a)         
    
    options=[opt[0] for opt in opts]
    for option in ['-i','-l','-u','-o']:
        if option not in options:
            print 'missing arguement'
            print 'Usage: fasta_length_sorter.py [-h] [-i in] [-l lower limit] [-u upper limit] [-o out]'
            sys.exit()      
    
    fasta_length_sorter(infile,lower,upper,outprefix)

if __name__=='__main__':
    main()
    sys.exit()
