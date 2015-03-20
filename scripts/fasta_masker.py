#!/usr/bin/python
import sys

def help_screen():
    print """
----------------------------------------------------------------
FASTA MASKER: Mask specified regions in your fasta file!
by D. Gonzales
               
This program takes a fasta file and masks given regions by
replacing them with 'N's. User should provide a mask
coordinate csv file with the format shown in the example
below [sequence_label*,start**,end]:

    seq1,0,100
    seq2,200,250
    seq3,150,300
       
Usage: fasta_masker.py [-h] [-m str] [-i in] [-o out]
[-h] = This helpful help screen.
[-m] = Mask coordinates file. 
[-i] = FASTA input file. 
[-o] = FASTA output file.

*Labels should not contain white spaces or commas.
**Index starts at 1 (not 0).
Requirements: python 2.7, sys, getopt, biopython.
----------------------------------------------------------------"""
    sys.exit()

def fasta_masker(masks,infile,outfile):
    mask_id=[]
    mask_coord=[]
    masks=open(masks,'r')
    for mask in masks:
        mask_id.append(mask.split(',')[0])
        mask_coord.append([mask.split(',')[1],mask.split(',')[2]])

    from Bio import SeqIO
    infile_parse=SeqIO.parse(infile,'fasta')
    outfile=open(outfile,'w')
    start=0
    for sequence in infile_parse:
        if start==0:
            label='>'
            start=1
        else:
            label='\n>'
        
        working_id=str(sequence.id)
        working_seq=str(sequence.seq)
        
        while working_id in mask_id:
            if int(mask_coord[mask_id.index(working_id)][0])<int(mask_coord[mask_id.index(working_id)][1]):
                start=int(mask_coord[mask_id.index(working_id)][0])
                end=int(mask_coord[mask_id.index(working_id)][1])
            else:
                end=int(mask_coord[mask_id.index(working_id)][0])
                start=int(mask_coord[mask_id.index(working_id)][1])   
                
            working_seq=working_seq[0:start-1]+''.join(['N' for n in range(start,end+1)])+working_seq[end:]
            mask_coord.remove(mask_coord[mask_id.index(working_id)])
            mask_id.remove(working_id)

        outfile.write(label+working_id)
        outfile.write('\n'+working_seq)            
            
    outfile.close()
    
def main():
    import getopt
    try:
        opts,args=getopt.getopt(sys.argv[1:],'m:i:o:h')
    except getopt.GetoptError as err:
        print str(err)
        print 'Usage: fasta_masker.py [-h] [-m str] [-i in] [-o out]'
        sys.exit()
    
    for o,a in opts:
        if o=='-h':
            help_screen()
        if o=='-m':
            masks=str(a)
        if o=='-i':
            infile=str(a)
        if o=='-o':
            outfile=str(a)         
    
    options=[opt[0] for opt in opts]
    for option in ['-m','-i','-o']:
        if option not in options:
            print 'missing arguement'
            print 'Usage: fasta_masker.py [-h] [-m str] [-i in] [-o out]'
            sys.exit()      
    
    fasta_masker(masks,infile,outfile)

if __name__=='__main__':
    main()
    sys.exit()