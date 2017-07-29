#!/usr/bin/python
import sys

def help_screen():
    print """
----------------------------------------------------------------
DNA CODE: Extract the DNA code from the protein code!
by D. Gonzales
               
This program takes a protein fasta file and then extracts the 
corresponding DNA sequences in its DNA fasta file. The input
protein fasta should have labels that specify the frame (1,2,3,
-1,-2,-3) of the sequence (i.e. >contig_1,1)
       
Usage: DNA_code.py [-h] [-p protein in] [-d dna in] [-o out]
[-h] = This helpful help screen.
[-p] = Protein FASTA input file. 
[-d] = DNA FASTA input file. 
[-o] = DNA output file.

*Labels should not contain white spaces.
*Requirements: python 2.7, sys, getopt, biopython. 
----------------------------------------------------------------"""
    sys.exit()
    

def dna_code(inprotein,indna,outfile):
    from Bio import SeqIO
    inprotein=SeqIO.index(inprotein,'fasta')
    indna=SeqIO.index(indna,'fasta')
    outfile=open(outfile,'w')        
    inprotein_keys=inprotein.keys()

    for key in inprotein_keys:
        frame=int(key.split(',')[1].strip())
        if frame>0:
            protein=str(indna[key.split(',')[0]].seq[(frame-1):].translate())
        elif frame<0:
            protein=str(indna[key.split(',')[0]].seq.reverse_complement()[(-1*frame)-1:].translate())
      
        for index in range(len(protein)):
            if protein[index:index+len(str(inprotein[key].seq))]==str(inprotein[key].seq) and frame>0:
                outfile.write('>'+str(key.split(',')[0])+','+str(frame)+'\n'+str(indna[key.split(',')[0]].seq[(index*3)+frame-1:((index+len(str(inprotein[key].seq)))*3)+frame-1]))
                break
            elif protein[index:index+len(str(inprotein[key].seq))]==str(inprotein[key].seq) and frame<0:
                outfile.write('>'+str(key.split(',')[0])+','+str(frame)+'\n'+str(indna[key.split(',')[0]].seq.reverse_complement()[(index*3)+(-1*frame)-1:((index+len(str(inprotein[key].seq)))*3)+(-1*frame)-1]))
                break
        
        if key!=inprotein_keys[len(inprotein_keys)-1]:
            outfile.write('\n')

    outfile.close()

def main():
    import getopt
    try:
        opts,args=getopt.getopt(sys.argv[1:],'p:d:o:h')
    except getopt.GetoptError as err:
        print str(err)
        print 'Usage: DNA_code.py [-h] [-p protein in] [-d dna in] [-o out]'
        sys.exit()
    
    for o,a in opts:
        if o=='-h':
            help_screen()
        if o=='-p':
            inprotein=str(a)
        if o=='-d':
            indna=str(a)
        if o=='-o':
            outfile=str(a)         
    
    options=[opt[0] for opt in opts]
    for option in ['-p','-d','-o']:
        if option not in options:
            print 'missing arguement'
            print 'Usage: DNA_code.py [-h] [-p protein in] [-d dna in] [-o out]'
            sys.exit()      
    
    dna_code(inprotein,indna,outfile)

if __name__=='__main__':
    main()
    sys.exit()
