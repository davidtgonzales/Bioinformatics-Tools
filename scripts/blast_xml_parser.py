#!/usr/bin/python
import sys

def help_screen():
    print """
----------------------------------------------------------------
BLAST XML PARSER: Convert blast xml formats into tables!
by D. Gonzales
               
This program takes a blast xml output format and converts it
into a csv file containing the f.f. information by column:
    
1 query_id
2 query_length
3 frame
4 start
5 end
6 hit_length
7 length_align
8 percent_identity
9 e-value
10 hit_id*
11 hit_def*
12 query_sequence
13 hit_sequence

*Replaces all commas with a dash to avoid csv errors.

Usage: blast_xml_parser.py [-h] [-i in] [-o out]
[-h] = This helpful help screen.
[-i] = xml input file. 
[-o] = csv output file.

*Requirements: python 2.7, sys, getopt, xml. 
----------------------------------------------------------------"""
    sys.exit()

def blast_xml_parser(infile,outfile):
    import xml.etree.cElementTree as ET
    
    infile=ET.parse(infile)
    infile=infile.getroot()
    
    for iter in infile[8].findall('Iteration'):    
        iterhit=iter[4].text
        if iterhit==None:
            infile[8].remove(iter)
       
    outfile=open(outfile,'w')
    outfile.write('query_id'+','+'query_length'+','+'frame'+','+'start'+','+'end'+','+'hit_length'+','+'length_align'+','+'percent_identity'+','+'e-value'+','+'hit_id'+','+'hit_def'+','+'query_sequence'+','+'hit_sequence'+'\n')

    for iter in infile[8].findall('Iteration'):
        for hit in iter[4].findall('Hit'):
            hit_id=hit[3].text.replace(',','-') #hit_id
            hit_def=hit[2].text.replace(',','-') #hit_def
            outfile.write(iter[2].text.replace(',','-')+','+iter[3].text+','+hit[5][0][8].text+','+hit[5][0][4].text+','+hit[5][0][5].text+','+hit[4].text+','+hit[5][0][13].text+','+str(float(hit[5][0][10].text)/float(hit[5][0][13].text)*100)+','+hit[5][0][3].text+','+hit_id+','+hit_def+','+hit[5][0][14].text+','+hit[5][0][15].text+'\n')

    outfile.close()

def main():
    import getopt
    try:
        opts,args=getopt.getopt(sys.argv[1:],'i:o:h')
    except getopt.GetoptError as err:
        print str(err)
        print 'Usage: blast_xml_parser.py [-h] [-i in] [-o out]'
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
            print 'Usage: blast_xml_parser.py [-h] [-i in] [-o out]'
            sys.exit() 
    
    blast_xml_parser(infile,outfile)
            
if __name__=='__main__':
    main()
    sys.exit()
