#!/usr/bin/env python3

import argparse
import os
import sys
import random
from Bio import SeqIO

def get_args():
    parser = argparse.ArgumentParser(description='Fasta random subsampler')
    parser.add_argument('file_list' , metavar ='FASTA_file_list', help='FASTA file list',
                        type=str, default='foo')
    parser.add_argument('-p', '--percentage', help='percentage of original reads to subsample',
                        metavar='int', type=int, default=75)
    parser.add_argument('-r', '--reps', help='repetitions for each subsampling percentage',
                        metavar='int', type=int, default=1)
    
    return parser.parse_args()

def main():
    args = get_args()
    list  = args.file_list
    percentage =args.percentage
    reps = args.reps
    
    per=percentage/100
    
    with open(list) as filelist:
        for line in filelist:
            file = str(line.rstrip())
            i = 1
            
    
            while i <= reps:
                outfile = file+'.subset'+str(i)
                i=i+1
                out_fh = open(outfile, 'w')
                             
                print (outfile)
                
                exec(open('sample_fastq.py'.read())
                
#                records = sum(1 for _ in open (file))/4
 #               N = records*per
  #              random_records = sorted([random.randint(0, records -1) for _ in xrange(N)])
   #             
    #            with open (file) as fh:
     #               rec_no = -1
      #              for rr in rand_records:
       #                 
        #                while rec_no <rr:
         #                   rec_no += 1
          #                  for i in range(4):
           #                     fh.readline()                
            #            for i in range(4):
             #               out_fh.write(fh.readline())
              #          rec_no += 1
               #     print ("wrote {} of the original {} records to {}.".format(N, records, outfile)     
          

if __name__ == '__main__':
    main()