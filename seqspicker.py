#!/usr/bin/env python
# -*- ecoding: utf-8 -*-

import pandas as pd
import sys


'''
@Author:  Eric Lee
@date:    2018-04-24
@fileName:seqspicker.py
@version: 0.2

This script only used for simplly pick up sequence blast results 
    from IDBA_UD scaffold fasta file.

Usage: python seqspicker.py scaffold.fa seqblastout
'''


def parse_blastn(blastn6):
   '''parse blastn tabular results'''
   with open(blastn6) as fi:
       blast_out_list = pd.read_csv(fi, sep='\t',header = None)
       scaffold_name = list(blast_out_list[0])
       index_name = list(blast_out_list[1])
       index_dict = dict(zip(scaffold_name,index_name))
       
   return index_dict


def read_fasta(scaffold_fa):
   '''read fasta file as dict'''
   fasta = {}
   with open(scaffold_fa) as file_one:
       for line in file_one:
           line = line.strip()
           if not line:
               continue
           if line.startswith(">"):
               active_sequence_name = line[1:]
               if active_sequence_name not in fasta:
                  fasta[active_sequence_name] = []
               continue
           sequence = line
           fasta[active_sequence_name].append(sequence)

   return fasta


def get_hits(scaffold, seq_blast):
   '''get seq targets'''
   index = parse_blastn(seq_blast)
   fasta = read_fasta(scaffold)
   hit = {}
   with open('hits.fas','w') as hits:
       for key, value in index.items():
           newkey = key+'_'+value
           hit[newkey] = ''.join(fasta[key])

       for newkey,hit_value in hit.items():
           hits.write('>{}\n{}\n'.format(newkey,hit_value))


if __name__ == '__main__':

   get_hits(sys.argv[1],sys.argv[2])
   print('\nCompleted!\n')
