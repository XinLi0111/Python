#!/usr/bin/python
#-*- coding: utf-8 -*-


'''
date: 2018-11-25
version: 0.1
history:
    0.1    create script
'''


import argparse
from Bio import SeqIO


parser = argparse.ArgumentParser(description = 'used for myself extract genbank features',
                                usage = 'python extract_features.py -g <sequence.gb> -m <cds>')
parser.add_argument('-g', '--genbank',
                    default = 'sequence.gb',
                    help = '<sequence.gb> file in genbank format as input')
parser.add_argument('-f', '--feature',
                    default = 'misc_feature',
                    required = True,
                    help = '<feature> specify the feature that will be extracted')
parser.add_argument('-o', '--outfile',
                    default = 'extract_results.fa',
                    help = '<extract_results.fa> file in fasta format as output')
parser.add_argument('-k', '--keywords',
                    action = 'store_true',
                    help = 'decide whether the keyword (RefSeq) is displayed')

args = parser.parse_args()


def organism_name(seq_record):
    'get organism name from seq_record and remove space in scientific name'
    name_string = seq_record.annotations["organism"]
    new_name = "_".join(name_string.split())
    return new_name


def keywords(seq_record):
    'get keyword frome seq_record class and transfer list into str'
    if args.keywords:
        keyword = "".join(seq_record.annotations["keywords"])
    else:
        keyword = ""
    return " " + keyword


def subseq(seq_record, seq_feature):
    'extract sub sequence from seq_record and features'
    sequence = seq_feature.location.extract(seq_record.seq)
    return sequence


for seq_record in SeqIO.parse(args.genbank, 'genbank'):
    #print(f'Dealing with GenBank record \n>{seq_record.id}')
    for seq_feature in seq_record.features:
        if seq_feature.type == args.feature:
            #with open(args.outfile, 'w') as fo:
            print(f'>{seq_record.name}|{organism_name(seq_record)}{keywords(seq_record)} {seq_feature.qualifiers["note"]}\
            \n{subseq(seq_record, seq_feature)}')
