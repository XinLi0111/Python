#!/usr/bin/env python
# -*-encoding:utf-8-*-

'''
This script used for idba NGS sequences assemble.
And two fastq file or two tar.gz file are neccessary!

# Author:   Eric Lee
# Date:     2019-05-28
# Version:  0.2.5

Usage: python idbassemble5.py or idbassemble5.py

baits file must have extention of <*.fa>

'''

'''

History:
0.1   create
0.1.1 add unzip function
      add choice of read length (150/250)
0.2   rewrite with def funtion, print color symbol notice
0.2.1 add prcessing state
0.2.2 add judge function to process the baits input filenames
0.2.3 fix the fq2fa parameters mistakes 'R1 R2 is not global variations
0.2.4 support extract more type compressed files
0.2.5 integrated seqspicker function in one script
'''


import os
import pandas as pd
from glob import glob
from multiprocessing import cpu_count


print(__doc__)


def plog(log):
    'print logs with green color'
    return print(f'\33[32m {log} \33[0m')


def perror(error):
    'print error with red color'
    return print(f'\33[31m {error} \33[0m')


######################
# seqspicker
def parse_blastn(blastn6):
    '''parse blastn tabular results'''
    with open(blastn6) as fi:
        blast_out_list = pd.read_csv(fi, sep='\t', header=None)
        scaffold_name = list(blast_out_list[0])
        index_name = list(blast_out_list[1])
        index_dict = dict(zip(scaffold_name, index_name))
        
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
    with open('hits.fas', 'w') as hits:
        for key, value in index.items():
            newkey = key+'_'+value
            hit[newkey] = ''.join(fasta[key])

        for newkey, hit_value in hit.items():
            hits.write('>{}\n{}\n'.format(newkey, hit_value))


##########################
# idba_ud parameters
cpu_num = cpu_count()   # get cpu numbers
if cpu_num < 4:
    perror(f'NOTICE: The CPU number is {cpu_num}, it is not recommanded to run following steps!')
    os._exit(0)
else:
    plog(f'The CPU number is {cpu_num}, and {cpu_num - 1} will be used in following steps.')


read_len = 0
while read_len not in (150, 250):    # get and check read length parameters
    read_len = int(input('PLEASE input proper sequence read length: [150 or 250] '))
else:
    plog(f'\n STEP 1/6: Read length has been assiged')


kmer_parameters = ''  #set kmers  
if read_len > 151:
    kmer_parameters = '--mink 80 --maxk 240'
    plog(f'STEP 2/6: Set kmer "{kmer_parameters}"')
else:
    kmer_parameters = '--mink 40 --maxk 140'
    plog(f'STEP 2/6: Set kmer "{kmer_parameters}"')
    

# unzip tar files
bz2_files = glob(r'*.tar.bz2')
tar_files = glob(r'*.tar.gz')
gz_files = glob(r'*.gz') 
for tarfile in tar_files:
    gz_files.remove(tarfile)


if bz2_files:
    plog('STEP 3/6: Start decompress files')
    for bz2_file in bz2_files:
        os.system(f'tar -jxvf {bz2_file}')
elif tar_files:
    plog('STEP 3/6: Start decompress files')
    for tar_file in tar_files:
        os.system(f'tar -zxvf {tar_file}')
elif gz_files:
    plog('STEP 3/6: Start decompress files')
    for gz_file in gz_files:
        os.system(f'gzip -d {gz_file}')
else:
    perror('NOTICE: Compressed files not exist! Maybe the decompressed files has existed.') 


fq_files = glob(r'*.fq')
fq_files.extend(glob(r'*.fastq'))    # assign fastq file name
#print(fq_files)

if not fq_files:
    perror('ERROR: fastq file not exist!') 
    os._exit(0)
#else:
#    R1 = fq_files[0]
#    R2 = fq_files[1]


###########################
# start idba_ud processes
plog('STEP 4/6: Start merge fastq files')
os.system(f'fq2fa --merge {fq_files[0]} {fq_files[1]} merged.fsa --filter')


plog('STEP 5/6: Start assemble, this process will take a lot of time')
os.system(f'idba_ud -r merged.fsa {kmer_parameters} --step 20 -o OutFolder --num_threads {cpu_num - 1} --min_contig 500 --similar 0.98')


plog('STEP 6/6: Start make blast index database')
baits = glob(r'*.fa')
while not baits:
    perror('ERROR: baits files not exist! Please provide your baits files')
    baits = list(input('please input your baits file names:').split(' '))

for bait in baits:
    print(f'\n\tStart make \33[32m {bait} \33[0m index blast database')
    os.system(f'makeblastdb -in {bait} -dbtype nucl -parse_seqids -out {bait}')
    print(f'\tStart blast using \33[32m {bait} \33[0m bait')
    os.system(f'blastn -query ./OutFolder/scaffold.fa -out {bait}.blastout -db {bait} -outfmt 6 -evalue 1e-5 -max_target_seqs 1 -num_threads {cpu_num - 1}')
        
    plog(f'#\n# STEP: Start pick up target sequences by {bait}\n#')
    blastoutfilename = f'{bait}.blastout'
    get_hits('./OutFolder/scaffold.fa', blastoutfilename)

print('''

　　＿＿＿＿＿＿＿
　 ｜Complete!｜
　　￣￣￣∨￣￣￣
　　　　∧＿∧　
　　　 (´・ω⊂ヽ゛
　　　 /　　 _ノ⌒⌒ヽ
、(￣⊂人　 ノシ⌒　　ノ
⊂ニニニニニニニニニニニ⊃

''')
