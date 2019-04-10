#!/usr/bin/env python
#-*-ecoding: utf-8-*-

'''
This script used for idba NGS sequences assemble.
And two fastq file or two tar.gz file are neccessary!

# Author:   Eric Lee
# Date:     2019-02-15
# Version:  0.2.3

Usage: python idbassemble-023.py or idbassemble-023.py

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
'''


import os
from glob import glob
from multiprocessing import cpu_count


print(__doc__)

def plog(log):
    'print logs with green color'
    return print(f'\33[32m {log} \33[0m')


def perror(error):
    'print error with red color'
    return print(f'\33[31m {error} \33[0m')


cpu_num = cpu_count() #get cpu numbers
if cpu_num < 4:
    perror(f'NOTICE: The CPU number is {cpu_num}, it is not recommanded to run following steps!')
    os._exit(0)
else:
    plog(f'The CPU number is {cpu_num}, and {cpu_num - 1} will be used in following steps.')


read_len = 0
while read_len not in (150, 250): #get and check read length parameters
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
    

tar_files = glob(r'*.tar.gz') # unzip tar files
if not tar_files:
    perror('NOTICE: zip file not exist! Maybe the unziped files has existed.') 
else:
    plog('STEP 3/6: Start unzip files')
    for tar_file in tar_files:
        os.system(f'tar -zxvf {tar_file}')


fq_files = glob(r'*.fq') # assign fastq file name
if not fq_files:
    perror('ERROR: fastq file not exist!') 
    os._exit(0)
#else:
#    R1 = fq_files[0]
#    R2 = fq_files[1]


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
    os.system(f'python seqspicker.py ./OutFolder/scaffold.fa {bait}.blastout')

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
