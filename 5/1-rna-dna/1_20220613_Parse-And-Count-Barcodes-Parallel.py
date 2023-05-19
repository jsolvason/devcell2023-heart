#!/usr/bin/env python
# coding: utf-8

# # Parse and count barcodes

# In[1]:


import pickle
import sys


# Sample dir
in_sample_dir='/home/solvason/scratch/ets-heart-mpra-rna-dna/'

in_sample=sys.argv[1]

# Get sample list
# in_sample_list=!ls {in_sample_dir}
# in_sample_list=[i for i in in_sample_list if i!='split-fastq.ipynb']
# in_sample_list

# Modify sample list
# in_sample_list= [
#     'head_109_DNA_S160_L001_R1_001.fastq'
# ]

out_path_parent='/home/solvason/projects/ets/heart/1-heart-mpra/1-rna-dna/data-processed/1-barcodes/'

done=[]

from Bio import SeqIO

# in_sample=in_sample_list[0]    

out_path=f'{out_path_parent}/'

# Real data
# in_fqgz_list=!ls {in_sample_dir}/{in_sample}

#     # Test
#     in_fqgz_list=!ls /home/solvason/projects/otxa/otxa-scrambled-constant/rna-dna/raw-data/subsample/split/{in_sample}/*.fastq

# print(f'Analyzing {in_sample}...')
# for in_fqgz in in_fqgz_list:

in_basename=in_sample.split('.fastq')[0]
out_fn=f'{out_path}/{in_basename}.Bc2ReadCount.pickle'

##################################################
# Read in barcodes and tabulate readcounts
##################################################

print(f'Reading {in_sample}...')

total_reads=0
Bc2ReadCount={}
with open(in_sample_dir+in_sample,'r') as handle:
    for record in SeqIO.parse(handle, "fastq"):
        total_reads+=1

        seq=str(record.seq)

        bc=seq[:25]
        if bc not in Bc2ReadCount:
            Bc2ReadCount[bc]=0
        Bc2ReadCount[bc]+=1

print(f'\tWriting {in_sample}...')
with open(out_fn, 'wb') as handle: pickle.dump(Bc2ReadCount, handle, protocol=pickle.HIGHEST_PROTOCOL)
print(f'\tTotal reads: {total_reads}')





# In[ ]:




