import numpy as np
import os
import json
from Bio import SeqIO
from settings import *

# human proteome of reviewed and canonical sequences
fastaPath = os.path.join(RAW_INPUT_PATH, 
    'uniprot-taxonomy%3A9606+AND+reviewed%3Ayes+AND+keyword%3A181.fasta')
record_dict = SeqIO.index(fastaPath, "fasta")

def get_all_prot():
    '''get ids from fasta file, Order is NOT guaranteed'''
    keys = record_dict.keys()
    return [(key, extract_pid(key)) for key in keys]

def get_seq(key):
    return record_dict.get(key, None)

def extract_pid(key):
    db, pid, gene_id = key.split('|')
    return pid

