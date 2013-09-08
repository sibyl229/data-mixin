import numpy as np
import os
import re
import json
from Bio import SeqIO
from settings import *

class MySequence(object):

    
    def __init__(self, species):
        seqName = species.context['seqName']
        seqPath = os.path.join(RAW_INPUT_PATH, 
                                    species.rel_path(seqName))
        self.record_dict = SeqIO.index(seqPath, "fasta")

    def get_all_prot(self):
        '''get ids from fasta file, Order is NOT guaranteed'''
        keys = self.record_dict.keys()
        return [(key, self.extract_pid(key)) for key in keys]

    def get_seq(self, key):
        return self.record_dict.get(key, None)

    @staticmethod
    def extract_pid(key):
        match = re.search(r'[A-Z0-9]{5,7}', key)
        if not match:
            print key
        pid = match.group()
        return pid

