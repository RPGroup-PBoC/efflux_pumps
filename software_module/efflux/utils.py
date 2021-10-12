import numpy as np
from Bio.Seq import Seq
from Bio.SeqIO import parse
import pandas as pd

import pathlib

def isint(var):
    """
    Check if variable is of type `int` or type `float`, but an integer.
    """
    
    if not (isinstance(var, int) or isinstance(var, np.int64)):
        if (isinstance(var, float) or isinstance(var, np.float64)):
            if not var.is_integer():
                return False
            else:
                return True
        else:
            return False
    else:
        return True
            
        
def seq2mat(seq,seq_dict):
    "Convert sequence to binary array."
    mat = np.zeros((len(seq_dict), len(seq)), dtype=int)
    for i,bp in enumerate(seq):
        mat[seq_dict[bp], i] = 1
    return mat


def choose_dict(dicttype, modeltype='MAT'):
    "Convert from bp to an index."
    if dicttype == 'dna':
        seq_dict = {'A':0,'C':1,'G':2,'T':3}
        inv_dict = {0:'A',1:'C',2:'G',3:'T'}
    elif dicttype == 'rna':
        seq_dict = {'A':0,'C':1,'G':2,'U':3}
        inv_dict = {0:'A',1:'C',2:'G',3:'U'}
    elif dicttype == 'protein':
        seq_dict = {
            '*':0,'A':1,'C':2,'D':3,'E':4,'F':5,'G':6,'H':7,'I':8,'K':9,'L':10,
            'M':11,'N':12,'P':13,'Q':14,'R':15,'S':16,'T':17,'V':18,'W':19,'Y':20}
        inv_dict = {v:k for k,v in seq_dict.items()}
    else:
        raise RuntimeError('Unkonwn dicttype: %s'%dicttype)

    if modeltype == 'NBR' or modeltype == 'PAIR':
        seq_dict = {
            ''.join([inv_dict[i],inv_dict[z]]):i*len(seq_dict)+z
            for i in range(len(seq_dict)) for z in range(len(seq_dict))}
        inv_dict = {seq_dict[i]:i for i in seq_dict.keys()}
    return seq_dict, inv_dict


def complement_seq(sequence, rev=False):
    com_dict = {
        "A": "T", 
        "T": "A", 
        "G": "C", 
        "C": "G", 
        "a": "t", 
        "t": "a",
        "c": "g",
        "g": "c"
    }
    rev_list = [com_dict[x] for x in sequence]
    if rev:
        com_sequence = "".join(reversed(rev_list))
    else:
        com_sequence = "".join(rev_list)
    return com_sequence


def import_primer_fwd(index):
   # Import primers       
    local_path = pathlib.Path(__file__).parent.absolute()
    primer_list = list(parse(str(local_path) + '/forward_finalprimers.fasta','fasta'))

    # Extract primer that is added
    primer = str(primer_list[index].seq)
    return primer
    


def import_primer_rev(index):
   # Import primers       
    local_path = pathlib.Path(__file__).parent.absolute()
    primer_list = list(parse(str(local_path) + '/reverse_finalprimers.fasta','fasta'))
    
    # Extract primers that is added
    primer = str(primer_list[index].seq)
    return primer

def _check_sequence_list(sequence_list):
    if type(sequence_list) not in [list, np.ndarray, pd.core.series.Series]:
        raise TypeError("sequence_list has to be list, numpy array or pandas series.")
    else:
        if any([type(seq) not in [str, Seq] for seq in sequence_list]):
            raise TypeError("entries in `sequence_list` have to be of type string or Bio.Seq.Seq.")
            
    for i,seq in enumerate(sequence_list):
        if type(seq) == str:
            sequence_list[i] = Seq(seq)
            
    return sequence_list
