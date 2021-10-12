import sys

import numpy as np
import pandas as pd

import copy

from Bio.Restriction import *
from Bio.Seq import Seq

from .seq_utils import _check_sequence_list
from .utils import isint,import_primer_fwd,import_primer_rev,_check_sequence_list,complement_seq

import warnings

def mutation_coverage(wildtype_seq, mutants_list, site_start=0, site_end=None):
    
    if not isint(site_start):
        raise TypeError("`site_start` is of type {} but has to be integer valued.".format(type(site_start)))
        
    if not (isint(site_end) or site_end == None):
        raise TypeError("`site_end` is of type {} but has to be integer valued.".format(type(site_end)))
    
    if site_end==None:
        wildtype_seq = wildtype_seq[site_start:]
        mutations = np.zeros([len(mutants_list), len(wildtype_seq)])
        wildtype_list = list(wildtype_seq)
        
        for i, sequence in enumerate(mutants_list):
            mutations[i, :] = [x != y for (x, y) in zip(list(sequence[site_start:]), wildtype_list)]
    else:
        wildtype_seq = wildtype_seq[site_start:site_end]
        mutations = np.zeros([len(mutants_list), len(wildtype_seq)])
        wildtype_list = list(wildtype_seq)
        for i, sequence in enumerate(mutants_list):
            mutations[i, :] = [x != y for (x, y) in zip(list(sequence[site_start:site_end]), wildtype_list)]
           
    return np.mean(mutations, axis=0)


def scan_enzymes(sequence_list, enzymes=[], ret_list=False):
    """Compute number if restriction sites in a list of sequences for list of enzymes.
    
    Parameters
    ----------
    sequence_list : array-like
    
    enzymes : array-like
        List of enzymes to check. If none is given, then all enzymes are tested. Not Recommended.
    ret_list : boolean, default `False`
        If `True`, returns list of enzymes.

    
    """
    
    # Check inputs
    if type(enzymes) not in [list, np.ndarray, pd.core.series.Series]:
        raise TypeError("enzymes has to be list, numpy array or pandas series.")
    else:
        if any([(not type(enz) == str) and (not type(enz).__bases__[-1] == Restriction.RestrictionType) for enz in enzymes]):
            raise TypeError("entries in `enzymes` have to be of type string or Bio.RestrictionType")
    
    # Check inputs
    if type(sequence_list) == pd.core.series.Series:
        sequence_list = copy.deepcopy(sequence_list.values)
    sequence_list = _check_sequence_list(sequence_list)
    
    
    # Choose all commercially available enzymes if none given
    if len(enzymes) == 0:
        enzymes = CommOnly
        ret_list = True
        
    num_sites = np.zeros(len(enzymes))
    for i, enz in enumerate(enzymes):
        sites = find_restriction_sites(enz, sequence_list)
        num_sites[i] = len([ele for sub in sites for ele in sub])
    
    if ret_list:
        return num_sites, enzymes
    else:
        return num_sites


def scan_enzymes_print(sequence_list, enzymes=[]):
    counts, enzs = scan_enzymes(sequence_list, enzymes, ret_list=True)
    ind_sorted = np.flip(np.argsort(counts))
    message = ""
    for count, enzyme in zip(counts[ind_sorted], np.array([enz for enz in enzs])[ind_sorted]):
        message += str(enzyme)+ ": "+str(count)+'\n'
    return message


def find_restriction_sites(enzyme, sequence_list):
    """Searches for restriction sites of a specific enzyme in a list of sequences
    
    Parameters
    ----------
    enzyme : Bio.Restriction or string
        Name of the enzyme.
    sequence_list: array-type
        
    """
    
    if type(enzyme) == str:
        try:
            getattr(sys.modules[__name__], enzyme)
        except AttributeError:
            raise ValueError("There is no restriction enzyme {}. Check spelling, naming is case-sensitive.".format(enzyme))
        enzyme = getattr(sys.modules[__name__], enzyme)
        
    # Enzyme classes in Biopython are weird but this works
    elif not (type(enzyme).__bases__[-1] == Restriction.RestrictionType): 
        raise TypeError("enzyme has to be either string of Bio.Restriction.Restriction.RestrictionType")
    
    # Transform sequence inputs
    sequence_list = _check_sequence_list(sequence_list)
   
    return [enzyme.search(sequence) for sequence in sequence_list]


def digital_PCR(sequence_list, fwd_primer, rev_primer):
    """Perfom a digital PCR on a sequence pool for a primer pair.

    Parameters
    ----------
    sequence_list : array-like
        List of sequences in the subpool. Sequences can be either strings or of Bio.Seq type.
    fwd_primer : string or int
        Forward primer for PCR. Can either be the sequence, or the index of the primer in the Kosuri list.
    rev_primer : string
        Reverse primer for PCR. Can either be the sequence, or the index of the primer in the Kosuri list.
    """
    if type(sequence_list) not in [list, np.ndarray, pd.core.series.Series]:
        raise TypeError("sequence_list has to be list, numpy array or pandas series.")
    
    if not (isint(fwd_primer) or type(fwd_primer) == str):
        raise TypeError("Primers have to be either the sequence (string), or the index (integer) of the primer in the Kosuri list")
    if not (isint(rev_primer) or type(rev_primer) == str):
        raise TypeError("Primers have to be either the sequence (string), or the index (integer) of the primer in the Kosuri list")
    
    # Import primers if given as indices
    if isint(fwd_primer):
        fwd_primer = import_primer_fwd(fwd_primer)

    if isint(rev_primer):
        rev_primer = import_primer_fwd(rev_primer)

    sequence_list = _check_sequence_list(sequence_list)

    forward_sites = [seq.find(fwd_primer) for seq in sequence_list]
    reverse_sites = [seq.find(rev_primer) for seq in sequence_list]

    amplified_seqs = [str(seq[i:j+20]) if (i != -1 and j != -1) else "x" for ((ind, seq), i,j) in zip(enumerate(sequence_list), forward_sites, reverse_sites) ]


    return amplified_seqs


def check_primers_pool_df(df_pool):
    """
    Check primer positions in DataFrame. Has to contain columns `seq` with sequences,
    and columns containing `forward_primer` and `reverse_primer`. Multiple of these columns are possible.
    Primer columns have to consist of tuples with index of primer in Kosuri list and position in the sequence.

    Parameters
    ----------
    df_pool : DataFrame
        DataFrame containing sequences as well as primer indices and positions.

    Returns
    -------
    no_warnings : boolean
        `True` if every sequence passes the test, `False` otherwise.
    """

    fwd_primer_columns = [x  for x in df_pool.columns if ("forward_primer" in x)]
    rev_primer_columns = [x  for x in df_pool.columns if ("reverse_primer" in x)]

    no_warnings = True
    for index, row in df_pool.iterrows():
        for fwd_column in fwd_primer_columns:
            ind, pos = row[fwd_column]
            if not ((ind == None) or (pos == None)):
                primer = import_primer_fwd(ind)
                primer_pos = Seq(row['seq'].upper()).find(primer)
                if primer_pos == -1:
                    warnings.resetwarnings()
                    warnings.warn("Primer {} is not in sequence {} as expected.".format(ind, index))
                    no_warnings = False
                elif pos != primer_pos:
                    warnings.resetwarnings()
                    warnings.warn("Primer {} is not at the right position in sequence {}. Found at position {}, supposed to be at {}.".format(ind, index, primer_pos, pos))
                    no_warnings = False

        for rev_column in rev_primer_columns:
            ind, pos = row[rev_column]
            if not ((ind == None) or (pos == None)):
                primer = complement_seq(import_primer_rev(ind), rev=True)
                primer_pos = Seq(row['seq'].upper()).find(primer)
                if primer_pos == -1:
                    warnings.resetwarnings()
                    warnings.warn("Primer {} is not in sequence {} as expected.".format(ind, index))
                    no_warnings = False
                elif pos != primer_pos:
                    warnings.resetwarnings()
                    warnings.warn("Primer {} is not at the right position in sequence {}. Found at position {}, supposed to be at {}.".format(ind, index, primer_pos, pos))
                    no_warnings = False

    return no_warnings



def check_primer_seq(seq_list, primer, position):
    """
    Check if primer is in every given sequence and at the right position.

    Parameters
    ----------
    seq_list : list 
        List of sequences
    primer : string
    position : position of the primer in the sequence

    Returns
    -------
    no_warnings : boolean
        `True` if every sequence passes the test, `False` otherwise.
    """
    no_warnings = True
    for i, seq in enumerate(seq_list):
        primer_pos = Seq(seq).find(primer)
        if primer_pos == -1:
            warnings.resetwarnings()
            warnings.warn("Primer is not in sequence {} as expected.".format(i))
            no_warnings = False
        elif position != primer_pos:
            warnings.resetwarnings()
            warnings.warn("Primer is not at the right position in sequence {}. Found at position {}, supposed to be at {}.".format(i, primer_pos, position))
            no_warnings = False
        return no_warnings