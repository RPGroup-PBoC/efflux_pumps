import sys
from contextlib import suppress
import numpy as np
import pandas as pd
import copy
import random
import numba
import itertools

from scipy.special import factorial

import pathlib

from Bio.Restriction import *
from Bio.Seq import Seq
from Bio.SeqIO import parse

from .utils import isint, choose_dict, import_primer_fwd, import_primer_rev,_check_sequence_list,complement_seq
from .quality_control import scan_enzymes
import warnings


def find_seq(TSS, strand, up, dn, genome):
    '''
    Find sequence upstream and downstream of TSS. If gene is reverse 
    transcribed, return reverse complement of sequence.
    
    Parameters
    ----------
    TSS : int
        Transcription start site. Refers to position in the genome.
    strand : string
        Either '+' for fwd strand or '-' for reverse strand.
    up : int
        Number of basepairs upstream of TSS to grab. Relative to TSS according to strand / direction.
    dn : int
        Number of basepairs downstream of TSS to grab. Relative to TSS according to strand / direction.
    genome : string
        Reference genome as a string 
        
    Returns
    -------
    output : string
        Genetic sequence upstream and downstream of TSS
    ''' 
        
    if strand == '-':
        
        gene = genome[TSS-dn:TSS+up]
        
        tempgene = Seq(gene)
        outgene = str(tempgene.reverse_complement())
        
        left_pos = TSS-dn
        right_pos = TSS+up
        
    elif strand == '+':
        outgene = genome[TSS-up:TSS+dn]
        
        left_pos = TSS-up
        right_pos = TSS+dn
        
    return (outgene, left_pos, right_pos)


def find_seq_df(row, TSS_col, strand_col, up, dn, genome):
    '''
    Find sequence upstream and downstream of TSS. If gene is reverse 
    transcribed, return reverse complement of sequence. Modified to be applied on a dataframe of many TSS and associated strands.
    
    Parameters
    ----------
    row : pandas dataframe row
        This row in a pandas dataframe will have columns TSS_col and strand_col. see find_seq()
    TSS_col : string
        The column in the dataframe that contains the TSS position. TSS values are ints.
    strand_col : string
        The column in the dataframe that contains the strand direction. Strand values are either '+' for fwd strand or '-' for reverse strand.
    up : int
        Number of basepairs upstream of TSS to grab. Relative to TSS according to strand / direction.
    dn : int
        Number of basepairs downstream of TSS to grab. Relative to TSS according to strand / direction.
    genome : string
        Reference genome as a string 
        
    Returns
    -------
    output : string
        Genetic sequence upstream and downstream of TSS
        
    example : df[['sequence','left_pos','right_pos']] = df.apply(find_seq_df, axis = 1, 
                                                        args = ('TSS_position','direction',115,45, genome), result_type = 'expand')
    '''
    TSS = row[TSS_col]
    strand = row[strand_col]
 
    return find_seq(TSS, strand, up, dn, genome)


def gen_rand_seq(length, letters=['A','T','C','G']):
    """
    Generate a random DNA sequence of defined length
    
    Parameters
    ----------
    length : int
        Length of generated sequence
    letters : list
        List of letters to choose from when generating the sequence.
    Returns
    ---------
    seq : string
         Random sequence
    """
    # Confirm argument types
    if type(letters) not in [list, np.ndarray, pd.core.series.Series]:
        raise TypeError("`letters` has to be a list-type.")
        
    if not isint(length):
        raise TypeError("`length` has to be an integer.")
    elif length < 0:
        raise ValueError("`length` has to be non-negative.")
    else:
        length = int(length)
            
    # Generate random sequence
    seq = ''.join(random.choice(letters) for i in range(length))
    return seq


def scramble(sequence, attempts, preserve_content=True):
    """
    Scramble letters of a given sequence. Returns most distant sequence.
    
    Parameters
    ----------
    sequence : string
        
    attempts : int
        Number of scrambles which are created. Most dissimilar one is chosen. Ignored if 
        content is not preserved.
    preserve_content : bool
        If True, shuffles the existing sequence. If False, a completely arbitrary sequence is created.
    
    Returns
    -------
    scramble_seq : string
        Most distant scrambled sequence
    """

    # Make multiple scrambles if content is preserved
    if preserve_content:
        scrambles = np.empty(attempts, dtype='<U{}'.format(len(sequence)))
        # Create scrambles
        if attempts == None:
            raise ValueError("`attempts` has to be given when site is shuffeled.")
        else:
            if not isint(attempts):
                raise ValueError("`attempts` has to be an integer.")
            else:
                # If type float, change to int
                attempts = int(attempts)

        for i in range(attempts):
            snip = "".join(random.sample(sequence, len(sequence)))
            scrambles[i] = snip

        # Compute number of mismatches
        distances = np.zeros(attempts)
        for i in range(attempts):
            distances[i] = np.sum([x != y for (x, y) in zip(sequence, scrambles[i])])
        scramble_seq = scrambles[np.argmax(distances)]
    else:
        scramble_seq = _make_distant_sequence(sequence)
    
    return scramble_seq
    

def _make_distant_sequence(sequence):
    seq_list = list(sequence)
    scramble_seq_list = []
    for x in seq_list:
        sub_alph = filter_mutation(x, np.array(["A", "C", "G", "T"]))
        scramble_seq_list.append(np.random.choice(sub_alph, 1)[0])
    return "".join(scramble_seq_list)


def scramble_list(
    sequence, 
    windowsize, 
    overlap,
    attempts=None,
    number=1,
    preserve_content=True, 
    ignore_imperfect_scrambling=False,
    check_enzymes=[]
):
    """
    Scramble letters in a window of given sequence. 
    
    Parameters
    ----------
    sequence : string
        
    windowsize : int
        Size of the window within letters are scrambled
    overlap : int
        Number of letters that each window overlaps with the neighbors
    attempts : int
        Number of scrambles which are created. Most dissimilar one is chosen.
    number : int
        Number of scrambles per window.
    preserve_content : bool, default True
        If True, shuffles the existing sequence. If Flase, a completely arbitrary sequence is created.
    ignore_imperfect_scrambling : bool, default False 
        If False, returns an error when scrambles are not created evenly (last scramble does not end at
        sequence end).
        
    Returns
    -------
    scrambled_sequences : list
        List of scrambled sequences
    """
    
    # Check argument types
    if not isint(windowsize):
        raise TypeError("`windowsize` has to be an integer.")
    else:
        # If type float, change to int
        windowsize = int(windowsize)
        
    if not isint(overlap):
        raise TypeError("`overlap` has to be an integer.")
    else:
        # If type float, change to int
        overlap = int(overlap)
        
    if not isint(attempts):
        raise TypeError("`attempts` has to be an integer.")
    else:
        # If type float, change to int
        attempts = int(attempts)
    
    if overlap >= windowsize:
        raise ValueError("overlap cannot be equal to or bigger than windowsize.")
        
    if not isint(number):
        raise TypeError("`windowsize` has to be an integer.")
    else:
        # If type float, change to int
        number = int(number)
        
    # Compute number of scrambles
    l_sequence = len(sequence)
    n_scrambles = (l_sequence - windowsize) / (windowsize - overlap) + 1
    
    # Test if scrambles can created throughout whole site
    if not n_scrambles.is_integer():
        if not ignore_imperfect_scrambling:
            raise ValueError("Cannot make scrambles with windowsize {} and overlap {}.".format(windowsize, overlap))
        else:
            n_scrambles = np.floor(n_scrambles)
            warnings.resetwarnings()
            warnings.warn("Imperfect scrambles. Last scramble is omitted.")
        
    
    # Create array to store sequences and add wild type
    scrambled_sequences = np.empty(number * int(n_scrambles) + 1, dtype='<U160')

    scrambled_sequences[0] = sequence
    
    # Check for restriction sites in wild type
    if len(check_enzymes) > 0:
        background_sites = scan_enzymes([sequence], enzymes=check_enzymes)
    
    # Function to find right indices for scrambles
    indices = lambda i: slice((windowsize-overlap) * i, (windowsize-overlap) * i + windowsize)
    
    # Take subsequence and get scramble
    for i in range(int(n_scrambles)):
        temp_sequence = list(sequence)
        for j in range(number):
            enzyme_check = True
            while enzyme_check:
                sc = scramble(sequence[indices(i)], attempts, preserve_content=preserve_content)
                temp_sequence[indices(i)] = sc.lower()
                if len(check_enzymes) > 0:
                    if sum(scan_enzymes(["".join(temp_sequence)], enzymes=check_enzymes)) == sum(background_sites):
                        enzyme_check = False
                else:
                    enzyme_check = False
            scrambled_sequences[i * number + j + 1] = "".join(temp_sequence)
    
    return scrambled_sequences



def create_scrambles(
    sequence, 
    windowsize, 
    overlap, 
    attempts=None,
    number=1,
    preserve_content=True, 
    ignore_imperfect_scrambling=False,
    check_enzymes=[]
):
    """
    Scramble letters in a window of given sequence. Returns pandas df with metadata.
    
    Parameters
    ----------
    sequence : string
        
    windowsize : int
        Size of the window within letters are scrambled
    overlap : int
        Number of letters that each window overlaps with the neighbors
    attempts : int
        Number of scrambles which are created. Most dissimilar one is chosen.
    number : int
        Number of scrambles per window.
    preserve_content : bool, default True
        If True, shuffles the existing sequence. If Flase, a completely arbitrary sequence is created.
    ignore_imperfect_scrambling : bool, default False 
        If False, returns an error when scrambles are not created evenly (last scramble does not end at
        sequence end).
    check_enzymes : list
        List of enzymes for which restriction sites the scrambles are checked for.
        
    Returns
    -------
    scrambled_sequences : Pandas.DataFrame
        DataFrame of scrambled sequences.
    """
    
    # Check argument types
    if not isint(windowsize):
        raise TypeError("`windowsize` has to be an integer.")
    else:
        # If type float, change to int
        windowsize = int(windowsize)
        
    if not isint(overlap):
        raise TypeError("`overlap` has to be an integer.")
    else:
        # If type float, change to int
        overlap = int(overlap)
        
    if not isint(attempts):
        raise TypeError("`attempts` has to be an integer.")
    else:
        # If type float, change to int
        attempts = int(attempts)
        
    if not isint(number):
        raise TypeError("`attempts` has to be an integer.")
    else:
        # If type float, change to int
        number = int(number)
    
    # Get scrambles
    scrambled_sequences = scramble_list(
        sequence, 
        windowsize, 
        overlap, 
        attempts,
        number,
        preserve_content,
        ignore_imperfect_scrambling,
        check_enzymes
    )
    
    # Get number of scrambles
    n_scrambles = int(np.floor((len(sequence) - windowsize) / (windowsize - overlap))) + 1
    
    # Compute start and end positions of scrambles
    start_pos = np.arange(0, int(n_scrambles), 1) * (windowsize - overlap)
    stop_pos = np.arange(0,int(n_scrambles), 1) * (windowsize - overlap) + windowsize
    
    # Account for multiple scrambles per window
    if number > 1:
        start_pos = np.array([start_pos[i] * np.ones([n_scrambles, number])[i, :] for i in range(n_scrambles)]).flatten()
        stop_pos = np.array([stop_pos[i] * np.ones([n_scrambles, number])[i, :] for i in range(n_scrambles)]).flatten()
    
    # Store sequences in data frame
    scramble_df = pd.DataFrame(
        {'start_pos':start_pos, 
         'stop_pos':stop_pos, 
         'sequence':scrambled_sequences[1:]}
    )
    
    # Compute center of scrambles
    scramble_df['center_pos'] = scramble_df[['start_pos','stop_pos']].mean(axis = 1)
    
    return scramble_df


def create_scrambles_df(
    df, 
    windowsize, 
    overlap, 
    attempts, 
    number=1, 
    preserve_content=True, 
    ignore_imperfect_scrambling=False,
    check_enzymes=[]
):
    """
    Function to apply create_scrambles to a dataframe of wt sequences.
    
    Parameters
    ----------
    df : pd.DataFrame
        Should contain columns 'wt_sequence', 'promoter', 'TSS', and 'direction', with the appropriate data.
    windowsize : int
        Size of the window within letters are scrambled. Passed to create_scramble_df.
    overlap : int
        Number of letters that each window overlaps with the neighbors. Passed to create_scramble_df.
    attempts : int
        Number of scrambles which are created. Most dissimilar one is chosen. Passed to create_scramble_df.
    number : int
        Number of scrambles per window. Passed to create_scramble_df.
    preserve_content : bool, default True
        If True, shuffles the existing sequence. If False, a completely arbitrary sequence is created. Passed to create_scramble_df.
    ignore_imperfect_scrambling : bool, default False 
        If False, returns an error when scrambles are not created evenly (last scramble does not end at
        sequence end). Passed to create_scramble_df.
    check_enzymes : list
        List of enzymes for which restriction sites the scrambles are checked for.
    Returns
    -------
    df_results : pd.DataFrame
        DataFrame of scrambled sequences for many different WT sequences.
    """
    
    df_tmp = pd.DataFrame()
    df_results = pd.DataFrame()
    
    for _,row in df.iterrows():
        wt_seq = str(row['wt_sequence'])
    
        df_tmp = create_scrambles(
            wt_seq, 
            windowsize, 
            overlap, 
            attempts, 
            number, 
            preserve_content, 
            ignore_imperfect_scrambling,
            check_enzymes
        )
        
        df_tmp['promoter'],df_tmp['TSS'],df_tmp['direction'] = row['promoter'],row['TSS'],row['direction']
        
    
        df_results = pd.concat([df_results,df_tmp])
    
    return df_results
    

def scramble_site(
    wt_seq, 
    rel_start_site, 
    rel_end_site, 
    attempts=100, 
    preserve_content=False, 
    minus_10_start=None, 
    minus_10_end=None, 
    minus_35_start=None, 
    minus_35_end=None
    ):
    
    """
    Scrambles a site within a sequence, defined by coordinates relative to the wt_seq string. 
    Will avoid scrambling -10/-35 regions if given. Assumes that site at least overlaps wt_seq. 
    Assumes minus 10/35 are contained entirely within wt_seq.
    
    Parameters
    ----------
    wt_seq : string
        The WT sequence to be scrambled. Indexes will be relative to this string, indexed from 0. Ideally string is upper case.
    rel_start_site : int
        The start position of the site to scramble. Relative to the beginning of the wt_seq string.
    rel_end_site : int
        The end position of the site to scramble. Relative to the beginning of the wt_seq string.
    attempts : int
        Number of attempts to draw the scramble from. Passed to scramble()
    preserve_content : bool
        Should the nucleotide content of the scramble site be preserved or should random nucleotides be used? Passed to scramble()
    minus_10_start : int
        The start position of the minus 10 to preserve. Relative to the beginning of the wt_seq string.
    minus_10_end : int
        The end position of the minus 10 to preserve. Relative to the beginning of the wt_seq string.
    minus_35_start : int  
        The start position of the minus 10 to preserve. Relative to the beginning of the wt_seq string.
    minus_35_end : int
        The start position of the minus 10 to preserve. Relative to the beginning of the wt_seq string.    
        
    Returns
    ---------
    scramble : string
        The wt sequence with the defined regions scrambled (or preserved). Scrambles shown in lowercase.
    """
    
    #Check that start and end of site are in bounds
    #set to 0 or max if out of bounds.
    if rel_start_site<0 :
        rel_start_site = 0
    
    if rel_end_site>len(wt_seq) : 
        rel_end_site = len(wt_seq)
        
    rel_start_site = int(rel_start_site)
    rel_end_site = int(rel_end_site)
    
    site = scramble(wt_seq[rel_start_site:rel_end_site], attempts, preserve_content)
        
    scramble_seq = wt_seq[0:rel_start_site] + site.lower() + wt_seq[rel_end_site:len(wt_seq)]

    minus_10 = (type(minus_10_start)==int) & (type(minus_10_end)==int)
    
    if type(minus_10_start)==float:
        minus_10 = (~np.isnan(minus_10_start)) & (~np.isnan(minus_10_end))
        
    minus_35 = (type(minus_35_start)==int) & (type(minus_35_end)==int)
    
    if type(minus_35_start)==float:
        minus_35 = (~np.isnan(minus_35_start)) & (~np.isnan(minus_35_end))
        
    if minus_10:
        
        minus_10_start = int(minus_10_start)
        minus_10_end = int(minus_10_end)
        
        scramble_seq = scramble_seq[0:minus_10_start] + wt_seq[minus_10_start:minus_10_end] + scramble_seq[minus_10_end:len(scramble_seq)]
    
    if minus_35:
        
        minus_35_start = int(minus_35_start)
        minus_35_end = int(minus_35_end)
        
        scramble_seq = scramble_seq[0:minus_35_start] + wt_seq[minus_35_start:minus_35_end] + scramble_seq[minus_35_end:len(scramble_seq)]
    
    return scramble_seq    


def add_restriction_site(sequence_list, enzyme):
    """
    Add a restriction site to the end of each sequence in a list.
    """
    # Enzyme classes in Biopython are weird but this works
    if not ((type(enzyme) == str) or (type(enzyme).__bases__[-1] == Restriction.RestrictionType)): 
        raise TypeError("enzyme has to be either string or of type of Bio.Restriction.Restriction.RestrictionType")
    
    if type(enzyme) == str:
        try:
            getattr(sys.modules[__name__], enzyme)
        except AttributeError:
            raise ValueError("There is no restriction enzyme {}. Check spelling, naming is case-sensitive.".format(enzyme))
        enzyme = getattr(sys.modules[__name__], enzyme)
        
    # Check inputs
    if type(sequence_list) == pd.core.series.Series:
        sequence_list = copy.deepcopy(sequence_list.values)
        
    sequence_list = _check_sequence_list(sequence_list)
    site = enzyme.site
    
    return np.array([str(seq + site) for seq in sequence_list], dtype=object)
    
    
def mutations_det(
    sequence, 
    num_mutants=None, 
    mut_per_seq=1, 
    site_start=0, 
    site_end=None, 
    alph_type="DNA",
    keep_wildtype=False
):
    """Create all possible mutants deterministically and choose subset or all.

    Choose this method for creating mutants from a sequence when the number of possible
    mutants is in the same order of magnitude as the number of chosen mutants to prevent
    duplicates.
    
    
    Parameters
    ----------
    - sequence : string
        DNA sequence that is going to be mutated.
    - num_mutants : int, default None
        Number of mutant sequences. If None, all possible mutatants are created.
    - mut_per_seq : int, default 1
        Number of mutations per sequence.
    - site_start : int, default 0
        Beginning of the site that is about to be mutated.
    - site_end : int, default -1
        End of the site that is about to be mutated.
    - alph_type : string, default "DNA"
        Can either be "DNA" for letter sequences, or "Numeric" for integer sequences.
    - keep_wildtype : bool, default False
        If True, adds wild type sequence as first sequence in the list.
    Returns
    -------
    """
    
    if not (isint(num_mutants) or num_mutants == None):
        raise TypeError("`num_mutants` is of type {} but has to be integer valued.".format(type(num_mutants)))
        
    if not isint(mut_per_seq):
        raise TypeError("`mut_per_seq` is of type {} but has to be integer valued.".format(type(mut_per_seq)))
        
    if not isint(site_start):
        raise TypeError("`site_start` is of type {} but has to be integer valued.".format(type(site_start)))
        
    if not (isint(site_end) or site_end == None):
        raise TypeError("`site_end` is of type {} but has to be integer valued.".format(type(site_end)))
        
    if alph_type not in ["DNA", "Numeric"]:
        raise ValueError("`alph_type` has to be either \"DNA\" or \"Numeric\" ")
        
    
    # Get site to mutate
    if site_end==None:
        mutation_window = sequence[site_start:]
    else:
        mutation_window = sequence[site_start:site_end]

    # Compute number of possible mutations
    poss_mutants = int(np.prod([3 * (len(mutation_window) - i) for i in range(mut_per_seq)]) / factorial(mut_per_seq))

    if num_mutants == None:
        num_mutants = poss_mutants
    elif num_mutants > poss_mutants:
        warnings.resetwarnings()
        warnings.warn("Trying to create more unique mutants than possible. Choosing maximal number of {} mutants.".format(poss_mutants))
        num_mutants = poss_mutants
    elif num_mutants/poss_mutants < 0.1:
        warnings.resetwarnings()
        warnings.warn("Sampling less than 10 percent of possible mutations. Consider generating mutations at random to save computational recources.")

    # Create list 
    if keep_wildtype:
        mutants = np.empty(num_mutants+1, dtype='U{}'.format(len(sequence)))
        mutants[0] = sequence
        i_0 = 1
    else:
        mutants = np.empty(num_mutants, dtype='U{}'.format(len(sequence)))
        i_0 = 0
        
    mutant_indeces = create_mutant_index(mutation_window, num_mutants, mut_per_seq)
    if alph_type == "DNA":
        letters = np.array(["A", "C", "G", "T"])
    elif alph_type == "Numeric":
        letters = np.arange(4)
    else:
        raise ValueError("Alphabet type has to be either \"DNA\" or \"Numeric\"")
    
    for i, x in enumerate(mutant_indeces):
        mutants[i + i_0] = sequence[0:site_start] + mutate_from_index(mutation_window, x, letters) + sequence[site_end:]
        
    return mutants



def create_mutant_index(sequence, num_mutants, mut_per_seq):
    """Create index for mutations. Each index is a tuple of a position in the
    sequence and an integer from 0 to 2 that resembles the base (excluding wild type
    base) for the mutation.
    
    Parameters
    ----------
    - sequence : string
        Wildtype sequence
    - num_mutants: int
        Number of mutations.
    - mut_per_seq : int
        Number of mutations per mutant sequence
        
    Returns
    -------
    - mutants : list of arrays of tuples
        List of mutations for each sequence.
        
    """
    mutants = []
    # Generate single mutants
    mutants = [(j, i) for i in range(3) for j in range(len(sequence))]
    if mut_per_seq > 1:
        somelists = mut_per_seq * [mutants]
        elements = np.array(list(itertools.product(*somelists)))
        mask = np.empty(len(elements), dtype=bool)
        for i, element in enumerate(elements):
            # No multiple mutations on the same site
            positions = [el[0] for el in element]
            mask[i] = (len(np.unique(positions)) == mut_per_seq)
            # Check for right order of mutants
            mask[i] = mask[i] * np.all(np.diff(positions) >= 0)
        mutants = elements[mask]
    else:
        mutants = [np.array([x]) for x in mutants]
    
        
    if num_mutants < len(mutants):
        mutants_inds = np.random.choice(np.arange(len(mutants), dtype=int), num_mutants, replace=False)
        temp_mutants = []
        for index in mutants_inds:
            temp_mutants.append(mutants[index])
            
        mutants = temp_mutants
    return mutants


def mutate_from_index(sequence, index, alph):
    """Generate mutated sequence from wild type sequence given a an index and 
    an alphabet.
    
    Parameters
    ----------
    - sequence : string
        Wild type seqence
    - index : array of tuples
        Each tuple is a position and the mutation at that position.
    - alph : array 
        Alphabet to pick mutation from
        
    Returns
    -------
    - str
        Sequence including mutation
    """
    seq_list = list(sequence)
    for loci, mutation in index:
        seq_list[loci] = filter_mutation(seq_list[loci], alph)[mutation].lower()
    return "".join(seq_list)


def mutations_rand(
    sequence, 
    num_mutants,
    rate,
    site_start=0, 
    site_end=None, 
    alph_type="DNA",
    number_fixed=False,
    keep_wildtype=False,
):
    """Creates single or double mutants.
    
    
    Parameters
    ----------
    - sequence : string
        DNA sequence that is going to be mutated.
    - num_mutants : int, default None
        Number of mutant sequences. If None, all possible mutatants are created.
    - mut_per_seq : int, default 1
        Number of mutations per sequence.
    - site_start : int, default 0
        Beginning of the site that is about to be mutated.
    - site_end : int, default -1
        End of the site that is about to be mutated.
    - alph_type : string, default "DNA"
        Can either be "DNA" for letter sequences, or "Numeric" for integer sequences.
    - number_fixed : bool
        If True, the number of mutations is fixed as the rate times length of the sequence.
    - keep_wildtype : bool, default False
        If True, adds wild type sequence as first sequence in the list.
        
    Returns
    -------
    - mutants : list
        List of mutant sequences. Each element is a string.
    """
    
    if not (isint(num_mutants) or num_mutants == None):
        raise TypeError("`num_mutants` is of type {} but has to be integer valued.".format(type(num_mutants)))
    
    if not type(rate) == float:
        raise TypeError("`rate` is of type {} but has to be a float.".format(type(rate)))
        
    if not isint(site_start):
        raise TypeError("`site_start` is of type {} but has to be integer valued.".format(type(site_start)))
        
    if not (isint(site_end) or site_end == None):
        raise TypeError("`site_end` is of type {} but has to be integer valued.".format(type(site_end)))
        
    if alph_type not in ["DNA", "Numeric"]:
        raise ValueError("`alph_type` has to be either \"DNA\" or \"Numeric\" ")
        
    
    # Get site to mutate
    
    if site_end==None:
        mutation_window = sequence[site_start:]
    else:
        mutation_window = sequence[site_start:site_end]

    
    # Create list
    if keep_wildtype:
        mutants = np.empty(num_mutants+1, dtype='U{}'.format(len(sequence)))
        mutants[0] = sequence
        i_0 = 1
    else:
        mutants = np.empty(num_mutants, dtype='U{}'.format(len(sequence)))
        i_0 = 0
        
    mutant_indeces = random_mutation_generator(mutation_window, rate, num_mutants, number_fixed)
    if alph_type == "DNA":
        letters = np.array(["A", "C", "G", "T"])
    elif alph_type == "Numeric":
        letters = np.arange(4)
    else:
        raise ValueError("Alphabet type has to be either \"DNA\" or \"Numeric\"")
    
    for i, x in enumerate(mutant_indeces):
        mutants[i + i_0] = sequence[0:site_start] + mutate_from_index(mutation_window, x, letters) + sequence[site_end:]
        
    return mutants


def random_mutation_generator(sequence, rate, num_mutants, number_fixed):
    mutant_list = np.empty(num_mutants, dtype=object)
    for i in range(num_mutants):
        mutant_list[i] = _random_mutation_generator(sequence, rate, number_fixed)
    
    return mutant_list
    

def add_primers(sequence_list, primer_index, autocomplete=False, len_to_complete=200, rev_only=False):
    """Add a orthogonal primers to sequences and possibily add random positions up to certain length.
    
    Parameters
    ----------
    - sequence_list : list of sequences
        List of sequences.
    - primer_index : int
        Index of the primer pair that is added to the sequences.
    - autocomplete : boolean
        If True, adds random base pairs to the end of the sequence until length matches `len_to_complete`.
    - len_to_complete : int
        Length of total sequence after primers are added.
        
    Returns
    - new_seq_list : list of sequences
        List of modified sequences with primer pairs.
    
    """
    
    # Check input
    if type(sequence_list) not in [list, np.ndarray, pd.core.series.Series]:
        raise TypeError("sequence_list has to be list, numpy array or pandas series.")
    else:
        if any([type(seq) not in [str, Seq, np.str_] for seq in sequence_list]):
            raise TypeError("entries in `sequence_list` have to be of type string or Bio.Seq.Seq.")
            
    if not isint(primer_index):
        raise TypeError("`primer_index` is of type {} but has to be integer valued.".format(type(primer_index)))
        
    if not isint(len_to_complete):
        raise TypeError("`len_to_complete` is of type {} but has to be integer valued.".format(type(len_to_complete)))
    elif len_to_complete < 0:
        raise ValueError("`len_to_complete` has to be non-negative.")
        
    if type(autocomplete) != bool:
        raise TypeError("`len_to_complete` is of type {} but has to be boolean.".format(type(autocomplete)))
        
    if type(rev_only) != bool:
        raise TypeError("`len_to_complete` is of type {} but has to be boolean.".format(type(rev_only)))
        
    # Import primers
    forward = import_primer_fwd(primer_index)
    reverse = complement_seq(import_primer_rev(primer_index), rev=True)
    
    # Copy sequence_list to not override input
    new_seq_list = copy.deepcopy(sequence_list)

    # Add primers to sequences and autocomplete if chosen to
    for i, seq in enumerate(sequence_list):
        if rev_only:
            new_seq_list[i] = seq + reverse
        else:
            new_seq_list[i] = forward + seq + reverse
        if autocomplete and len(sequence_list[i]) < len_to_complete:
            new_seq_list[i] = new_seq_list[i] + gen_rand_seq(len_to_complete - len(new_seq_list[i]))
        elif len(new_seq_list[i]) > len_to_complete:
            warnings.resetwarnings()
            warnings.warn("Sequence {} is longer than {}bp.".format(i, len_to_complete))
      
    return new_seq_list


    
@numba.njit
def _random_mutation_generator(sequence, rate, number_fixed):
    if number_fixed:
        num_mutations = int(rate*len(sequence))
    else:
        num_mutations = np.random.poisson(len(sequence) * rate)
    positions = np.random.choice(np.arange(len(sequence)), num_mutations, replace=False)
    mutants = np.random.choice(np.arange(3), num_mutations)
    return  [(x, y) for (x, y) in zip(positions, mutants)]
    

@numba.njit
def filter_letter(x, letter):
    return x != letter


@numba.njit
def filter_mutation(letter, alph):
    j = 0
    for i in range(alph.size):
        if filter_letter(alph[i], letter):
            j += 1
    result = np.empty(j, dtype=alph.dtype)
    j = 0
    for i in range(alph.size):
        if filter_letter(alph[i], letter):
            result[j] = alph[i]
            j += 1
    return result




