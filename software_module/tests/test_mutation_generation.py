import pytest
import wgregseq
import numpy as np

def test_rand_seq_generetor():
    # Test input types
    with pytest.raises(TypeError):
        wgregseq.gen_rand_seq(1.5)
    with pytest.raises(TypeError):
        wgregseq.gen_rand_seq("A")

    # Test sequence length   
    seq = wgregseq.gen_rand_seq(10)
    assert len(seq) == 10
    

def test_single_mutants():
    seq = wgregseq.gen_rand_seq(10)
    
    # Error for asking for too many mutants
    with pytest.warns(UserWarning):
        wgregseq.mutations_det(seq, mut_per_seq=1, num_mutants=100)
    
    # Single mutants
    mutants = wgregseq.mutations_det(seq, mut_per_seq=1)
    ## Test number of mutants
    assert len(mutants) == 30
    ## Test number of mutations
    assert sum([sum([x != y for (x, y) in zip (list(seq), list(mutant))]) != 1 for mutant in mutants]) == 0
    ## Test for duplicates
    assert len(np.unique(mutants)) == len(mutants)
    
    
def test_double_mutants():
    seq = wgregseq.gen_rand_seq(10)
    # Double mutants
    mutants = wgregseq.mutations_det(seq, mut_per_seq=2)
    ## Test number of mutants
    assert len(mutants) == (30 * 27) / 2
    ## Test number of mutations
    assert sum([sum([x != y for (x, y) in zip (list(seq), list(mutant))]) != 2 for mutant in mutants]) == 0
    ## Test for duplicates
    assert len(np.unique(mutants)) == len(mutants)
    

def test_triple_mutants():
    seq = wgregseq.gen_rand_seq(10)
    # Triple mutants
    mutants = wgregseq.mutations_det(seq, mut_per_seq=3)
    ## Test number of mutants
    assert len(mutants) == (30 * 27 * 24) / 6
    ## Test number of mutations
    assert sum([sum([x != y for (x, y) in zip (list(seq), list(mutant))]) != 3 for mutant in mutants]) == 0
    ## Test for duplicates
    assert len(np.unique(mutants)) == len(mutants)
    

    
    