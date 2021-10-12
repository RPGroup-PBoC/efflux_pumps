import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api
from .utils import isint
import numba

def gen_emat_rand(site_size, mean=1, sd=1):
    """
    Generate a random energy matrix for a defined sequence length. Arbitrary values for each possible base, normally distributed around mean 1 with standard deviation 1.
    
    Parameters
    ----------
    site_size : int
        Length of the sequence to generate the energy matrix for, in bp.
    mean : float
        Mean of entries in energy matrix.
    sd : float
        Standard deviation of entries in energy matrix.
    Returns
    ----------
    energy_matrix : np.array
    """
    # Check argument types
    if not isint(site_size):
        raise ValueError("`site_size` has to be an integer.")
    else:
        # If type float, change to int
        site_size = int(site_size)
        
    energy_matrix = np.random.normal(mean, sd, (site_size, 4))
    return energy_matrix


def gen_emat_single_site(
    seq, 
    site_start, 
    site_size, 
    site_mean=1, 
    site_sd=1, 
    background_mean=0, 
    background_sd=0):
    """
    Generate energy matrix for sequence with one site. 
    Mean and sd values can be set for site and non site positions.
    WT sequence is set to zero.
    
    Parameters
    ----------
    seq : string
        Sequence. Used to set entries for wild type binding site to zero.
    site_start : int
        First base of binding site
    site_size : int
        Length of binding site.
    site_mean: float
        mean energy for site mutations, for np.random.normal
    site_sd: float
        standard deviation of energy for site mutations, for np.random.normal
    background_mean: float
        mean energy for non site mutations, for np.random.normal
    background_sd: float
        standard deviation of energy for non site mutations, for np.random.normal
    
    Returns
    ---------
    seq_emat : pd.DataFrame
        generated energy matrix
    """
    
    # Check argument types
    if not isint(site_start):
        raise ValueError("`site_start` has to be an integer.")
    else:
        # If type float, change to int
        site_start = int(site_start)
        
    if not isint(site_size):
        raise ValueError("`site_size` has to be an integer.")
    else:
        # If type float, change to int
        site_size = int(site_size)
        
    if not isinstance(site_mean, int) or isinstance(site_mean, float):
        raise ValueError("`site_mean` has to be an integer or float.")
    
    if not isinstance(site_sd, int) or isinstance(site_sd, float):
        raise ValueError("`site_sd` has to be an integer or float.")
        
    if not isinstance(background_mean, int) or isinstance(background_mean, float):
        raise ValueError("`background_mean` has to be an integer or float.")
        
    if not isinstance(background_sd, int) or isinstance(background_sd, float):
        raise ValueError("`background_sd` has to be an integer or float.")
        
        
    # Set background values
    seq_emat = np.random.normal(background_mean, background_sd, (len(seq), 4))
    
    # Set site values
    seq_emat[site_start:(site_start + site_size), :] = np.random.normal(site_mean, site_sd,(site_size, 4))
    
    # Convert np.array to pd.DataFrame
    seq_emat = pd.DataFrame(data=seq_emat, columns=('A','T','C','G'))
    
    # Set WT values = 0
    for ind,char in enumerate(seq, start=0):
        seq_emat.iloc[ind][char] = 0
    
    return seq_emat 


def sum_emat(seq, emat):
    """
    Retrieve and sum the energy matrix values for a given sequence variant and matrix.
    
    Parameters
    ----------
    seq : string
    
    emat : pd.DataFrame with columns A T C G
    
    Returns
    ---------
    sum : float
    """
    
    mat_vals = np.zeros(len(emat.index))
    
    for ind,char in enumerate(seq, start = 0):
        mat_vals[ind] = (emat.iloc[ind][char])
        
    return np.sum(mat_vals) 


@numba.njit()
def sum_emat_arr(seq, emat):
    """
    Retrieve and sum the energy matrix values for a given sequence variant and matrix.
    Uses numba to speed up computation.
    
    Parameters
    ----------
    seq : string
    
    emat : numpy array with columns A C G T
    
    Returns
    ---------
    sum : float
    """
    
    mat_vals = np.zeros(len(seq))
    letter_to_int = {"A": 0, "C": 1, "G": 2, "T": 3}
    for ind in range(len(seq)):
        mat_vals[ind] = emat[ind, letter_to_int[seq[ind]]]
        
    return np.sum(mat_vals) 


def sum_emat_df(scrambles_df, emat):
    """
    Sums energy matrices for a dataframe of scrambles with `sum_emat()`.
    
    Parameters
    ------------
    scrambles_df : pd.DataFrame 
        Output by `create_scrambles_df()`
    emat : pd.DataFrame
        Energy matrix output by `gen_emat_single_site()`
        
    Returns
    ----------
    scrambles_df : pd.DataFrame
        Including the additive 'effect' column next to each scramble
    """
    
    scrambles_df['effect'] = np.nan
    emat_arr = emat[["A", "C", "G", "T"]].to_numpy()
    
    for ind, scr_seq in enumerate(scrambles_df['sequence'], start = 0):
        scrambles_df.at[ind, 'effect'] = sum_emat_arr(seq = scr_seq, emat = emat_arr)
        
    return(scrambles_df)


def gen_barcode_effects(barcode_num, barcode_noise, df):
    """
    Generate barcode effects for each scramble. Effects are drawn from a normal around the defined effect.
    Wildtype effects are sampled from normal centered at 0.
    
    Parameters
    ------------
    barcode_num : int
        Number of barcodes to generate for each scramble
    barcode_noise : float
        standard deviation of normal distribution to draw effects from for each barcode
    df : pd.DataFrame
        Dataframe with scrambles and effects output by `sum_emat_df()`
    """
    
    #Generate wildtype barcode effects from normal(0, barcode_noise)
    wt_bc_effects = np.random.normal(loc = 0, scale = barcode_noise, size = barcode_num)
    
    #Initialize new columns
    df['barcode_effects'] = ''
    df['wt_barcode_effects'] = ''
    df['p_val'] = ''
    
    # Iterate through the scrambles in the dataframe
    for i in range(len(df)):
        
        #Generate barcode effects from normal(effect, barcode_noise) and calculate barcode mean
        barcode_effects = np.random.normal(loc=df.iloc[i]['effect'], scale=barcode_noise, size = barcode_num)
        bc_mean = np.mean(barcode_effects)
        
        #Add vals to dataframe
        df.at[i,'barcode_effects'] = barcode_effects
        df.at[i,'bc_mean'] = bc_mean
        df.at[i,'wt_barcode_effects'] = wt_bc_effects
        
        #Perform t-test on scramble barcode effects vs. wt barcode effects 
        df.at[i,'p_val'] = stats.ttest_ind(df.iloc[i]['barcode_effects'], df.iloc[i]['wt_barcode_effects'], equal_var = False)[1]
    
    #Correct for multiple significance tests. Gives adjusted pval and significance call.
    stats_corrected = statsmodels.stats.multitest.multipletests(df['p_val'])
    
    df['adj_p_val'] = stats_corrected[1]
    df['sig'] = stats_corrected[0]

    return(df)

def gen_scramble_dataset(
    seq_length = 50,
    replicates = 100,
    windowsize = 10, 
    overlap = 5, 
    attempts = 100, 
    preserve_content = True, 
    site_start = 20, 
    site_size = 10, 
    site_mean = 1, 
    site_sd = 1, 
    background_mean = 0, 
    background_sd = 0,
    barcode_num = 10,
    barcode_noise = 1):
    """
    Generate a scramble dataset with replicate sequences drawn from the same parameters. 
    Wraps gen_rand_seq(), gen_emat_single_site(), create_scrambles_df(), sum_emat_df(), and gen_barcode_effects().
    
    Parameters
    -------------
    seq_length : int
        Length of sequence to generate in bp
    replicates : int
        Number of individual sequences to generate and scramble
    windowsize : int
        Size of scramble in bp
    overlap : int
        Overlap of scrambles in bp
    attempts : int
        Number of scrambles which are created. Most dissimilar one is chosen.
    preserve_content : bool
        If True, shuffles the existing sequence. If False, a completely arbitrary sequence is created.
    site_start : int
        First base of binding site
    site_size : int
        Length of binding site.
    site_mean: float
        mean energy for site mutations, for np.random.normal
    site_sd: float
        standard deviation of energy for site mutations, for np.random.normal
    background_mean: float
        mean energy for non site mutations, for np.random.normal
    background_sd: float
        standard deviation of energy for non site mutations, for np.random.normal
    barcode_num : int
        Number of barcodes to generate for each scramble
    barcode_noise : float
        standard deviation of normal distribution to draw effects from for each barcode
        
    Returns
    ---------
    results : pd.DataFrame
        Dataframe containing generated scrambles for each sequence, WT and barcode effects, and significance test results.
    """
    
    results = pd.DataFrame()
    
    for i in range(replicates):
        
        #if (i % 10)==0: print(i)
        # Generate WT sequence
        seq = wgregseq.gen_rand_seq(seq_length)
        
        # Generate energy matrix
        emat = wgregseq.gen_emat_single_site(seq = seq, 
                                             site_start = site_start, 
                                             site_size = site_size, 
                                             site_mean = site_mean, 
                                             site_sd = site_sd, 
                                             background_mean = background_mean, 
                                             background_sd = background_sd)
        
        # Generate scrambles
        scrambles = wgregseq.create_scrambles_df(sequence = seq, 
                                                 windowsize = windowsize, 
                                                 overlap = overlap, 
                                                 attempts = attempts, 
                                                 preserve_content = True)
        
        # Sum effects for each scramble
        scramble_effects = wgregseq.sum_emat_df(scrambles_df = scrambles, emat = emat)
        scramble_effects['rep'] = i

        barcode_effects = wgregseq.gen_barcode_effects(barcode_num = barcode_num, barcode_noise = barcode_noise, df = scramble_effects)
        
        results = results.append(barcode_effects)
    
    return(results)
        
def merge_sig_scrambles(df):
    """
    Merge significant adjacent scrambles with same sign effect (positive or negative) to find regulatory sites.
    
    Parameters
    -----------
    df : pd.DataFrame
        Dataframe of generated scrambles and effects, probably from gen_scramble_dataset(). 
        Must contain columns sig : bool, and bc_mean : float
        
    Returns
    ----------
    site_positions : pd.DataFrame
        An aggregated dataframe of all unique site positions from merged scrambles. 
    """
    sig_pos_bool = (df['sig'] == True) & (df['bc_mean']>0)
    sig_neg_bool = (df['sig'] == True) & (df['bc_mean']<0)
    
    df_sig_pos = df[sig_pos_bool]
    df_sig_neg = df[sig_neg_bool]
    
    df_sig_pos['site_id']=df_sig_pos.groupby((~sig_pos_bool).cumsum()).grouper.group_info[0]
    df_sig_neg['site_id']=df_sig_neg.groupby((~sig_neg_bool).cumsum()).grouper.group_info[0]
    
    site_positions_pos = df_sig_pos.groupby(['rep','site_id']).agg({'start_pos':'min', 'stop_pos':'max'})
    site_positions_pos['effect_sign'] = '+'

    site_positions_neg = df_sig_neg.groupby(['rep','site_id']).agg({'start_pos':'min', 'stop_pos':'max'})
    site_positions_neg['effect_sign'] = '-'

    site_positions = site_positions_pos.append(site_positions_neg)
    site_positions['center_pos'] = site_positions[['start_pos','stop_pos']].mean(axis = 1)
    site_positions['site_size'] = site_positions['stop_pos'] - site_positions['start_pos']

    return(site_positions)