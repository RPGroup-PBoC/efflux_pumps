#%%
import os
import itertools
import subprocess
import glob
import git

#%% 
# Date for sequencing run
DATE = 20221209
# Description for project. Must be something that shortly summarizes
# the reason for this sequencing run. It MUST end with /
DESCRIPTION = '_mapping/'

# Conda environment name where qiime2 is installed
CONDA_ENV = 'qiime2'

# Find project parental directory
repo = git.Repo('./', search_parent_directories=True)
homedir = repo.working_dir

# Path containing all raw input information
INPUT_DIR = f'{homedir}/data/sequencing/{DATE}{DESCRIPTION}' 

# Path containing raw input sequences
RAW_DIR = INPUT_DIR + '221209_M05340_0384_000000000-GFLHM/Data/Intensities/BaseCalls/'

# Define requirements for qiime2 demultiplexing
# Path to barcode TSV file
BARCODE_LIST = INPUT_DIR + 'sequencing_barcodes_qiime2.tsv'

# Path to output folder
OUTPUT_DIR = f'{homedir}/data/demux_sequencing/' +\
             f'{DATE}{DESCRIPTION}'

# Define temporary directory
TMP_DIR = OUTPUT_DIR + 'tmp/'
# Define diretory to save all of qiime2 files
QIIME2_DIR = OUTPUT_DIR + 'qiime2_output/'

# Generate output directory if it doesn't exist
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# Generate tmp directory if it doesn't exist
if not os.path.exists(TMP_DIR):
    os.makedirs(TMP_DIR)

# Generate qiime2_output directory if it doesn't exist
if not os.path.exists(QIIME2_DIR):
    os.makedirs(QIIME2_DIR)

#%%

# Pattern on each of the fastq files
FORWARD_STR = 'Undetermined*R1*.fastq.gz'
REVERSE_STR = 'Undetermined*R2*.fastq.gz'
BARCODE_STR = 'Undetermined*I1*fastq.gz'
# Path to each file
FORWARD_SEQ = glob.glob(RAW_DIR + FORWARD_STR)[0]
REVERSE_SEQ = glob.glob(RAW_DIR + REVERSE_STR)[0]
BARCODE_SEQ = glob.glob(RAW_DIR + BARCODE_STR)[0]

# Define requirements for qiime2
sequence_type = 'EMPPairedEndSequences'  # Paired end illumina sesquences
artifact_file = f'{QIIME2_DIR}{DATE}_sequence_artifact.qza'  # output artifact

# Check if this process hasn't been done before
print('checking if qiime2 sequence artifact has been generated before.')
#if not os.path.exists(artifact_file):
# Temporarily rename files so that qiime2 can recognize them
print('temporarily moving fastq files...')
os.rename(FORWARD_SEQ, TMP_DIR + 'forward.fastq.gz')
os.rename(REVERSE_SEQ, TMP_DIR + 'reverse.fastq.gz')
os.rename(BARCODE_SEQ, TMP_DIR + 'barcodes.fastq.gz')

# Define bash command to generate artifact for qiime to register seqeunces
bash_artifact = f'''qiime tools import \
--type {sequence_type} \
--input-path  {TMP_DIR} \
--output-path {artifact_file}
'''

# Run subprocess activating qiime2 conda environmentß
print(f'generating qiime2 artifact on conda {CONDA_ENV} environment.')
subprocess.run(f'source activate {CONDA_ENV} && {bash_artifact}',
                shell=True)

# Rename files back to original name
print('renaming fastq files back to original name...')
os.rename(TMP_DIR + 'forward.fastq.gz', FORWARD_SEQ)
os.rename(TMP_DIR + 'reverse.fastq.gz', REVERSE_SEQ)
os.rename(TMP_DIR + 'barcodes.fastq.gz', BARCODE_SEQ)
print('DONE!')

#%%
# Define requirements for qiime2 demultiplexing
# Column in metadata containing barcodes
barcode_col = 'barcode-sequence'
# File to save the demultiplexed qiime2 artifact
demux_file = f'{QIIME2_DIR}{DATE}_demux_artifact.qza'
# File to save the details for the demultiplexed sequences
demux_details = f'{QIIME2_DIR}{DATE}_demux_details.qza'

# Define shell command to demultiplex sequences
# NOTE: the --p-no-golay-error-correction option is added because we have 8nt
# barcodes rather than 12nt barcodes, so the so-called Golay correction 
# cannot be applied
bash_demux = f'''qiime demux emp-paired \
    --m-barcodes-file {BARCODE_LIST} \
    --m-barcodes-column {barcode_col} \
    --i-seqs {artifact_file} \
    --o-per-sample-sequences {demux_file} \
    --o-error-correction-details {demux_details} \
    --p-no-golay-error-correction
  '''

# Check if this process hasn't been done before
print('checking if sequences have been demultiplexed already.')
#if not os.path.exists(demux_file):
# Run subprocess activating qiime2 conda environment
print(f'demultiplexing sequences.')
subprocess.run(f'source activate {CONDA_ENV} && {bash_demux}', shell=True)
print('DONE!')
#%%
# Path to save interactive plots output for the analysis
demux_viz = f'{QIIME2_DIR}{DATE}_demix_viz.qzv'  

# bash commnand to generate interactive plots and export sequences
bash_viz = f'''qiime demux summarize \
  --i-data {demux_file} \
  --o-visualization {demux_viz} \
'''
bash_exp = f'''qiime tools export \
  --input-path {demux_file} \
  --output-path {OUTPUT_DIR} 
'''
# Check if this process hasn't been done before
print('checking if data has already been exported.')
if not os.path.exists(demux_viz):
    print(f'generating interactive plots and exporting data.')
    subprocess.run(f'source activate {CONDA_ENV} && {bash_viz} && {bash_exp}',
                   shell=True)

print('DONE!')