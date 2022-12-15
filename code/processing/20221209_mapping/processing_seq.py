#%%
import os
import glob
import pandas as pd
import git

#%%
# Date for sequencing run of the library to map
DATE = 20221209
# Description to be attached to folder names
DESCRIPTION = '_mapping/'

# Find project parental folder
repo = git.Repo('./', search_parent_directories=True)
homedir = repo.working_dir

# Path to input folder
INPUT_DIR = f'{homedir}/data/demux_sequencing/{DATE}{DESCRIPTION}'

# Path to output folder
OUTPUT_DIR = f'{homedir}/data/processed_sequencing/{DATE}{DESCRIPTION}'

# Generate output directory if it doesn't exist
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

#%%
# Read list of demultiplexed sequences
seq_table = pd.read_csv(INPUT_DIR + 'MANIFEST')

# Group sequences by index to feed them into the fastp
seq_group = seq_table.groupby('sample-id')

for group, seq in seq_group:
    # Extract file names
    forward = seq[seq.direction == 'forward'].filename.values[0]
    reverse = seq[seq.direction == 'reverse'].filename.values[0]
    # Define inputs for fastp
    in1 = f'{INPUT_DIR}{forward}'  # forward input read
    in2 = f'{INPUT_DIR}{reverse}'  # reverse input read

    # Define outputs
    out1 = f'{OUTPUT_DIR}{group}_R1.fastq.gz'  
    out2 = f'{OUTPUT_DIR}{group}_R2.fastq.gz'  
    merged_out = f'{OUTPUT_DIR}{group}_merged.fastq.gz'
    html_report = f'{OUTPUT_DIR}{DATE}_{group}_fastp_report.html'
    json_report = f'{OUTPUT_DIR}{DATE}_{group}_fastp_report.json'
    report_title = f'{DATE}{DESCRIPTION} fastp report'
    
    # Define string to be ran on the terminal
    fastp = f'''fastp \
        --in1 {in1} \
        --in2 {in2} \
        --out1 {out1} \
        --out2 {out2} \
        --merge \
        --merged_out {merged_out} \
        --verbose \
        --disable_length_filtering \
        --correction \
        --overlap_len_require {10} \
        --html {html_report} \
        --json {json_report} \
        --report_title "{html_report}" \
        --thread 6
    '''

    # Run program
    os.system(fastp)