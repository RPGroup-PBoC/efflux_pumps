import git
import skbio
import pandas as pd
import numpy as np
import os


# Find home directory for repo
repo = git.Repo("./", search_parent_directories=True)
homedir = repo.working_dir

# Define data directory
datadir = f"{homedir}/data/processed_sequencing/"

# List all fastq.gz files
fastq_file = datadir + "efflux_merged.fastq.gz"

# Use skbio to have a generator to iterate over fastq
print("Loading sequences...\n")

seqs = skbio.io.read(
    fastq_file,
    format="fastq",
    verify="false",
    variant="illumina1.8",
)

# Initialize list to save sequence objects
seq_list = list()
# Iterate over sequences
print("Extracting sequence information...\n")
for seq in seqs:
    # Extract sequence information
    seq_id = seq.metadata["id"]
    quality = seq._positional_metadata["quality"].values
    sequence = str(skbio.DNA(sequence=seq, validate=True))
    # Append to list
    seq_list.append([sequence, quality, np.sum(quality)])

# Initialize dataframe to save sequences
names = ["sequence", "quality", "total_quality"]
df_seq = pd.DataFrame.from_records(seq_list, columns=names)

# Add index and sequence length to dataframe
df_seq["seq_len"] = df_seq.sequence.apply(len)


# Filter sequences and find barcode/promoter
print("Finding promoters and barcodes...\n")
df_filt = df_seq[df_seq.seq_len == 220].reset_index(drop=True)
df_filt.insert(4, "barcode", [seq[0:20] for seq in df_filt.sequence])
df_filt.insert(4, "promoter", [str(skbio.DNA(sequence=seq[-160:]).complement(reverse=True)) for seq in df_filt.sequence])

# Load ordered sequences
df_list = pd.read_csv(f"{homedir}/data/efflux_pumps_twist.csv")
df_list.insert(2, "promoter", [seq[20:180] for seq in df_list.sequence])

# Find promoters
print("Identifying promoters...\n")
name_list = df_list.name.values
list_promoters = df_list.promoter.values
sequenced_promoters = df_filt.promoter.values
ident_list = []

for promoter in sequenced_promoters:
    index = np.where(promoter == list_promoters)
    if len(index[0]) != 0:
        ident_list.append(name_list[index[0][0]])
    else:
        ident_list.append("None")
        
df_filt.insert(5, "identified_promoter", ident_list)

# Return found sequences
print("Saving results...\n")
df_return = df_filt[["identified_promoter", "promoter", "barcode"]]
df_return[df_return.identified_promoter != 'None']


# Define output dir
outputdir = f"{homedir}/data/barcodes/"

if not os.path.exists(outputdir):
    os.mkdir(outputdir)

df_return.to_csv(f"{outputdir}/barcode_mapping.csv", index=False)