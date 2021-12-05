# `processing`

This folder contains all code executed to transform or generate data. Within
this directory, some subdirectories contain summaries of specific experiments,
the code used to process and transform the data, and when indicated the output
of any processing functions.

# Sequencing data
The following steps are taken to transform raw `fastq` files into tidy data
formats that can be easily handled in Python for further analysis.

## Sequencing demultiplexing

All scripts named `demultiplex_seq.py` split `fastq` files given a series of
index reads. To demultiplex our sequencing runs we used the
[`qiime2`](https://qiime2.org) platform. The `qiime2` developing team strongly
suggests installing the platform in a separate `conda` environment. To do so,
there must be already a functional [`Anaconda
distribution`](https://www.anaconda.com/distribution/) on your computer. 
Instructions on installing Qiime can be found [here](https://docs.qiime2.org/2021.4/install/native/).
Here we show the installation instructions on  for the version of `qiime2` that is used for analysis of data
in this project. There might be newer versions of `qiime2` available, but we cannot guarantee that
the analysis pipeline works for these versions. Note that the only difference between the official instructions
and the ones shown here is the name of the conda environment that is created.
**macOS/OS X (64-bit)**
```
conda update conda
conda install wget

wget https://data.qiime2.org/distro/core/qiime2-2019.7-py36-osx-conda.yml
conda env create -n qiime2 --file qiime2-2019.7-py36-osx-conda.yml

rm qiime2-2019.7-py36-osx-conda.yml
```

**Linux (64-bit)**
```
conda update conda

wget https://data.qiime2.org/distro/core/qiime2-2019.7-py36-linux-conda.yml
conda env create -n qiime2 --file qiime2-2019.7-py36-linux-conda.yml

rm qiime2-2019.7-py36-linux-conda.yml
```
This will create a new conda environment and install the necessary libraries to run `qiime2` basic functions.

In general, to run `qiime2` the environment needs to be activated. This is done by writing
in the command line (replace `conda` with `source` in Windows)
```
conda activate qiime2
```

The to switch back the base conda environment, deactivate the current environment by
```
conda deactivate
```

To run our customized python scripts the environment does not need to be activated.
It is done automatically by the `subprocess` module.

The demultiplexing scripts assume the file structure for this project that
looks like:
```
+---code
|   +---processing
|        +---demultiplex_seq.py
|
+---data
    +---raw_sequencing
    |    +---DATE_SEQUENCING-INFORMATION
    |        +---sequencing_barcodes_qiime2.tsv (barcodes list)
    |        +---MISEQ (files as given by Illumina)
    |            +---Data
    |                +---Intensities
    |                    +---BaseCalls
    |                        +---R1.fastq.gz (forward read)
    |                        +---R2.fastq.gz (reverse read)
    |                        +---I1.fastq.gz (index read)
```
Where all names written in capital letters are given as variables in the `demultiplex_seq.py` script. The folder structure for
data is quite strict, but can be modified, as long as the `demultiplex_seq.py` is adapted appropriately.

To run the script, simply execute the `demultiplex_seq.py` file from the terminal.

The index barcodes needed for demultplexing are kept in `sequencing_barcodes_qiime2.tsv`. For this file, make sure the
columns are separated by tabs, and not spaces, otherwise `qiime` complains.

## Demultiplexed sequencing processing

After the raw sequences have been split into individual `fastq` files for each
of the indexes, the sequences need to be processed based on their quality,
length, and in the case of paired-end reads, the forward and reverse read must
be stuck together. To perform these tasks, we use
[`fastp`](https://github.com/OpenGene/fastp). Follow the installation instructions
in the repository. Make sure to install a newer version than 0.19.1, since merging was
introduced in this version.
After that, all scripts named `processing_seq.py` can be run to process the
short-reads. These scripts assume that the data has been demultiplexed already
as it takes the resulting `fastq` files from the `demux_sequencing` folder
indicated above. The output of these scripts is saved in
```
+---data
    +---processed_sequencing
        +---DATE_experiment_description
```
Specifically `fastp` generates individual `fastq` files with all the reads
(merged for paired-end runs) for each index. It also generates `HTML` and
`JSON` summaries of the sequencing processing, listing average read quality,
base composition per position, among other useful quantities.
