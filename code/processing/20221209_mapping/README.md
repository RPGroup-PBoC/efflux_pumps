---
status: Accepted
---

# 2022-12-09 Barcode Mapping run

## Purpose

In this sequencing run we map the promoters two barcodes for two libraries, one for site scrambles (for TF identification) and the other is a set of efflux pumps.
## Platform

MiSeq (Caltech)
Run ID: 221209_M05340_0384_000000000-GFLHM

## Template

## Primers used

Sequencing Read 1 Primer (site scramble):

CCCTGAAactagaTCTAAACAGTTAGGCCCAGG (SC247)

Sequencing Index Primer (site scramble):

ctcgacCATCGGGTGGGATTTAGCTACAGGTAT (SC251)

Sequencing Read 2 Primer (site scramble):

ATACCTGTAGCTAAATCCCACCCGATGgtcgag (SC248)

Sequencing Read 1 Primer (efflux pumps):

Standard Illumina primers

Sequencing Index Primer (efflux pumps):

GTCGACTTCAGGGATAACATGGCACTATGCACG (SC128)

Sequencing Read 2 Primer (efflux pumps):

CGTGCATAGTGCCATGTTATCCCTGAAGTCGAC (SC127)


## Sequencing kit

NovaSeq S2 flowcell, SE reads, 26 cycles on read1

## Materials

| **id** | **barcode-sequence** | **description** |
| :--: | :--: | :--: | :--:|
efflux pumps	|	ATGGCT | Set of  *E. coli* efflux pumps
scramble	|	ACACTG | scrambled TF binding sites for identification of 


## Processing
Store the sequencing files in the following format:

```
project
│   README.md  
│
└───data
│   └───sequencing
|       └───20221209_mapping
│           │   221209_M05340_0384_000000000-GFLHM/
│           │   sequencing_barcodes_qiime2.tsv

```

The files are processed using `fastp`. The first step is to make the bash scripts in this folder executable. This can be done by using 

```
chmod +x *.sh
```

Next, we filter the sequences for quality scores and trim the trailing six bases from each read, since the first 20 bases of the read are the barcode, and the following bases are conserved. This is done by running the `processing_seq.sh` script in this folder. Run 

```
./processing_seq.sh
```

Then, we extract the barcodes and promoters from the sequencing data. Run

```
./extract_barcodes.sh
```

The script creates files containing each barcode, as well as their counts. The results will be stored in a list (and `.fasta` file), which will be used to map the sequences to promoters.
**Insert here part about barcode qc**