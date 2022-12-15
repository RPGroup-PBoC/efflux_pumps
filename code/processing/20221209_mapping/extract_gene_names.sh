#!/bin/bash

# Find working directiory

PARENT_PATH=$(dirname $(greadlink -f $0))
RESULT="${PARENT_PATH##*/}"

FILE='efflux_merged.fastq.gz'

FILENAME=${FILE%.fastq.gz}
FILENAME=${FILENAME#*/}

# Go back path
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}
PARENT_PATH=${PARENT_PATH%/*}

# Find data directory
PROCESSING_FOLDER=$PARENT_PATH'/data/ident_promoters/'$RESULT
SAM_FILE=$PROCESSING_FOLDER'/'$FILENAME'_collapsed.sam'


BBMAP_PATH=$(find $PARENT_PATH -name "bbmap.sh")

# Make directories for stored data
OUT_DIR=$PARENT_PATH'/data/barcodes/'$RESULT'/'$FILENAME'_per_gene'
if [ -d $OUT_DIR ] 
then
    echo $OUT_DIR" exists"
else 
    mkdir $OUT_DIR
fi

echo $BBMAP_PATH

$BBMAP_PATH ref=$PARENT_PATH'/data/efflux_promoters.fasta' path=$PROCESSING_FOLDER

$BBMAP_PATH ambiguous='best' indelfilter='0' nfilter='0' minid='0.85' trimreaddescriptions='t' in=$PROCESSING_FOLDER/$FILENAME'.fasta' out=$SAM_FILE t='8' path=$PROCESSING_FOLDER


# Start the process
start=$SECONDS

# Input file name
echo $OUT_DIR
echo "Assigning gene names to promoter-barcode pairs..."

# Extract promoter-bc pairs and corresponding gene names
awk -v o="$OUT_DIR" 'BEGIN{FS="\t";OFS = ","} !(NR%500000){print NR " Promoters Processed"}; NF>10{gsub(/_/, ",", $1); print $10,$1,$3 >> (o"/"$3"_barcodes.txt")}' "$SAM_FILE"

# terminal output message
#echo "All Promoters Processed, now adding headers..."
#
# Loop through output directory and add a header to each file
#cd $out_dir
#echo "promoter,barcode,counts,name" > headerfile
#for file in *barcodes.txt; do cat headerfile $file > tmpfile2; mv tmpfile2 "$file"; done

#rm headerfile

# terminal output message
echo "done! Output files written to " "$OUT_DIR"
end=$SECONDS
duration=$(( end - start ))
echo
echo "time elapsed: $duration seconds"
