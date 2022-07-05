# 💾SAM_demultiplex.py
---
## Description
Used to demultiplex an Illumina Hi-seq run using either one or two indexing primers
Written for Python 3.5+
___
## Inputs:
- Multiplexed SAM file with i7 barcodes stored in "BC:Z:" tags and i5 barcodes stored in "B2:Z" tags.
- Tab-separated table with three columns: sample name, i7 sequence, i5 sequence

## Usage:
```
python3 SAM_demultiplex.py -S <SAMfile> -T <indexTable> -O <outputDir>
```
___
## Arguments:
```
argument       options           description
-S, --sam      <str> (required)  filepath to multiplexed SAM file
-T, --table    <str> (required)  reference TSV of indices
-I, --indices  i5, i7, both      which index feature(s) to demultiplex on (default: both)
-M, --mismatch 0, 1, 2, 3        number of mismatches to allow in index sequence (default: 0)
-O, --output   <str>             output directory (default: .)
--paired                         input is paired-end (default: False)
```
___
## Output:
- One FASTQ file (or two files if paired-end) for each sample listed in the index table, written to the specified output directory.

## Testing:
Clone this repository and test SAM_demultiplex.py with the following commands:
```
git clone https://github.com/maschon0/SAM_demultiplex
cd SAM_demultiplex
python3 SAM_demultiplex.py -S test/test_file.sam -T test/test_table.tsv -O test/test_output
```

You should see a folder with 12 files identical to the file found in ['test/test_output_validation'](test/test_output_validation). The 10 reads in the SAM file are assigned to sample_1 through sample_10, while sample_11.fastq and unassigned.fastq receive no reads.
