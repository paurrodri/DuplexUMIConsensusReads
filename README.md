# DuplexUMIConsensusReads

calls consensus reads from paired-end, Duplex-UMI reads, preserving mapping information



## Description




## Installation

Download the script `DuplexUMIConsensusReads.py`

Or clone the repository:
```
git clone https://github.com/paurrodri/DuplexUMIConsensusReads.git
```


## Dependencies

DuplexUMIConsensusReads is dependent on Python >= 3.6. It requires the following libraries:

* sys
* argparse
* math
* random
* numpy
* pysam



## Usage

Run:
```
$ python3 DuplexUMIConsensusReads.py -i <inputfile.bam>
``` 

### Required arguments

The only required argument is the input file, given as `-i` or `--input_file`.
  
The input file is assumed to:
* be a binary Sequence Alignment/Map file (`.bam` format)
* have reads generated by paired-end sequencing, from a library with Duplex Unique Molecular Identifiers (Duplex-UMI) attached to the DNA fragments
* have the reads sorted by position
* have the Duplex-UMI sequence specified in every read, with the tag `RX`
* have the family code specified in every read, with the tag `MI`. This family code is unique for each combination of Duplex UMI sequences

The current version of `DuplexUMIConsensusReads` is recommended to be used for input files that had been processed by [fulcrumgenomics GroupReadsByUMI tool](http://fulcrumgenomics.github.io/fgbio/tools/latest/GroupReadsByUmi.html)


### Optional arguments


Optional arguments in `DuplexUMIConsensusReads` are:

| argument | flag | long option | type | description | default |
|-------------------------------------|------|----------------------------|------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------|
| output file | -o | --output_file | str | name of the output .bam file with consensus reads | <input_filename>_cons.bam |
| verbose | -v | --verbose |  | add this option to display how the program is running |  |
| minimum  mapping  quality | -q | --min_map_quality | int | minimum mapping quality (given as Phred score) of a read to be used to generate the consensus read | 20                                                                                                                                                                                                      |
| minimum base quality |  | --min_base_quality | int | minimum base quality (given as Phred score) for a base not to be masked | 20                                                                                                                                                                                                      |
| minimum read number |  | --min_reads | int | minimum number of reads from the same sub-family to generate a single-strand consensus read | 1         |
| maximum read number |  | --max_reads | int | maximum number of reads from the same sub-family to generate a single-strand consensus read | 100                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| maximum base quality |  | --max_base_quality | int | Recommended to use in case of overestimated base qualities. All bases with a quality greather than this value (given as a Phred score) will be capped before calling consensus | 60                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| base quality shift |  | --base_quality_shift | int | Recommended to use in case of overestimated base qualities. All base qualities will be substracted this value (given as a Phred score) before calling consensus | 0 |
| error rate post labeling |  | --error_rate_post_labeling | int | Recommended to use in case of having prior knowledge about this error. Error rate for an error after the UMIs have been integrated to the fragment, but prior to sequencing. Given as a Phred score. | 0 |
| error rate pre labeling |  | --error_rate_pre_labeling | int | Recommended to use in case of having prior knowledge about this error. Error rate for an errror before the UMIs have been integrated to the fragment. Given as a Phred score. | 0 |
| deletion score |  | --deletion_score | int | quality score (Phred score) for a deletion | 30                                                                                                                                                                                                                                                                                                          |
| score for the absence of  insertion |  | --no_insertion_score | int | quality score (Phred score) for the "absence of an insertion".  This refers to reads not having any insertion in a position where other reads from the same family have an insertion. | 30                                                                                                                                                                                                                                                                                                          |


## Help

To get inline help about the usage and arguments, run: 
```
python3 DuplexUMIConsensusReads.py -h
```


## Documentation

For full documentation, see [insert true link](https://github.com/paurrodri/DuplexUMIConsensusReads.git)


## Authors

* Paula Rodríguez García (<p.rodriguezgar@hotmail.com>)



_created for the Centre for Genomic Medicine - Rigshospitalet, Copenhagen, Denmark_
_2020_