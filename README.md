# methclone: Detect the dynamic evolution of clonal epialleles in DNA methylation sequencing data

## Installation 
make 

## usage
Usage: methClone [OPTION]

methclone: a program to detect the dynamic evolution of clonal epialleles in DNA
methylation sequencing data. 

Options:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --f1=STRING           First methylation bam file for methclone
  --f2=STRING           Second methylation bam file for methclone
  -o STRING             Name of the output file for methclone
  -s STRING             Sample name
  -m INT                Cutoff for methylation difference between two samples.
Default: 0
  -d INT                Distance cutoff to consider methylated values as one
                        value. Default: 72
  -c INT                Minimum read coverage. Default: 60



