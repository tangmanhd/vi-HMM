# vi-HMM

## Overview

Our methid is coded in MATLAB. It is a novel method for finding SNPs adn Indels by Hidden Markov Model (HMM).  The model takes read alignments (SAM files with Phred+33 encoded quality scores) and a reference genome (FASTA file) as input to determine the most likely mutation state at each position in the reference. The it reports variants based on that state information in TXT format. We also provide a code to convert a TXT format to variant call file (VCF) format.


## Obtaining and Usage

Unzip the package. Change the current directory in Matlab to the 'vi-HMM' folder containing the code organized into subfolders. In order to run the programs, 'vi-HMM' and its sudirectories have to be added to the path. This can be achieved with the following command at the Matlab prompt

```
>> addpath(genpath(pwd))

```

The data to be analyzed with vi-HMM has to be placed in the folder 'data'. vi-HMM comes already with example data, ref.fa, example.sam and truevar.txt.

vi-HMM runs with default transition matrix and heterozygous rate for our example.  For instance:

```
>> runviHMM
```


