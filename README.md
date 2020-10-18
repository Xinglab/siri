# SIRI: Systematic Investigation of Retained Introns

## <a name="siri"></a>What is SIRI?
SIRI is a general tool for quantifying the percent of intron inclusion in RNA-seq data. It adopts parts of functions from [`SQUID`](https://github.com/Xinglab/SQUID) and [`rMATS-turbo`](https://github.com/Xinglab/rmats-turbo). 

<p>
  <figure class="figure1" data-title="HOMER motif"><img alt="" src="docs/intron_type.png" />
  <figcaption>
  </figcaption>
  </figure>
</p>

## Table of Contents

- [What is SIRI?](#siri)
- [Installation](#install)
  - [Dependencies](#depen)
  - [Install isoCirc from source](#src)
- [Getting started](#start)
- [Input and output](#input_output)
  - [Input files](#input_file)
  - [Output files](#output_file)
- [Contact](#contact)

## <a name="install"></a>Installation
### <a name="depen"></a>Dependencies
SIRI is dependent on [`pysam`](https://pypi.org/project/pysam/0.8.4/).
Please make sure that SIRI is running under python 2.7.x environment and pysam is installed before running SIRI. 

### <a name="src"></a>Install SIRI from source
You can install **SIRI** from source:
```
git clone https://github.com/Xinglab/siri.git
```
 
## <a name="start"></a>Getting started with toy example in `test`
```
../bin/SIRI --gtf test.gtf --bam test_R1.bam,test_R2.bam --anchor 8 --length 100 --lib first --read P -o SIRI_Output 
```

Detailed arguments:
```
usage: SIRI
      --gtf:
        gtf files provided for PI estimation.
      --bam/bam_files:
        s1.bam/s1.sam[,s2.bam/s2.sam]. STAR mapping results for all of samples in bam/sam format. Different samples are separated by commas.
      --anchor:
        The anchor length in nucleotide. The program will only count reads spanning junctions with at least this anchor length on each side. The default is 8.
      --length:
      The length of reads specified by user.
      --lib:
        The library type with choices of unstrand/first/second. The default is unstrand.
      -o:
        The output directory of results.
      --read:
        The sequencing strategy of producing reads with choices of P(paired end) or S (single end). The default is P
```

## <a name="input_output"></a>Input and output
### <a name="input_file"></a>Input files
SIRI takes a bam file or multiple bam files as input.

It also requires a reference gene annotation.

### <a name="output_file"></a>Output files
The folder contains the final output of PI quantification (intron.PI.txt).
  
    The detailed explanation of column names for Intron_PI.txt.
    
      Intron_id:                         Intron Id representing the chromosome position, start and end (1-index based)
      Gene_id:                           Gene id of intron residing genes
      Strand:                            Strand of intron residing genes
      Chr:                               Chromosome name of introns
      Start:                             Start coordinate of introns
      End:                               End coordinate of introns
      Annotated:                         Whether this intron was annotated in the gtf file as retained intron event
      Attributes:                        One of the four intron types, U/E/I/EI
      Inclusion_counts:                  Inclusion junction counts separated by commas
      Skipping_counts:                   Skipping counts separated by commas
      Inclusion_counts_with_intron_body: Inclusion junction counts plus counts from intron body.
      Inclusion_length:                  Effective inclusion junction length
      Skipping_length:                   Effective skipping length
      Intron_body_length:                Effective inclusion and intron body length
      PI_Junction:                       PI_Junction value separated by commas
      PI_JunctionIntron:                 PI_JunctionIntron separated by commas
      
#### gtfs
An intermediate folder contains different types of gtf files to run the program. Use test.gtf as examples.
  
    test.gtf: the ensembl gtf files. This file should be provided by user.
    Exon_test.gtf: the gtf file contains exons only
    Intron_test.gtf: the gtf file contains intron only
    Intron_Annotated_test.gtf: the gtf file contains the attributes whether the intron was annotated as retained introns in the original gtf files
    Intron_attri_test.gtf: the gtf file contains the attributes whether the intron was overlapped with Exon and whether the intron is overlapped with other introns.

#### counts
An intermediate folder contains all of the count files

count_intron.txt: a file contains the counts for all of the introns   
  
    column 1:            Intron Id representing the chromosome position, start and end.
    column 2:            Gene id
    column 3:            Strand
    column 4:            Comma separated logical values to denote Whether this intron was intron (E) or intron (I) based on read info
    column 5:            Chromosome name
    column 6:            Start coordinate
    column 7:            End coordinate    
    column 8:            Inclusion counts at 5' splice sites for sample 1
    column 9:            Skipping counts at 5' splice sites for sample 1
    column 10:           Inclusion counts at 3' splice sites for sample 1
    column 11:           Skipping counts at 3' splice sites for sample 1
    column 12:           Skipping counts of the intron for sample 1
    column 13:           counts lying in the intron for sample 1
    column 14~6*(n+1)+1: more counts for samples 2~n
    
count_exon.txt: a file contains the counts for all of the exon in each gene
    
    column 1:            Gene id
    column 2~ n+1:       Gene counts in samples 2~n
    column n+2:          Gene strand
    column n+3:          The chromosome of the gene

## Copyright and License Information
Copyright (C) 2020 Children’s Hospital of Philadelphia, University of Pennsylvania
Yi Xing

Contributors: Zhicheng Pan, Shaofang Li, Yi Xing

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see http://www.gnu.org/licenses/. 

## <a name="FAQ"></a>FAQ
## <a name="contact"></a>Contact

Zhicheng Pan zcpan1016@gmail.com

Yi Xing yi.xing@pennmedicine.upenn.edu

[github issues](https://github.com/Xinglab/siri/issues)

