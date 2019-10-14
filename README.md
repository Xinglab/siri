## SQUID: Stringent Quantitation of Unspliced Intron by Deep-sequencing

Introduction
----------
SQUID is a general tool for quantifying the percent of intron inclusions in RNA-seq data. 
<p><figure class="figure1" data-title="HOMER motif"><img alt="" src="docs/intron_type.png" /><figcaption></figcaption></figure></p>

#### Types of PI (Percent of Introns) Calculation

  PI_Junction: 
    Inclusion counts divided by the sum of inclusion and skipping junction counts
  PI_Density:
    The observed counts divided by the expected counts of the intron
        
#### Types of intron retention

  U: 
    Introns that are not partly overlapped with exons or overlapped with other introns.
        
  E:
          Introns that are partly overlapped with exons but not overlapped with other introns.
        
    I:  
    Introns that are overlapped with other introns but not overlapped with exons.
    
  EI: 
    Introns that are both overlapped with exons and other introns.

Requirements
------------
1. Install [Python 2.7.x](https://www.python.org/downloads)
2. Install [pysam](https://pypi.python.org/pypi/pysam/0.8.4)
3. Install [STAR](https://github.com/alexdobin/STAR) for the run without alignment files provided. The star index of human and mouse can be downloaded through the following link
http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/. Install [kallisto](https://pachterlab.github.io/kallisto/download) for the run without FPKM files provided.
5. Install [cufflinks](http://cole-trapnell-lab.github.io/cufflinks/install/) for the fun without FPKM and fastq files provided.  
6. Install [DEXSeq](https://bioconductor.org/packages/devel/bioc/html/DEXSeq.html) to run differential spliced intron analysis

Installation
------------
The source code can be directly called from Python.

Usage
--------------------------------
Run SQUID with provided fastq files
  
  python ../SQUID.py --GTF ./test.gtf  --fastq ./test_R1_1.fq:./test_R1_2.fq,./test_R2_1.fq:./test_R2_2.fq,./control_R1_1.fq:./control_R1_2.fq,./control_R2_1.fq:./control_R2_2.fq --check_len true --index_kallisto ./kallisto/test --index_star ./star --anchor 8 --length 100 --lib first --read P --Cal All  --c1 0.05  --p 1 --Comparison ./Comparison --analysis U -o ./bam1 --resume true

Please note that this run require high RAM due to STAR alignment  

Run SQUID with provided alignment files and fastq  files

  python ../SQUID.py --GTF ./test.gtf  --fastq ./test_R1_1.fq:./test_R1_2.fq,./test_R2_1.fq:./test_R2_2.fq,./control_R1_1.fq:./control_R1_2.fq,./control_R2_1.fq:./control_R2_2.fq --align ./test_R1.bam,./test_R2.bam,./control_R1.bam,./control_R2.bam --check_len true --index_kallisto ./kallisto/test  --anchor 8 --length 100 --lib first --read P --Cal All  --c1 0.05  --p 1 --Comparison ./Comparison --analysis U -o ./bam2

Run SQUID with provided alignment files and FPKM files

  python ../SQUID.py --GTF ./test.gtf --align ./test_R1.bam,./test_R2.bam,./control_R1.bam,./control_R2.bam --FPKM transcript_exp.txt --anchor 8 --length 100 --lib unstrand --read P --Cal All  --c1 0.05  --p 1 --Comparison ./Comparison --analysis U -o ./bam3

Required Parameters
------------
  --align:
    s1.bam/s1.sam[,s2.bam/s2.sam]. Mapping results for all of samples in bam/sam format. Different samples are separated by commas
  --GTF:
    The gtf file
  --fastq: 
    s1_1.fq[:s1_2.fq][,s1_1.fq[:s2_2.fq],...]. The raw sequencing reads in fastq or fastq format that is required to call kallisto to calculate FPKM value
  --index_star:
    The path to the star index that is required to do the alignment using STAR
   
Optional Parameters
------------  
  --FPKM:
                A file providing the FPKM value for each sample, the first column is transcript ID with the following column being the FPKM value for each sample. If it is not provided, kallisto will be called to calculate FPKM value
  --index_kallisto:
                The path to the kallisto index that is required to run kallisto from raw reads
  --o/--output:
    The output directory. The default is current directory
  --check_len: 
    Whether to generate new fastq files to with equal read length. The default value is false
  --l:
    Estimated average fragment length. The parameter to run kallisto with default value of 200
  --s:
    Estimated standard deviation of fragment length. The parameter to run kallisto with default value of 100
  --update:
    Whether to update the attributes of introns using spliced reads. The default is false
  --lib:
    The library type with choices of unstrand/first/second. The details are explained in the parameter of library-type in tophat2. The default is unstrand
  --read: 
    The sequencing strategy of producing reads with choices of P(paired end) or S (single end). The default is P
  --length: 
    The read length of sequencing reads in nucleotide. The default length is 100
  --anchor: 
    The anchor length in nucleotide. The program will only count reads spanning junctions with at least this anchor length on each side. The default is 8
  --Cal: 
    Which  part of the program user choose to run, the choices are All/count. All means run the whole program including count and differential analysis of spliced introns, count means only run the PI value calculation part. The default is All.
  --Comparison: 
    A file providing the sample pairs to calculate the differential RI level.The format should be column 1(name of comparions), column 2 (sample 1 order in the align files,replicates separated by commas), column 3 (sample 2 order in the align files,replicates separated by commas), column 4 (optional, if present as 'pool', the replicates are combined together in rMATS calculation). If absent, the step of calculation of differential spliced introns  will be skipped
  --analysis: 
    Type of rMATS analysis to perform. analysisType is either P or U. P is for paired analysis and U is for unpaired analysis. Default is U
  --c1: 
    The cutoff of splicing difference of rMATS run using Junction method. The cutoff used in the null hypothesis test for differential splicing. The default is 0.0001
  --p: 
    The number of threads used to run rMATS. The default is 1;
  --F_deltaPI: 
    The cutoff of delta to output differential spliced introns.The default is 0.05;
  --F_FDR: 
    The cutoff of combined FDR to output differential spliced introns.The default is 0.1;
  --resume:
    Whether to resume previous run. The default is false.

Output list
------------

Notes: $n denotes the number of samples provided in SQUID run

#### Result
The folder contains the final output of four types of files: 
    <code>
    intron_PI.txt,
      Diff_$comparison_intron_PI.txt,
      Decrease_$comparison_intron_PI.txt,
      Increase_$comparison_intron_PI.txt.
    </code>
$comparison denotes the label of comparison performed

    intron_PI.txt is the file containing the info of annotated introns
    Diff_$comparison_intron_PI.txt containing the differential spliced introns info of all annotated introns
    Decrease_$comparison_intron_PI.txt containing the differential spliced introns info of the introns with decreased PI in sample2 
    Increase_$comparison_intron_PI.txt containing the differential spliced introns info of the introns with increased PI in sample2

Common columns in all files:

    Intron_id:         Intron Id representing the chromosome position, start and end
    Gene_id:           Gene id of intron residing genes
    Strand:            Strand of intron residing genes
    Chr:               Chromosome name of introns
    Start:             Start coordinate of introns
    End:               End coordinate of introns
    Annotated:         Whether this intron was annotated in the gtf file as retained intron event
    Attributes:        One of the four intron types, U/E/I/EI
    Inclusion_counts:  Inclusion counts separated by commas
    Skipping_counts:   Skipping counts separated by commas
    Inclusion_length:  Effective inclusion length
    Skipping_length:   Effective skipping length
    PI_Junction:       PI_Junction value separated by commas
    Observed_counts:   Observed counts separated by commas
    Expected_counts:   Expected counts separated by commas
    PI_Density:        PI_Density separated by commas


Extra columns in Diff_$comparison_intron_PI.txt, Decrease_$comparison_intron_PI.txt,Increase_$comparison_intron_PI.txt.
    
    PValue_rMATS:      p-value from rMATS
    FDR_rMATS:         FDR from rMATS
    PValue_DEXSeq:     p-value from DEXSeq
    FDR_DEXSeq:        FDR from DEXSeq
    Combined_FDR:    FDR of rank product test of rMATS and DEXSeq

#### test
A folder contains test files to run the program

#### gtf_files
An intermediate folder contains different types of gtf files to run the program. Use mouse genome as examples.

  Mus_musculus.Ensembl.GRCm38.78.gtf: the ensemble gtf files. This file should be provided by user. 
  Exon_Mus_musculus.Ensembl.GRCm38.78.gtf: the gtf file contains exons only
  Intron_Mus_musculus.Ensembl.GRCm38.78.gtf: the gtf file contains intron only
  Intron_Annotated_Mus_musculus.Ensembl.GRCm38.78.gtf: the gtf file contains the attributes whether the intron was annotated as retained introns in the original gtf files
  Intron_attri_Mus_musculus.Ensembl.GRCm38.78.gtf: the gtf file contains the attributes whether the intron was overlapped with Exon and whether the intron is overlapped with other intron. 

#### fq
An intermediate folder contains all of the fastq with equal read length files

#### align
An intermediate folder contains all of the alignment files

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
    column n+3:          The chromosome of the gene residing 
    
count_all_Density.txt: a file contains the observed counts and expected counts for all of the introns

    column 1:           Intron id representing the chromosome position, start and end.
    column 2:           The length of introns
    column 3~n+2:       The observed counts
    column n+3~2n+2:    The expected counts
    
Total.txt: a file contains total number of unique reads in each sample

    column 1~n:        Total number of unique reads in samples 1~n

#### FPKM   
An intermediate optional folder contains the result of FPKM result and gene expression files for the squid run without gene expression file provided. 

  kallisto_$n
    The result of kallisto of each sample
  transcript_exp.txt
    The file is FPKM file that contains FPKM value for each gene.
    column 1:        Gene ID
    column 2~n+1:    FPKM value for samples 
    
#### rMATS_files
An intermediate folder contains all of the rMATS input and output files

#### DEXSeq_files
An intermediate folder contains all of the DEXSeq input and output files


Contacts and bug reports
------------------------
Yi Xing
xingyi@email.chop.edu

If you found a bug or mistake in this project, we would like to know about it.
Before you send us the bug report though, please check the following:

1. Are you using the latest version? The bug you found may already have been
   fixed.
2. Check that your input is in the correct format and you have selected the
   correct options.
3. Please reduce your input to the smallest possible size that still produces
   the bug; we will need your input data to reproduce the problem, and the
   smaller you can make it, the easier it will be.


Copyright and License Information
---------------------------------
Copyright (C) 2015 University of California, Los Angeles (UCLA)
Yi Xing

Authors: Zhicheng Pan, Shaofang Li, Yi Xing

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see http://www.gnu.org/licenses/.


Frequently Asked Questions
---------------------------------
Q: How much RAM does SQUID require?

A: It depends on the input you used in the SQUID. If the input is alignment files and FPKM files, 4G is enough for the human genome. If FPKM is not provided, FPKM will be generated by Kallisto which may takes 1G extra RAM. If alignment files is not provided, STAR will generate alignment files that may takes up to 40G to for the human genome.

Q: Are the observation counts always larger than the expected counts?

A: There is a small portion of introns that have larger observation counts and PI_Density is assigned as 1 in these introns.

Q: Does the index in the comparison file 0-based or 1-based?

A: The index in the comparison file is 1-based. 

Q: Can SQUID resume previous run?

A: Yes, SQUID will resume previous run by setting parameter resume as true. Please delete the intermediate files that you’d like to regenerate. 

Q: Do the reads in the fastq files have to be the same length?

A: No. If the reads in the fastq files do not have equal length, please set parameter check_len to true to generate new fastq files with equal length.

Q: Can I run SQUID to output intron retention values without doing differential analysis? 

A: Yes. You can run SQUID by ignoring --Comparison parameters. 