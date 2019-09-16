#!/bin/sh


##run SQUID with provided fastq files
python ../SQUID.py --GTF ./test.gtf  --fastq ./test_R1_1.fq:./test_R1_2.fq,./test_R2_1.fq:./test_R2_2.fq,./control_R1_1.fq:./control_R1_2.fq,./control_R2_1.fq:./control_R2_2.fq --check_len true --index_kallisto ./kallisto/test --index_star ./star --anchor 8 --length 100 --lib first --read P --Cal All  --c1 0.05  --p 1 --Comparison ./Comparison --analysis U -o ./bam1

##run SQUID with provided alignment files and fastq  files
python ../SQUID.py --GTF ./test.gtf  --fastq ./test_R1_1.fq:./test_R1_2.fq,./test_R2_1.fq:./test_R2_2.fq,./control_R1_1.fq:./control_R1_2.fq,./control_R2_1.fq:./control_R2_2.fq --align ./test_R1.bam,./test_R2.bam,./control_R1.bam,./control_R2.bam --check_len true --index_kallisto ./kallisto/test  --anchor 8 --length 100 --lib first --read P --Cal All  --c1 0.05  --p 1 --Comparison ./Comparison --analysis U -o ./bam2

##run SQUID with provided alignment files and FPKM files
python ../SQUID.py --GTF ./test.gtf --align ./test_R1.bam,./test_R2.bam,./control_R1.bam,./control_R2.bam --FPKM transcript_exp.txt --anchor 8 --length 100 --lib unstrand --read P --Cal All  --c1 0.05  --p 1 --Comparison ./Comparison --analysis U -o ./bam3

