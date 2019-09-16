#!/bin/sh

gtf_to_fasta ../test.gtf ./test.fa test_transcripts.fa
cat test_transcripts.fa | perl -ne 'if ($_ =~/^\>\d+\s+\w+\s+(ERCC\S+)[\+\-]/){print ">$1\n"}elsif($_ =~ /\d+\s+(ENST\d+)/){print ">$1\n"}else{print $_}' > test_transcripts.clean.fa
cat test_transcripts.clean.fa | grep ">" | perl -ne '$_ =~ s/\>//; print $_' | sort | uniq > transcript_id_list.txt
kallisto index --index=test test_transcripts.clean.fa

