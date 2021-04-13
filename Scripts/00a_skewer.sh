# ==============================================================================
# 00a_skewer.sh
# adapter + low base quality trimming, CentOS Linux 7
# ==============================================================================


#!/bin/bash

# sequencing run 1
reads_dir=/home/liucu/MJ_RNAseq/Run1
output_dir=/trimmed_reads 
skewer=/Software/skewer-0.2.2-linux-x86_64

mkdir $output_dir

for fq in $reads_dir/*001.fastq.gz; do
  bn=$(basename $fq .fastq.gz)
  $skewer --compress -q 10 -l 31 --output $output_dir/$bn $fq
done


# sequencing run 2
reads_dir_two=/Run2
output_dir=/trimmed_reads 
skewer=/Software/skewer-0.2.2-linux-x86_64
for fq in $reads_dir_two/*001.fastq.gz; do
  bn=$(basename $fq .fastq.gz)
  $skewer --compress -q 10 -l 31 --output $output_dir/$bn $fq
done
