# ==============================================================================
# 00b_kallisto.sh
# .fastq --> pseudoaligned reads
# ==============================================================================

# modified by a script originally written by Dr. Pamela Russell (thanks!)

# requirements:
# * Homo sapiens, Ensembl Transcriptome v96 retrieved from
#  https://github.com/pachterlab/kallisto-transcriptome-indices/releases
# * kallisto v0.46.0  (used Mac binary)






# core kallisto function
# ------------------------------------------------------------------------------

cd /home/liucu/MJ_RNAseq/

cat > runkallisto.sh
#!/bin/bash

if [ "$#" -ne 1 ]; then
  echo "Usage: runkallisto <sample id>"
  exit -1
fi

sample=$1
reads_dir=trimmed_reads/
output_dir=kallisto_output
index=Software/kallisto/homo_sapiens/transcriptome.idx
kallisto_path=Software/kallisto/kallisto

fqs=$(ls -l $reads_dir/*trimmed.fastq.gz | grep "/${sample}_" | awk 'BEGIN {x = ""} {x = x " " $9} END {print x}')

$kallisto_path quant \
--index=$index \
--output-dir=$output_dir/$sample \
--single \
--rf-stranded \
--fragment-length=200 \
--sd=20 \
$fqs




# extract sample names from .fastq.gz to use runkallisto.sh
# ------------------------------------------------------------------------------


# extract sample names from .fastq.gz files -->
# apply "runkallisto.sh" to each
ls trimmed_reads/*trimmed.fastq.gz > run_all_kallisto.sh
sed -i '' "s/trimmed_reads\\//.\\/runkallisto.sh /" run_all_kallisto.sh
sed -i '' 's/_.*/; sleep 30; /' run_all_kallisto.sh
sed -i '' -n 'g;n;p' run_all_kallisto.sh






chmod u+x runkallisto.sh
chmod u+x run_all_kallisto.sh
chmod u+x Software/kallisto/kallisto # $kallisto_path

./run_all_kallisto.sh >> run_all_kallisto.log
