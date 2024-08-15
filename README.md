# count_indels_integrations

Originally written by Shengdar Tsai, modified by Yichao Li.

# Usage

```

module load parallel
module load bwa
module load samtools/1.7
module load trimmomatic/0.36
module load flash
module load fastqc/0.11.5


mkdir {{jid}}/${COL3}_results

fastqc ${COL1} -o {{jid}}/log_files/
fastqc ${COL2} -o {{jid}}/log_files/

trimmomatic PE ${COL1} ${COL2} -threads 4 -baseout {{jid}}/${COL3}_results/${COL3} ILLUMINACLIP:/hpcf/apps/trimmomatic/install/0.36/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3


# Run flash to combine reads 

cd {{jid}}/${COL3}_results

flash ${COL3}_1P ${COL3}_2P --max-overlap=120 --min-overlap={{FLASH_min_overlap}} --output-prefix=${COL3}


bwa mem -t 6 /research_jude/rgs01_jude/groups/tsaigrp/projects/Genomics/common/genomes/hg38/hg38.fa ${COL3}.extendedFrags.fastq > ${COL3}.sam
samtools sort -o ${COL3}.bam ${COL3}.sam
samtools index ${COL3}.bam
samtools flagstat ${COL3}.bam > ${COL3}.flagstat
samtools stats ${COL3}.bam > ${COL3}.stats


module load conda3/202011

source activate /home/yli11/.conda/envs/cutadaptenv/

cd ../..

#bed_fix.py ${COL4} output.bed

/home/yli11/.conda/envs/cutadaptenv/bin/python /research_jude/rgs01_jude/groups/tsaigrp/projects/Genomics/common/src/cas9engineering/count_indels_integrations.py --bed ${COL4} --ref /research_jude/rgs01_jude/groups/tsaigrp/projects/Genomics/common/genomes/hg38/hg38.fa --bam {{jid}}/${COL3}_results/${COL3}.bam --out {{jid}}/${COL3}_results/${COL3}.counts.tsv --site ${COL5}


```
