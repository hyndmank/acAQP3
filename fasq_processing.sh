### DOWNLOAD REFERENCE GENOME MUS MUSCULUS GRCm39
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-107/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz .
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-107/gtf/mus_musculus/ .


### UNPACK FILES IN THE DIRECTORY
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz 
gunzip Mus_musculus.GRCm39.107.gtf.gz


### GENERATE GENOME FILES USING STAR PACKAGE VER. 2.7.10a 
cd /data/user/vanhuynh/RefGenome/Mouse107/
#!/bin/bash
STAR --runThreadN 6 --runMode genomeGenerate --genomeDir /data/user/vanhuynh/RefGenome/Mouse107/ --genomeFastaFiles Mus_musculus.GRCm39.dna.primary_assembly.fa --sjdbGTFfile Mus_musculus.GRCm39.107.gtf --sjdbOverhang 99


### TRIM ILLUMINA ADAPTOR SEQUENCE USING TRIM-GALORE PACKAGE
#!/bin/bash 
for f1 in *_R1_001.fastq.gz 
do
        f2=${f1%%_R1_001.fastq.gz}"_R2_001.fastq.gz"
        trim_galore --illumina --paired --fastqc -o /data/user/vanhuynh/BulkRNASeq/QRMutants/AdLib/Trim $f1 $f2 
done


### SEQUENCE QUALITY CONTROL
cd /data/user/vanhuynh/BulkRNASeq/QRMutants/AdLib
fastqc *.fastq.gz # CAN RUN AFTER TRIMING 


### GENERATE BAM AND BAI FILES USING STAR PACKAGE VER. 2.7.10a 
#!/bin/bash
index=/data/user/vanhuynh/RefGenome/Mouse107/
prj=/data/user/vanhuynh/BulkRNASeq/QRMutants/AdLib/Trim
fq1=(${prj}/*_R1_001_val_1.fq.gz)
fq2=(${prj}/*_R2_001_val_2.fq.gz)

#### MAPPING TO GENOME
for ((i=0;i<"${#fq1[@]}";i++)); do
  sample="${fq1[$i]%%_R*}"
  rgline="ID:${sample}    PU:${sample}    SM:${sample}    PL:ILLUMINA LB:${sample}"
  outprefix="${sample}"

  STAR \
  --runThreadN 6 \
  --genomeDir $index \
  --readFilesIn "${fq1[$i]}" "${fq2[$i]}" \
  --readFilesCommand gunzip -c \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMunmapped Within \
  --outFileNamePrefix $outprefix \
  --outSAMattrRGline $rgline \
  --quantMode GeneCounts
done

#### GENERATE .BAI FILES FOR EACH BAM 
ml SAMtools
cd /data/user/vanhuynh/BulkRNASeq/QRMutants/AdLib/Trim/
#!/bin/bash
for i in *.bam
do
echo "Indexing: "$i
samtools index $i $i".bai"
done


### FEATURECOUNTS (PART OF SUBREAD PACKAGE)
featureCounts -p -t exon -g gene_id -a /data/user/vanhuynh/RefGenome/Mouse107/Mus_musculus.GRCm39.107.gtf -o counts.txt *.bam


### TRIM THE TABLE TO CUT THE FEATURE INFO, KEEPING FIRST COLUMN ENSEMBL GENE AND COUNTS FOR EACH SAMPLE 
cut -f 1,7-36 featureCounts_counts.txt > featureCounts_counts_clean.txt


