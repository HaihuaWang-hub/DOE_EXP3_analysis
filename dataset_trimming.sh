####################################################
#This is for metatranscriptomic triming analysis of plant_emf_sap_interaction project
####################################################
#!/bin/bash
dir=/Volumes/T7/plant_emf_sap_interaction

cd $dir
mkdir $dir/fastqc_rawdata
mkdir $dir/cleandata
mkdir $dir/fastqc_cleandata

#  section 1: quality trimming
####################################################
#1. MD5 check
#conda install perl-digest-md5
#brew install md5sha1sum

cd $dir/rawdata
md5 *.gz > md5.txt
md5sum -c md5.txt

# 2. fastqc quality control
cd $dir/rawdata

ls $dir/rawdata/*.gz|cut -d "_" -f 1,2,3,4,5,6,7 |sort -u |while read id;do
   fastqc ${id}_001.fastq.gz \
   -t 3 \
   -o $dir/fastqc_rawdata 
   
done
   


cd $dir/fastqc_rawdata
multiqc *.zip

##3. triming with trim_galore
cd $dir/rawdata



ls *.gz|cut -d"_" -f 1,2,3 |sort -u |while read id;do
trim_galore -q 25 --phred33 --stringency 3 --length 100 \
--paired ${id}_R1_001.fastq.gz    ${id}_R2_001.fastq.gz \
--gzip \
--cores 10 \
-o ../cleandata_trimglore  

done 




ls *.gz|cut -d"_" -f 1 |sort -u |while read id;do
   echo $id
   if [ -f cleandata_rmdup/${id}_rmdup_val_2.fq.gz ]; then
      echo "${id} has been analyzed"
   else
      echo "Start to process ${id}"
      pigz -d ${id}_val_1.fq.gz 
      pigz -d ${id}_val_2.fq.gz
      cd-hit-dup  \
          -i  ${id}_val_1.fq \
          -i2 ${id}_val_2.fq \
          -o cleandata_rmdup/${id}_rmdup_val_1.fq \
          -o2 cleandata_rmdup/${id}_rmdup_val_2.fq \
          -e 0
      pigz ${id}_val_1.fq
      pigz ${id}_val_2.fq
      pigz cleandata_rmdup/${id}_rmdup_val_1.fq
      pigz cleandata_rmdup/${id}_rmdup_val_2.fq
   fi      
done


trim_galore -q 25 --phred33 --stringency 3 --length 100 \
--paired ${id}_R1_001.fastq.gz    ${id}_R2_001.fastq.gz \
--gzip \
--cores 10 \
-o ../cleandata_trimglore  






2>$dir/3_suillus_alignment/1_bowtie2_log_file/${id}_bowtie2.log


ls *.gz|cut -d"_" -f 1,2,3 |sort -u |while read id;do
trimmomatic PE -threads 3 \
${id}_R1_001.fastq.gz ${id}_R2_001.fastq.gz \
../cleandata/${id}_paired_R1.fq.gz \
../cleandata/${id}_unpaired_R1.fq.gz \
../cleandata/${id}_paired_R2.fq.gz  \
../cleandata/${id}_unpaired_R2.fq.gz  \
ILLUMINACLIP:/opt/miniconda3/envs/RNASeq/share/trimmomatic-0.39-1/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:5 TRAILING:5 -phred33 \
SLIDINGWINDOW:5:15  MINLEN:25 
done



cd $dir/cleandata
fastqc *.gz -t 4 -o $dir/fastqc_cleandata

cd $dir/fastqc_cleandata
multiqc *.zip


