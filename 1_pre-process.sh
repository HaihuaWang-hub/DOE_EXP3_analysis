#############################################################################################
##part 1: pre-process  the sequences from metatranscriptomic and metagenomic dataset
for i in $(ls ./) ; do
    echo $i
    mv $i/Raw_Data/*.gz ./
done

splite_paired_seq(){
i=$1
base=$(basename "$i" .fastq.gz)
if [ ! -f "${base}_R1.fastq.gz" ]; then
  pigz -d $i
  cat $base.fastq | paste - - - - - - - - \
    | tee >(cut -f 1-4 | tr "\t" "\n" > ${base}_R1.fastq) \
    | cut -f 5-8 | tr "\t" "\n" > ${base}_R2.fastq
 fi
pigz $base.fastq  ${base}_R1.fastq ${base}_R2.fastq
}

export -f splite_paired_seq

time parallel -j 2 --eta --load 99% --noswap  splite_paired_seq ::: $(ls *.fastq.gz)





#quality control for the fastq sequences
##############################################################################################
fastqc rawdata/*.gz -t 20 -o fastqc_rawdata


for i in rawdata/*_R1.fastq.gz; do
   basename=$(basename "$i" _R1.fastq.gz)
   echo $basename
   if [ ! -f cleandata/${basename}_R2_val_2.fq.gz ]; then
     trim_galore -q 25 --phred33 --stringency 3 --length 110 \
               --paired rawdata/${basename}_R1.fastq.gz    rawdata/${basename}_R2.fastq.gz \
               --gzip \
               --basename $basename \
               --cores 24 \
               -o cleandata 
   fi
done



#extract the rRNA sequences https://github.com/biocore/sortmerna
############################################################################################
mkdir rRNA_removed_data
remove_rRNA(){
 nohup for name in $(ls cleandata/*_val_1.fq.gz |cut -d "/" -f 2); do
  base=$(basename $name _val_1.fq.gz)
  if [ ! -f rRNA_removed_data/$base.non.rRNA_fwd.fq.gz ]; then
      mkdir rRNA_removed_data/$base
     sortmerna \
       --workdir rRNA_removed_data/$base \
       --ref /home/microbiome/tools/sortmerna/data/rRNA_databases/rfam-5.8s-database-id98.fasta \
       --ref /home/microbiome/tools/sortmerna/data/rRNA_databases/rfam-5s-database-id98.fasta \
       --ref /home/microbiome/tools/sortmerna/data/rRNA_databases/silva-arc-16s-id95.fasta \
       --ref /home/microbiome/tools/sortmerna/data/rRNA_databases/silva-arc-23s-id98.fasta \
       --ref /home/microbiome/tools/sortmerna/data/rRNA_databases/silva-bac-16s-id90.fasta \
       --ref /home/microbiome/tools/sortmerna/data/rRNA_databases/silva-bac-23s-id98.fasta \
       --ref /home/microbiome/tools/sortmerna/data/rRNA_databases/silva-euk-18s-id95.fasta \
       --ref /home/microbiome/tools/sortmerna/data/rRNA_databases/silva-euk-28s-id98.fasta \
       --paired_out --out2  --zip-out \
       --reads cleandata/${base}_val_1.fq.gz \
       --reads cleandata/${base}_val_2.fq.gz \
       --fastx  --aligned rRNA_removed_data/$base/$base.rRNA \
       --other rRNA_removed_data/$base.non.rRNA \
       --threads 24
     echo "processing of ${base} is done"
   fi
 done &
}

export -f remove_rRNA


#nohup time parallel -j 2 --eta --load 100% --noswap  remove_rRNA ::: $(ls cleandata/*_val_1.fq.gz |cut -d "/" -f 2) &


#sub-sample the reads （seqtk）https://cloud.tencent.com/developer/article/1674827
#############################################################################
  seqtk sample -s 100 read1.fq.gz 10000 | gzip > sub1.fq.gz
  seqtk sample -s 100 read2.fq.gz 10000 | gzip > sub2.fq.gz
