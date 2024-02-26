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
     trim_galore -q 25 --phred33 --stringency 3 --length 40 \
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
 for name in $(ls cleandata/*_val_1.fq.gz |cut -d "/" -f 2); do
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
#

#calculate the sequence number
#############################################################################
rm -f reads_num
for i in cleandata/*.gz; do
   num=$(cat $i |pigz -d | grep -c '^+$') 
   echo -e "${i}\t${num}" >> reads_num
done

#sub-sample the reads （seqtk）https://cloud.tencent.com/developer/article/1674827
#############################################################################
for i in cleandata/*_val_1.fq.gz; do
  base=$(basename $i _val_1.fq.gz)
  seed=$(echo $RANDOM % 100 + 1 | bc)
  seqtk sample -s $seed cleandata/${base}_val_1.fq.gz 6000000 | pigz > subsample_data/${base}.subsample_R1.fq.gz
  seqtk sample -s $seed cleandata/${base}_val_2.fq.gz 6000000 | pigz > subsample_data/${base}.subsample_R2.fq.gz
  echo "${base} finished"
done



#reference preparation
bowtie2-build --threads 4  Suicot1_AssemblyScaffolds_Repeatmasked.fasta Suicot_genome

#mapping
##################################################################

function run_alignment () {

file_folder=$(echo $work_dir/$file_folder)
dir=$work_dir
Suicot_genome="/home/microbiome/data_storage/SATA2/RNA_data/genome_reference/Suicot_genome"
gtf_suicot="/home/microbiome/data_storage/SATA2/RNA_data/genome_reference/Suicot1_GeneCatalog_20171209.gtf"

fwd=$1
rev=$2

mkdir $dir/3_suillus_alignment
mkdir $dir/3_suillus_alignment/1_suillus_aligned_fastq
mkdir $dir/3_suillus_alignment/1_bowtie2_met_file
mkdir $dir/3_suillus_alignment/1_bowtie2_log_file
mkdir $dir/3_suillus_alignment/2_bam_flagstat_file
mkdir $dir/3_suillus_alignment/3_suillus_count_file
mkdir $dir/3_suillus_alignment/2_suillus_bam_file
mkdir $dir/3_suillus_alignment/1_suillus_unaligned_fastq

  if [ -f "$dir/3_suillus_alignment/3_suillus_count_file/${fwd}_suillus_gene_id_count.txt" ]; then
      echo "${fwd} has been analyzed"
  else
     mkdir $dir/3_suillus_alignment/${fwd}.temp
     time bowtie2 -p 24 -x $Suicot_genome \
                  -1 $file_folder/$fwd \
                  -2 $file_folder/$rev \
                  -S $dir/3_suillus_alignment/${fwd}.temp/${fwd}_suillus.sam \
                  --al-conc-gz $dir/3_suillus_alignment/1_suillus_aligned_fastq/${fwd}_aligned.fastq.gz \
                  --un-conc-gz $dir/3_suillus_alignment/1_suillus_unaligned_fastq/${fwd}_unaligned.fastq.gz \
                  --met-file $dir/3_suillus_alignment/1_bowtie2_met_file/${fwd}_met.txt \
                  2>$dir/3_suillus_alignment/1_bowtie2_log_file/${fwd}_bowtie2.log

     samtools sort -o bam -@ 3 -o $dir/3_suillus_alignment/2_suillus_bam_file/${fwd}_suillus.bam $dir/3_suillus_alignment/${fwd}.temp/${fwd}_suillus.sam
     samtools flagstat -@ 3 $dir/3_suillus_alignment/${fwd}.temp/${fwd}_suillus.sam > $dir/3_suillus_alignment/2_bam_flagstat_file/${fwd}.flagstat
     featureCounts -t exon -F GTF -g gene_id -T 4 -a $gtf_suicot -o $dir/3_suillus_alignment/3_suillus_count_file/${fwd}_suillus_gene_id_count.txt \
                    $dir/3_suillus_alignment/2_suillus_bam_file/${fwd}_suillus.bam
     rm -rf $dir/3_suillus_alignment/${fwd}.temp
     mv $dir/3_suillus_alignment/1_suillus_aligned_fastq/${fwd}_aligned.fastq.1.gz $dir/3_suillus_alignment/1_suillus_aligned_fastq/${fwd}_aligned_R1.fastq.gz
     mv $dir/3_suillus_alignment/1_suillus_aligned_fastq/${fwd}_aligned.fastq.2.gz $dir/3_suillus_alignment/1_suillus_aligned_fastq/${fwd}_aligned_R2.fastq.gz
  fi
}
export -f run_alignment

work_dir="/home/microbiome/data_storage/SATA3/RNA_data/FISH_RNA/SuiPinnscriptome_6_download/SuiPinnscriptome_6"
file_folder="rRNA_removed_data"

ls *.fq.gz | cut -d "." -f 1,2,3,4 |sort -u |while read id
 do
   echo $id
   fwd=$(echo ${id}.non.rRNA_fwd.fq.gz)
   rev=$(echo ${id}.non.rRNA_rev.fq.gz)
   echo $fwd $rev
   run_alignment $fwd $rev
done
  



dir="/home/microbiome/data_storage/SATA2/RNA_data/DOE_EXP3"
gtf_suicot="/home/microbiome/data_storage/SATA2/RNA_data/genome_reference/Suicot1_GeneCatalog_20171209.gtf"
gtf_suicot="$dir/reference/Suicot1_all_genes_20171209.gff"
#featureCounts -t exon -F GFF -f -p -g name -M -T 4 -a $gtf_suicot -o $dir/3_suillus_alignment/3_suillus_count_file/suillus_all_gene_id_count.txt *.bam  1>$dir/3_suillus_alignment/3_suillus_count_file/counts.all_gene_id.log 2>&1
#featureCounts -t exon -F GFF -f -p -O -g name -M -T 4 -a $gtf_suicot -o $dir/3_suillus_alignment/3_suillus_count_file/suillus_all_gene_id_count_overlap.txt *.bam  1>$dir/3_suillus_alignment/3_suillus_count_file/counts.all_gene_id_overlap.log 2>&1
featureCounts -t exon -F GTF -g gene_id -p -M -T 24 -a $gtf_suicot -o $dir/3_suillus_alignment/3_suillus_count_file/suillus_catalog_gene_id_count.txt *.bam  1>$dir/3_suillus_alignment/3_suillus_count_file/counts.catalog_gene_id.log 2>&1
#featureCounts -t exon -F GTF -g gene_id -p -O -M -T 4 -a $gtf_suicot -o $dir/3_suillus_alignment/3_suillus_count_file/suillus_catalog_gene_id_count_overlap.txt *.bam  1>$dir/3_suillus_alignment/3_suillus_count_file/counts.catalog_gene_id_overlap.log 2>&1


























