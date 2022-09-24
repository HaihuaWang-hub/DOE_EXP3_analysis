####################################################
#This is for metatranscriptomic filter analysis of plant_emf_sap_interaction project
####################################################
#!/bin/bash
dir=/Volumes/T7/plant_emf_sap_interaction

cd $dir
mkdir $dir/rRNA_db
mkdir $dir/rRNA_filter
mkdir $dir/filtered_data
mkdir $dir/fungi_rRNA
mkdir $dir/bacteria_rRNA
mkdir $dir/archaea_rRNA

#  section 1: prepare database
####################################################
cd $dir/rRNA_db
#download the database
wget -c https://rdp.cme.msu.edu/download/current_Bacteria_unaligned.fa.gz
wget -c https://rdp.cme.msu.edu/download/current_Archaea_unaligned.fa.gz
wget -c https://rdp.cme.msu.edu/download/current_Fungi_unaligned.fa.gz








################################################
#####extract fungal rRNA LSU D1D2 regin
################################################
dir=/Volumes/T7/plant_emf_sap_interaction

mkdir $dir/fungal_rRNA_D1D2
bowtie2-build --threads 4 /Volumes/T7/plant_emf_sap_interaction/rRNA_db/Step5_D1D2_29342_plus_PMI93_82_plus_Scoth_Tr210_database29539.fasta /Volumes/T7/plant_emf_sap_interaction/rRNA_db/Fungi_rRNA_D1D2
cd $dir/cleandata


fungi_db="/Volumes/T7/plant_emf_sap_interaction/rRNA_db/Fungi_rRNA_D1D2"

ls *.gz|cut -d"_" -f 1,2,3 |sort -u |while read id;do
time bowtie2 -p 4 -x $fungi_db \
-1 ${id}_paired_R1.fq.gz \
-2 ${id}_paired_R2.fq.gz \
-S $dir/fungal_rRNA_D1D2/${id}_fungi_rRNA.sam \
--al-conc-gz $dir/fungal_rRNA_D1D2/${id}_fungi_rRNA.fastq.gz \
--met-file $dir/fungal_rRNA_D1D2/${id}_met.txt \
2>$dir/fungal_rRNA_D1D2/${id}_bowtie2.log
rm $dir/fungal_rRNA_D1D2/${id}_fungi_rRNA.sam
mv $dir/fungal_rRNA_D1D2/${id}_fungi_rRNA.fastq.1.gz $dir/fungal_rRNA_D1D2/${id}_fungi_rRNA_R1.fastq.gz
mv $dir/fungal_rRNA_D1D2/${id}_fungi_rRNA.fastq.2.gz $dir/fungal_rRNA_D1D2/${id}_fungi_rRNA_R2.fastq.gz
done





################################################
#####extract fungal rRNA ITS regin
################################################
dir=/Volumes/T7/plant_emf_sap_interaction

mkdir $dir/2_rRNA_fungal_ITS
bowtie2-build unite_its_plus_suillus_cothurnatus.fasta unite --threads 4
cd $dir/1_cleandata


fungi_db="/Volumes/T7/plant_emf_sap_interaction/database_rRNA_db/unite"

ls *.gz|cut -d"_" -f 1,2,3 |sort -u |tail -n 53|while read id;do
time bowtie2 -p 4 -x $fungi_db \
-1 ${id}_paired_R1.fq.gz \
-2 ${id}_paired_R2.fq.gz \
-S $dir/2_rRNA_fungal_ITS/${id}_fungi_rRNA.sam \
--al-conc-gz $dir/2_rRNA_fungal_ITS/${id}_fungi_rRNA.fastq.gz \
--met-file $dir/2_rRNA_fungal_ITS/${id}_met.txt \
2>$dir/2_rRNA_fungal_ITS/${id}_bowtie2.log
rm $dir/2_rRNA_fungal_ITS/${id}_fungi_rRNA.sam
mv $dir/2_rRNA_fungal_ITS/${id}_fungi_rRNA.fastq.1.gz $dir/2_rRNA_fungal_ITS/${id}_fungi_rRNA_R1.fastq.gz
mv $dir/2_rRNA_fungal_ITS/${id}_fungi_rRNA.fastq.2.gz $dir/2_rRNA_fungal_ITS/${id}_fungi_rRNA_R2.fastq.gz
done




time bowtie2 -p 4 -x $fungi_db \
-1 merged_R1.fastq.gz \
-2 merged_R2.fastq.gz \
-S merged.sam \
2>merged_bowtie2.log
ls *.sam|while read id ;do (samtools sort -o bam -@ 8 -o $(basename ${id} ".sam").bam ${id});done
rm merged.sam





################################################
#####extract bacteria 16S rRNA
################################################
dir=/Volumes/T7/plant_emf_sap_interaction
mkdir $dir/2_bacteria+Archae
bowtie2-build --threads 4 $dir/database_rRNA_db/Bacteria_unaligned.fa $dir/database_rRNA_db/Bacteria_unaligned


bacteria_db="$dir/database_rRNA_db/Bacteria_unaligned"
ls *.gz|cut -d"_" -f 1,2,3,4,5,6 |sort -u |while read id;do

bowtie2 -p 4 -x $bacteria_db \
-1 ${id}_R1.fastq.gz \
-2 ${id}_R2.fastq.gz \
-S $dir/2_bacteria+Archae/${id}_bacteria_rRNA.sam \
--al-conc-gz $dir/2_bacteria+Archae/${id}_bacteria+Archae_rRNA.fastq.gz \
--met-file $dir/2_bacteria+Archae/${id}_met.txt \
2>$dir/2_bacteria+Archae/${id}_bowtie2.log
rm $dir/2_bacteria+Archae/${id}_bacteria_rRNA.sam
done


















#extract the rRNA sequences https://github.com/biocore/sortmerna
############################################################################################
mkdir 1_pinus_unaligned_fastq_nonrRNA
extract_rRNA (){
  name=$1
  base=$(basename $name _R1.fastq.gz)
  if [ ! -f 1_pinus_unaligned_fastq_nonrRNA/$base.rRNA_fwd.fq.gz ]; then
      mkdir 1_pinus_unaligned_fastq_nonrRNA/$base
      echo "start to process ${base}" 
     sortmerna \
       --workdir 1_pinus_unaligned_fastq_nonrRNA/$base \
       --ref /home/microbiome/tools/sortmerna/data/rRNA_databases/rfam-5.8s-database-id98.fasta \
       --ref /home/microbiome/tools/sortmerna/data/rRNA_databases/rfam-5s-database-id98.fasta \
       --ref /home/microbiome/tools/sortmerna/data/rRNA_databases/silva-arc-16s-id95.fasta \
       --ref /home/microbiome/tools/sortmerna/data/rRNA_databases/silva-arc-23s-id98.fasta \
       --ref /home/microbiome/tools/sortmerna/data/rRNA_databases/silva-bac-16s-id90.fasta \
       --ref /home/microbiome/tools/sortmerna/data/rRNA_databases/silva-bac-23s-id98.fasta \
       --ref /home/microbiome/tools/sortmerna/data/rRNA_databases/silva-euk-18s-id95.fasta \
       --ref /home/microbiome/tools/sortmerna/data/rRNA_databases/silva-euk-28s-id98.fasta \
       --paired_out --out2  --zip-out \
       --reads 1_pinus_unaligned_fastq/${base}_R1.fastq.gz \
       --reads 1_pinus_unaligned_fastq/${base}_R2.fastq.gz \
       --fastx  --other 1_pinus_unaligned_fastq_nonrRNA/$base.nonrRNA \
       --threads 24
     echo "processing of ${base} is done"
   else
     echo "${base} has been anlyzed"
   fi
}

export -f extract_rRNA


nohup time parallel -j 2 --eta --load 100% --noswap  extract_rRNA ::: $(cat 1_pinus_unaligned_fastq/file_list) &










































































































file_list.txt*
AA232R_S61_L003_paired_met.txt*              file_list_1.5.txt*
AA233R_S62_L003_paired_bowtie2.log*          file_list_1.txt*
AA233R_S62_L003_paired_fungi_rRNA_R1.fastq*  file_list_2.5.txt*
AA233R_S62_L003_paired_fungi_rRNA_R2.fastq*  file_list_2.txt*
AA233R_S62_L003_paired_met.txt*              file_list_3.5.txt*
AA240R_S63_L003_paired_bowtie2.log*          file_list_3.txt


time fastuniq -i file_list_1.5.txt  -t q -o fungi_rRNA_1_R1.fastq -p fungi_rRNA_1_R2.fastq
time fastuniq -i file_list_1.txt  -t q -o fungi_rRNA_2_R1.fastq -p fungi_rRNA_2_R2.fastq
time fastuniq -i file_list_2.5.txt  -t q -o fungi_rRNA_3_R1.fastq -p fungi_rRNA_3_R2.fastq
time fastuniq -i file_list_2.txt  -t q -o fungi_rRNA_4_R1.fastq -p fungi_rRNA_4_R2.fastq
time fastuniq -i file_list_3.5.txt  -t q -o fungi_rRNA_5_R1.fastq -p fungi_rRNA_5_R2.fastq
time fastuniq -i file_list_3.txt  -t q -o fungi_rRNA_6_R1.fastq -p fungi_rRNA_6_R2.fastq










time trinity --CPU 2 \
    --seqType fq \
    --left B_1_R1_val_1.fq.gz,B_2_R1_val_1.fq.gz,control_1_R1_val_1.fq.gz,control_2_R1_val_1.fq.gz,D_1_R1_val_1.fq.gz,D_2_R1_val_1.fq.gz   \
    --right   B_1_R2_val_2.fq.gz,B_2_R2_val_2.fq.gz,control_1_R2_val_2.fq.gz,control_2_R2_val_2.fq.gz,D_1_R2_val_2.fq.gz,D_2_R2_val_2.fq.gz \
    --max_memory 200G \
    --output ../trinity_output

















