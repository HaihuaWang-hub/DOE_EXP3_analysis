
#Apply MicroFisher to classify the fungal community
for file in $(ls *_1.fq.gz); do 
   i=$(echo $file |cut -d "_" -f 1)
   echo $i
   workspace=$(echo microfisher/$i)
   mkdir $workspace
   MicroFisher preset --preset_db ITS+LSU \
   --db_path /home/microbiome/data_storage/SATA3/Fisher_test/short_DBs/MicroFisher_DBs     \
   --workspace ./  \
   --paired ${i}_val_1.fq.gz ${i}_val_2.fq.gz \
   --out_dir microfisher \
   --out_prefix $i \
   --threads 24; 
done






















#######################################################
##fastqc to show the duplications
#######################################################
fastqc 
multiqc

#######################################################
##fastuniq remove duplicates 
#######################################################

#merge paired end fastq with seqtk
ls *.gz |cut -d"_" -f 1,2,3,4,5,6|sort -u |while read id; do
/home/wang/genome_tools/merge-paired-reads.sh ${id}_R1.fastq.gz ${id}_R2.fastq.gz  ${id}.fastq.gz
done

 ls *.gz |cut -d"_" -f 1,2,3,4,5,6 |while read id; do
 seqkit rmdup --by-seq  --threads 5 ${id}.fastq.gz > dup_removed_data/${id}_dup_removed_fastq
 done


 seqkit rmdup --by-seq --threads 5 AA376R_S97_L003_paired_fungi_rRNA.fastq.gz > AA376R_S97_L003_paired_fungi_rRNA_dup_removed.fastq.gz


ls *.gz |cut -d"_" -f 1,2,3,4,5,6 |while read id; do
mv ${id}_dup_removed_fastq.gz ${id}_dup_removed.fastq
done

ls *.fastq |cut -d"." -f 1 |while read id; do
/home/wang/genome_tools/unmerge-paired-reads.sh ${id}.fastq \
${id}_R1.fastq ${id}_R2.fastq
done


fastuniq -i file_list -t q -o fungi_rRNA_R1.fastq -p fungi_rRNA_R2.fastq -c 0





Trinity --CPU 6 \
    --seqType fq \
    --left  dup_removedmerged_paried_R1.fastq   \
    --right dup_removedmerged_paried_R1.fastq  \
    --max_memory 15G \
    --output Trinity_output_ref_free_2 &
    
    
    Trinity --CPU 6 \
    --seqType fq \
    --samples_file file_list  \
    --max_memory 15G \
    --output ../trinity_output
    
    
    
/home/wang/miniconda3/envs/genomic/opt/trinity-2.12.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl \
  Trinity.fasta > trinity_unigene.fasta
  
  




#estimate abundance with fungal_D1D2 files
######################################################
mkdir /home/wang/usb/fungal_rRNA_D1D2/alignment_dir/
align_and_estimate_abundance.pl \
--transcripts /home/wang/usb/fungal_rRNA_D1D2/trinity_output/trinity_unigene.fasta \
--seqType fq \
--samples_file file_list \
--est_method RSEM \
--aln_method bowtie2 \
--trinity_mode \
--prep_reference  \
--output_dir /home/wang/usb/fungal_rRNA_D1D2/alignment_dir/ \
--thread_count 6


find * -name '*.isoforms.results'> quant.file 
abundance_estimates_to_matrix.pl --est_method RSEM \
--gene_trans_map /home/wang/usb/fungal_rRNA_D1D2/trinity_output/Trinity.fasta.gene_trans_map \
--name_sample_by_basedir \
--quant_files quant.file



#estimate abundance with fungal files
######################################################
align_and_estimate_abundance.pl \
--transcripts /home/wang/usb/fungi_rRNA/fungal_trinity_unigene.fasta \
--seqType fq \
--samples_file file_list \
--est_method RSEM \
--aln_method bowtie2 \
--trinity_mode \
--prep_reference  \
--output_dir /home/wang/usb/fungi_rRNA/alignment_dir/ \
--thread_count 3

find * -name '*.isoforms.results'> quant.file 
abundance_estimates_to_matrix.pl --est_method RSEM \
--gene_trans_map /home/wang/usb/fungi_rRNA/trinity_output/Trinity.fasta.gene_trans_map \
--name_sample_by_basedir \
--quant_files quant.file





#estimate abundance with bacteria files
######################################################
ls *.gz |cut -d"_" -f 1,2,3,4,5,6,7|sort -u > 1
ls *_R1.fastq.gz > 2
ls *_R2.fastq.gz > 3
paste 1 1 2 3 > file_list
align_and_estimate_abundance.pl \
--transcripts /home/wang/usb/bacteria_rRNA/bacteria_trinity_unigene.fasta \
--seqType fq \
--samples_file file_list \
--est_method RSEM \
--aln_method bowtie2 \
--trinity_mode \
--prep_reference  \
--output_dir /home/wang/usb/bacteria_rRNA/alignment_dir/ \
--thread_count 3

find * -name '*.isoforms.results'> quant.file 
abundance_estimates_to_matrix.pl --est_method RSEM \
--gene_trans_map /home/wang/usb/bacteria_rRNA/trinity_output/Trinity.fasta.gene_trans_map \
--name_sample_by_basedir \
--quant_files quant.file


#estimate abundance with Archae files
######################################################
ls *.gz |cut -d"_" -f 1,2,3,4,5,6,7,8|sort -u > 1
ls *_R1.fastq.gz > 2
ls *_R2.fastq.gz > 3
paste 1 1 2 3 > file_list
align_and_estimate_abundance.pl \
--transcripts /home/wang/usb/Archae_rRNA/Archae_trinity_unigene.fasta \
--seqType fq \
--samples_file file_list \
--est_method RSEM \
--aln_method bowtie2 \
--trinity_mode \
--prep_reference  \
--output_dir /home/wang/usb/Archae_rRNA/alignment_dir/ \
--thread_count 3

find * -name '*.isoforms.results'> quant.file 
abundance_estimates_to_matrix.pl --est_method RSEM \
--gene_trans_map /home/wang/usb/Archae_rRNA/trinity_output/Trinity.fasta.gene_trans_map \
--name_sample_by_basedir \
--quant_files quant.file









grep -r "Boletales" fixrank_fungal_LSU_D1D2.fasta_classified.txt > Boletales.txt
grep -r "Hysterangiales" fixrank_fungal_LSU_D1D2.fasta_classified.txt > Hysterangiales.txt
grep -r "Atractiellales" fixrank_fungal_LSU_D1D2.fasta_classified.txt > Atractiellales.txt
cat Boletales.txt |cut -d ";" -f 1 > Boletales_id.txt
cat Hysterangiales.txt |cut -d ";" -f 1 > Hysterangiales_id.txt
cat Atractiellales.txt |cut -d ";" -f 1 > Atractiellales_id.txt
/home/wang/genome_tools/faSomeRecords trinity_unigene.fasta Boletales_id.txt Boletales_seq.fasta
/home/wang/genome_tools/faSomeRecords trinity_unigene.fasta Hysterangiales_id.txt Hysterangiales_seq.fasta
/home/wang/genome_tools/faSomeRecords trinity_unigene.fasta Atractiellales_id.txt Atractiellales_seq.fasta




#re-assembly
###############################################################
dir="/home/wang/usb/fungal_rRNA_D1D2/first_3_order"
mkdir $dir/database
mkdir $dir/Boletales_assembly_files
mkdir $dir/Atractiellales_assembly_files
mkdir $dir/Hysterangiales_assembly_files
mkdir $dir/Boletales_assemblies
mkdir $dir/Atractiellales_assemblies
mkdir $dir/Hysterangiales_assemblies


bowtie2-build --threads 6  $dir/Hysterangiales_seq.fasta $dir/database/Hysterangiales_seq
bowtie2-build --threads 6  $dir/Boletales_seq.fasta $dir/database/Boletales_seq
bowtie2-build --threads 6  $dir/Atractiellales_seq.fasta $dir/database/Atractiellales_seq


cd /home/wang/usb/fungal_rRNA_D1D2/rawdata
#########Boletales
ls *.gz |cut -d"_" -f 1,2,3 |sort -u |while read id; do
mkdir $dir/Boletales_assembly_files/${id}
bowtie2 -p 5 -x $dir/database/Boletales_seq \
-1 ${id}_fungi_rRNA_R1.fastq.gz \
-2 ${id}_fungi_rRNA_R2.fastq.gz \
-S $dir/Boletales_assembly_files/${id}/${id}.sam \
--al-conc-gz $dir/Boletales_assembly_files/${id}/${id}.fastq.gz \
2>$dir/Boletales_assembly_files/${id}/${id}.bowtie2.log

mv $dir/Boletales_assembly_files/${id}/${id}.fastq.1.gz $dir/Boletales_assembly_files/${id}/${id}_1.fastq.gz
mv $dir/Boletales_assembly_files/${id}/${id}.fastq.2.gz $dir/Boletales_assembly_files/${id}/${id}_2.fastq.gz
rm $dir/Boletales_assembly_files/${id}/${id}.sam

Trinity --CPU 6 \
    --seqType fq \
    --left  $dir/Boletales_assembly_files/${id}/${id}_1.fastq.gz  \
    --right $dir/Boletales_assembly_files/${id}/${id}_2.fastq.gz  \
    --max_memory 15G \
    --output $dir/Boletales_assembly_files/${id}/Trinity_output

/home/wang/miniconda3/envs/genomic/opt/trinity-2.12.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl \
  $dir/Boletales_assembly_files/${id}/Trinity_output/Trinity.fasta > $dir/Boletales_assemblies/${id}_trinity_unigene.fasta
  
mv $dir/Boletales_assembly_files/${id}/Trinity_output/Trinity.fasta $dir/Boletales_assemblies/${id}_Trinity.fasta

done



######Atractiellales
dir="/home/wang/usb/fungal_rRNA_D1D2/first_3_order"
cd /home/wang/usb/fungal_rRNA_D1D2/rawdata
ls *.gz |cut -d"_" -f 1,2,3 |sort -u |while read id; do
mkdir $dir/Atractiellales_assembly_files/${id}

bowtie2 -p 5 -x $dir/database/Atractiellales_seq \
-1 ${id}_fungi_rRNA_R1.fastq.gz \
-2 ${id}_fungi_rRNA_R2.fastq.gz \
-S $dir/Atractiellales_assembly_files/${id}/${id}.sam \
--al-conc-gz $dir/Atractiellales_assembly_files/${id}/${id}.fastq.gz \
2>$dir/Atractiellales_assembly_files/${id}/${id}.bowtie2.log

mv $dir/Atractiellales_assembly_files/${id}/${id}.fastq.1.gz $dir/Atractiellales_assembly_files/${id}/${id}_1.fastq.gz
mv $dir/Atractiellales_assembly_files/${id}/${id}.fastq.2.gz $dir/Atractiellales_assembly_files/${id}/${id}_2.fastq.gz
rm $dir/Atractiellales_assembly_files/${id}/${id}.sam

Trinity --CPU 3 \
    --seqType fq \
    --left  $dir/Atractiellales_assembly_files/${id}/${id}_1.fastq.gz  \
    --right $dir/Atractiellales_assembly_files/${id}/${id}_2.fastq.gz  \
    --max_memory 15G \
    --output $dir/Atractiellales_assembly_files/${id}/Trinity_output

/home/wang/miniconda3/envs/genomic/opt/trinity-2.12.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl \
  $dir/Atractiellales_assembly_files/${id}/Trinity_output/Trinity.fasta > $dir/Atractiellales_assemblies/${id}_trinity_unigene.fasta
  
mv $dir/Atractiellales_assembly_files/${id}/Trinity_output/Trinity.fasta $dir/Atractiellales_assemblies/${id}_Trinity.fasta

done






######Hysterangiales
dir="/home/wang/usb/fungal_rRNA_D1D2/first_3_order"
cd /home/wang/usb/fungal_rRNA_D1D2/rawdata
ls *.gz |cut -d"_" -f 1,2,3 |sort -u |while read id; do
mkdir $dir/Hysterangiales_assembly_files/${id}

bowtie2 -p 5 -x $dir/database/Hysterangiales_seq \
-1 ${id}_fungi_rRNA_R1.fastq.gz \
-2 ${id}_fungi_rRNA_R2.fastq.gz \
-S $dir/Hysterangiales_assembly_files/${id}/${id}.sam \
--al-conc-gz $dir/Hysterangiales_assembly_files/${id}/${id}.fastq.gz \
2>$dir/Hysterangiales_assembly_files/${id}/${id}.bowtie2.log

mv $dir/Hysterangiales_assembly_files/${id}/${id}.fastq.1.gz $dir/Hysterangiales_assembly_files/${id}/${id}_1.fastq.gz
mv $dir/Hysterangiales_assembly_files/${id}/${id}.fastq.2.gz $dir/Hysterangiales_assembly_files/${id}/${id}_2.fastq.gz
rm $dir/Hysterangiales_assembly_files/${id}/${id}.sam

Trinity --CPU 6 \
    --seqType fq \
    --left  $dir/Hysterangiales_assembly_files/${id}/${id}_1.fastq.gz  \
    --right $dir/Hysterangiales_assembly_files/${id}/${id}_2.fastq.gz  \
    --max_memory 15G \
    --output $dir/Hysterangiales_assembly_files/${id}/Trinity_output

/home/wang/miniconda3/envs/genomic/opt/trinity-2.12.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl \
  $dir/Hysterangiales_assembly_files/${id}/Trinity_output/Trinity.fasta > $dir/Hysterangiales_assemblies/${id}_trinity_unigene.fasta
  
mv $dir/Hysterangiales_assembly_files/${id}/Trinity_output/Trinity.fasta $dir/Hysterangiales_assemblies/${id}_Trinity.fasta

done


#################
#rename seq name with file name
##################
ls *_unigene.fasta|while read id;do
awk '/^>/{print ">" substr(FILENAME,1,length(FILENAME)-6); next} 1' ${id} > ${id}.renamed.fasta
done 

cat *.renamed.fasta > @@@_assembly_seqs.fasta

seqkit rename --by-name 00_assembly_seqs.fasta > assembly_seqs.fasta

grep -r "Suillus;100%" allrank_assembly_seqs.fasta_classified.txt |cut -d ";" -f 1

grep -r "Suillus;100%" allrank_assembly_seqs.fasta_classified.txt |cut -d ";" -f 1 > selected_suillus_gene.txt


/home/wang/genome_tools/faSomeRecords assembly_seqs.fasta selected_suillus_gene.txt selected_suillus_gene.fasta























/home/wang/genome_tools/merge-paired-reads.sh merged_R1.fastq merged_R2.fastq  merged_paried.fastq

 seqkit rmdup --by-seq  --threads 8 merged_paried.fastq > dup_removedmerged_paried.fastq

/home/wang/genome_tools/unmerge-paired-reads.sh dup_removedmerged_paried.fastq dup_removedmerged_paried_R1.fastq dup_removedmerged_paried_R2.fastq

bowtie2-build --threads 6 Step5_D1D2_29342_plus_PMI93_82_plus_Scoth_Tr210_database29539.fasta fungal_rRNA_D1D2

ls *.gz |cut -d"_" -f 1,2,3 |while read id; do
bowtie2 -p 10 -x /home/wang/usb/fungal_rRNA_D1D2/database/fungal_rRNA_D1D2 \
-1 dup_removedmerged_paried_R1.fastq.gz \
-2 dup_removedmerged_paried_R2.fastq.gz \
-S dup_removed_merged.sam 
done


ls *.sam|while read id ;do (samtools sort -o bam -@ 8 -o $(basename ${id} ".sam").bam ${id});done

Trinity --genome_guided_bam dup_removed_merged.bam \
        --max_memory 50G \
        --genome_guided_max_intron 10000 \
        --genome_guided_min_coverage  300 \
        --CPU 8 --output Trinity_output_merged_00

/home/wang/miniconda3/envs/genomic/opt/trinity-2.12.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl Trinity-GG.fasta > Trinity-GG_unigene.fasta




mkdir /home/wang/usb/fungal_rRNA_D1D2/alignment_dir_ref_guide_assembly/


align_and_estimate_abundance.pl \
--transcripts /home/wang/usb/fungal_rRNA_D1D2/ref_guide_assembly/Trinity_output_merged/Trinity-GG_unigene.fasta \
--seqType fq \
--samples_file file_list \
--est_method RSEM \
--aln_method bowtie2 \
--trinity_mode \
--prep_reference  \
--output_dir /home/wang/usb/fungal_rRNA_D1D2/alignment_dir_ref_guide_assembly/ \
--thread_count 3





find * -name '*.isoforms.results'> quant.file 

abundance_estimates_to_matrix.pl --est_method RSEM \
--gene_trans_map /home/wang/usb/fungal_rRNA_D1D2/ref_guide_assembly/Trinity_output_merged/Trinity-GG.fasta.gene_trans_map \
--name_sample_by_basedir \
--quant_files quant.file





##########################################################################################
# de novo assembly for each sample
##########################################################################################
dir=/Volumes/T7/plant_emf_sap_interaction

mkdir $dir/2_rRNA_fungal_D1D2/denovo_assembly_for_each_sample

cd $dir/2_rRNA_fungal_D1D2/rawdata
ls $dir/2_rRNA_fungal_D1D2/rawdata/*.gz |cut -d"/" -f 7 |cut -d"_" -f 1,2,3,4,5 |sort -u|while read id; do
/usr/local/bin/trinityrnaseq/Trinity --CPU 6 \
    --seqType fq \
    --left  $dir/2_rRNA_fungal_D1D2/rawdata/${id}_R1.fastq.gz  \
    --right $dir/2_rRNA_fungal_D1D2/rawdata/${id}_R2.fastq.gz  \
    --max_memory 1G \
    --output $dir/2_rRNA_fungal_D1D2/denovo_assembly_for_each_sample/${id}_Trinity_output
    done
/home/wang/miniconda3/envs/genomic/opt/trinity-2.12.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl \
$dir/2_rRNA_fungal_D1D2/denovo_assembly_for_each_sample/${id}_Trinity_output/Trinity.fasta > $dir/2_rRNA_fungal_D1D2/denovo_assembly_for_each_sample/${id}_Trinity_unigene.fasta

done












































