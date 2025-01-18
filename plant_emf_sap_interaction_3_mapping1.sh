####################################################
#This is for metatranscriptomic mapping analysis of plant_emf_sap_interaction project
####################################################
#!/bin/bash
conda activate RNASeq
dir=


#  section 1: reference preparation
####################################################
cd $dir/reference

gunzip *.gz
gffread Suicot1_GeneCatalog_20171209.gff3 -T -o Suicot1_GeneCatalog_20171209.gtf

bowtie2-build --threads 4  Suicot1_AssemblyScaffolds_Repeatmasked.fasta Suicot_genome
bowtie2-build --threads 4 TrPtA_269336_P_taeda_mRNAdatabase_328662.fasta Pintaeda_genome


#3.mapping
##################################################################
Suicot_genome="/data_storage/SATA2/haihua/RNA_data/genome_reference/Suicot_genome"
Pintaeda_genome="/home/microbiome/data_storage/SATA2/plant_genome/pita/index/pita"
gtf_suicot="/data_storage/SATA2/haihua/RNA_data/genome_reference/Suicot1_GeneCatalog_20171209.gtf"
gtf_pintaeda="/home/microbiome/data_storage/SATA2/plant_genome/pita/Pita.2_01.gtf"



#  section2: mapping to suillus
################################################################
dir=#Path to file
mkdir $dir/3_suillus_alignment
mkdir $dir/3_suillus_alignment/1_suillus_aligned_fastq
mkdir $dir/3_suillus_alignment/1_bowtie2_met_file
mkdir $dir/3_suillus_alignment/1_bowtie2_log_file
mkdir $dir/3_suillus_alignment/2_bam_flagstat_file
mkdir $dir/3_suillus_alignment/3_suillus_count_file
mkdir $dir/3_suillus_alignment/2_suillus_bam_file
mkdir $dir/3_suillus_alignment/1_suillus_unaligned_fastq

#bowtie2-build --threads 4  $dir/reference/Suicot1_AssemblyScaffolds_Repeatmasked.fasta $dir/reference/Suicot_genome
#gffread $dir/reference/Suicot1_GeneCatalog_20171209.gff3 -T -o $dir/reference/Suicot1_GeneCatalog_20171209.gtf
#Suicot_genome="$dir/reference/Suicot_genome"
#gtf_suicot="$dir/reference/Suicot1_GeneCatalog_20171209.gtf"



cd $dir/cleandata
ls *.gz|cut -d"_" -f 1 |sort -u  |while read id;do
  echo $id
  fwd="cleandata_rmdup/${id}_rmdup_val_1.fq.gz"
  rev="cleandata_rmdup/${id}_rmdup_val_2.fq.gz"
  if [ -f "$dir/3_suillus_alignment/3_suillus_count_file/${id}_suillus_gene_id_count.txt" ]; then
      echo "${id} has been analyzed"
      else
      mkdir $dir/3_suillus_alignment/${id}.temp
      time bowtie2 -p 24 -x $Suicot_genome \
           -1 $fwd \
           -2 $rev \
           -S $dir/3_suillus_alignment/${id}.temp/${id}_suillus.sam \
           --al-conc-gz $dir/3_suillus_alignment/1_suillus_aligned_fastq/${id}_aligned.fastq.gz \
           --un-conc-gz $dir/3_suillus_alignment/1_suillus_unaligned_fastq/${id}_unaligned.fastq.gz \
           --met-file $dir/3_suillus_alignment/1_bowtie2_met_file/${id}_met.txt \
           2>$dir/3_suillus_alignment/1_bowtie2_log_file/${id}_bowtie2.log

      samtools sort -o bam -@ 3 -o $dir/3_suillus_alignment/2_suillus_bam_file/${id}_suillus.bam $dir/3_suillus_alignment/${id}.temp/${id}_suillus.sam
      samtools flagstat -@ 3 $dir/3_suillus_alignment/${id}.temp/${id}_suillus.sam > $dir/3_suillus_alignment/2_bam_flagstat_file/${id}.flagstat
      featureCounts -t exon -F GTF -g gene_id -T 4 -a $gtf_suicot -o $dir/3_suillus_alignment/3_suillus_count_file/${id}_suillus_gene_id_count.txt $dir/3_suillus_alignment/2_suillus_bam_file/${id}_suillus.bam  1>$dir/3_suillus_alignment/3_suillus_count_file/${id}_suillus.log 2>&1
      rm -rf $dir/3_suillus_alignment/${id}.temp
      mv -f $dir/3_suillus_alignment/1_suillus_aligned_fastq/${id}_aligned.fastq.1.gz $dir/3_suillus_alignment/1_suillus_aligned_fastq/${id}_aligned_R1.fastq.gz
      mv -f $dir/3_suillus_alignment/1_suillus_aligned_fastq/${id}_aligned.fastq.2.gz $dir/3_suillus_alignment/1_suillus_aligned_fastq/${id}_aligned_R2.fastq.gz
   fi
done


gtf_suicot="$dir/reference/Suicot1_GeneCatalog_20171209.gtf"
gtf_suicot="$dir/reference/Suicot1_all_genes_20171209.gff"
featureCounts -t exon -F GFF -f -p -g name -M -T 24 -a $gtf_suicot -o $dir/3_suillus_alignment/3_suillus_count_file/suillus_all_gene_id_count.txt *.bam  1>$dir/3_suillus_alignment/3_suillus_count_file/counts.all_gene_id.log 2>&1
featureCounts -t exon -F GFF -f -p -O -g name -M -T 24 -a $gtf_suicot -o $dir/3_suillus_alignment/3_suillus_count_file/suillus_all_gene_id_count_overlap.txt *.bam  1>$dir/3_suillus_alignment/3_suillus_count_file/counts.all_gene_id_overlap.log 2>&1
featureCounts -t exon -F GTF -g gene_id -p -M -T 24 -a $gtf_suicot -o $dir/3_suillus_alignment/3_suillus_count_file/suillus_catalog_gene_id_count.txt *.bam  1>$dir/3_suillus_alignment/3_suillus_count_file/counts.catalog_gene_id.log 2>&1
featureCounts -t exon -F GTF -g gene_id -p -O -M -T 24 -a $gtf_suicot -o $dir/3_suillus_alignment/3_suillus_count_file/suillus_catalog_gene_id_count_overlap.txt *.bam  1>$dir/3_suillus_alignment/3_suillus_count_file/counts.catalog_gene_id_overlap.log 2>&1


#extract mapping rate
echo -e "fileName\trate(%)" > Suillus_mappingRate.txt
for file in $(ls *.log); do
  i=$(basename $file _bowtie2.log)
  rate=$(tac $file |head -n 1 |cut -d"%" -f 1)
  echo -e ${i}"\t"${rate} >> Suillus_mappingRate.txt
done
  
  
  

done




#  section 3: mapping rhizopogon
####################################################
bowtie2-build --threads 4 Rhisa1_AssemblyScaffolds_Repeatmasked.fasta Rhisa1_genome
bowtie2-build --threads 4 Rhitru1_AssemblyScaffolds_Repeatmasked.fasta Rhitru1_genome
bowtie2-build --threads 4 Rhives1_AssemblyScaffolds_Repeatmasked.fasta Rhives1_genome
bowtie2-build --threads 4 Rhivi1_AssemblyScaffolds_Repeatmasked.fasta Rhivi1_genome
bowtie2-build --threads 4 Rhivul1_AssemblyScaffolds_Repeatmasked.fasta Rhivul1_genome


Rhisa1_genome="/Volumes/T7/plant_emf_sap_interaction/rhizopogon_ref/Rhisa1_genome"
Rhitru2_genome="/Volumes/T7/plant_emf_sap_interaction/rhizopogon_ref/Rhitru1_genome"
Rhives3_genome="/Volumes/T7/plant_emf_sap_interaction/rhizopogon_ref/Rhives1_genome"
Rhivi4_genome="/Volumes/T7/plant_emf_sap_interaction/rhizopogon_ref/Rhivi1_genome"
Rhivul5_genome="/Volumes/T7/plant_emf_sap_interaction/rhizopogon_ref/Rhivul1_genome"


ls *.gz|cut -d"_" -f 1,2,3 |sort -u  |while read id;do
bowtie2 -p 4 -x $Rhisa1_genome \
   -1 ${id}_paired_R1.fq.gz \
   -2 ${id}_paired_R2.fq.gz \
   -S /Volumes/T7/plant_emf_sap_interaction/5_rhizoppgon_mapping/${id}_Rhisa1.sam \
   2>/Volumes/T7/plant_emf_sap_interaction/5_rhizoppgon_mapping/${id}_bowtie2_Rhisa1.log
rm /Volumes/T7/plant_emf_sap_interaction/5_rhizoppgon_mapping/${id}_Rhisa1.sam
done

ls *.gz|cut -d"_" -f 1,2,3 |sort -u  |while read id;do
bowtie2 -p 4 -x $Rhitru2_genome \
   -1 ${id}_paired_R1.fq.gz \
   -2 ${id}_paired_R2.fq.gz \
   -S /Volumes/T7/plant_emf_sap_interaction/5_rhizoppgon_mapping/${id}_Rhitru2.sam \
   2>/Volumes/T7/plant_emf_sap_interaction/5_rhizoppgon_mapping/${id}_bowtie2_Rhitru2.log
rm /Volumes/T7/plant_emf_sap_interaction/5_rhizoppgon_mapping/${id}_Rhitru2.sam
done


ls *.gz|cut -d"_" -f 1,2,3 |sort -u  |while read id;do
bowtie2 -p 4 -x $Rhives3_genome \
   -1 ${id}_paired_R1.fq.gz \
   -2 ${id}_paired_R2.fq.gz \
   -S /Volumes/T7/plant_emf_sap_interaction/5_rhizoppgon_mapping/${id}_Rhives3.sam \
   2>/Volumes/T7/plant_emf_sap_interaction/5_rhizoppgon_mapping/${id}_bowtie2_Rhives3.log
rm /Volumes/T7/plant_emf_sap_interaction/5_rhizoppgon_mapping/${id}_Rhives3.sam
done


ls *.gz|cut -d"_" -f 1,2,3 |sort -u  |while read id;do
bowtie2 -p 4 -x $Rhivi4_genome \
   -1 ${id}_paired_R1.fq.gz \
   -2 ${id}_paired_R2.fq.gz \
   -S /Volumes/T7/plant_emf_sap_interaction/5_rhizoppgon_mapping/${id}_Rhivi4.sam \
   2>/Volumes/T7/plant_emf_sap_interaction/5_rhizoppgon_mapping/${id}_bowtie2_Rhivi4.log
rm /Volumes/T7/plant_emf_sap_interaction/5_rhizoppgon_mapping/${id}_Rhivi4.sam
done



ls *.gz|cut -d"_" -f 1,2,3 |sort -u  |while read id;do
bowtie2 -p 4 -x $Rhivul5_genome \
   -1 ${id}_paired_R1.fq.gz \
   -2 ${id}_paired_R2.fq.gz \
   -S /Volumes/T7/plant_emf_sap_interaction/5_rhizoppgon_mapping/${id}_Rhivul5.sam \
   2>/Volumes/T7/plant_emf_sap_interaction/5_rhizoppgon_mapping/${id}_bowtie2_Rhivul5.log
rm /Volumes/T7/plant_emf_sap_interaction/5_rhizoppgon_mapping/${id}_Rhivul5.sam
done




#  section 3: mapping pinus #trinity
####################################################
#3.1 mapping for pinus teada

Suicot_genome="$dir/reference/Suicot_genome"
Pintaeda_genome="/home/microbiome/data_storage/SATA2/plant_genome/pita/index/pita"
gtf_suicot="$dir/reference/Suicot1_GeneCatalog_20171209.gtf"
gtf_pintaeda="/home/microbiome/data_storage/SATA2/plant_genome/pita/Pita.2_01.gtf"
#gtf_pintaeda="$dir/reference/TrPtA_269336_P_taeda_mRNAdatabase_328662.gtf"
Pintaeda_transcript="/home/microbiome/data_storage/SATA2/plant_genome/pita/Pita.2_01.fa"

/home/microbiome/miniconda3/envs/RNASeq/bin/align_and_estimate_abundance.pl \
  --transcripts $Pintaeda_transcript \
  --seqType fq \
  --samples_file file_list \
  --est_method RSEM \
  --aln_method bowtie2 \
  --trinity_mode \
  --prep_reference \
  --output_dir 4_pinus_rsem_estimate_outdir \
  --thread_count 24


/home/microbiome/miniconda3/envs/RNASeq/bin/abundance_estimates_to_matrix.pl \
  --est_method RSEM \
  --cross_sample_norm TMM \
  --quant_files file.listing_target_files.txt \
  --gene_trans_map /home/microbiome/data_storage/SATA2/RNA_data/Liu_RNA/Trinity_output/Trinity.fasta.gene_trans_map \
  --name_sample_by_basedir \
  --basedir_index basedir_index.txt \
  --out_prefix  Quantification_Result



#  section 3: mapping pinus #bowtie2
####################################################
#3.1 mapping for pinus teada
dir=
mkdir $dir/4_alignment_pinus
mkdir $dir/4_alignment_pinus/1_pinus_aligned_fastq
mkdir $dir/4_alignment_pinus/1_pinus_unaligned_fastq
mkdir $dir/4_alignment_pinus/1_bowtie2_met_file
mkdir $dir/4_alignment_pinus/1_bowtie2_log_file
mkdir $dir/4_alignment_pinus/2_bam_flagstat_file
mkdir $dir/4_alignment_pinus/2_bam_file
mkdir $dir/4_alignment_pinus/3_pinus_count_express
mkdir $dir/4_alignment_pinus/3_pinus_count_featurecount

bowtie2-build --threads 4 $dir/reference/TrPtA_269336_P_taeda_mRNAdatabase_328662.fasta $dir/reference/Pintaeda_genome

Suicot_genome="$dir/reference/Suicot_genome"
Pintaeda_genome="/home/microbiome/data_storage/SATA2/plant_genome/pita/index/pita"
gtf_suicot="$dir/reference/Suicot1_GeneCatalog_20171209.gtf"
gtf_pintaeda="/home/microbiome/data_storage/SATA2/plant_genome/pita/Pita.2_01.gtf"
#gtf_pintaeda="$dir/reference/TrPtA_269336_P_taeda_mRNAdatabase_328662.gtf"
Pintaeda_transcript="/home/microbiome/data_storage/SATA2/plant_genome/pita/Pita.2_01.fa"

dir=/Volumes/T7/plant_emf_sap_interaction
cd $dir/1_cleandata

ls *.gz|cut -d"_" -f 1,2,3 |sort -u  |while read id;do
  echo $id
  fwd="cleandata_rmdup/${id}_rmdup_val_1.fq.gz"
  rev="cleandata_rmdup/${id}_rmdup_val_2.fq.gz"
  if [ -f "$dir/4_alignment_pinus/1_pinus_unaligned_fastq/${id}_unaligned_R2.fastq.gz" ]; then
    echo "${id} has been analyzed"
  else
    mkdir $dir/4_alignment_pinus/${id}.temp
    time bowtie2 -p 24 -x $Pintaeda_genome \
      -1  $fwd \
      -2  $rev \
      -S $dir/4_alignment_pinus/${id}.temp/${id}_pinus.sam \
      --al-conc-gz $dir/4_alignment_pinus/1_pinus_aligned_fastq/${id}_aligned.fastq.gz \
      --un-conc-gz $dir/4_alignment_pinus/1_pinus_unaligned_fastq/${id}_unaligned.fastq.gz \
      --met-file $dir/4_alignment_pinus/1_bowtie2_met_file/${id}_met.txt \
      2>$dir/4_alignment_pinus/1_bowtie2_log_file/${id}_bowtie2.log

  samtools sort -o bam -@ 24 -o $dir/4_alignment_pinus/2_bam_file/${id}_pinus.bam $dir/4_alignment_pinus/${id}.temp/${id}_pinus.sam
  samtools flagstat -@ 24 $dir/4_alignment_pinus/${id}.temp/${id}_pinus.sam > $dir/4_alignment_pinus/2_bam_flagstat_file/${id}.flagstat

  featureCounts -t gene -F GTF -g gene_id -T 24 -a $gtf_pintaeda \
    -o $dir/4_alignment_pinus/3_pinus_count_featurecount/${id}_pinus_gene_id_count.txt \
    $dir/4_alignment_pinus/2_bam_file/${id}_pinus.bam  \
    1>$dir/4_alignment_pinus/3_pinus_count_featurecount/${id}_pinus.log 2>&1

# samtools sort -n -@ 24 -o $dir/4_alignment_pinus/${id}.temp/${id}_pinus.sort.bam $dir/4_alignment_pinus/2_bam_file/${id}_pinus.bam
# express -o $dir/4_alignment_pinus/3_pinus_count_express $Pintaeda_transcript $dir/4_alignment_pinus/${id}.temp/${id}_pinus.sort.bam
# mv $dir/4_alignment_pinus/3_pinus_count_express/results.xprs $dir/4_alignment_pinus/3_pinus_count_express/${id}.bam.tab
# mv $dir/4_alignment_pinus/3_pinus_count_express/params.xprs  $dir/4_alignment_pinus/3_pinus_count_express/${id}.bam.params.xprs

  rm -rf $dir/4_alignment_pinus/${id}.temp
  mv $dir/4_alignment_pinus/1_pinus_aligned_fastq/${id}_aligned.fastq.1.gz $dir/4_alignment_pinus/1_pinus_aligned_fastq/${id}_aligned_R1.fastq.gz
  mv $dir/4_alignment_pinus/1_pinus_aligned_fastq/${id}_aligned.fastq.2.gz $dir/4_alignment_pinus/1_pinus_aligned_fastq/${id}_aligned_R2.fastq.gz
  mv $dir/4_alignment_pinus/1_pinus_unaligned_fastq/${id}_unaligned.fastq.1.gz $dir/4_alignment_pinus/1_pinus_unaligned_fastq/${id}_unaligned_R1.fastq.gz
  mv $dir/4_alignment_pinus/1_pinus_unaligned_fastq/${id}_unaligned.fastq.2.gz $dir/4_alignment_pinus/1_pinus_unaligned_fastq/${id}_unaligned_R2.fastq.gz
fi
done


featureCounts -t gene -F GTF -g gene_id -T 24 -a $gtf_pintaeda \
   -o $dir/4_alignment_pinus/3_pinus_count_featurecount/pinus_gene_id_count.txt \
   $dir/4_alignment_pinus/2_bam_file/*_pinus.bam  \
   1>$dir/4_alignment_pinus/3_pinus_count_featurecount/pinus_gene_id_count.log 2>&1


featureCounts -t gene -F GTF -g gene_id -p -M -T 24 -a $gtf_pintaeda \
    -o $dir/4_alignment_pinus/3_pinus_count_featurecount/pinus_gene_id_count_overlap.txt \
   $dir/4_alignment_pinus/2_bam_file/*_pinus.bam  \
   1>$dir/4_alignment_pinus/3_pinus_count_featurecount/pinus_gene_id_count_overlap.log 2>&1





















