####################################################
#This is for metatranscriptomic mapping analysis of plant_emf_sap_interaction project
####################################################
#!/bin/bash
conda activate RNASeq
dir=/Volumes/T7/plant_emf_sap_interaction


#  section 1: reference preparation
####################################################
cd $dir/reference

gunzip *.gz
gffread Suicot1_GeneCatalog_20171209.gff3 -T -o Suicot1_GeneCatalog_20171209.gtf

bowtie2-build --threads 4  Suicot1_AssemblyScaffolds_Repeatmasked.fasta Suicot_genome
bowtie2-build --threads 4 TrPtA_269336_P_taeda_mRNAdatabase_328662.fasta Pintaeda_genome


#3.mapping
##################################################################
Suicot_genome="$dir/reference/Suicot_genome"
Pintaeda_genome="/home/microbiome/data_storage/SATA2/plant_genome/pita/index/pita"
gtf_suicot="$dir/reference/Suicot1_GeneCatalog_20171209.gtf"
gtf_pintaeda="/home/microbiome/data_storage/SATA2/plant_genome/pita/Pita.2_01.gtf"



#  section2: mapping to suillus
################################################################
dir=/Volumes/T7/plant_emf_sap_interaction
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


dir=/Volumes/T7/plant_emf_sap_interaction
cd $dir/cleandata
ls *.gz|cut -d"_" -f 1,2,3 |sort -u |while read id;do
if [ -f "$dir/3_suillus_alignment/3_suillus_count_file/${id}_suillus_gene_id_count.txt" ]; then
    echo "${id} has been analyzed"
    else
mkdir $dir/3_suillus_alignment/${id}.temp
time bowtie2 -p 4 -x $Suicot_genome \
-1 ${id}_paired_R1.fq.gz \
-2 ${id}_paired_R2.fq.gz \
-S $dir/3_suillus_alignment/${id}.temp/${id}_suillus.sam \
--al-conc-gz $dir/3_suillus_alignment/1_suillus_aligned_fastq/${id}_aligned.fastq.gz \
--un-conc-gz $dir/3_suillus_alignment/1_suillus_unaligned_fastq/${id}_unaligned.fastq.gz \
--met-file $dir/3_suillus_alignment/1_bowtie2_met_file/${id}_met.txt \
2>$dir/3_suillus_alignment/1_bowtie2_log_file/${id}_bowtie2.log

samtools sort -o bam -@ 3 -o $dir/3_suillus_alignment/2_suillus_bam_file/${id}_suillus.bam $dir/3_suillus_alignment/${id}.temp/${id}_suillus.sam
samtools flagstat -@ 3 $dir/3_suillus_alignment/${id}.temp/${id}_suillus.sam > $dir/3_suillus_alignment/2_bam_flagstat_file/${id}.flagstat
featureCounts -t exon -F GTF -g gene_id -T 4 -a $gtf_suicot -o $dir/3_suillus_alignment/3_suillus_count_file/${id}_suillus_gene_id_count.txt $dir/3_suillus_alignment/2_suillus_bam_file/${id}_suillus.bam  1>$dir/suillus_alignment/3_suillus_count_file/${id}_suillus.log 2>&1
rm -rf $dir/3_suillus_alignment/${id}.temp
mv $dir/3_suillus_alignment/1_suillus_aligned_fastq/${id}_aligned.fastq.1.gz $dir/3_suillus_alignment/1_suillus_aligned_fastq/${id}_aligned_R1.fastq.gz
mv $dir/3_suillus_alignment/1_suillus_aligned_fastq/${id}_aligned.fastq.2.gz $dir/3_suillus_alignment/1_suillus_aligned_fastq/${id}_aligned_R2.fastq.gz
fi
done


gtf_suicot="$dir/reference/Suicot1_GeneCatalog_20171209.gtf"
gtf_suicot="$dir/reference/Suicot1_all_genes_20171209.gff"
featureCounts -t exon -F GFF -f -p -g name -M -T 4 -a $gtf_suicot -o $dir/3_suillus_alignment/3_suillus_count_file/suillus_all_gene_id_count.txt *.bam  1>$dir/3_suillus_alignment/3_suillus_count_file/counts.all_gene_id.log 2>&1
featureCounts -t exon -F GFF -f -p -O -g name -M -T 4 -a $gtf_suicot -o $dir/3_suillus_alignment/3_suillus_count_file/suillus_all_gene_id_count_overlap.txt *.bam  1>$dir/3_suillus_alignment/3_suillus_count_file/counts.all_gene_id_overlap.log 2>&1
featureCounts -t exon -F GTF -g gene_id -p -M -T 4 -a $gtf_suicot -o $dir/3_suillus_alignment/3_suillus_count_file/suillus_catalog_gene_id_count.txt *.bam  1>$dir/3_suillus_alignment/3_suillus_count_file/counts.catalog_gene_id.log 2>&1
featureCounts -t exon -F GTF -g gene_id -p -O -M -T 4 -a $gtf_suicot -o $dir/3_suillus_alignment/3_suillus_count_file/suillus_catalog_gene_id_count_overlap.txt *.bam  1>$dir/3_suillus_alignment/3_suillus_count_file/counts.catalog_gene_id_overlap.log 2>&1





time bowtie2 -p 5 -x $Suicot_genome \
-1 10424.5.160042.GGCTAC.anqrpht_R1.fastq.gz \
-2 10424.5.160042.GGCTAC.anqrpht_R2.fastq.gz \
-S suillus_conthurnatus_culture_RNA.sam \
--al-conc-gz suillus_conthurnatus_culture_RNA_aligned.fastq.gz \
--un-conc-gz suillus_conthurnatus_culture_RNA_unaligned.fastq.gz \
--met-file suillus_conthurnatus_culture_RNA_met.txt \
2>suillus_conthurnatus_culture_RNA_bowtie2.log

samtools sort -o bam -@ 3 -o suillus_conthurnatus_culture_RNA.bam suillus_conthurnatus_culture_RNA.sam
samtools flagstat -@ 3 suillus_conthurnatus_culture_RNA.sam > suillus_conthurnatus_culture_RNA.flagstat

featureCounts -t exon -F GTF -g gene_id -p -O -M -T 4 -a $gtf_suicot -o suillus_catalog_gene_id_count_overlap.txt *.bam  1>counts.catalog_gene_id_overlap.log 2>&1



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







#  section 4: mapping rhizopogon suillus 
####################################################
dir=/Volumes/T7/plant_emf_sap_interaction


ls *.gz |while read id;do
trim_galore -q 25 --phred33 --stringency 3 --length 25 \
${id} \
--gzip \
--cores 3 \
-o $dir/1_suillus_culture_RNA
done 


bowtie2-build --threads 4 Rhivul1_AssemblyScaffolds_Repeatmasked.fasta Rhivul_genome
bowtie2-build --threads 4 Suilu4_AssemblyScaffolds_Repeatmasked.fasta Suilut_genome

Suicot_genome="$dir/reference/Suicot_genome"
Suilut_genome="/Volumes/T7/plant_emf_sap_interaction/1_suillus_culture_RNA/suillus_luteus/Suilut_genome"
Rhivul_genome="/Volumes/T7/plant_emf_sap_interaction/1_suillus_culture_RNA/Rhizopogon_vulgaris/Rhivul_genome"


ls *_trimmed.fq.gz |cut -d"_" -f 1,2,3 |while read id;do

time bowtie2 -p 4 -x $Suilut_genome \
   -U ${id}_trimmed.fq.gz \
   -S ${id}_Suilut.sam \
   2>${id}_Suilut_bowtie2.log
rm ${id}_Suilut.sam


time bowtie2 -p 4 -x $Rhivul_genome \
   -U ${id}_trimmed.fq.gz \
   -S ${id}_Rhivul.sam \
   2>${id}_Rhivul_bowtie2.log
rm ${id}_Rhivul.sam
done




bowtie2 -p 4 -x $Suicot_genome \
   -U ${id} \
   -S ${id}_Suicot.sam \
   2>${id}_Suicot_bowtie2.log
rm ${id}_Suicot.sam








#  section 3: mapping pinus
####################################################
#3.1 mapping for pinus teada
dir=/Volumes/T7/plant_emf_sap_interaction
mkdir $dir/4_alignment_pinus
mkdir $dir/4_alignment_pinus/1_pinus_aligned_fastq
mkdir $dir/4_alignment_pinus/1_bowtie2_met_file
mkdir $dir/4_alignment_pinus/1_bowtie2_log_file
mkdir $dir/4_alignment_pinus/2_bam_flagstat_file
mkdir $dir/4_alignment_pinus/2_bam_file
mkdir $dir/4_alignment_pinus/3_pinus_count_express
mkdir $dir/4_alignment_pinus/3_pinus_count_featurecount

bowtie2-build --threads 4 $dir/reference/TrPtA_269336_P_taeda_mRNAdatabase_328662.fasta $dir/reference/Pintaeda_genome

Pintaeda_genome="$dir/reference/Pintaeda_genome"
gtf_pintaeda="$dir/reference/TrPtA_269336_P_taeda_mRNAdatabase_328662.gtf"
Pintaeda_transcript="$dir/reference/TrPtA_269336_P_taeda_mRNAdatabase_328662.fasta"

dir=/Volumes/T7/plant_emf_sap_interaction
cd $dir/1_cleandata
ls *.gz|cut -d"_" -f 1,2,3 |sort -u  |while read id;do
if [ -f "$dir/4_alignment_pinus/3_pinus_count_express/${id}.bam.tab" ]; then
    echo "${id} has been analyzed"
    else
    
mkdir $dir/4_alignment_pinus/${id}.temp
time bowtie2 -p 4 -x $Pintaeda_genome \
   -1 ${id}_paired_R1.fq.gz \
   -2 ${id}_paired_R2.fq.gz \
   -S $dir/4_alignment_pinus/${id}.temp/${id}_pinus.sam \
   --al-conc-gz $dir/4_alignment_pinus/1_pinus_aligned_fastq/${id}_aligned.fastq.gz \
   --un-conc-gz $dir/4_alignment_pinus/${id}.temp/${id}_unaligned.fastq.gz \
   --met-file $dir/4_alignment_pinus/1_bowtie2_met_file/${id}_met.txt \
   2>$dir/4_alignment_pinus/1_bowtie2_log_file/${id}_bowtie2.log

samtools sort -o bam -@ 24 -o $dir/4_alignment_pinus/2_bam_file/${id}_pinus.bam $dir/4_alignment_pinus/${id}.temp/${id}_pinus.sam
samtools flagstat -@ 24 $dir/4_alignment_pinus/${id}.temp/${id}_pinus.sam > $dir/4_alignment_pinus/2_bam_flagstat_file/${id}.flagstat

featureCounts -t exon -F GTF -g gene_id -T 24 -a $gtf_pintaeda \
   -o $dir/4_alignment_pinus/3_pinus_count_featurecount/${id}_pinus_gene_id_count.txt \
   $dir/4_alignment_pinus/2_bam_file/${id}_pinus.bam  \
   1>$dir/4_alignment_pinus/3_pinus_count_featurecount/${id}_pinus.log 2>&1

samtools sort -n -@ 4 -o $dir/4_alignment_pinus/${id}.temp/${id}_pinus.sort.bam $dir/4_alignment_pinus/2_bam_file/${id}_pinus.bam
express -o $dir/4_alignment_pinus/3_pinus_count_express $Pintaeda_transcript $dir/4_alignment_pinus/${id}.temp/${id}_pinus.sort.bam
mv $dir/4_alignment_pinus/3_pinus_count_express/results.xprs $dir/4_alignment_pinus/3_pinus_count_express/${id}.bam.tab
mv $dir/4_alignment_pinus/3_pinus_count_express/params.xprs  $dir/4_alignment_pinus/3_pinus_count_express/${id}.bam.params.xprs

rm -rf $dir/4_alignment_pinus/${id}.temp
mv $dir/4_alignment_pinus/1_pinus_aligned_fastq/${id}_aligned.fastq.1.gz $dir/4_alignment_pinus/1_pinus_aligned_fastq/${id}_aligned_R1.fastq.gz
mv $dir/4_alignment_pinus/1_pinus_aligned_fastq/${id}_aligned.fastq.2.gz $dir/4_alignment_pinus/1_pinus_aligned_fastq/${id}_aligned_R2.fastq.gz
fi
done



























