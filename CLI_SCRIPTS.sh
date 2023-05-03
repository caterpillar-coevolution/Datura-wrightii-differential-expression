# Datura-wrightii-genome: UNIX scripts associated with assembly, annotation, and analyses
# Broken down into 8 steps: 1-assembly (hifiasm), 3-blobtools, 4-RepeatMasker, 5-structural annotation, 6-functional annotation, 7-BUSCO, 8-Gene counting (STAR), 9-Paleopolyploidy analyses (DupPipe and MCScanX)

#Step 1: Assembly
cat RAW_DATA_FILES > merged.fastq.gz
hifiasm -o dwrightii_hifi_1.asm -t 32 merged.fastq.gz
#Note: Assembly stats were obtained using the Bandage GUI and no script was used

#Step 3: blobtools (contamination check)

#blast against nt to get Taxa IDs
module load blast/2.13.0  
blastn \
-task megablast \
-query ./dwrightii_hifi_1.asm.bp.p_ctg.fa \
-db nt \
-outfmt '6 qseqid staxids bitscore std' \
-max_target_seqs 1 \
-max_hsps 1 \
-num_threads 32 \
-evalue 1e-25 \
-out Dwri.redo.vs.nt.redownloaded.mts1.hsp1.1e25.megablast.out

#align raw reads against finished assembly
minimap2 -ax map-hifi ./dwrightii_hifi_1.asm.bp.p_ctg.fa ./Raw_Dwri_HiFi_reads.fasta.gz > Dwri_coverage.sam
module load samtools
samtools view -bS Dwri_coverage.sam > Dwri_coverage.bam
samtools sort Dwri_coverage.bam -o Dwri_coverage_sorted.bam
samtools index Dwri_coverage_sorted.bam

#run the blobtools pipeline

blobtools create \
 -i ./dwrightii_hifi_1.asm.bp.p_ctg.fa \
 -b ./Dwri_coverage_sorted.bam \
 -t ./Dwri.redo.vs.nt.redownloaded.mts1.hsp1.1e25.megablast.out \
 -o ./Dwri_blobplot \
 --db nodesDB.txt \
 -x bestsumorder

/xdisk/judieb/jaykgold/blob_tools/blobtools/blobtools view \
 -i ./Dwri_blobplot.blobDB.json \
 -x bestsumorder \
 -r family \
 -o ./
 
 grep '^##' example/Dwri_blobplot.blobDB.table.txt ; \
 grep -v '^##' example/Dwri_blobplot.blobDB.table.txt | \
 column -t -s $'\t'
 
/xdisk/judieb/jaykgold/blob_tools/blobtools/blobtools plot \
 -i ./Dwri_blobplot.blobDB.json \
 -x bestsumorder \
 -r family \
 -o ./

#Step 4: Repeat Masking
awk '/^S/{print ">"$2"\n"$3}' dwrightii_hifi_1.asm.bp.hap1.p_ctg.gfa | fold > dwrightii_hifi_1.asm.bp.hap1.p_ctg.fa

singularity exec -B ~/trf409.linux64:/opt/trf:ro tetools_1_1.sif BuildDatabase -name Dwri_genomic.DB -engine rmblast dwrightii_hifi_1.asm.bp.p_ctg.fa

singularity exec -B ~/trf409.linux64:/opt/trf:ro tetools_1_1.sif RepeatModeler -pa 32 -database Dwri_genomic.DB -LTRStruct

singularity exec -B ~/trf409.linux64:/opt/trf:ro tetools_1_1.sif RepeatMasker -lib Dwri_genomic.DB-families.fa -xsmall -pa 32 -gff -e ncbi dwrightii_hifi_1.asm.bp.p_ctg.fa

#Step 5: Structural annotation

singularity exec /xdisk/judieb/jaykgold/Funannotate/funannotate.sif funannotate train \
    --cpus 32 \
    -i dwrightii_hifi_1.asm.bp.p_ctg.fa.masked \
    -o Dwri_corrected_annotation \
    -l Dwri_root_1.fq.gz Dwri_IB_1.fq.gz \
    -r Dwri_root_2.fq.gz Dwri_IB_2.fq.gz \
    -s flowcell152_lane8_ATCACGLateField.fastq.gz \
        flowcell152_lane8_TTAGGCEarlyField.fastq.gz \
        flowcell147_lane2_ACTTGA2CLC3.fastq.gz \
        flowcell147_lane2_ATCACG2CLA1.fastq.gz \
        flowcell147_lane2_GATCAG2LA3.fastq.gz \
        flowcell147_lane2_GGCTAC2LC2.fastq.gz \
        flowcell147_lane2_TAGCTT2LB3.fastq.gz \
        flowcell147_lane2_TTAGGC2CLB2.fastq.gz \
        flowcell203_lane7_ACAGTG_3D0.fastq.gz \
        flowcell203_lane7_ATCACG_1C0.fastq.gz \
        flowcell203_lane7_CGATGT_3C0.fastq.gz \
        flowcell203_lane7_GCCAAT4D0.fastq.gz \
        flowcell203_lane7_TGACCA1D0.fastq.gz \
        flowcell203_lane7_TTAGGC4C0.fastq.gz \
        flowcell206_lane5_ACTTGA3c12.fastq.gz \
        flowcell206_lane5_CAGATC1c12.fastq.gz \
        flowcell206_lane5_CTTGTA4d12.fastq.gz \
        flowcell206_lane5_GATCAG4c12.fastq.gz \
        flowcell206_lane5_GGCTAC3d12.fastq.gz \
        flowcell206_lane5_TAGCTT1d12.fastq.gz \
        flowcell210_lane6_ACTTGA3c24.fastq.gz \
        flowcell210_lane6_GATCAG4c24.fastq.gz \
        flowcell210_lane6_GGCTAC3d24.fastq.gz \
        flowcell220_lane8_ATCACG1c24c.fastq.gz \
        flowcell215_lane2_ACAGTG3D.fastq.gz \
        flowcell215_lane2_ATCACG1C.fastq.gz \
        flowcell215_lane2_CGATGT3C.fastq.gz \
        flowcell215_lane2_GCCAAT4D.fastq.gz \
        flowcell215_lane2_TGACCA1D.fastq.gz \
        flowcell215_lane2_TTAGGC4C.fastq.gz \
        flowcell215_lane3_ACTTGA3C.fastq.gz \
        flowcell215_lane3_GATCAG4C.fastq.gz \
        flowcell215_lane3_GGCTAC3D.fastq.gz \
        flowcell215_lane3_TAGCTT1D.fastq.gz \
        flowcell223_lane1_CAGATC1C.fastq.gz \
        flowcell223_lane1_CTTGTA4D.fastq.gz \

singularity exec /xdisk/judieb/jaykgold/Funannotate/funannotate.sif funannotate predict \
    -i dwrightii_hifi_1.asm.bp.p_ctg.fa.masked \
    -o Dwri_corrected_annotation  \
    -s "Dwri_corrected_annotation" \
    --cpus 32 \
    --max_intronlen 10000 \
    --organism other \
    --busco_db embryophyta \
    -d /xdisk/judieb/jaykgold/Funannotate/funannotate_databases \
    --repeats2evm \
    --protein_evidence dwri_peptides.fa

#Step 6: Functional annotation (InterProScan and funannotate)
#Note: InterProScan was run on a private server due to version conflicts with key dependencies
nohup my_interproscan/interproscan-5.57-90.0/interproscan.sh -i Dwri_corrected_annotation.proteins.fa &
singularity exec funannotate.sif funannotate annotate \
    -i ./Dwri_corrected_annotation/ \
    --cpus 32 \
    --iprscan Dwri_corrected_annotation.proteins.fa.xml

#Step 7: BUSCO analysis and figure making
busco_5.1.3.sif busco -i dwrightii_hifi_1.asm.bp.p_ctg.fa -l solanales_odb10 -o Dwri_assembly_BUSCO -m genome -f
busco_5.1.3.sif busco -i ./Dwri_corrected_annotation.proteins.fa -l solanales_odb10 -o Dwri_transcriptome_busco_analysis -m prot -f
module load python/3.9/3.9.10
singularity exec busco_5.1.3.sif python3 generate_plot.py --working_directory /path/to/WD

#Step 8: generate gene counts for differential expression study
module load star
module load cufflinks
module load bcftools
module load samtools
module load perl
module load sratoolkit

gffread Dwri_corrected_annotation.gff3 -T -o Dwri_corrected_annotation.gtf

STAR --runThreadN 15 \
--runMode genomeGenerate \
--genomeDir Dwri/ \
--genomeFastaFiles dwrightii_hifi_1.asm.bp.p_ctg.fa.masked  \
--sjdbGTFfile Dwri_corrected_annotation.gtf \
--genomeSAindexNbases 13 \
--sjdbOverhang 74 

#Run trimmomatic
java -jar /xdisk/judieb/jaykgold/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 \
fasta/replace.fastq.gz \
fasta/replace_trimmed.fa \
ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

STAR --genomeDir Dwri/ \
--runThreadN 15 \
--readFilesIn fasta/replace_trimmed.fa \
--outFileNamePrefix results/replace \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
--limitSjdbInsertNsj=3000000 \
--limitOutSJcollapsed=3000000


#index the bam file
samtools index results/replaceAligned.sortedByCoord.out.bam

#Step 9: paleopolyploidy analyses (DupPipe and MCScanX)
#Note: these analyses were performed by my coauthor, Michael McKibben so exact scripts are not available to me at the time of this writing.
