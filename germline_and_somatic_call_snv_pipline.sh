#!/bin/sh

set -eu

# Update with the fullpath location of your sample fastq
TUMOR_SM=$1 #sample name
TUMOR_RGID="rg_$TUMOR_SM" #read group ID
# PL="DNB" #or other sequencing platform
PL="illumina" #or other sequencing platform
# FASTQ_FOLDER="/home/pipeline/samples"
TUMOR_FASTQ_1=$2
TUMOR_FASTQ_2=$3

GENOME_VERSION="hg19"

SAMPLE_INFO="/home/myhone/sample_infomation_table.txt"

# Update with the location of the reference data files
FASTA_DIR="/home/myhone/Reference/GRCh37_Reference/Sequence/BWAIndex/bwa/"
FASTA="/home/myhone/Reference/GRCh37_Reference/Sequence/BWAIndex/bwa/Homo_sapiens_assembly19.fasta"
# MOUSE_FASTA="/share/Workspace/4.Database/Gencode/GRCm39/GRCm39.primary_assembly.genome.fa"
INTERVAL_FILE="/home/myhone/Reference/GRCh37_Reference/Annotation/intervals/NEXome_Core_Panel.target.hg19_HumanG1Kv37.nochr.noPAR.interval_list"
PROBE_FILE="/home/myhone/Reference/GRCh37_Reference/Annotation/intervals/NEXome_Core_Panel.probe_merged.hg19_HumanG1Kv37.nochr.noPAR.interval_list"
TARGET_BED="/home/myhone/Reference/GRCh37_Reference/Annotation/intervals/NEXome_Core_Panel.target.hg19_HumanG1Kv37.nochr.noPAR.bed"
KNOWN_DBSNP="/home/myhone/Reference/GRCh37_Reference/Annotation/GATKBundle/Homo_sapiens_assembly19.dbsnp138.vcf.gz"
KNOWN_INDELS="/home/myhone/Reference/GRCh37_Reference/Annotation/GATKBundle/1000G_phase1.indels.b37.vcf.gz"
KNOWN_MILLS="/home/myhone/Reference/GRCh37_Reference/Annotation/GATKBundle/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
# Update with the location of the panel of normal and CosmicDB vcf files
# We recommend that you create the panel of normal file with the corresponding algorithm that you plan to use for the somatic mutation calling. 
# PANEL_OF_NORMAL="/share/Workspace/4.Database/hg19/Annotation/GATKBundle/Mutect2-exome-panel.vcf"
GERMLINE_VCF="/home/myhone/Reference/GRCh37_Reference/Annotation/GermlineResource/gnomAD.r2.1.1.GRCh37.PASS.AC.AF.only.vcf.gz"

# Other settings
NT=$(nproc); if ((NT > 12)); then NT=12; fi; echo "使用核心数："$NT #number of threads to use in computation, set to number of cores in the server
START_DIR=$4 #Determine where the output files will be stored

QC_script="/share/result/WES_DX/20230808_NFYKDX_PXH_Human_Liver/scripts/Statistical_QC.py"
GATK=/home/myhone/Softwares/GATK4/gatk-4.2.6.1/gatk
qualimap=/home/myhone/Softwares/QualiMap/qualimap_v2.3/qualimap
vcf2maf=/home/myhone/Softwares/vcf2maf/vcf2maf.pl
vep_path=/home/myhone/miniconda3/envs/wes_tumor/bin/
vep_data=/home/myhone/Reference/VEP/VEP_106_GRCh37
annovar_path=/home/myhone/Softwares/ANNOVAR/annovar
annovar_data=/home/myhone/Reference/ANNOVAR/hg19


# You do not need to modify any of the lines below unless you want to tweak the pipeline

# ************************************************************************************************************************************************************************

# ******************************************
# 0. Setup
# ******************************************
WORKDIR="$START_DIR"
mkdir -p $WORKDIR
cd $WORKDIR
mkdir -p 1.raw_fq 2.clean_fq 3.qc/ 4.align/{flagstat,stats,mosdepth,picard} 5.gatk/{flagstat,stats,mosdepth,qualimap,picard} 6.gvcf 6.mutect  7.annotation/{vep,annovar} 8.cnv 9.downstream tmp log tmp_QC
LOGFILE=$WORKDIR/log/${TUMOR_SM}_run.log
exec >$LOGFILE 2>&1

# ******************************************
# 1a. Fastp Quanlity Control
# ******************************************
dir=2.clean_fq
biosoft=fastp
if [  ! -f  ./${dir}/ok.${TUMOR_SM}_${biosoft}.status ]; 
then

    fastp \
    --in1 ${TUMOR_FASTQ_1} \
    --in2 ${TUMOR_FASTQ_2} \
    --out1 ./${dir}/${TUMOR_SM}_R1.clean.fastq.gz \
    --out2 ./${dir}/${TUMOR_SM}_R2.clean.fastq.gz \
    --json ./${dir}/${TUMOR_SM}.fastp.json \
    --html ./${dir}/${TUMOR_SM}.fastp.html \
    --thread $NT \
    --detect_adapter_for_pe \
    -q 20 -u 50 -n 10 -l 40 \
    1> ./${dir}/${TUMOR_SM}.fastp.log 2>&1

    if [ $? -eq 0 ]
    then
    touch ./${dir}/ok.${TUMOR_SM}_${biosoft}.status 
    else
    touch ./${dir}/fail.${TUMOR_SM}_${biosoft}.status
    fi
    # sleep 60s
fi

pre_dir=2.clean_fq
dir=3.qc
biosoft=fastp_plot
if [ ! -f ./${dir}/ok.${TUMOR_SM}_${biosoft}.status ];
then
    conda \
    run \
    -p /home/myhone/miniconda3/envs/mustseq_pipe \
    Rscript \
    /home/myhone/Documents/RNAseq_DEG_pipeline/fastp_qc.R \
    -i ./${pre_dir}/${TUMOR_SM}.fastp.json \
    -o ./${dir}/${TUMOR_SM} \
    -s ${TUMOR_SM}

    if [ $? -eq 0 ]
    then
    touch ./${dir}/ok.${TUMOR_SM}_${biosoft}.status
    else
    touch ./${dir}/fail.${TUMOR_SM}_${biosoft}.status
    fi
fi

# ******************************************
# 1b. Mapping reads with BWA-MEM, sorting for tumor sample
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000 
dir=4.align
biosoft=bwa_${GENOME_VERSION}
if [  ! -f  ./${dir}/ok.${TUMOR_SM}_${biosoft}.status ]; 
then



    (bwa mem -R "@RG\tID:$TUMOR_RGID\tSM:$TUMOR_SM\tPL:$PL" \
    -t $NT -K 10000000 \
    $FASTA \
    ./2.clean_fq/${TUMOR_SM}_R1.clean.fastq.gz \
    ./2.clean_fq/${TUMOR_SM}_R2.clean.fastq.gz || \
    { echo -n 'BWA error'; touch ./${dir}/fail.${TUMOR_SM}_${biosoft}.status; exit 1;} ) | \
    samtools sort -@ $NT -o ./${dir}/${TUMOR_SM}.sort.bam - || \
    { echo "Alignment failed"; touch ./${dir}/fail.${TUMOR_SM}_${biosoft}.status; exit 1; }
    samtools index -@ $NT ./${dir}/${TUMOR_SM}.sort.bam

    if [ $? -eq 0 ]
    then
    touch ./${dir}/ok.${TUMOR_SM}_${biosoft}.status 
    else
    touch ./${dir}/fail.${TUMOR_SM}_${biosoft}.status
    fi
    # sleep 60s
fi

# ******************************************
# 2a. Metrics for tumor sample
# ******************************************
pre_dir=4.align
dir=4.align/flagstat
biosoft=samtools_flagstat
if [  ! -f  ./${dir}/ok.${TUMOR_SM}_${biosoft}.status ]; 
then

    samtools flagstat -@ $NT                   \
    ./${pre_dir}/${TUMOR_SM}.sort.bam          \
    > ./${dir}/${TUMOR_SM}_sort_flagstat.txt

    if [ $? -eq 0 ]
    then
    touch ./${dir}/ok.${TUMOR_SM}_${biosoft}.status 
    else
    touch ./${dir}/fail.${TUMOR_SM}_${biosoft}.status
    fi
    # sleep 60s
fi

pre_dir=4.align
dir=4.align/stats
biosoft=samtools_stats
if [  ! -f  ./${dir}/ok.${TUMOR_SM}_${biosoft}.status ]; 
then

    samtools stats -@ $NT                      \
    ./${pre_dir}/${TUMOR_SM}.sort.bam          \
    > ./${dir}/${TUMOR_SM}_sort_stats.txt

    if [ $? -eq 0 ]
    then
    touch ./${dir}/ok.${TUMOR_SM}_${biosoft}.status 
    else
    touch ./${dir}/fail.${TUMOR_SM}_${biosoft}.status
    fi
    # sleep 60s
fi

pre_dir=4.align
dir=4.align/mosdepth
biosoft=mosdepth
if [  ! -f  ./${dir}/ok.${TUMOR_SM}_${biosoft}.status ]; 
then

    mosdepth                                         \
    --threads $NT                                    \
    --by ${TARGET_BED}                               \
    --no-per-base                                    \
    --fast-mode                                      \
    --thresholds 0,1,2,3,4,5,10,20,30,50,100,200,300 \
    ${dir}/${TUMOR_SM}_dup                                  \
    ./${pre_dir}/${TUMOR_SM}.sort.bam


    if [ $? -eq 0 ]
    then
    touch ./${dir}/ok.${TUMOR_SM}_${biosoft}.status 
    else
    touch ./${dir}/fail.${TUMOR_SM}_${biosoft}.status
    fi
    # sleep 60s
fi

pre_dir=4.align
dir=4.align/picard
biosoft=picard
if [  ! -f  ./${dir}/ok.${TUMOR_SM}_${biosoft}.status ]; 
then

    ${GATK} --java-options "-Xmx30G"  CollectHsMetrics \
        -R $FASTA                                      \
        -I ./${pre_dir}/${TUMOR_SM}.sort.bam           \
        -O ./${dir}/${TUMOR_SM}.hs_metrics.txt         \
        -BI ${PROBE_FILE}                              \
        -TI ${INTERVAL_FILE}                           \
        --CLIP_OVERLAPPING_READS false                 \
        --TMP_DIR ./tmp                                \
        1>./${dir}/${TUMOR_SM}_CollectHsMetrics.log 2>&1

    ${GATK} --java-options "-Xmx30G"  CollectMultipleMetrics \
        -R $FASTA                                            \
        -I ./${pre_dir}/${TUMOR_SM}.sort.bam                 \
        -O ./${dir}/${TUMOR_SM}.multiple_metrics             \
        --TMP_DIR ./tmp                                      \
        1>./${dir}/${TUMOR_SM}_CollectMultipleMetrics.log 2>&1 

    if [ $? -eq 0 ]
    then
    touch ./${dir}/ok.${TUMOR_SM}_${biosoft}.status 
    else
    touch ./${dir}/fail.${TUMOR_SM}_${biosoft}.status
    fi
    # sleep 60s
fi

# ******************************************
# 3a. Remove Duplicate Reads for tumor
# sample. It is possible
# to mark instead of remove duplicates
# by ommiting the -REMOVE_DUPLICATES option in Dedup
# ******************************************
dir=5.gatk
biosoft=MarkDuplicates
if [  ! -f  ./${dir}/ok.${TUMOR_SM}_${biosoft}.status ]; 
then

    ${GATK} --java-options "-Xmx30G" MarkDuplicates \
        --REFERENCE_SEQUENCE $FASTA \
        --INPUT ./4.align/${TUMOR_SM}.sort.bam \
        --OUTPUT ./${dir}/${TUMOR_SM}.marked.bam \
        --METRICS_FILE ./${dir}/${TUMOR_SM}.metircs \
        --TMP_DIR ./tmp \
        --REMOVE_DUPLICATES false --VALIDATION_STRINGENCY LENIENT --CREATE_INDEX true \
        1>./${dir}/${TUMOR_SM}_mark.log 2>&1

    if [ $? -eq 0 ]
    then
    touch ./${dir}/ok.${TUMOR_SM}_${biosoft}.status 
    else
    touch ./${dir}/fail.${TUMOR_SM}_${biosoft}.status
    fi
    # sleep 60s
fi

# ******************************************
# 4a. Base recalibration for tumor sample
# ******************************************
dir=5.gatk
biosoft=BQSR
if [  ! -f  ./${dir}/ok.${TUMOR_SM}_${biosoft}.status ]; 
then

    ${GATK} --java-options "-Xmx30G"  BaseRecalibrator \
        -R $FASTA \
        -I ./${dir}/${TUMOR_SM}.marked.bam \
        -O ./${dir}/${TUMOR_SM}_recal.table \
        --known-sites ${KNOWN_DBSNP}        \
        --known-sites ${KNOWN_INDELS}        \
        --known-sites ${KNOWN_MILLS}        \
        --tmp-dir ./tmp \
        1>./${dir}/${TUMOR_SM}_recal.log 2>&1

    ${GATK} --java-options "-Xmx30G"  ApplyBQSR \
        -R $FASTA \
        -I ./${dir}/${TUMOR_SM}.marked.bam     \
        -bqsr ./${dir}/${TUMOR_SM}_recal.table \
        -O ./${dir}/${TUMOR_SM}.BQSR.bam       \
        --tmp-dir ./tmp \
        1>./${dir}/${TUMOR_SM}_bqsr.log 2>&1 

    if [ $? -eq 0 ]
    then
    touch ./${dir}/ok.${TUMOR_SM}_${biosoft}.status 
    else
    touch ./${dir}/fail.${TUMOR_SM}_${biosoft}.status
    fi
    # sleep 60s
fi

# ******************************************
# 4b. Metrics for tumor sample
# ******************************************
pre_dir=5.gatk
dir=5.gatk/flagstat
biosoft=samtools_flagstat
if [  ! -f  ./${dir}/ok.${TUMOR_SM}_${biosoft}.status ]; 
then

    samtools flagstat -@ $NT                   \
    ./${pre_dir}/${TUMOR_SM}.BQSR.bam          \
    > ./${dir}/${TUMOR_SM}_BQSR_flagstat.txt

    if [ $? -eq 0 ]
    then
    touch ./${dir}/ok.${TUMOR_SM}_${biosoft}.status 
    else
    touch ./${dir}/fail.${TUMOR_SM}_${biosoft}.status
    fi
    # sleep 60s
fi

pre_dir=5.gatk
dir=5.gatk/stats
biosoft=samtools_stats
if [  ! -f  ./${dir}/ok.${TUMOR_SM}_${biosoft}.status ]; 
then

    samtools stats -@ $NT                      \
    ./${pre_dir}/${TUMOR_SM}.BQSR.bam          \
    > ./${dir}/${TUMOR_SM}_BQSR_stats.txt

    if [ $? -eq 0 ]
    then
    touch ./${dir}/ok.${TUMOR_SM}_${biosoft}.status 
    else
    touch ./${dir}/fail.${TUMOR_SM}_${biosoft}.status
    fi
    # sleep 60s
fi

pre_dir=5.gatk
dir=5.gatk/mosdepth
biosoft=mosdepth
if [  ! -f  ./${dir}/ok.${TUMOR_SM}_${biosoft}.status ]; 
then

    mosdepth                                         \
    --threads $NT                                    \
    --no-per-base                                    \
    --fast-mode                                      \
    --by ${TARGET_BED}                               \
    --thresholds 0,1,2,3,4,5,10,20,30,50,100,200,300 \
    ${dir}/${TUMOR_SM}_bqsr                                 \
    ./${pre_dir}/${TUMOR_SM}.BQSR.bam

    if [ $? -eq 0 ]
    then
    touch ./${dir}/ok.${TUMOR_SM}_${biosoft}.status 
    else
    touch ./${dir}/fail.${TUMOR_SM}_${biosoft}.status
    fi
    # sleep 60s
fi

pre_dir=5.gatk
dir=5.gatk/qualimap
biosoft=qualimap
if [ ! -f ./${dir}/ok.${TUMOR_SM}_${biosoft}.status ];
then

    ${qualimap}                            \
    bamqc                                  \
    -bam ./${pre_dir}/${TUMOR_SM}.BQSR.bam \
    -cl 100                                \
    -gd ${GENOME_VERSION}                  \
    -gff ${TARGET_BED}                     \
    --skip-duplicated                      \
    --java-mem-size=30G                    \
    -outdir ${dir}/${TUMOR_SM}             \
    -outfile ${TUMOR_SM}.report.pdf        \
    -outformat PDF:HTML

    if [ $? -eq 0 ]
    then
    touch ./${dir}/ok.${TUMOR_SM}_${biosoft}.status 
    else
    touch ./${dir}/fail.${TUMOR_SM}_${biosoft}.status
    fi
fi

pre_dir=5.gatk
dir=5.gatk/picard
biosoft=picard
if [  ! -f  ./${dir}/ok.${TUMOR_SM}_${biosoft}.status ]; 
then

    ${GATK} --java-options "-Xmx30G"  CollectHsMetrics \
        -R $FASTA                                      \
        -I ./${pre_dir}/${TUMOR_SM}.BQSR.bam           \
        -O ./${dir}/${TUMOR_SM}.hs_metrics.txt         \
        -BI ${PROBE_FILE}                              \
        -TI ${INTERVAL_FILE}                           \
        --CLIP_OVERLAPPING_READS false                 \
        --TMP_DIR ./tmp                                \
        1>./${dir}/${TUMOR_SM}_CollectHsMetrics.log 2>&1

    ${GATK} --java-options "-Xmx30G"  CollectMultipleMetrics \
        -R $FASTA                                            \
        -I ./${pre_dir}/${TUMOR_SM}.BQSR.bam                 \
        -O ./${dir}/${TUMOR_SM}.multiple_metrics             \
        --TMP_DIR ./tmp                                      \
        1>./${dir}/${TUMOR_SM}_CollectMultipleMetrics.log 2>&1 

    if [ $? -eq 0 ]
    then
    touch ./${dir}/ok.${TUMOR_SM}_${biosoft}.status 
    else
    touch ./${dir}/fail.${TUMOR_SM}_${biosoft}.status
    fi
    # sleep 60s
fi

dir=tmp_QC
biosoft=Statistical_QC
if [  ! -f  ./${dir}/ok.${TUMOR_SM}_${biosoft}.status ]; 
then

    python3                                                             \
    ${QC_script}                                                        \
    -j ./2.clean_fq/${TUMOR_SM}.fastp.json                              \
    -m ./5.gatk/picard/${TUMOR_SM}.multiple_metrics.insert_size_metrics \
    -r ./5.gatk/${TUMOR_SM}.metircs                                     \
    -s ${SAMPLE_INFO}                                                   \
    -o ${dir}                                                           \
    -b ./5.gatk/${TUMOR_SM}.BQSR.bam                                    \
    -u ./4.align/stats/${TUMOR_SM}_sort_stats.txt                       \
    --dup_metrics ./4.align/picard/${TUMOR_SM}.hs_metrics.txt           \
    --rmdup_metrics ./5.gatk/picard/${TUMOR_SM}.hs_metrics.txt

    if [ $? -eq 0 ]
    then
    touch ./${dir}/ok.${TUMOR_SM}_${biosoft}.status 
    else
    touch ./${dir}/fail.${TUMOR_SM}_${biosoft}.status
    fi
    # sleep 60s
fi

# ******************************************
# 5a. Somatic Variant Calling - Mutect2
# ******************************************
dir=6.mutect
biosoft=Mutect2
if [  ! -f  ./${dir}/ok.${TUMOR_SM}_${biosoft}.status ]; 
then

    ${GATK} --java-options "-Xmx30G" Mutect2 \
    -R $FASTA                                \
    -I ./5.gatk/${TUMOR_SM}.BQSR.bam         \
    -L ${INTERVAL_FILE}                      \
    --germline-resource ${GERMLINE_VCF}      \
    -O ./${dir}/${TUMOR_SM}_mutect2.vcf.gz   \
    --tmp-dir ./tmp                          \
    1>./${dir}/${TUMOR_SM}_mutect.log 2>&1
    
    ${GATK} --java-options "-Xmx30G" FilterMutectCalls \
    -R $FASTA                                          \
    -V ./${dir}/${TUMOR_SM}_mutect2.vcf.gz             \
    -O ./${dir}/${TUMOR_SM}_somatic.vcf.gz             \
    --tmp-dir ./tmp                                    \
    1>>./${dir}/${TUMOR_SM}_mutect.log 2>&1
    
    ${GATK} --java-options "-Xmx30G"  SelectVariants \
    -R $FASTA                                        \
    -V ${dir}/${TUMOR_SM}_somatic.vcf.gz             \
    -select "vc.isNotFiltered()"                     \
    -O ${dir}/${TUMOR_SM}_filter.vcf                 \
    --tmp-dir ./tmp                                  \
    1>>./${dir}/${TUMOR_SM}_mutect.log 2>&1

    if [ $? -eq 0 ]
    then
    touch ./${dir}/ok.${TUMOR_SM}_${biosoft}.status 
    else
    touch ./${dir}/fail.${TUMOR_SM}_${biosoft}.status
    fi
    # sleep 60s
fi

# ******************************************
# 5b. Germline Variant Calling - HaplotypeCaller
# ******************************************
dir=6.gvcf
biosoft=HaplotypeCaller
if [  ! -f  ./${dir}/ok.${TUMOR_SM}_${biosoft}.status ]; 
then
    
    ${GATK} --java-options "-Xmx30G" HaplotypeCaller  \
      -R $FASTA                        \
      -I ./5.gatk/${TUMOR_SM}.BQSR.bam \
      -O ./${dir}/${TUMOR_SM}_0raw.vcf.gz   \
      --dbsnp ${KNOWN_DBSNP}           \
      --intervals ${INTERVAL_FILE}     \
      --tmp-dir ./tmp                  \
      1>./${dir}/${TUMOR_SM}_HC.log 2>&1

    # SNP\MNP

    ${GATK} --java-options "-Xmx30G"  SelectVariants \
    -R $FASTA                                         \
    -V ./${dir}/${TUMOR_SM}_0raw.vcf.gz                    \
    --select-type-to-include SNP                      \
    --select-type-to-include MNP                      \
    -O ${dir}/${TUMOR_SM}_1raw_SNP.vcf                 \
    --tmp-dir ./tmp                                  \
    1>>./${dir}/${TUMOR_SM}_HC.log 2>&1

    ${GATK} --java-options "-Xmx30G"  VariantFiltration \
    -R $FASTA                                           \
    -V ${dir}/${TUMOR_SM}_1raw_SNP.vcf                   \
    -filter "QD < 2.0" --filter-name "SNP_QD"           \
    -filter "MQ < 40.0" --filter-name "SNP_MQ"          \
    -filter "FS > 60.0" --filter-name "SNP_FS"          \
    -filter "SOR > 3.0" --filter-name "SNP_SOR"         \
    -filter "MQRankSum < -12.5" --filter-name "SNP_MQRankSum"          \
    -filter "ReadPosRankSum < -8.0" --filter-name "SNP_ReadPosRankSum" \
    -filter "QUAL < 30.0" --filter-name "SNP_QUAL"    \
    -filter "DP < 10.0" --filter-name "SNP_DP"        \
    -O ${dir}/${TUMOR_SM}_2filter_SNP.vcf              \
    --tmp-dir ./tmp                                  \
    1>>./${dir}/${TUMOR_SM}_HC.log 2>&1

    ${GATK} --java-options "-Xmx30G"  SelectVariants \
    -R $FASTA                                        \
    -V ${dir}/${TUMOR_SM}_2filter_SNP.vcf            \
    -select "vc.isNotFiltered()"                     \
    -O ${dir}/${TUMOR_SM}_SNP_filtered.vcf           \
    --tmp-dir ./tmp                                  \
    1>>./${dir}/${TUMOR_SM}_HC.log 2>&1

    # INDEL

    ${GATK} --java-options "-Xmx30G"  SelectVariants \
    -R $FASTA                                        \
    -V ./${dir}/${TUMOR_SM}_0raw.vcf.gz \
    --select-type-to-include INDEL                     \
    -O ${dir}/${TUMOR_SM}_1raw_INDEL.vcf                 \
    --tmp-dir ./tmp                                  \
    1>>./${dir}/${TUMOR_SM}_HC.log 2>&1

    ${GATK} --java-options "-Xmx30G"  VariantFiltration \
    -R $FASTA                                           \
    -V ${dir}/${TUMOR_SM}_1raw_INDEL.vcf                   \
    -filter "QD < 2.0" --filter-name "INDEL_QD" \
    -filter "ReadPosRankSum < -20.0" --filter-name "INDEL_ReadPosRankSum" \
    -filter "InbreedingCoeff < -0.8" --filter-name "INDEL_InbreedingCoeff" \
    -filter "FS > 200.0" --filter-name "INDEL_FS" \
    -filter "SOR > 10.0" --filter-name "INDEL_SOR" \
    -filter "QUAL < 30.0" --filter-name "INDEL_QUAL" \
    -filter "DP < 10.0" --filter-name "INDEL_DP"       \
    -O ${dir}/${TUMOR_SM}_2filter_INDEL.vcf              \
    --tmp-dir ./tmp                                     \
    1>>./${dir}/${TUMOR_SM}_HC.log 2>&1

    ${GATK} --java-options "-Xmx30G"  SelectVariants \
    -R $FASTA                                        \
    -V ${dir}/${TUMOR_SM}_2filter_INDEL.vcf            \
    -select "vc.isNotFiltered()"                     \
    -O ${dir}/${TUMOR_SM}_INDEL_filtered.vcf           \
    --tmp-dir ./tmp                                  \
    1>>./${dir}/${TUMOR_SM}_HC.log 2>&1

    if [ $? -eq 0 ]
    then
    touch ./${dir}/ok.${TUMOR_SM}_${biosoft}.status 
    else
    touch ./${dir}/fail.${TUMOR_SM}_${biosoft}.status
    fi
    # sleep 60s
fi

# ******************************************
# 6a. Somatic Variant Annotation - VEP/vcf2maf
# ******************************************
dir=7.annotation/vep
biosoft=vep
if [  ! -f  ./${dir}/ok.${TUMOR_SM}_${biosoft}.status ]; 
then

    perl ${vcf2maf}             \
    --input-vcf  6.mutect/${TUMOR_SM}_filter.vcf  \
    --output-maf ./${dir}/${TUMOR_SM}_vep.maf \
    --ref-fasta  ${FASTA}      \
    --tumor-id   ${TUMOR_SM}       \
    --vep-path   ${vep_path}    \
    --vep-data   ${vep_data}    \
    --species    homo_sapiens   \
    --ncbi-build GRCh37         \
    1>./${dir}/${TUMOR_SM}_vep.log 2>&1 

    if [ $? -eq 0 ]
    then
    touch ./${dir}/ok.${TUMOR_SM}_${biosoft}.status 
    else
    touch ./${dir}/fail.${TUMOR_SM}_${biosoft}.status
    fi
    # sleep 60s
fi

# ******************************************
# 6b. Somatic Variant Annotation - ANNOVAR
# ******************************************
dir=7.annotation/annovar
biosoft=annovar
if [  ! -f  ./${dir}/ok.${TUMOR_SM}_${biosoft}.status ]; 
then

    perl ${annovar_path}/table_annovar.pl   \
    ./6.mutect/${TUMOR_SM}_filter.vcf       \
    ${annovar_data}/                     \
    -buildver  ${GENOME_VERSION}         \
    -otherinfo -remove                   \
    -nastring .                          \
    -vcfinput                            \
    -protocol refGene,genomicSuperDups,avsnp150,1000g2015aug_all,1000g2015aug_eas,cg69,esp6500siv2_all,exac03nonpsych,gnomad211_exome,gnomad211_genome,cosmic96_coding,clinvar_20220528,intervar_20180118,VarCards,revel,dbscsnv11,dbnsfp42c \
    -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f                   \
    -out ./7.annotation/annovar/${TUMOR_SM}  \
    1>7.annotation/annovar/${TUMOR_SM}_annovar.log 2>&1 

    if [ $? -eq 0 ]
    then
    touch ./${dir}/ok.${TUMOR_SM}_${biosoft}.status 
    else
    touch ./${dir}/fail.${TUMOR_SM}_${biosoft}.status
    fi
    # sleep 60s
fi

# ******************************************
# 6c. Germline Variant Annotation - ANNOVAR
# ******************************************
predir=6.gvcf
dir=7.annotation/annovar
biosoft=annovar_germline
if [  ! -f  ./${dir}/ok.${TUMOR_SM}_${biosoft}.status ]; 
then

    bcftools \
    norm -m-both \
    -o ${dir}/${TUMOR_SM}_germline_SNP.step1.vcf \
    ${predir}/${TUMOR_SM}_SNP_filtered.vcf \
    1>${dir}/${TUMOR_SM}_annovar.log 2>&1 

    bcftools \
    norm -f ${FASTA} \
    -o ${dir}/${TUMOR_SM}_germline_SNP.step2.vcf \
    ${dir}/${TUMOR_SM}_germline_SNP.step1.vcf \
    1>>${dir}/${TUMOR_SM}_annovar.log 2>&1 

    perl ${annovar_path}/table_annovar.pl   \
    ${dir}/${TUMOR_SM}_germline_SNP.step2.vcf  \
    ${annovar_data}/                     \
    -buildver  ${GENOME_VERSION}         \
    -otherinfo -remove                   \
    -nastring .                          \
    -vcfinput                            \
    -protocol refGene,genomicSuperDups,avsnp150,1000g2015aug_all,1000g2015aug_eas,esp6500siv2_all,exac03nonpsych,gnomad211_exome,clinvar_20230813,intervar_20180118,dbscsnv11,dbnsfp42a,gene4denovo201907 \
    -operation g,r,f,f,f,f,f,f,f,f,f,f,f                  \
    -out ${dir}/${TUMOR_SM}_germline_SNP  \
    1>>${dir}/${TUMOR_SM}_annovar.log 2>&1

    bcftools \
    norm -m-both \
    -o ${dir}/${TUMOR_SM}_germline_INDEL.step1.vcf \
    ${predir}/${TUMOR_SM}_INDEL_filtered.vcf \
    1>>${dir}/${TUMOR_SM}_annovar.log 2>&1 

    bcftools \
    norm -f ${FASTA} \
    -o ${dir}/${TUMOR_SM}_germline_INDEL.step2.vcf \
    ${dir}/${TUMOR_SM}_germline_INDEL.step1.vcf \
    1>>${dir}/${TUMOR_SM}_annovar.log 2>&1 

    perl ${annovar_path}/table_annovar.pl   \
    ${dir}/${TUMOR_SM}_germline_INDEL.step2.vcf  \
    ${annovar_data}/                     \
    -buildver  ${GENOME_VERSION}         \
    -otherinfo -remove                   \
    -nastring .                          \
    -vcfinput                            \
    -protocol refGene,genomicSuperDups,avsnp150,1000g2015aug_all,1000g2015aug_eas,esp6500siv2_all,exac03nonpsych,gnomad211_exome,clinvar_20230813,intervar_20180118,dbscsnv11,dbnsfp42a,gene4denovo201907 \
    -operation g,r,f,f,f,f,f,f,f,f,f,f,f                  \
    -out ${dir}/${TUMOR_SM}_germline_INDEL  \
    1>>${dir}/${TUMOR_SM}_annovar.log 2>&1

    if [ $? -eq 0 ]
    then
    touch ./${dir}/ok.${TUMOR_SM}_${biosoft}.status 
    else
    touch ./${dir}/fail.${TUMOR_SM}_${biosoft}.status
    fi
    # sleep 60s
fi