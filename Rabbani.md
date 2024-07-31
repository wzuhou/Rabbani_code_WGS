# **Code Availability**

**Author**: MAG Rabbani

**Affilliation**: The Roslin institute


Data analyses were completed using standard bioinformatic tools running on the Scientific Linux 7 system. The version and code/parameters of the main software tools are as described below:

**_Step_1:_** **Quality control of FASTQ files (FastQC-v0.11.7)**

```sh
fastqc -t 1 ${READs}.fastq.gz -o ${READs}
```

**_Step_2:_** **Mapping/Alignment of raw reads against reference genome**
**a.** **Mapping reads with BWA (bwa-v0.7.15)**

```sh
$bwa mem -t 4 -M -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:Illumina\tLB:${SAMPLE}\tPU:unkn-0.0" ${REF_GENOME} ${READS_1} ${READS_2} > ${SAMPLE}.sam
```

**b.** **Sort SAM into coordinated order and save as BAM (picard-v2.25.4)**

```sh
picard SortSam \
    I=${SAMPLE}.sam \
    O=${SAMPLE}_sorted.bam \
    SORT_ORDER=coordinate \
    TMP_DIR= tmp_${SAMPLE}
```

**c.** **Mark duplicates and create bam index (picard-v2.25.4)**

```sh
picard MarkDuplicates \
    I=${SAMPLE}_sorted.bam \
    O=${SAMPLE}_mdup.bam \
    CREATE_INDEX=true \
    M=metrics/${SAMPLE}_mdup_metrics.txt \
    TMP_DIR=tmp_${SAMPLE} \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=4000 \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500
```

**d.** **Validate BAM file (picard-v2.25.4)**

```sh
picard ValidateSamFile \
    I= ${SAMPLE}_mdup.bam \
    MODE=SUMMARY \
    MAX_OPEN_TEMP_FILES=8000 \
    TMP_DIR=tmp_${SAMPLE}
```

**e.** **Flagstat of BAM file (samtools-v1.13)**

```sh
$samtools flagstat ${SAMPLE}_sorted.bam
$samtools flagstat ${SAMPLE}_mdup.bam
```

**_Step_3:_** **BQSR using GATK and Picard tools**
**a. Analyze of the patterns of covariation in the sequence dataset (GATK-v4.0.10.1)**

```sh
gatk BaseRecalibrator \
    -R ${REF_GENOME} \
    -I ${SAMPLE}_mdup.bam \
    --known-sites ${VCF} \
    --output ${SAMPLE}_recal_data.table
```

**b. A second pass to analyze covariation post-recalibration (GATK-v4.0.10.1)**

```sh
gatk BaseRecalibrator \
    -R ${REF_GENOME} \
    -I {SAMPLE}_mdup.bam \
    --known-sites ${VCF} \
    --output ${SAMPLE}_post_recal_data.table
```

**c. Generate before/after plots (GATK-v4.0.10.1)**

```sh
gatk AnalyzeCovariates \
    -before ${SAMPLE}_recal_data.table \
    -after ${SAMPLE}_post_recal_data.table \
    -plots ${SAMPLE}_recalibration_plots.pdf
```

**d. Apply the recalibration to your sequence data (GATK-v4.0.10.1)**

```sh
gatk PrintReads \
    -R ${REF_GENOME} \
    -I ${SAMPLE}_mdup.bam \
    -O ${SAMPLE}_recal.bam
```

**e. Validate the recalibrated BAM file** **(picard-v2.25.4)**

```sh
picard ValidateSamFile \
    I= ${SAMPLE}_recal.bam \
    MODE=SUMMARY \
    TMP_DIR=tmp_${SAMPLE} \
    MAX_OPEN_TEMP_FILES=4000
```

**f. Flagstat of recalibrated BAM file** **(samtools-v1.13)**

```sh
$samtools flagstat ${SAMPLE}_recal.bam
```

**g. Additional step: Insert Size metrics** **(picard-v2.25.4)**

```sh
picard CollectInsertSizeMetrics \
    I=${SAMPLE}_mdup.bam \
    O=${SAMPLE}_mdup_insertSize_metrics.txt \
    HISTOGRAM_FILE=${SAMPLE}_mdup_insertSize_metrics.pdf
```

**_Step_4:_** **Variant calling using GATK**
**Call variants in gVCF mode for cohort analysis (GATK-v4.0.10.1)**

```sh
gatk HaplotypeCaller \
    -R ${REF_GENOME} \
    -I ${SAMPLE}_recal.bam \
    -O ${SAMPLE}.g.vcf.gz \
    -ERC GVCF
```

**_Step_5:_** **Joint genotyping of a cohort of samples:**
**a.** **GenomicsDBImport (GATK-v4.0.10.1):**

```sh
gatk GenomicsDBImport \
    --genomicsdb-workspace-path my_database \
    --intervals interval.bed \
    -V ${SAMPLE}1.g.vcf.gz \
    -V ${SAMPLE}3.g.vcf.gz \
    -V ${SAMPLE}3.g.vcf.gz
```

**b.** **GenotypeGVCFs (GATK-v4.0.10.1)**

```sh
gatk GenotypeGVCFs \
    -R ${REF_GENOME} \
    -V gendb://my_database \
    -O ${SAMPLE}.vcf.gz
```

**_Step_6:_** **VQSR**
**a.** **VariantRecalibrator (GATK-v4.0.10.1)**

```sh
gatk VariantRecalibrator \
    -R ${REF_GENOME} \
    -V ${SAMPLE}.vcf.gz \
    --resource GRCg7b_dbSNP,known=true,training=false,truth=false,prior=2.0:${KNOWNVAR} \
    --resource GCRg6a_validated_snp,known=false,training=true,truth=true,prior=12.0:${TRUEVAR} \
    -an DP -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode SNP \
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
    -O ${SAMPLE}.SNPs.recal.gz \
    --tranches-file ${SAMPLE}.SNPs.tranches \
    --rscript-file ${SAMPLE}_recalSNPS.plots.R
```

**b.** **ApplyVQSR (GATK-v4.0.10.1)**

```sh
gatk ApplyVQSR \
    -R ${REF_GENOME} \
    -V ${SAMPLE}.vcf.gz \
    --mode SNP \
    --truth-sensitivity-filter-level 99.0 \
    --recal-file ${SAMPLE}.SNPs.recal.gz \
    --tranches-file${SAMPLE}.SNPs.tranches \
    -O ${SAMPLE}_recalSNPs.vcf.gz
```

**_Step_7:_** **Select Variants** **(GATK-v4.0.10.1)**

```sh
gatk SelectVariants \
    -R ${REF_GENOME}\
    -V ${SAMPLE}_recalSNPs.vcf.gz \
    --select-type SNP \
    --restrict-alleles-to BIALLELIC \
    -O ${SAMPLE}.vcf.gz
```

**_Step_8:_** **Filtration of variants (VCFtools-v0.1.13)**

```sh
vcftools \
    --gzvcf ${SAMPLE}.vcf.gz \
    --hwe 0.00001 \
    --max-missing 0.9 \
    --minGQ 20.0 \
    --minDP 3 \
    --min-alleles 2 \
    --max-alleles 2 \
    --recode \
    --recode-INFO-all \
    --out ${SAMPLE}_fil.vcf.gz
```
