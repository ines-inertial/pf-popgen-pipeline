#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ==================== PROCESSES =======================
// ------------------------------------------------------
// Step 1: FASTQC + MULTIQC on RAW reads
// ------------------------------------------------------
process FASTQC_RAW {
    tag "$sample_id"
    cpus 8
    conda "/hps/software/users/jlees/ines/miniforge3/envs/fastqc"
    publishDir "${params.outdir}/fastqc_raw", mode: 'copy'
    
    input:
    tuple val(sample_id), path(r1), path(r2)
    output:
    path "*.{zip,html}", emit: reports
    tuple val(sample_id), path(r1), path(r2), emit: reads
    script:
    """
    fastqc -t $task.cpus "$r1" "$r2"
    """
}

process MULTIQC_RAW {
    cpus 8
    conda "/hps/software/users/jlees/ines/miniforge3/envs/fastqc"
    publishDir "${params.outdir}/multiqc_raw", mode: 'copy'
    
    input:
    path fastqc_reports
    output:
    path "multiqc_raw_report.html", emit: report
    script:
    """
    multiqc . -n multiqc_raw_report.html
    """
}

// ------------------------------------------------------
// Step 2: FASTP TRIMMING on RAW reads
// ------------------------------------------------------
process FASTP_TRIM {
    tag "$sample_id"
    cpus 8
    conda "/hps/software/users/jlees/ines/miniforge3/envs/trim"
    publishDir "${params.outdir}/fastp_output", mode: 'copy'
    
    input:
    tuple val(sample_id), path(r1), path(r2)
    output:
    tuple val(sample_id), path("${sample_id}_1_trimmed.fastq.gz"), path("${sample_id}_2_trimmed.fastq.gz"), emit: trimmed_reads
    path "${sample_id}_fastp.html", emit: fastp_html
    path "${sample_id}_fastp.json", emit: fastp_json
    script:
    """
    fastp \
        -i "$r1" \
        -I "$r2" \
        -o "${sample_id}_1_trimmed.fastq.gz" \
        -O "${sample_id}_2_trimmed.fastq.gz" \
        -h "${sample_id}_fastp.html" \
        -j "${sample_id}_fastp.json" \
        --thread $task.cpus \
        --detect_adapter_for_pe \
        --correction \
        --cut_tail \
        --cut_front \
        --qualified_quality_phred 20 \
        --length_required 50
    """
}

// ------------------------------------------------------
// Step 3: FASTQC + MULTIQC on TRIMMED reads
// ------------------------------------------------------
process FASTQC_TRIM {
    tag "$sample_id"
    cpus 8
    conda "/hps/software/users/jlees/ines/miniforge3/envs/fastqc"
    publishDir "${params.outdir}/fastqc_trim", mode: 'copy'
    
    input:
    tuple val(sample_id), path(r1), path(r2)
    output:
    path "*.{zip,html}", emit: reports
    tuple val(sample_id), path(r1), path(r2), emit: reads
    script:
    """
    fastqc -t $task.cpus "$r1" "$r2"
    """
}

process MULTIQC_TRIM {
    cpus 8
    conda "/hps/software/users/jlees/ines/miniforge3/envs/fastqc"
    publishDir "${params.outdir}/multiqc_trim", mode: 'copy'
    
    input:
    path fastqc_reports
    output:
    path "multiqc_trim_report.html", emit: report
    script:
    """
    multiqc . -n multiqc_trim_report.html
    """
}

// ------------------------------------------------------
// Step 4: KRAKEN2 classification (host reads filtering)
// ------------------------------------------------------
process KRAKEN2_CLASSIFY {
    tag "$sample_id"
    cpus 8
    conda "/hps/software/users/jlees/ines/miniforge3/envs/kraken2"
    publishDir "${params.outdir}/kraken2_output", mode: 'copy'
    
    input:
    tuple val(sample_id), path(r1), path(r2)
    output:
    tuple val(sample_id), path("${sample_id}_classified_1.fastq.gz"), path("${sample_id}_classified_2.fastq.gz"), emit: classified_reads
    tuple val(sample_id), path("${sample_id}_unclassified_1.fastq.gz"), path("${sample_id}_unclassified_2.fastq.gz"), emit: unclassified_reads
    path "${sample_id}_kraken2_output.txt", emit: outputs
    path "${sample_id}_kraken2_report.txt", emit: reports
    script:
    """
    kraken2 \
        --db "${params.kraken_db}" \
        --threads $task.cpus \
        --paired \
        --classified-out "${sample_id}_classified#.fastq.gz" \
        --unclassified-out "${sample_id}_unclassified#.fastq.gz" \
        --output "${sample_id}_kraken2_output.txt" \
        --report "${sample_id}_kraken2_report.txt" \
        "$r1" "$r2"
    """
}

// ------------------------------------------------------
// Step 5: BWA-MEM Alignment (Add Read groups) + Sort & Index on classified reads
// ------------------------------------------------------
process BWA_MEM_MAP {
    tag "$sample_id"
    cpus 8
    conda "/hps/software/users/jlees/ines/miniforge3/envs/assembly"
    publishDir "${params.outdir}/mapping_output", mode: 'copy'
    
    input:
    tuple val(sample_id), path(r1), path(r2)
    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: sorted_bam
    script:
    def rg = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tLB:${sample_id}\\tPL:${params.pl}\\tPM:${params.pm}"
    """
    bwa mem -t $task.cpus -R "$rg" "${params.reference_map}" "$r1" "$r2" | \
    samtools sort -@ $task.cpus -o "${sample_id}.sorted.bam" -
    samtools index "${sample_id}.sorted.bam"
    """
}

// ------------------------------------------------------
// Step 6: Samtools QC + MultiQC on Mapped reads
// ------------------------------------------------------
process ALIGNMENT_QC_MAP {
    tag "$sample_id"
    cpus 8
    conda "/hps/software/users/jlees/ines/miniforge3/envs/assembly"
    publishDir "${params.outdir}/alignment_qc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    output:
    path "${sample_id}.flagstat.txt", emit: flagstat
    path "${sample_id}.coverage.txt", emit: coverage
    path "${sample_id}_qualimap", emit: qualimap_dir
    script:
    """
    samtools flagstat "$bam" > "${sample_id}.flagstat.txt"
    samtools coverage "$bam" > "${sample_id}.coverage.txt"
    qualimap bamqc -bam "$bam" -outdir "${sample_id}_qualimap" -nt $task.cpus
    """
}

process MULTIQC_MAP {
    cpus 8
    conda "/hps/software/users/jlees/ines/miniforge3/envs/fastqc"
    publishDir "${params.outdir}/multiqc_map", mode: 'copy'
    input:
    path alignmentqc_reports
    output:
    path "multiqc_map_report.html", emit: report
    script:
    """
    multiqc . -n multiqc_map_report.html
    """
}

// ------------------------------------------------------
// Step 7: Mark & Remove Duplicates
// ------------------------------------------------------
process MARK_DUPLICATES {
    tag "$sample_id"
    cpus 8
    conda "/hps/software/users/jlees/ines/miniforge3/envs/variantcall"
    publishDir "${params.outdir}/markdedup_output", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam"), path("${sample_id}.dedup*.bai"), emit: dedup_bam
    path "${sample_id}.dedup_metrics.txt", emit: metrics
    script:
    """
    picard MarkDuplicates \
        INPUT="$bam" \
        OUTPUT="${sample_id}.dedup.bam" \
        METRICS_FILE="${sample_id}.dedup_metrics.txt" \
        REMOVE_DUPLICATES=true \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=SILENT \
        TMP_DIR="${params.tmpdir}"
    """
}

// ------------------------------------------------------
// Step 8: GATK HaplotypeCaller
// ------------------------------------------------------
process HAPLOTYPE_CALLER {
    tag "$sample_id"
    cpus 8
    conda "/hps/software/users/jlees/ines/miniforge3/envs/variantcall"
    publishDir "${params.outdir}/variants/snpindels", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    output:
    tuple val(sample_id), path("${sample_id}_snpindels.g.vcf.gz"), path("${sample_id}_snpindels.g.vcf.gz.tbi"), emit: gvcf
    script:
    """
    gatk HaplotypeCaller \
        -R "${params.reference_variant}" \
        -I "$bam" \
        -O "${sample_id}_snpindels.g.vcf.gz" \
        --native-pair-hmm-threads $task.cpus \
        -ploidy 1 \
        -ERC GVCF \
        --tmp-dir "${params.tmpdir}"
    """
}

// ------------------------------------------------------
// Step 9: Genotype GVCF files
// ------------------------------------------------------
process GENOTYPE_GVCF {
    tag "$sample_id"
    cpus 8
    conda "/hps/software/users/jlees/ines/miniforge3/envs/variantcall"
    publishDir "${params.outdir}/variants/genotyped", mode: 'copy'
    
    input:
    tuple val(sample_id), path(gvcf), path(gvcf_idx)
    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi"), emit: vcf
    script:
    """
    gatk GenotypeGVCFs \
        -R "${params.reference_variant}" \
        -V "$gvcf" \
        -O "${sample_id}.vcf.gz"
    """
}

// ------------------------------------------------------
// Step 10: Select SNPs & INDELs (from GATK Variants)
// ------------------------------------------------------
process SELECT_VARIANTS {
    tag "$sample_id"
    cpus 8
    conda "/hps/software/users/jlees/ines/miniforge3/envs/variantcall"
    publishDir "${params.outdir}/selected", mode: 'copy'
    
    input:
    tuple val(sample_id), path(vcf), path(vcf_idx)
    output:
    tuple val(sample_id), path("${sample_id}_raw_snps.vcf.gz"), path("${sample_id}_raw_snps.vcf.gz.tbi"), emit: raw_snps
    tuple val(sample_id), path("${sample_id}_raw_indels.vcf.gz"), path("${sample_id}_raw_indels.vcf.gz.tbi"), emit: raw_indels
    script:
    """
    gatk SelectVariants \
        -R "${params.reference_variant}" \
        -V "$vcf" \
        --select-type-to-include SNP \
        -O "${sample_id}_raw_snps.vcf.gz"
    gatk SelectVariants \
        -R "${params.reference_variant}" \
        -V "$vcf" \
        --select-type-to-include INDEL \
        -O "${sample_id}_raw_indels.vcf.gz"
    """
}

// ------------------------------------------------------
// Step 11: Filter SNPs & INDELs (from GATK Variants)
// ------------------------------------------------------
process FILTER_SNPs {
    tag "$sample_id"
    cpus 8
    conda "/hps/software/users/jlees/ines/miniforge3/envs/variantcall"
    publishDir "${params.outdir}/variants/filtered_snps", mode: 'copy'
    
    input:
    tuple val(sample_id), path(raw_snps), path(raw_snps_idx)
    output:
    tuple val(sample_id), path("${sample_id}_filtered_snps.vcf.gz"), path("${sample_id}_filtered_snps.vcf.gz.tbi"), emit: filtered_snps
    script:
    """
    gatk VariantFiltration \
        -R "${params.reference_variant}" \
        -V "$raw_snps" \
        -O "${sample_id}_filtered_snps.vcf.gz" \
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
        --filter-name "LOWQUAL"
    gatk IndexFeatureFile -I "${sample_id}_filtered_snps.vcf.gz"
    """
}

process FILTER_INDELs {
    tag "$sample_id"
    cpus 8
    conda "/hps/software/users/jlees/ines/miniforge3/envs/variantcall"
    publishDir "${params.outdir}/variants/filtered_indels", mode: 'copy'
    
    input:
    tuple val(sample_id), path(raw_indels), path(raw_indels_idx)
    output:
    tuple val(sample_id), path("${sample_id}_filtered_indels.vcf.gz"), path("${sample_id}_filtered_indels.vcf.gz.tbi"), emit: filtered_indels
    script:
    """
    gatk VariantFiltration \
        -R "${params.reference_variant}" \
        -V "$raw_indels" \
        -O "${sample_id}_filtered_indels.vcf.gz" \
        --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
        --filter-name "LOWQUAL"
    gatk IndexFeatureFile -I "${sample_id}_filtered_indels.vcf.gz"
    """
}

// ------------------------------------------------------
// Step 12: Mask Hypervariable regions from Variants
// ------------------------------------------------------
process MASK_SNP_VARIANTS {
    tag "$sample_id"
    cpus 4
    conda "/hps/software/users/jlees/ines/miniforge3/envs/variantcall"
    publishDir "${params.outdir}/masked_snps", mode: 'copy'
    
    input:
    tuple val(sample_id), path(filtered_snps), path(filtered_snps_idx)
    output:
    tuple val(sample_id), path("${sample_id}_masked_snps.vcf.gz"), path("${sample_id}_masked_snps.vcf.gz.tbi"), emit: masked_snps
    script:
    """
    bcftools view -T "${params.mask_bed}" -Oz -o ${sample_id}_masked_snps.vcf.gz ${filtered_snps}
    bcftools index -t ${sample_id}_masked_snps.vcf.gz
    """
}

process MASK_INDEL_VARIANTS {
    tag "$sample_id"
    cpus 4
    conda "/hps/software/users/jlees/ines/miniforge3/envs/variantcall"
    publishDir "${params.outdir}/masked_indels", mode: 'copy'
    
    input:
    tuple val(sample_id), path(filtered_indels), path(filtered_indels_idx)
    output:
    tuple val(sample_id), path("${sample_id}_masked_indels.vcf.gz"), path("${sample_id}_masked_indels.vcf.gz.tbi"), emit: masked_indels
    script:
    """
    bcftools view -T "${params.mask_bed}" -Oz -o ${sample_id}_masked_indels.vcf.gz ${filtered_indels}
    bcftools index -t ${sample_id}_masked_indels.vcf.gz
    """
}

// ------------------------------------------------------
// Step 13: Annotate Variants (from GATK)
// ------------------------------------------------------
process SNPEFF_Annotate_SNPs {
    tag "$sample_id"
    cpus 8
    conda "/hps/software/users/jlees/ines/miniforge3/envs/snpeff"
    publishDir "${params.outdir}/annotated_snps", mode: 'copy'
    
    input:
    tuple val(sample_id), path(masked_snps), path(masked_snps_idx)
    output:
    tuple val(sample_id), path("${sample_id}_snps.ann.vcf.gz"), path("${sample_id}_snps.ann.vcf.gz.tbi"), emit: ann_snps
    script:
    """
    # Run snpEff annotation (using built-in Plasmodium falciparum database)
    snpEff -v "${params.snpeff_db}" "$masked_snps" > ${sample_id}_snps.ann.vcf
    # Compress & index
    bgzip -f ${sample_id}_snps.ann.vcf
    tabix -p vcf ${sample_id}_snps.ann.vcf.gz
    """
}

process SNPEFF_Annotate_INDELs {
    tag "$sample_id"
    cpus 8
    conda "/hps/software/users/jlees/ines/miniforge3/envs/snpeff"
    publishDir "${params.outdir}/annotated_indels", mode: 'copy'
    
    input:
    tuple val(sample_id), path(masked_indels), path(masked_indels_idx)
    output:
    tuple val(sample_id), path("${sample_id}_indels.ann.vcf.gz"), path("${sample_id}_indels.ann.vcf.gz.tbi"), emit: ann_indels
    script:
    """
    # Run snpEff annotation (using built-in Plasmodium falciparum database)
    snpEff -v "${params.snpeff_db}" "$masked_indels" > ${sample_id}_indels.ann.vcf
    # Compress & index
    bgzip -f ${sample_id}_indels.ann.vcf
    tabix -p vcf ${sample_id}_indels.ann.vcf.gz
    """
}

// ----------------------------------------------------------
// Step 14: Split annotated SNPs and INDELs by chromosome
// ----------------------------------------------------------
process SPLIT_ANNOTATED_SNPS_BY_CHROM {
    tag "${sample_id} - ${chrom}"
    cpus 8
    conda "/hps/software/users/jlees/ines/miniforge3/envs/variantcall"
    publishDir "${params.outdir}/per_chrom_snps", mode: 'copy'

    input:
    tuple val(sample_id), path(ann_snps), path(ann_snps_idx), val(chrom)
    output:
    path "${sample_id}_${chrom}_snps.ann.vcf.gz"
    script:
    """
    bcftools view -r ${chrom} ${ann_snps} -Oz -o ${sample_id}_${chrom}_snps.ann.vcf.gz
    tabix -p vcf ${sample_id}_${chrom}_snps.ann.vcf.gz
    """
}

process SPLIT_ANNOTATED_INDELS_BY_CHROM {
    tag "${sample_id} - ${chrom}"
    cpus 8
    conda "/hps/software/users/jlees/ines/miniforge3/envs/variantcall"
    publishDir "${params.outdir}/per_chrom_indels", mode: 'copy'

    input:
    tuple val(sample_id), path(ann_indels), path(ann_indels_idx), val(chrom)
    output:
    path "${sample_id}_${chrom}_indels.ann.vcf.gz"
    script:
    """
    bcftools view -r ${chrom} ${ann_indels} | \
    if grep -q "^#" ; then
        bcftools view -r ${chrom} -Oz -o ${sample_id}_${chrom}_indels.ann.vcf.gz ${ann_indels}
        tabix -p vcf ${sample_id}_${chrom}_indels.ann.vcf.gz
    else
        touch ${sample_id}_${chrom}_indels.ann.vcf.gz.empty_placeholder
    fi
    """
}

// ==================== WORKFLOW ====================
workflow {
    // ==================== CHANNELS ====================
    Channel
    .fromFilePairs(params.reads, size: 2)
    .map { sample_id, reads -> tuple(sample_id, reads[0], reads[1]) }
    .set { read_pairs_ch }

    chrom_ch = Channel.fromList(params.chroms)

    // Step 1: Raw QC
    FASTQC_RAW(read_pairs_ch)
    MULTIQC_RAW(FASTQC_RAW.out.reports.collect())

    // Step 2: Trimming (use raw reads pass-through)
    FASTP_TRIM(FASTQC_RAW.out.reads)

    // Step 3: Trimmed QC
    FASTQC_TRIM(FASTP_TRIM.out.trimmed_reads)
    MULTIQC_TRIM(FASTQC_TRIM.out.reports.collect())

    // Step 4: Classification (use trimmed reads pasis-through)
    KRAKEN2_CLASSIFY(FASTQC_TRIM.out.reads)

    // Step 5: Mapping
    BWA_MEM_MAP(KRAKEN2_CLASSIFY.out.classified_reads)

    // Step 6: Alignment QC
    ALIGNMENT_QC_MAP(BWA_MEM_MAP.out.sorted_bam)
    // Collect all QC files for MultiQC (flagstat + coverage + qualimap dirs)
    alignment_qc_files = ALIGNMENT_QC_MAP.out.flagstat
        .mix(ALIGNMENT_QC_MAP.out.coverage, ALIGNMENT_QC_MAP.out.qualimap_dir)
        .collect()
    MULTIQC_MAP(alignment_qc_files)

    // Step 7: Deduplication
    MARK_DUPLICATES(BWA_MEM_MAP.out.sorted_bam)

    // Step 8: Variant Calling
    HAPLOTYPE_CALLER(MARK_DUPLICATES.out.dedup_bam)

    // Step 9: Genotype Calls
    GENOTYPE_GVCF(HAPLOTYPE_CALLER.out.gvcf)

    // Step 10: Select Variants
    SELECT_VARIANTS(GENOTYPE_GVCF.out.vcf)

    // Step 11: Filter Variants
    FILTER_SNPs(SELECT_VARIANTS.out.raw_snps)
    FILTER_INDELs(SELECT_VARIANTS.out.raw_indels)

    // Step 12: Masked Variants
    MASK_SNP_VARIANTS(FILTER_SNPs.out.filtered_snps)
    MASK_INDEL_VARIANTS(FILTER_INDELs.out.filtered_indels)

    // Step 13: Annotate Variants
    SNPEFF_Annotate_SNPs(MASK_SNP_VARIANTS.out.masked_snps)
    SNPEFF_Annotate_INDELs(MASK_INDEL_VARIANTS.out.masked_indels)

    // Step 14: Split annotated SNPs by chromosome
    split_snps_input = SNPEFF_Annotate_SNPs.out.ann_snps.combine(chrom_ch)
    SPLIT_ANNOTATED_SNPS_BY_CHROM(split_snps_input)
    split_indels_input = SNPEFF_Annotate_INDELs.out.ann_indels.combine(chrom_ch)
    SPLIT_ANNOTATED_INDELS_BY_CHROM(split_indels_input)
}
