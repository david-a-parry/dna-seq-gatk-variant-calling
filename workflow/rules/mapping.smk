rule trim_reads_se:
    input:
        unpack(get_fastq),
    output:
        temp("results/trimmed/{sample}-{unit}.fastq.gz"),
    params:
        **config["params"]["trimmomatic"]["se"],
        extra="",
    resources:
        mem_mb=8096,
    log:
        "logs/trimmomatic/{sample}-{unit}.log",
    wrapper:
        "0.84.0/bio/trimmomatic/se"


rule trim_reads_pe:
    input:
        unpack(get_fastq),
    output:
        r1=temp("results/trimmed/{sample}-{unit}.1.fastq.gz"),
        r2=temp("results/trimmed/{sample}-{unit}.2.fastq.gz"),
        r1_unpaired=temp("results/trimmed/{sample}-{unit}.1.unpaired.fastq.gz"),
        r2_unpaired=temp("results/trimmed/{sample}-{unit}.2.unpaired.fastq.gz"),
    params:
        **config["params"]["trimmomatic"]["pe"],
        extra="",
    resources:
        mem_mb=8096,
    log:
        "logs/trimmomatic/{sample}-{unit}.log",
    wrapper:
        "0.84.0/bio/trimmomatic/pe"


rule map_reads:
    input:
        reads=get_trimmed_reads,
        idx=rules.bwa_index.output,
    output:
        temp("results/mapped/{sample}-{unit}.sorted.bam"),
    log:
        "logs/bwa_mem/{sample}-{unit}.log",
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=get_read_group,
        sorting="samtools",
        sort_order="coordinate",
    threads: 8
    wrapper:
        "v1.5.0/bio/bwa/mem"


rule mark_duplicates:
    input:
        "results/mapped/{sample}-{unit}.sorted.bam",
    output:
        bam=temp("results/dedup/{sample}-{unit}.bam"),
        metrics="results/qc/dedup/{sample}-{unit}.metrics.txt",
    log:
        "logs/picard/dedup/{sample}-{unit}.log",
    params:
        config["params"]["picard"]["MarkDuplicates"],
    resources:
        mem_mb=4096,
    wrapper:
        "0.84.0/bio/picard/markduplicates"


rule recalibrate_base_qualities:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref="resources/genome.fasta",
        dict="resources/genome.dict",
        known=get_variation_vcf(),
        known_idx=get_variation_vcf() + ".tbi",
    output:
        recal_table="results/recal/{sample}-{unit}.grp",
    log:
        "logs/gatk/bqsr/{sample}-{unit}.log",
    params:
        extra=get_regions_param() + config["params"]["gatk"]["BaseRecalibrator"],
    resources:
        mem_mb=1024,
    wrapper:
        "0.84.0/bio/gatk/baserecalibrator"


rule apply_base_quality_recalibration:
    input:
        bam=get_recal_input(),
        bai=get_recal_input(bai=True),
        ref="resources/genome.fasta",
        dict="resources/genome.dict",
        recal_table="results/recal/{sample}-{unit}.grp",
    output:
        bam=protected("results/recal/{sample}-{unit}.bam"),
    log:
        "logs/gatk/apply-bqsr/{sample}-{unit}.log",
    params:
        extra=get_regions_param(),
    resources:
        mem_mb=1024,
    wrapper:
        "0.84.0/bio/gatk/applybqsr"


rule samtools_index:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    log:
        "logs/samtools/index/{prefix}.log",
    wrapper:
        "0.84.0/bio/samtools/index"
