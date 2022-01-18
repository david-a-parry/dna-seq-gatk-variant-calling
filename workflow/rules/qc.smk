rule fastqc:
    input:
        unpack(get_fastq),
    output:
        html="results/qc/fastqc/{sample}-{unit}.html",
        zip="results/qc/fastqc/{sample}-{unit}.zip",
    log:
        "logs/fastqc/{sample}-{unit}.log",
    wrapper:
        "0.84.0/bio/fastqc"


rule samtools_stats:
    input:
        "results/recal/{sample}-{unit}.bam",
    output:
        "results/qc/samtools-stats/{sample}-{unit}.txt",
    log:
        "logs/samtools-stats/{sample}-{unit}.log",
    wrapper:
        "0.84.0/bio/samtools/stats"


if config["processing"].get("restrict-regions"):
    rule mosdepth_bed:
        input:
            bam=get_sample_bams,
            bed=config["processing"].get("restrict-regions", ""),
        output:
            "results/qc/mosdepth/{sample}.mosdepth.global.dist.txt",
            "results/qc/mosdepth/{sample}.mosdepth.region.dist.txt",
            "results/qc/mosdepth/{sample}.regions.bed.gz",
            "results/qc/mosdepth/{sample}.thresholds.bed.gz",
            summary="results/qc/mosdepth/{sample}.mosdepth.summary.txt",
        log:
            "logs/mosdepth_by/{sample}.log",
        params:
            thresholds="1,5,10,20,30",
            extra="--fast-mode --no-per-base --mapq 20",
        threads:
            4
        wrapper:
            "v0.85.1/bio/mosdepth"
else:
    rule mosdepth_global:
        input:
            bam=get_sample_bams,
        output:
            "results/qc/mosdepth/{sample}.mosdepth.global.dist.txt",
            "results/qc/mosdepth/{sample}.mosdepth.region.dist.txt",
            "results/qc/mosdepth/{sample}.regions.bed.gz",
            "results/qc/mosdepth/{sample}.thresholds.bed.gz",
            summary="results/qc/mosdepth/{sample}.mosdepth.summary.txt",
        log:
            "logs/mosdepth_by/{sample}.log",
        params:
            by=500,
            thresholds="1,5,10,20,30",
            extra="--fast-mode --no-per-base --mapq 20",
        threads:
            4
        wrapper:
            "v0.84.0/bio/mosdepth"


rule gatk_varianteval:
    input:
        vcf="results/genotyped/all.vcf.gz",
        ref="resources/genome.fasta",
        known=get_variation_vcf(),
    output:
        vcf="results/qc/varianteval.grp"
    log:
        "logs/gatk/varianteval.log"
    params:
        extra=get_variant_eval_extra(),
        java_opts="",
    resources:
        mem_gb=48
    wrapper:
        "v0.85.1/bio/gatk/varianteval"


rule multiqc:
    input:
        expand(
            "results/qc/samtools-stats/{u.sample}-{u.unit}.txt",
            u=units.itertuples()),
        expand(
            "results/qc/fastqc/{u.sample}-{u.unit}.zip",
            u=units.itertuples()),
        expand(
            "results/qc/dedup/{u.sample}-{u.unit}.metrics.txt",
            u=units.itertuples()),
        expand(
            "results/qc/mosdepth/{u.sample}.mosdepth.summary.txt",
            u=units.itertuples()),
        "results/annotated/all.vcf.gz",
        "results/qc/varianteval.grp",
    output:
        report(
            "results/qc/multiqc.html",
            caption="../report/multiqc.rst",
            category="Quality control",
        ),
    conda:
        "../envs/multiqc.yaml"
    log:
        "logs/multiqc.log",
    shell:
        "multiqc --force -o results/qc -n multiqc results/qc results/recal "
        "results/annotated results/stats logs/trimmomatic/ 2>&1 > {log}"
