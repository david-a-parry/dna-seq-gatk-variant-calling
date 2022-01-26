if "restrict-regions" in config["processing"]:

    rule compose_regions:
        input:
            config["processing"]["restrict-regions"],
        output:
            "results/called/{contig}.regions.bed",
        log:
            "logs/bedops/extract_{contig}.log",
        shell:
            "grep -w '^{wildcards.contig}' {input} > {output} || true"


rule call_variants:
    input:
        bam=get_sample_bams,
        ref="resources/genome.fasta",
        idx="resources/genome.dict",
        known=get_variation_vcf(),
        tbi=get_variation_vcf() + ".tbi",
        regions=(
            "results/called/{contig}.regions.bed"
            if config["processing"].get("restrict-regions")
            else []
        ),
    output:
        gvcf=protected("results/called/{sample}.{contig}.g.vcf.gz"),
    log:
        "logs/gatk/haplotypecaller/{sample}.{contig}.log",
    params:
        extra=get_call_variants_params,
        java_opts="",
    resources:
        mem_mb=4096
    wrapper:
        "0.84.0/bio/gatk/haplotypecaller"


rule combine_calls:
    input:
        ref="resources/genome.fasta",
        gvcfs=expand(
            "results/called/{sample}.{{contig}}.g.vcf.gz", sample=samples.index
        ),
    output:
        gvcf="results/called/all.{contig}.g.vcf.gz",
    log:
        "logs/gatk/combinegvcfs.{contig}.log",
    wrapper:
        "0.84.0/bio/gatk/combinegvcfs"


rule genotype_variants:
    input:
        ref="resources/genome.fasta",
        gvcf="results/called/all.{contig}.g.vcf.gz",
        known=get_variation_vcf(),
    output:
        vcf=temp("results/genotyped/all.{contig}.vcf.gz"),
    params:
        extra=get_genotype_variants_params,
    resources:
        mem_mb=4096
    log:
        "logs/gatk/genotypegvcfs.{contig}.log",
    wrapper:
        "0.84.0/bio/gatk/genotypegvcfs"


rule merge_variants:
    input:
        vcfs=lambda w: expand(
            "results/genotyped/all.{contig}.vcf.gz", contig=get_contigs()
        ),
    output:
        vcf="results/genotyped/all.vcf.gz",
    log:
        "logs/picard/merge-genotyped.log",
    wrapper:
        "0.84.0/bio/picard/mergevcfs"
