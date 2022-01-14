rule select_calls:
    input:
        ref="resources/genome.fasta",
        vcf="results/genotyped/all.vcf.gz",
    output:
        vcf=temp("results/filtered/all.{vartype}.vcf.gz"),
    params:
        extra=get_select_vartype_arg,
    log:
        "logs/gatk/selectvariants/{vartype}.log",
    wrapper:
        "0.84.0/bio/gatk/selectvariants"


rule hard_filter_calls:
    input:
        ref="resources/genome.fasta",
        vcf="results/filtered/all.{vartype}.vcf.gz",
    output:
        vcf=temp("results/filtered/all.{vartype}.hardfiltered.vcf.gz"),
    params:
        filters=get_filter,
    log:
        "logs/gatk/variantfiltration/{vartype}.log",
    wrapper:
        "0.84.0/bio/gatk/variantfiltration"


rule recalibrate_snvs:
    input:
        vcf="results/filtered/all.snvs.vcf.gz",
        ref="resources/genome.fasta",
        hapmap=config['filtering']['recal']['resources']['hapmap']['file'],
        omni=config['filtering']['recal']['resources']['omni']['file'],
        g1k=config['filtering']['recal']['resources']['g1k']['file'],
        dbsnp=config['filtering']['recal']['resources']['dbsnp']['file'],
    output:
        vcf="results/filtered/all.snvs.recal",
        tranches="results/filtered/all.snvs.tranches",
    params:
        extra=config["params"]["gatk"]["VariantRecalibrator"]["snvs"],
        mode="SNP",
        resources=get_vqsr_resources("snvs"),
        annotation=get_vqsr_annotations("snvs"),
        java_opts=""
    resources:
        mem_gb=48
    log:
        "logs/gatk/variantrecalibrator/snvs.log",
    wrapper:
        "0.84.0/bio/gatk/variantrecalibrator"


rule recalibrate_indels:
    input:
        vcf="results/filtered/all.indels.vcf.gz",
        ref="resources/genome.fasta",
        mills=config['filtering']['recal']['resources']['mills']['file'],
        axiom=config['filtering']['recal']['resources']['axiom']['file'],
        dbsnp=config['filtering']['recal']['resources']['dbsnp']['file'],
    output:
        vcf="results/filtered/all.indels.recal",
        tranches="results/filtered/all.indels.tranches",
    params:
        extra=config["params"]["gatk"]["VariantRecalibrator"]["indels"],
        mode="INDEL",
        resources=get_vqsr_resources("indels"),
        annotation=get_vqsr_annotations("indels"),
        java_opts=""
    resources:
        mem_gb=48
    log:
        "logs/gatk/variantrecalibrator/indels.log",
    wrapper:
        "0.84.0/bio/gatk/variantrecalibrator"


rule apply_vqsr:
    input:
        vcf="results/filtered/all.{vartype}.vcf.gz",
        recal="results/filtered/all.{vartype}.recal",
        tranches="results/filtered/all.{vartype}.tranches",
        ref="resources/genome.fasta",
    output:
        vcf=temp("results/filtered/all.{vartype}.recalibrated.vcf.gz"),
    params:
        mode=get_recal_mode,
        extra=get_vqsr_sensitivity,
    log:
        "logs/gatk/variantrecalibrator/apply_recal_{vartype}.log",
    wrapper:
        "0.84.0/bio/gatk/applyvqsr"


rule merge_calls:
    input:
        vcfs=expand(
            "results/filtered/all.{vartype}.{filtertype}.vcf.gz",
            vartype=["snvs", "indels"],
            filtertype="recalibrated"
            if config["filtering"]["vqsr"]
            else "hardfiltered",
        ),
    output:
        vcf="results/filtered/all.vcf.gz",
    log:
        "logs/picard/merge-filtered.log",
    wrapper:
        "0.84.0/bio/picard/mergevcfs"
