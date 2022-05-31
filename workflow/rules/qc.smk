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
            bam="results/recal/{sample}-{unit}.bam",
            bed=config["processing"].get("restrict-regions", ""),
        output:
            "results/qc/mosdepth/{sample}-{unit}.mosdepth.global.dist.txt",
            "results/qc/mosdepth/{sample}-{unit}.mosdepth.region.dist.txt",
            "results/qc/mosdepth/{sample}-{unit}.regions.bed.gz",
            "results/qc/mosdepth/{sample}-{unit}.thresholds.bed.gz",
            summary="results/qc/mosdepth/{sample}-{unit}.mosdepth.summary.txt",
        log:
            "logs/mosdepth_by/{sample}-{unit}.log",
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
            bam="results/recal/{sample}-{unit}.bam",
        output:
            "results/qc/mosdepth/{sample}-{unit}.mosdepth.global.dist.txt",
            "results/qc/mosdepth/{sample}-{unit}.mosdepth.region.dist.txt",
            "results/qc/mosdepth/{sample}-{unit}.regions.bed.gz",
            "results/qc/mosdepth/{sample}-{unit}.thresholds.bed.gz",
            summary="results/qc/mosdepth/{sample}-{unit}.mosdepth.summary.txt",
        log:
            "logs/mosdepth_by/{sample}-{unit}.log",
        params:
            by=500,
            thresholds="1,5,10,20,30",
            extra="--fast-mode --no-per-base --mapq 20",
        threads:
            4
        wrapper:
            "v0.85.1/bio/mosdepth"


rule gatk_varianteval:
    input:
        vcf="results/annotated/all.vcf.gz",
        tbi="results/annotated/all.vcf.gz.tbi",
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
        mem_mb=8192
    wrapper:
        "v0.85.1/bio/gatk/varianteval"


if config.get('ped') and config['ref']['species'].lower() == 'homo_sapiens':
    rule peddy:
        input:
            vcf="results/genotyped/all.vcf.gz",
            ped=config['ped']
        output:
            "results/qc/peddy/all.peddy.ped"
        log:
            "logs/qc/peddy.log"
        params:
            extra="--sites hg38" if config['ref']['build'].endswith('38') else ""
        conda:
            "../envs/peddy.yaml"
        threads:
            8
        shell:
            "python -m peddy {params.extra} -p {threads} --plot " +
            "--prefix results/qc/peddy/all {input.vcf} {input.ped}"


rule multiqc:
    input:
        get_multiqc_input
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
