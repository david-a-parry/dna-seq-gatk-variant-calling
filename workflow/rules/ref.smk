from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

if (config["ref"].get("fasta_url") is None or
        config["ref"].get("variation_url") is None):
    for x in ["species", "release", "build"]:
        if x not in config["ref"]:
            raise ValueError(
                    "You must specify species, release and build for " +
                    "retrieval of ensembl genome/variation data in " +
                    "config/config.yaml or else provide values for fasta_url" +
                    "and variation_url.")

if config["ref"].get("fasta_url") is not None:
    # genome will be retrieved from URL
    HTTP = HTTPRemoteProvider()
    rule get_genome_remote:
        input:
            HTTP.remote(config["ref"]["fasta_url"])
        output:
            "resources/genome.fasta",
        params:
            genome_tmp="resources/genome_tmp.fa.gz"
        log:
            "logs/get-genome.log",
        cache: True
        run:
            if config["ref"]["fasta_url"].endswith(".gz"):
                shell("mv -v {input} {params.genome_tmp} 2> {log}")
                shell("gzip -dc {params.genome_tmp} > {output} 2>> {log}")
            else:
                shell("mv -v {input} {output} 2> {log}")
else:
    # genome will be retrieved from ensembl using species/version
    rule get_genome_ensembl:
        output:
            "resources/genome.fasta",
        log:
            "logs/get-genome.log",
        params:
            species=config["ref"]["species"],
            datatype="dna",
            build=config["ref"]["build"],
            release=config["ref"]["release"],
        cache: True
        wrapper:
            "0.74.0/bio/reference/ensembl-sequence"


checkpoint genome_faidx:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.fasta.fai",
    log:
        "logs/genome-faidx.log",
    cache: True
    wrapper:
        "0.74.0/bio/samtools/faidx"


rule genome_dict:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.dict",
    log:
        "logs/samtools/create_dict.log",
    conda:
        "../envs/samtools.yaml"
    cache: True
    shell:
        "samtools dict {input} > {output} 2> {log} "


if config["ref"].get("variation_url") is not None:
    HTTP = HTTPRemoteProvider()
    rule get_known_variation_remote:
        input:
            vcf=HTTP.remote(config["ref"]["variation_url"]),
            fai="resources/genome.fasta.fai",
        output:
            "resources/variation.vcf.gz"
        log:
            "logs/get-known-variants.log",
        cache: True
        shell:
            "mv -v {input.vcf} {output} 2> {log}"
else:
    rule get_known_variation_ensembl:
        input:
            # use fai to annotate contig lengths for GATK BQSR
            fai="resources/genome.fasta.fai",
        output:
            vcf="resources/variation.vcf.gz",
        log:
            "logs/get-known-variants.log",
        params:
            species=config["ref"]["species"],
            build=config["ref"]["build"],
            release=config["ref"]["release"],
            type="all",
        cache: True
        wrapper:
            "0.74.0/bio/reference/ensembl-variation"


rule remove_iupac_codes:
    input:
        "resources/variation.vcf.gz",
    output:
        "resources/variation.noiupac.vcf.gz",
    log:
        "logs/fix-iupac-alleles.log",
    conda:
        "../envs/rbt.yaml"
    cache: True
    shell:
        "rbt vcf-fix-iupac-alleles < {input} | bcftools view -Oz > {output}"


rule tabix_known_variants:
    input:
        "resources/variation.noiupac.vcf.gz",
    output:
        "resources/variation.noiupac.vcf.gz.tbi",
    log:
        "logs/tabix/variation.log",
    params:
        "-p vcf",
    cache: True
    wrapper:
        "0.74.0/bio/tabix"


rule bwa_index:
    input:
        "resources/genome.fasta",
    output:
        multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index.log",
    resources:
        mem_mb=369000,
    cache: True
    wrapper:
        "0.74.0/bio/bwa/index"


rule get_vep_cache:
    output:
        directory("resources/vep/cache"),
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    log:
        "logs/vep/cache.log",
    wrapper:
        "0.74.0/bio/vep/cache"


rule get_vep_plugins:
    output:
        directory("resources/vep/plugins"),
    log:
        "logs/vep/plugins.log",
    params:
        release=config["ref"]["release"],
    wrapper:
        "0.74.0/bio/vep/plugins"
