import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.18.0")


report: "../report/workflow.rst"


container: "continuumio/miniconda3:4.8.2"


###### Config file and sample sheets #####
configfile: "config/config.yaml"


validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(
    ["sample", "unit"], drop=False
)
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels]
)  # enforce str in index
validate(units, schema="../schemas/units.schema.yaml")


##### Wildcard constraints #####
wildcard_constraints:
    vartype="snvs|indels",
    sample="|".join(samples.index),
    unit="|".join(units["unit"]),


##### Helper functions #####

# contigs in reference genome/target regions
def get_contigs():
    with checkpoints.genome_faidx.get().output[0].open() as fai:
        if "restrict-regions" in config["processing"]:
            bed = pd.read_csv(config["processing"]["restrict-regions"],
                              header=None, sep='\t', usecols=[0], dtype=str,
                              squeeze=True)
            return bed.unique()
        return pd.read_csv(fai, header=None, usecols=[0], squeeze=True,
                           sep='\t', dtype=str)


def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}


def is_single_end(sample, unit):
    """Return True if sample-unit is single end."""
    return pd.isnull(units.loc[(sample, unit), "fq2"])


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:{platform}'".format(
        sample=wildcards.sample,
        platform=units.loc[(wildcards.sample, wildcards.unit), "platform"],
    )


def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand(
            "results/trimmed/{sample}-{unit}.{group}.fastq.gz",
            group=[1, 2],
            **wildcards
        )
    # single end sample
    return "results/trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)


def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand(
        "results/recal/{sample}-{unit}.bam",
        sample=wildcards.sample,
        unit=units.loc[wildcards.sample].unit,
    )


def get_regions_param(regions=config["processing"].get("restrict-regions"), default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        padding = config["processing"].get("region-padding")
        if padding:
            params += "--interval-padding {}".format(padding)
        return params
    return default


def get_call_variants_params(wildcards, input):
    return (
        get_regions_param(
            regions=input.regions, default="--intervals {}".format(wildcards.contig)
        )
        + config["params"]["gatk"]["HaplotypeCaller"]
    )


def get_genotype_variants_params(wildcards):
    args = config["params"]["gatk"]["GenotypeGVCFs"]
    if config.get('ped'):
        args += ' --pedigree {}'.format(config['ped'])
    return args


def get_recal_input(bai=False):
    # case 1: no duplicate removal
    f = "results/mapped/{sample}-{unit}.sorted.bam"
    if config["processing"]["remove-duplicates"]:
        # case 2: remove duplicates
        f = "results/dedup/{sample}-{unit}.bam"
    if bai:
        if config["processing"].get("restrict-regions"):
            # case 3: need an index because random access is required
            f += ".bai"
            return f
        else:
            # case 4: no index needed
            return []
    else:
        return f


def get_snpeff_reference():
    return "{}.{}".format(config["ref"]["build"], config["ref"]["snpeff_release"])


def get_select_vartype_arg(wildcards):
    if wildcards.vartype == "snvs":
        return "--select-type-to-include SNP " + "--select-type-to-include MNP"
    return "--select-type-to-include INDEL " + "--select-type-to-include MIXED"


def get_recal_mode(wildcards):
    return "SNP" if wildcards.vartype == "snvs" else "INDEL"


def get_filter(wildcards):
    return {"snv-hard-filter": config["filtering"]["hard"][wildcards.vartype]}


def get_variation_vcf():
    if config["ref"].get('fix_non_iupac'):
        return "resources/variation.noiupac.vcf.gz"
    return "resources/variation.vcf.gz"


def get_vqsr_annotations(vartype):
    return config["filtering"]["recal"]["annotations"][vartype]


def get_vqsr_resources(vartype):
    resource_defaults = {
      "hapmap": {"known": False, "training": True, "truth": True, "prior": 15.0},
      "omni": {"known": False, "training": True, "truth": False, "prior": 12.0},
      "g1k": {"known": False, "training": True, "truth": False, "prior": 10.0},
      "mills": {"known": False, "training": True, "truth": True, "prior": 12.0},
      "axiom": {"known": False, "training": True, "truth": False, "prior": 10.0},
      "dbsnp": {"known": True, "training": False, "truth": False, "prior": 2.0}
    }
    if vartype == "snvs":
        datasets = [x for x in ["hapmap", "omni", "g1k", "dbsnp"] if x in
                    config['filtering']['recal']['resources']]
    else:
        datasets = [x for x in ["mills", "axiom", "dbsnp"] if x in
                    config['filtering']['recal']['resources']]
    res = dict()
    for x in datasets:
        res[x] = resource_defaults[x]
        if x in config['filtering']['recal']['resources']:
            res[x]['prior'] = config['filtering']['recal']['resources'].get(
                   'prior',
                   resource_defaults[x]['prior'])
    return res


def get_vqsr_sensitivity(wildcards):
    return "--truth-sensitivity-filter-level {}".format(
          config['filtering']['recal']['truth_sensitivity'][wildcards.vartype])


def get_variant_eval_extra():
    args = ''
    if config['processing'].get('restrict-regions'):
        args += ' --intervals {}'.format(
                config['processing']['restrict-regions'])
    if config.get('ped'):
        args += ' --pedigree {}'.format(config['ped'])
        args = ' -EV MendelianViolationEvaluator'
    return args


def get_multiqc_input(wildcards):
    inputs = expand(
        "results/qc/samtools-stats/{u.sample}-{u.unit}.txt",
        u=units.itertuples())
    inputs.extend(
            expand("results/qc/fastqc/{u.sample}-{u.unit}.zip",
                u=units.itertuples()))
    inputs.extend(
            expand("results/qc/dedup/{u.sample}-{u.unit}.metrics.txt",
                u=units.itertuples()))
    inputs.extend(
            expand(
                "results/qc/mosdepth/{u.sample}.mosdepth.summary.txt",
                u=units.itertuples()))
    inputs.append("results/annotated/all.vcf.gz")
    inputs.append("results/qc/varianteval.grp")
    if config.get('ped') and config['ref']['species'].lower() == 'homo_sapiens':
       inputs.append("results/qc/peddy/all.peddy.ped")
    return inputs
