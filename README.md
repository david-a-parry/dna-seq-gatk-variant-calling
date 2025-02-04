### _This is a fork of the workflow described below, modified to allow the user to specify their own reference genome and variant data._

## Additions/Improvements

* Implement VQSR correctly
* Allow specification of non-Ensembl reference genome
* Allow specification of non-Ensembl variant data
* Adjust rules so that FASTQC actually runs
* Avoid frequent segfaults in `compose_regions` rule when using a target file
* Update wrapper versions
* Add variant evaluation and mosdepth to QC rules
* Add more inputs to MultiQC
* Add SGE cluster configuration

# Snakemake workflow: dna-seq-gatk-variant-calling

[![DOI](https://zenodo.org/badge/139045164.svg)](https://zenodo.org/badge/latestdoi/139045164)
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.1.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling/workflows/Tests/badge.svg?branch=main)](https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling/actions?query=branch%3Amain+workflow%3ATests)

This Snakemake pipeline implements the [GATK best-practices workflow](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-) for calling small germline variants.

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=snakemake-workflows%2Fdna-seq-gatk-variant-calling).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and its DOI (see above).
