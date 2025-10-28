# nf-core/sammyseq: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/sammyseq/usage](https://nf-co.re/sammyseq/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

sammyseq is a workflow designed for the analysis of Sequential Analysis of MacroMolecules accessibilitY sequencing (SAMMY-seq) data, a cheap and effective methodology to analyze chromatin state in cells. SAMMY-seq is an innovative technique based on the separation of chromatin in fractions, each progressively based on their solubility and accessibility, and extraction and sequencing of the DNA present in each of them.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/sammyseq -r dev \
    -profile docker \
    --fasta ./reference_genome.fa\
    --input ./samplesheet.csv \
    --outdir ./results
```

This will launch the pipeline with the `docker` configuration profile. See [below](#profile) for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/sammyseq -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
fasta: './genome.fa'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

## Samplesheet input

Before running the pipeline, you will need to create a samplesheet with information about the samples you would like to analyze. Use this parameter to specify its location:

```bash
--input '[full path to samplesheet file]'
```

It has to be a comma-separated file with 5 columns, and a header row as shown in the examples below.

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample fraction sequenced across 2 lanes:

```console
sample,fastq_1,fastq_2,experimentalID,fraction,sample_group
CONTROL_REP1_S2,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,CONTROL_REP1,S2,CONTROL
CONTROL_REP1_S2,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz,CONTROL_REP1,S2,CONTROL
```

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can contain a mixture of single- and paired-end but in case of multiple runs of the same `sample` they have to be of the same type to be correctly merged. There can be additional columns but the first 5 have to match those defined in the table below.

```console
sample,fastq_1,fastq_2,experimentalID,fraction,sample_group
CTRL004_S2,/home/sammy/test_data/CTRL004_S2_chr22only.fq.gz,,CTRL004,S2,CTRL
CTRL004_S3,/home/sammy/test_data/CTRL004_S3_chr22only.fq.gz,,CTRL004,S3,CTRL
CTRL004_S4,/home/sammy/test_data/CTRL004_S4_chr22only.fq.gz,,CTRL004,S4,CTRL
```

| Column           | Description                                                                                                                                                                            |
| ---------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`         | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1`        | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2`        | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `experimentalID` | Experimental sample identifier. This represents the biological specimen of interest and will be the same for all fractions exctracted.                                                 |
| `fraction`       | Fraction derived from SAMMY protocol, e.g. depending on the protocol it can be S2, S2L, S2S, S3, S4.                                                                                   |
| `sample_group`   | Identifier used to group samples that belong to the same biological condition condition.                                                                                               |
|                  |

### Pairwise comparisons

The pipeline offers two different methods for generating these comparisons, selected with the `--comparison_maker` parameter.

The `spp` method (default) smooths fraction read density profiles using a Gaussian kernel, calculates differences between fractions, and outputs results in bigwig format, following the approach described in Kharchenko PK, Tolstorukov MY, Park PJ "Design and analysis of ChIP-seq experiments for DNA-binding proteins" Nat Biotech [doi](https://doi.org/10.1038/nbt.1508).

```bash
--comparison_maker spp
```

The `bigwigcompare` method partitions the genome into bins of equal size defined by the `--bw_resolution` parameter (defaults to 1bp), counts reads per bin, and calculates the log2 ratio between samples (other operations can be selected by changing the `--bigwigcompare_operation` parameter, please see [the software documentation](https://deeptools.readthedocs.io/en/latest/content/tools/bigwigCompare.html)). This method is required for paired-end data as `spp` does not support this data type.

```
--comparison_maker bigwigcompare
```

It is possible to generate one or more pairwise comparisons between fractions from the same experimental replicate by providing the `--comparison` parameter. You can specify a single comparison or multiple comparisons separated by commas:

**Single comparison:**

```bash
--comparison S2SvsS3
```

**Multiple comparisons:**

```
--comparison S2SvsS3,S2SvsS4,S4vsS3
```

For 4f-SAMMYseq protocols (S2S, S2L, S3, S4), valid comparisons are:

    S2SvsS3  - Compare S2S fraction vs S3 fraction
    S2LvsS3  - Compare S2L fraction vs S3 fraction
    S2SvsS4  - Compare S2S fraction vs S4 fraction
    S2LvsS4  - Compare S2L fraction vs S4 fraction
    S4vsS3   - Compare S4 fraction vs S3 fraction

> [!NOTE]
> For 3f-SAMMYseq protocols (S2, S3, S4), valid comparisons are:

    S2vsS3   - Compare S2 fraction vs S3 fraction
    S2vsS4   - Compare S2 fraction vs S4 fraction
    S4vsS3   - Compare S4 fraction vs S3 fraction

The pipeline will automatically create comparisons only between fractions from the same `experimentalID` (biological replicate), ensuring that comparisons are made within the same experimental condition rather than across different replicates.

Alternatively, it is possible to generate any pairwise comparisons between any fraction by providing a list with the parameter `--comparison_file` to indicate the full path to a comma-separated file with 2 columns:

`comparisons.csv`:

```csv
sample1,sample2
CTRL004_S2,CTRL004_S3
CTRL004_S2,CTRL004_S4
```

It can contain any combination of sample identifiers, they have to correspond to identifiers present in the `sample` column in the input file.

### Differential Solubility Analysis

The Differential Solubility Analysis has been developed to investigate differences in solubility patterns between experimental conditions and is enabled by setting the `--differential_solubility` parameter. This analysis uses the comparison tracks requested by the `--comparison` parameter (e.g., S2SvsS3, at least one comparison _has_ to be selected) to identify genomic regions with significantly different accessibility patterns:

The `--compare_groups` parameter specifies which sample groups (defined in the `sample_group` column of the samplesheet) to compare for differential analysis. In the format "GroupBvsGroupA", GroupA serves as the reference group against which GroupB is compared.

**Single group comparison:**

```
--compare_groups "GroupBvsGroupA"
```

**Multiple group comparisons:**

```
--compare_groups "GroupBvsGroupA,GroupCvsGroupA"
```

The `--gtf` parameter (optional): When provided, generates gene lists by identifying protein-coding genes whose promoter regions overlap with significantly different bins

### Combine fractions

Optionally, the fractions extracted from the same `experimentalID` can be combined together for later use by setting the parameter `--combine_fractions`.

## Reference files

### Genome

The minimum reference genome requirements is the FASTA file, provided with the mandatory parameter `--fasta`, the aligner index will be generated by the pipeline and can be saved for later reuse if the `--save_reference` parameter is passed. The index building step can be quite a time-consuming process and it permits their reuse for future runs of the pipeline to save disk space, if already present it can be passed using the `--bwa_index '/path/to/bwa/index/'` or `--bowtie2_index '/path/to/bowtie2/index/'` parameter, depending on the chosen algorithm. Also the `--fai` fasta index and the `--chrom_sizes` chromosome sizes file can be passed if available, otherwise will be generated. The genome coordinates is binned into windows of the size defined by the `--binsize` parameter for downstream analysis.

### Blacklist bed file

A blacklist of regions that will be excluded by signal tracks can be provided using the optional parameter `--blacklist` with full path to a coordinate file in bed format. Blacklist files for several genome builds can be found in the [ENCODE Blacklist Project](https://github.com/Boyle-Lab/Blacklist).

### Keep regions bed file

A list of regions that will be kept in the output after filtering the alignment with samtools using -L option (in addition to the flag and quality threshold filters) can be passed with the optional parameter `--keep_regions_bed` with full path to a coordinate file in bed format.

### TSS bed file

Path to BED file containing TSS regions provided using the optional parameter `--tss_bed`, it will be used to a file for each sample with fraction signal profiles across the TSS coordinates.

### GTF file

Path to GTF file containing genes coordinates provided using the optional parameter `--gtf`. When used with `--differential_solubility`, it enables gene-level annotation of results.

## Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/sammyseq
```

## Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/sammyseq releases page](https://github.com/nf-core/sammyseq/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
