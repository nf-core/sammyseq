# nf-core/sammyseq: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [FastQC](#fastqc)
- [Trim adapters](#trim-adapters)
- [Align to Reference](#align-to-reference)
- [Mark Duplicates](#mark-duplicates)
- [`Signal track generation`](#signal-track-generation)
- [`Comparisons`](#comparisons)
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about the sequenced reads. It provides information about the quality score distribution across reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

#### Trim adapters

[`Trimmomatic`](http://www.usadellab.org/cms/?page=trimmomatic) is a software for read trimming, which means it trims all reads in the front or the tail. This function is useful since sometimes you want to drop some cycles of a sequencing run. In the current implementation in Sarek
`--detect_adapter_for_pe` is set by default which enables auto-detection of adapter sequences. For more information on how to fine-tune adapter trimming, take a look into the parameter docs.

The resulting files are intermediate and by default not kept in the final files delivered to users. Set `--save_trimmed` to enable publishing of the files in:

<details markdown="1">
<summary>Output files for all samples</summary>

**Output directory: `{outdir}/preprocessing/fastp/<sample>`**

- `<sample>_<lane>_{1,2}.fastp.fastq.gz>`
  - Bgzipped FastQ file

</details>

:::note
The FastQC plots displayed in the MultiQC report shows both _untrimmed_ and _trimmed_ reads. They can be compared to check removal of adapter sequence and potentially regions with low quality after the trimming step.
:::

### Align to Reference

[BWA](https://github.com/lh3/bwa) is a software package for mapping low-divergent sequences against a large reference genome. The aligned reads are then sorted with [samtools](https://www.htslib.org/doc/samtools.html).

<details markdown="1">
<summary>Output files</summary>

- `alignment/bwa/`
  - `<sample>.bam` and `<sample>.bam.bai`

</details>

### Mark Duplicates

<details markdown="1">
<summary>Output files</summary>

- `alignment/markduplicates/`
  - `<sample>.md.bam` and `<sample>.md.bam.bai`
- `reports/markduplicates/`
  - `<sample>.md.MarkDuplicates.metrics.txt`

</details>

Read pairs that are likely to have originated from duplicates of the same original DNA fragments through some artificial processes are identified. These are considered to be non-independent observations, so all but a single read pair within each set of duplicates are marked, not removed from the bam file.

### Signal track generation

<details markdown="1">
<summary>Output files</summary>

- `single_tracks/deeptools/`
  - `<sample>.bigWig`

</details>

[deepTools](https://deeptools.readthedocs.io/en/develop/content/list_of_tools.html) is used to generate single fraction signals in [bigWig](https://genome.ucsc.edu/goldenpath/help/bigWig.html) format, an indexed binary format useful for displaying dense, continuous data in Genome Browsers such as the [UCSC](https://genome.ucsc.edu/cgi-bin/hgTracks) and [IGV](http://software.broadinstitute.org/software/igv/). The bigWig format is also supported by various bioinformatics software for downstream processing such as meta-profile plotting.

### Comparisons

When `--comparisonFile` is set, the difference between sample1 and sample2 read density profile, smoothed by the Gaussian kernel, is calculated and saved in bigwig format, as described in Kharchenko PK, Tolstorukov MY, Park PJ "Design and analysis of ChIP-seq experiments for DNA-binding proteins" Nat. Biotech. doi:10.1038/nbt.1508

<details markdown="1">
<summary>Output files</summary>

- `comparisons/spp_mle/`
  - `<sample1>.md_VS_<sample2>.md.bw`

</details>

### MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

### Reference genome files

<details markdown="1">
<summary>Output files</summary>

- `genome/`

  - A number of genome-specific files are generated by the pipeline in order to aid in the filtering of the data, and because they are required by standard tools such as BEDTools. These can be found in this directory.
  - `bwa/`: Directory containing BWA indices.

  - If the `--save_reference` parameter is provided then the alignment indices generated by the pipeline will be saved in this directory. This can be quite a time-consuming process so it permits their reuse for future runs of the pipeline or for other purposes.

</details>

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
