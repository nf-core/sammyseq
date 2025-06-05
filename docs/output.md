# nf-core/sammyseq: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [FastQC](#fastqc)
- [Trim reads](#trim-reads)
- [Alignment on Reference](#alignment-on-reference)
- [Mark Duplicate reads](#mark-duplicate-reads)
- [Samtools reads filtering](#samtools-reads-filtering)
- [Signal track generation](#signal-track-generation)
- [DeepTools based QC](#deeptools-based-qc)
- [Comparisons](#comparisons)
- [MultiQC](#multiqc)
- [Pipeline information](#pipeline-information)

### Read quality check

#### FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about the sequenced reads. It provides information about the quality score distribution across reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

#### Trim reads

The task of trim adapter sequences and low quality bases from the end can be performed using either [`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore) or [`Trimmomatic`](http://www.usadellab.org/cms/?page=trimmomatic) and quality check after this step is performed again with Fastqc.

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.
  - `*_trim_fastqc.html`: FastQC report containing quality metrics for trimmed reads.
  - `*_trim_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images for trimmed reads.

</details>

:::note
The FastQC plots displayed in the MultiQC report shows both _untrimmed_ and _trimmed_ reads so they can be directly compared.
:::

### Alignment on Reference

The alignment will be performed using [BWA-MEM](https://github.com/lh3/bwa) as default algorithm, alternatives BWA-ALN or [Bowtie2](https://github.com/BenLangmead/bowtie2) can be chosen with the parameter `--aligner`.
The aligned reads are then sorted by chromosome coordinates with [samtools](https://www.htslib.org/doc/samtools.html).

<details markdown="1">
<summary>Parameters options</summary>

- `--aligner bwamem` (default)
- `--aligner bwaaln`
- `--aligner bowtie2`
</details>

### Mark Duplicate reads

Read pairs that are likely to have originated from duplicates of the same original DNA fragments through some artificial processes are identified. These are considered to be non-independent observations, so all but a single read pair within each set of duplicates are marked, not removed from the bam file.

<details markdown="1">
<summary>Output files</summary>

- `alignment/markduplicates/`
  - `<sample>.md.bam` and `<sample>.md.bam.bai`
- `reports/markduplicates/`
  - `<sample>.md.MarkDuplicates.metrics.txt`

</details>

### Samtools reads filtering

The BAM files generated are further processed with SAMtools for filtering (based on samtools flags and quality score) and indexing, as well as to generate read mapping statistics.

<details markdown="1">
<summary>Output files</summary>

- `alignment/filtered/`
  - `<sample>.<q_score>.bam` and `<sample>.<q_score>.bam.bai`
- `/reports/samtools_stats/<sample>/filtered/`
  - `<sample>/filtered.idxstats`
  - `<sample>/filtered.flagstat`
  - `<sample>/filtered.stats`
  </details>

### Signal track generation

[deepTools](https://deeptools.readthedocs.io/en/develop/content/list_of_tools.html) is used to generate single fraction signals in [bigWig](https://genome.ucsc.edu/goldenpath/help/bigWig.html) format, an indexed binary format useful for displaying dense, continuous data in Genome Browsers such as the [UCSC](https://genome.ucsc.edu/cgi-bin/hgTracks) and [IGV](http://software.broadinstitute.org/software/igv/). The bigWig format is also supported by various bioinformatics software for downstream processing such as meta-profile plotting.
The generated signal tracks represent read coverage and can be normalized using different methods: RPKM (default option),CPM, BPM and RPGC.

<details markdown="1">
<summary>Output files</summary>

- `single_tracks/deeptools/`
  - `<sample>.<q_score>.<normalizeUsing>.bw`

</details>

### DeepTools based QC

DeepTools is used to perform quality control analysis at the aligned fraction level. The pipeline uses several DeepTools commands to generate comprehensive QC metrics and visualizations.

### MultiBAMSummary

The process starts with multiBamSummary, which computes the read coverage over the entire genome (or a specified region) for multiple BAM files. This creates a matrix of read counts that serves as input for the subsequent analyses. By default, the bin size is set to 50000, but you can adjust this using the --bam_binsize parameter.

<details markdown="1"> <summary>Output files</summary>

    deeptools/quality_control/multibamsummary/
        ${meta.id}.npz: Binary file containing the read coverage matrix
        outRawCounts.txt: Text file with raw read counts

</details>

### PCA (Principal Component Analysis)

PCA is used to analyze and visualize variability in high-dimensional datasets. In the context of sequencing data analysis, PCA helps to determine if samples show greater variability between experimental conditions than between replicates of the same treatment.

<details markdown="1"> <summary>Output files</summary>

    deeptools/quality_control/plotpca/
        ${meta.id}.pdf: PCA plot
        ${meta.id}.tab: Table with PCA coordinates

</details>

### Correlation Heatmap

The correlation analysis computes the overall similarity between samples based on read coverage. The result is visualized as a heatmap of correlation coefficients, indicating the strength of the relationship between samples. You can specify the correlation method (e.g., 'spearman', 'pearson') if the parameter qc_corr_method is provided (default is pearson)

<details markdown="1"> <summary>Output files</summary>

    deeptools/quality_control/plotcorrelation/
        ${meta.id}.pdf: Correlation heatmap
        ${meta.id}.tab: Table with correlation coefficients

</details>

### Fingerprint Plot

This fingerprint plot is particularly useful for assessing the strength of the experiment for factors with enrichment in well-defined and relatively narrow regions.

Two types of fingerprint plots are generated:

    Global Fingerprint Plot: Covers the entire genome

<details markdown="1"> <summary>Output files</summary>

    deeptools/quality_control/plotfingerprint/global/
        ${meta.id}_global.pdf: Global fingerprint plot
        ${meta.id}_global.raw.txt: Raw data for the global fingerprint plot

</details>

    Region-specific Fingerprint Plot: Focuses on a user-specified genomic region (if --region parameter is provided (e.g., 'chr1', 'chr2:1000000-2000000'))

<details markdown="1"> <summary>Output files</summary>

    deeptools/quality_control/plotfingerprint/${params.region}/
        ${meta.id}_region_${params.region}.pdf: Region-specific fingerprint plot
        ${meta.id}_region_${params.region}.raw.txt: Raw data for the region-specific fingerprint plot

</details>

### Comparisons

When `--comparisonFile` is set, the difference between sample1 and sample2 read density profile smoothed by the Gaussian kernel is calculated and saved in bigwig format, as described in Kharchenko PK, Tolstorukov MY, Park PJ "Design and analysis of ChIP-seq experiments for DNA-binding proteins" Nat. Biotech. doi:10.1038/nbt.1508

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

A number of genome-specific files if required by some of the analysis steps. If the `--save_reference` parameter is provided then the alignment indices generated by the pipeline will be saved in this directory.

<details markdown="1">
<summary>Output files</summary>

- `genome/`
  - `bwa/`: Directory containing BWA indices.

</details>

### Pipeline information

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>
