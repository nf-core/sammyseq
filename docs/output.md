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
- [Differential Solubility Analysis](#differential-solubility-analysis)
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

### MultiBigwigSummary

The process starts with multiBigwigSummary, which computes the read coverage for multiple BigWig files. This creates a matrix of read counts that serves as input for correlation and PCA analyses.

#### PCA (Principal Component Analysis)

PCA is used to analyze and visualize variability in high-dimensional datasets. In the context of sequencing data analysis, PCA helps to determine if samples show greater variability between experimental conditions than between replicates of the same treatment.

<details markdown="1">
<summary>Output files</summary>

- `deeptools/quality_control/plotpca/`
  - `${meta.id}.pdf: PCA plot`
  - `${meta.id}.tab: Table with PCA coordinates`

</details>

#### Correlation Heatmap

The correlation analysis computes the overall similarity between samples based on read coverage. The result is visualized as a heatmap of correlation coefficients, indicating the strength of the relationship between samples. You can specify the correlation method (e.g., 'spearman', 'pearson') if the parameter qc_corr_method is provided (default is spearman)

<details markdown="1"> <summary>Output files</summary>

- `deeptools/quality_control/plotcorrelation/`
  - `${meta.id}.pdf: Correlation heatmap`
  - `${meta.id}.tab: Table with correlation coefficients`

</details>

#### Fingerprint Plot (Optional)

This fingerprint plot is particularly useful for assessing the strength of the experiment for factors with enrichment in well-defined and relatively narrow regions. This analysis uses BAM files and can be disabled using --plotfingerprint false.

#### Global Fingerprint Plot (Optional)

Covers the entire genome to assess overall signal enrichment patterns.

<details markdown="1"> <summary>Output files</summary>

- `deeptools/quality_control/plotfingerprint/global/`
  - `${meta.id}_global.pdf: Global fingerprint plot`
  - `${meta.id}_global.raw.txt: Raw data for the global fingerprint plot`

</details>

#### Region-specific Fingerprint Plot (Optional)

Focuses on a user-specified genomic region (if --region parameter is provided, e.g., 'chr1', 'chr2:1000000-2000000').

<details markdown="1"> <summary>Output files</summary>

- `deeptools/quality_control/plotfingerprint/${params.region}/`
  - `${meta.id}_region_${params.region}.pdf: Region-specific fingerprint plot`
  - `${meta.id}_region_${params.region}.raw.txt: Raw data for the region-specific fingerprint plot`

</details>

### TSS Profile Plot (Optional)

This analysis uses computeMatrix in reference-point mode to generate coverage profiles centered around transcription start sites (TSS). By supplying a BED file with genomic TSS coordinates through the `--tss_bed` parameter, the pipeline produces a line plot for each fraction of every sample, showing signal distribution around the TSS. This helps assess patterns of accessibility or enrichment near promoters.

<details markdown="1"> <summary>Output files</summary>

- `deeptools/quality_control/plotprofile/`
  - `${meta.id}_tss.pdf: TSS centered coverage profile`
  - `${meta.id}_tss.tab: Matrix of values used for the plot`

</details>

### Comparisons

Pairwise comparisons can be generated by setting the parameter `--comparison_maker spp` (default) either by `--comparison` or `--comparison_file`. The difference between each fraction read density profile, smoothed by the Gaussian kernel is calculated and saved in bigwig format, as described in Kharchenko PK, Tolstorukov MY, Park PJ "Design and analysis of ChIP-seq experiments for DNA-binding proteins" Nat. Biotech. doi:10.1038/nbt.1508

Alternatively, `--comparison_maker bigwigcompare` uses [deepTools](https://deeptools.readthedocs.io/en/develop/content/tools/bigwigCompare.html) and must be used for paired-end data as SPP does not handle them. This tool compares two bigWig files based on the number of mapped reads: the genome is partitioned into bins of equal size defined by the `--bw_resolution` parameter (defaults to 1bp), then the number of reads found in each bigWig file are counted per bin and the log2 ratio is calculated and reported (other operations can be selected by changing the `--bigwigcompare_operation` parameter).

<details markdown="1">
<summary>Output files</summary>

- `comparisons/spp_mle/`
  - `<sample1>.md_VS_<sample2>.md.bw`

- `comparisons/deeptools_bigwigcompare_bigwig/`
  - `<sample1>_VS_<sample2>.q<q_score>.<normalizeUsing>.bigWig`

</details>

### Differential Solubility Analysis

When the `--differential_solubility` parameter is enabled, the pipeline performs a differential solubility analysis between two sample groups. For each sample associated to the specified groups, the comparison tracks requested with the `--comparison` parameter (e.g., S2SvsS3, at least one comparison _has_ to be selected) are used to identify genomic regions with significantly different solubility patterns between the groups. So the relevant parameters are:

- `--compare_groups`: specifies which sample groups to compare (format: "GroupAvsGroupB" or "GroupAvsGroupB,GroupAvsGroupC") and they must be present in the `sample_group` column of the input samplesheet
- `--comparison`: the ratios used for analysis (e.g., "S2SvsS3")
- `--keep_regions_bed`: genomic regions to include in the analysis

The genome is partitioned into bins of equal size (default: `--binsize 50000`) within the regions specified by `--keep_regions_bed`. For each bin, the pipeline calculates the average value per sample, normalizes them using quantile algorithm, filters out the values below a configurable threshold (`--solubility_threshold`, default: 0.1), and performs statistical comparisons between the groups calculating Cohen's d to identifiy bins with significantly different solubility patterns.

When a GTF file is provided using `--gtf`, the pipeline processes protein-coding gene annotations and performs coordinate-based overlap between significantly different bins and gene promoter regions (2,500 bp upstream to 500 bp downstream of transcription start sites), generating lists of genes whose regulatory elements intersect with differential accessibility patterns.

<details markdown="1">
<summary>Output files</summary>

- `differential_solubility/csv/`
  - `*_all_bins_complete.csv`: Complete analysis results for all genomic bins
  - `*_selected_bins_filtered.csv`: Filtered results containing only significantly different bins

- `differential_solubility/rdata/`
  - `*_analysis.rds`: Complete R analysis objects for further downstream analysis

- `differential_solubility/regions/`
  - `*_regions.bed`: BED files containing coordinates of significantly different regions by category (up/down regulation per fraction)

- `differential_solubility/genes/`
  - `*_genes.txt`: Gene lists overlapping with significantly different regions (when GTF annotation is provided)

- `differential_solubility/reports/`
  - `*_analysis_summary.txt`: Summary report of the differential analysis results

</details>

### MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory. This step can be skipped using the parameter `--skip_multiqc`

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
