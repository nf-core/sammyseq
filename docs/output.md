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
- [Samtools bam filtering](#samtools-bam-filtering)
- [Signal track generation](#signal-track-generation)
- [DeepTools based QC](#deeptools-based-qc)
- [Compartments Analysis](#compartments-analysis)
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

Read duplicate marking is carried out on aligned BAM using the [Picard](https://github.com/broadinstitute/picard) MarkDuplicates command. Read pairs that are likely to have originated from duplicates of the same original DNA fragments through some artificial processes are identified. These are considered to be non-independent observations, so all but a single read pair within each set of duplicates are marked, not removed from the BAM file.

<details markdown="1">
<summary>Output files</summary>

- `alignment/markduplicates/`
  - `<sample>.md.bam` and `<sample>.md.bam.bai`
- `reports/markduplicates/`
  - `<sample>.md.MarkDuplicates.metrics.txt`

</details>

### Samtools bam filtering

BAM files generated after alignment and duplicate marking are further processed with [samtools](https://www.htslib.org/doc/samtools.html) to apply filtering based on mapping quality (default q_score > 1) and SAM flags (default 1540). This step also includes indexing the filtered BAM files and generating various alignment statistics, such as read counts per chromosome, overall alignment rate, and flag summaries.

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

#### Correlation Heatmap

The [deepTools](https://deeptools.readthedocs.io/en/develop/content/list_of_tools.html) plotCorrelation command is used to compute the overall similarity between samples based on genome-wide read coverage. The result is visualized as a heatmap of correlation coefficients, indicating the strength of the relationship between samples. You can specify the correlation method (e.g., 'spearman', 'pearson') by setting the `--qc_corr_method` parameter (default is 'pearson').

<details markdown="1"> <summary>Output files</summary>

    reports/deeptools/plotcorrelation/
        ${meta.id}.pdf: Correlation heatmap
        ${meta.id}.tab: Table with correlation coefficients

</details>

#### Fingerprint Plot

The [deepTools](https://deeptools.readthedocs.io/en/develop/content/list_of_tools.html) plotFingerprint command is useful for assessing the strength of the experiment for factors with enrichment in well-defined and relatively narrow regions.

Two types of fingerprint plots are generated:

Global Fingerprint Plot: Covers the entire genome

<details markdown="1"> <summary>Output files</summary>

    reports/deeptools/plotfingerprint/global/
        ${meta.id}_global.pdf: Global fingerprint plot
        ${meta.id}_global.raw.txt: Raw data for the global fingerprint plot

</details>

Region-specific Fingerprint Plot: Focuses on a user-specified genomic region (if `--region` parameter is provided (e.g., 'chr1', 'chr2:1000000-2000000'))

<details markdown="1"> <summary>Output files</summary>

    reports/deeptools/plotfingerprint/${params.region}/
        ${meta.id}_region_${params.region}.pdf: Region-specific fingerprint plot
        ${meta.id}_region_${params.region}.raw.txt: Raw data for the region-specific fingerprint plot

</details>

#### PCA (Principal Component Analysis)

The [deepTools](https://deeptools.readthedocs.io/en/develop/content/list_of_tools.html) plotPCA command is used to determine whether samples vary more between experimental conditions than between replicates. The output PDF includes both the PCA plot and the corresponding scree plot, which displays the proportion of variance explained by each principal component.

<details markdown="1"> <summary>Output files</summary>

    deeptools/quality_control/plotpca/
        ${meta.id}.pdf: PCA plot (including scree plot)
        ${meta.id}.tab: Table with PCA coordinates

</details>

#### Plot Profile

If the `--tss_bed` parameter is provided, the [deepTools](https://deeptools.readthedocs.io/en/develop/content/list_of_tools.html) plotProfile command will generate TSS-centered signal profile plots, which help visualize the average distribution of sequencing signal (e.g. coverage or enrichment) around transcription start sites (TSS). All fractions belonging to the same sample are grouped and their signal tracks aggregated to produce a single profile per sample.

<details markdown="1"> <summary>Output files</summary>

    reports/deeptools/plotprofile/{$params.tss_bed}/
        ${meta.id}.${params.tss_bed}.plotProfile.pdf: Line plot showing the average signal across TSS for all fractions of a given sample.
        ${meta.id}.${params.tss_bed}.plotProfile.tab: Tabular file with the raw values used in the plot.

</details>

### Compartments Analysis

When `--compartmentsAnalysis` is enabled, a module is triggered to infer A/B chromatin compartments from SAMMY-seq signal tracks. The analysis is based on fixed size genomic binning, performed per chromosome using bedtools makewindows, which divides the genome into windows of a defined size (default 50000) set by the `--binsize` parameter. To restrict the analysis to specific chromosomes, a BED file with only the chromosomes to include can be provided via the `--keep_regions_bed` parameter. The `--gtf` parameter is also required, as gene annotations are used in downstream steps.
Each fraction is identified using the `experimentalID` column in the samplesheet allowing fractions of the same sample to be analyzed together. The Compartments calling is based on [CALDER2](https://github.com/CSOgroup/CALDER2) algorithm, which builds a correlation matrix across genomic bins and applies eigenvector decomposition to classify each bin as either A (open) or B (closed) compartment. After compartment calling, results from all analyzed chromosomes are combined into a single compartment BED and a single BedGraph file with eigenvalues for each sample.

Consensus profiles are then computed by merging compartment calls from all replicates within the same sample group. For each genomic bin—defined according to the resolution set by the `--binsize` parameter, the consensus label is assigned based on the majority vote among replicates. A bin is labeled as A or B if it receives more compartment calls for that label than for the other; otherwise, it is marked as NA. The input for this step is the set of combined compartment BED files, and the consensus output is sorted and formatted for IGV visualization.

<details markdown="1"><summary>Output files</summary>

    compartments/

        <sample>_combined_compartments.bed: BED file with genomic bins annotated as A or B compartments, combined from all analyzed chromosomes.

        <sample>_combined_compartments.bedGraph: BedGraph file with PC1 eigenvector values for each bin.

    compartments/consensus/

        <sample_group>_compartments_consensus.bed: consensus BED file.

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
