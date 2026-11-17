<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-methylseq_logo_dark.png">
    <img alt="nf-core/methylseq" src="docs/images/nf-core-methylseq_logo_light.png">
  </picture>
</h1>

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/nf-core/methylseq)
[![GitHub Actions CI Status](https://github.com/nf-core/methylseq/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-core/methylseq/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/methylseq/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/methylseq/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/methylseq/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.1343417-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.1343417)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.1)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/methylseq)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23methylseq-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/methylseq)[![Follow on Bluesky](https://img.shields.io/badge/bluesky-%40nf__core-1185fe?labelColor=000000&logo=bluesky)](https://bsky.app/profile/nf-co.re)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

This extended methylation analysis pipeline is developed at the [National Facility for Data Handling and Analysis, Human Technopole](https://humantechnopole.it/en/facilities/national-facility-for-data-handling-analysis/), building on **nf-core/methylseq** to bring standard processing, differential methylation analysis, and genomic annotation into one Nextflow workflow.

**nf-core/methylseq** is a bioinformatics analysis pipeline used for Methylation (Bisulfite) sequencing data. It pre-processes raw data from FastQ inputs, aligns the reads and performs extensive quality-control on the results.

This repository extends the cloned **nf-core/methylseq** pipeline with downstream differential methylation analysis and genomic annotation. Standard processing and advanced analysis can now run within the same Nextflow pipeline: start from FASTQ files and continue through differential analysis and annotation, or run advanced analysis independently on existing Bismark coverage files.

The original standard-analysis description is retained below. The additions specific to this repository are described in [Advanced analysis features](#advanced-analysis-features) and [Selecting standard and advanced analysis](#selecting-standard-and-advanced-analysis). Commands using `nf-core/methylseq` refer to the upstream pipeline; use `nextflow run .` from this cloned repository to access the added features.

![nf-core/methylseq metro map](docs/images/4.2.0_metromap.png)

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker / Singularity / Podman / Charliecloud / Apptainer containers making installation trivial and results highly reproducible.

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources.The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/methylseq/results).

> Read more about **Bisulfite Sequencing & Three-Base Aligners** used in this pipeline [here](docs/usage/bs-seq-primer.md)

## Pipeline Summary

### Standard analysis features

The pipeline allows you to choose between running either [Bismark](https://github.com/FelixKrueger/Bismark), [bwa-meth](https://github.com/brentp/bwa-meth) / [MethylDackel](https://github.com/dpryan79/methyldackel) or [BWA-Mem](https://github.com/lh3/bwa) plus [rastair](https://bitbucket.org/bsblabludwig/rastair/src/master/) for TAPS data processing. rastair can also be used with bwa-meth aligned reads by setting the aligner to `--aligner bwameth` and adding the flag `--taps`.

Choose between workflows by using `--aligner bismark` (default, uses bowtie2 for alignment), `--aligner bismark_hisat`, `--aligner bwameth` or `--aligner bwamem`. For higher performance, the pipeline can leverage the [Parabricks implementation of bwa-meth (fq2bammeth)](https://docs.nvidia.com/clara/parabricks/latest/documentation/tooldocs/man_fq2bam_meth.html) and the [Parabricks implementation of bwa-mem (fq2bammemh)](https://docs.nvidia.com/clara/parabricks/latest/documentation/tooldocs/man_fq2bam.html), which implement the baseline tools `bwa-meth` and `bwa-mem`. To use this option, include the `gpu` profile along with `--aligner bwameth` or `--aligner bwamem`.

Note: For faster CPU runs with BWA-Meth, enable the BWA-MEM2 algorithm using `--use_mem2`. The GPU pathway (Parabricks) requires `-profile gpu` and a container runtime (Docker, Singularity, or Podman); Conda/Mamba are not supported for the GPU module.

| Step                                         | Bismark workflow         | bwa-meth workflow     | bwa-mem + TAPS workflow    |
| -------------------------------------------- | ------------------------ | --------------------- | -------------------------- |
| Generate Reference Genome Index _(optional)_ | Bismark                  | bwa-meth              | bwa index                  |
| Merge re-sequenced FastQ files               | cat                      | cat                   | cat                        |
| Raw data QC                                  | FastQC                   | FastQC                | FastQC                     |
| Adapter sequence trimming                    | Trim Galore!             | Trim Galore!          | Trim Galore!               |
| Align Reads                                  | Bismark (bowtie2/hisat2) | bwa-meth              | bwa mem                    |
| Deduplicate Alignments                       | Bismark                  | Picard MarkDuplicates | Picard MarkDuplicates      |
| Extract methylation calls                    | Bismark                  | MethylDackel          | TAPS subworkflow (rastair) |
| Sample report                                | Bismark                  | -                     | -                          |
| Summary Report                               | Bismark                  | -                     | -                          |
| Alignment QC                                 | Qualimap _(optional)_    | Qualimap _(optional)_ | Qualimap _(optional)_      |
| Sample complexity                            | Preseq _(optional)_      | Preseq _(optional)_   | Preseq _(optional)_        |
| Project Report                               | MultiQC                  | MultiQC               | MultiQC                    |

Optional targeted sequencing analysis is available via `--run_targeted_sequencing` and `--target_regions_file`; see the [usage documentation](https://nf-co.re/methylseq/usage) for details.

### Advanced analysis features

The advanced workflow adds the following steps to standard methylation processing:

- **Methylation preprocessing and QC with methylKit:** coverage filtering and normalization, shared CpG tables, PCA, sample clustering, and coverage and methylation summaries.
- **Differential methylation analysis:** choose `methylkit` for differentially methylated cytosines (DMCs), `dss` for locus-level results (DMLs) and differentially methylated regions (DMRs), or `dmr_detector` for region detection from paired sample methylation profiles. Multiple methods can be selected in one run.
- **Comparison design:** generate comparisons from sample metadata, optionally within strata, or supply an explicit comparison table.
- **Genomic annotation of methylKit and DSS results:** gene parts, distance to transcription start sites (TSSs), CpG islands and shores, and GTF-based promoter, gene-body, and nearest-TSS annotations. Outputs include annotation tables, plots, and enrichment summaries where the required classes are available.

The advanced workflow currently reads **Bismark coverage files**. The standard aligner options described above remain available for primary analysis. `dmr_detector` currently supports one case and one control sample per comparison; its results do not pass through the shared annotation step.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,genome
SRR389222_sub1,https://github.com/nf-core/test-datasets/raw/methylseq/testdata/SRR389222_sub1.fastq.gz,,
SRR389222_sub2,https://github.com/nf-core/test-datasets/raw/methylseq/testdata/SRR389222_sub2.fastq.gz,,
SRR389222_sub3,https://github.com/nf-core/test-datasets/raw/methylseq/testdata/SRR389222_sub3.fastq.gz,,
Ecoli_10K_methylated,https://github.com/nf-core/test-datasets/raw/methylseq/testdata/Ecoli_10K_methylated_R1.fastq.gz,https://github.com/nf-core/test-datasets/raw/methylseq/testdata/Ecoli_10K_methylated_R2.fastq.gz,
```

> Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

Now, you can run the pipeline using default parameters as:

```bash
nextflow run nf-core/methylseq --input samplesheet.csv --outdir <OUTDIR> --genome GRCh37 -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/methylseq/usage) and the [parameter documentation](https://nf-co.re/methylseq/parameters).

### Selecting standard and advanced analysis

Two parameters control which parts of this repository's pipeline run. **Standard analysis** corresponds to `primary_analysis`; **advanced analysis** corresponds to `advanced_analysis`.

| Mode | `--primary_analysis` | `--advanced_analysis` | Input |
| --- | --- | --- | --- |
| Standard only (default) | `true` | `false` | FASTQ samplesheet described above |
| Standard followed by advanced | `true` | `true` | FASTQ samplesheet and metadata CSV |
| Advanced only | `false` | `true` | Coverage samplesheet and metadata CSV |

Setting both parameters to `false` is invalid. Run the following examples from the root of this cloned repository.

**Standard only:** retains the original preprocessing, alignment, methylation calling, and QC workflow.

```bash
nextflow run . -profile docker \
    --input samplesheet.csv \
    --outdir results/standard \
    --genome mm10 \
    --primary_analysis true \
    --advanced_analysis false
```

**Standard followed by advanced:** passes the Bismark coverage files directly into differential analysis and annotation.

```bash
nextflow run . -profile docker \
    --input samplesheet.csv \
    --metadata metadata.csv \
    --outdir results/combined \
    --genome mm10 \
    --aligner bismark \
    --primary_analysis true \
    --advanced_analysis true \
    --differential_analysis_methods methylkit,dss \
    --sample_suffix '<COVERAGE_FILENAME_SUFFIX>'
```

**Advanced only:** reuses coverage files without repeating FASTQ processing or alignment.

```bash
nextflow run . -profile docker \
    --input coverage_samplesheet.csv \
    --metadata metadata.csv \
    --outdir results/advanced \
    --genome mm10 \
    --primary_analysis false \
    --advanced_analysis true \
    --differential_analysis_methods methylkit,dss \
    --sample_suffix '<COVERAGE_FILENAME_SUFFIX>'
```

Replace `<COVERAGE_FILENAME_SUFFIX>` with the exact suffix to remove from coverage filenames to obtain the sample IDs in metadata. For example, `case_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz` requires the suffix `_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz` to recover the ID `case`. The suffix depends on trimming, alignment, paired-end status, and deduplication; check the actual filenames for your run.

#### Advanced-analysis inputs

For advanced-only runs, `coverage_samplesheet.csv` contains `sample` and `meth_cov` columns. Use absolute paths for coverage files:

```csv
sample,meth_cov
control,/data/control_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
case,/data/case_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz
```

Metadata is a CSV whose first column contains sample IDs matching the coverage filenames after suffix removal. The default grouping column is `genotype`, and the default case label is `Mut`:

```csv
sample,genotype
control,WT
case,Mut
```

Use `--group_column` and `--group_case` to change these defaults. Automatic comparisons use the case label against all other labels within each stratum. Set `--comparison_stratify_by tissue,sex`, for example, to generate separate comparisons within combinations of those metadata columns. Each stratum needs both case and control samples to generate a comparison.

Alternatively, provide `--comparisons comparisons.tsv` with these **tab-separated** columns; sample lists within a cell are comma-separated:

```text
comparison_id	case_samples	control_samples
Mut_vs_WT	case	control
```

Choose methods with `--differential_analysis_methods methylkit`, `dss`, `dmr_detector`, a comma-separated combination, or `all`. The default is `methylkit`. DSS supports either one sample per group or at least two samples in both groups; selecting `all` also requires the one-case/one-control design supported by `dmr_detector`.

#### Filtering and annotation settings

The main methylKit settings are `--lo_count`, `--lo_perc`, `--hi_count`, and `--hi_perc` for coverage filtering; `--destrand` and `--min_per_group` for merging CpGs; and `--diff_cutoff`, `--qvalue_cutoff`, `--overdispersion`, and `--adjust` for differential testing. The default methylKit effect-size and q-value cutoffs are 25 percentage points and 0.01. DSS uses its own calling defaults; these methylKit testing options do not configure the DSS caller.

Annotation for methylKit and DSS requires `--refseq_bed`, `--cpg_bed`, and a GTF. The GTF is resolved from `--genome` when available, or can be supplied explicitly with `--gtf`. Promoter windows are controlled by `--annotation_prom_up` (default 1000 bp) and `--annotation_prom_down` (default 200 bp).

The bundled annotation BED files and default `--assembly` label are **mm10**. For another reference assembly, supply matching BED files and `--assembly`, and ensure the GTF matches the reference used for alignment.

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/methylseq/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the [output documentation](https://nf-co.re/methylseq/output).

Advanced analysis adds the following directories under `--outdir`, according to the selected methods:

| Directory | Contents |
| --- | --- |
| `methylkit/preprocessing/` | Preprocessed methylation objects, methylKit database, shared CpG matrix, QC tables, and plots; also used by `dmr_detector` |
| `methylkit/dma/<comparison_id>/` | Raw differential results, significant hyper/hypomethylated DMC tables, and plots |
| `dss/dma/<comparison_id>/` | DML and DMR tables, comparison summaries, and saved R objects |
| `dmr_detector/<comparison_id>/` | Hyper/hypomethylated region and window tables |
| `<caller>/annotation/<comparison_id>/` | Annotated methylKit or DSS results, gene/TSS and CpG-context tables, plots, and enrichment summaries |

Standard reports remain part of primary analysis; advanced-only runs do not generate the primary-analysis MultiQC report.

## Credits

The extensions in this repository were developed at the [National Facility for Data Handling and Analysis, Human Technopole](https://humantechnopole.it/en/facilities/national-facility-for-data-handling-analysis/) by:

- Youssef Abili ([@YussAb](https://github.com/YussAb))
- Giulia Scotti ([@laPepatencia](https://github.com/laPepatencia))

Their contributions extend the upstream pipeline with integrated advanced analysis and genomic annotation.

nf-core/methylseq was originally written by Phil Ewels ([@ewels](https://github.com/ewels)), and Sateesh Peri ([@sateeshperi](https://github.com/sateeshperi)) is its active maintainer.

We thank the following people for their extensive assistance in the development of this pipeline:

- Felix Krueger ([@FelixKrueger](https://github.com/FelixKrueger))
- Edmund Miller ([@EMiller88](https://github.com/emiller88))
- Rickard Hammarén ([@Hammarn](https://github.com/Hammarn/))
- Alexander Peltzer ([@apeltzer](https://github.com/apeltzer/))
- Patrick Hüther ([@phue](https://github.com/phue/))
- Maxime U Garcia ([@maxulysse](https://github.com/maxulysse/))

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#methylseq` channel](https://nfcore.slack.com/channels/methylseq) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

When using this extended pipeline, please acknowledge the [National Facility for Data Handling and Analysis, Human Technopole](https://humantechnopole.it/en/facilities/national-facility-for-data-handling-analysis/) and the development contributions of Youssef Abili and Giulia Scotti, alongside the upstream pipeline and relevant analysis tools. Suggested acknowledgement:

> Methylation analysis was performed using the extended nf-core/methylseq pipeline developed by Youssef Abili and Giulia Scotti at the National Facility for Data Handling and Analysis, Human Technopole, integrating standard processing, differential methylation analysis, and genomic annotation.

If you use nf-core/methylseq for your analysis, please cite it using the following doi: [10.5281/zenodo.1343417](https://doi.org/10.5281/zenodo.1343417)

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
