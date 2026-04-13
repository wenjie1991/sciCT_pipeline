# CoCUT&Tag Preprocessing Pipeline

This repository contains a Nextflow pipeline for sequencing preprocessing steps used in CoCUT&Tag analysis.

The workflow supports two input modes:

- `paired_fastq`: start from already demultiplexed paired FASTQs
- `demux`: start from raw sciCUT&Tag-style `R1/R2/I1/I2` FASTQs plus primer and Tn5 annotation tables

The pipeline covers:

- Demultiplexing of FASTQs using the provided index annotation tables
- barcode correction / FASTQ header rewriting using a Well-ID matrix
- adapter trimming with Cutadapt
- paired-end alignment with Bowtie2
- SAM to BAM conversion
- BAM to compressed BED fragment generation
- normalized BigWig generation

## Main Files

- `main.nf`: main Nextflow workflow
- `nextflow.config`: runtime profiles and resource defaults
- `bin/rewrite_fastq_barcodes`: wrapper that prefers the compiled barcode-rewrite binary
- `bin/rewrite_fastq_barcodes.py`: barcode rewrite helper extracted from the notebook
- `bin/modify_scict_header.sh`: generic header normalizer used for sciCT demultiplexing mode
- `src/rewrite_fastq_barcodes.cpp`: fast C++ implementation of barcode rewriting
- `tools/build_rewrite_fastq_barcodes.sh`: build script for the C++ binary
- `envs/cuttag-preprocess.yml`: Conda environment definition
- `envs/cuttag-preprocess-container.yml`: lighter Conda environment used inside the Singularity/Apptainer image
- `containers/cuttag-preprocess.def`: Singularity/Apptainer definition file
- `containers/cuttag-preprocess.sif`: default Singularity image path expected by the config
- `METHODS.md`: extended workflow notes and examples

## Faster Barcode Rewriting

The barcode rewrite step now prefers a compiled C++ implementation because this stage is dominated by gzip I/O and header string processing, which is much faster in native code than in pure Python.

Build the fast binary once:

```bash
./tools/build_rewrite_fastq_barcodes.sh
```

This creates:

```bash
bin/rewrite_fastq_barcodes_cpp
```

The wrapper `bin/rewrite_fastq_barcodes` will:

- use the compiled binary if it already exists
- try to build it automatically with `g++` if available
- fall back to the original Python implementation if compilation is not possible

## Inputs

Required parameters:

- `--input_mode`: `paired_fastq` or `demux`
- `--input_dir`: directory containing paired `*_R1.fq.gz` and `*_R2.fq.gz`
- `--ref`: Bowtie2 index basename
- `--chrom_sizes`: chromosome sizes file for BigWig generation

For `demux` mode, use these instead of `--input_dir`:

- `--primer_annot`: Primer annotation file.
- `--tn5_annot`: Tn5 barcode annotation file.
- `--fastq1`: The R1 FASTQ file
- `--fastq2`: The R2 FASTQ file
- `--umi1`: The I1 FASTQ file
- `--umi2`: The I2 FASTQ file
- `--demux_min_reads`: Minimum total reads for a sample to be retained after demultiplexing, default `10000`.
- `--demux_swap_index_ends`: swap the first and last 8 bases of `I1/I2` before demultiplexing, default `true`.

Both annotation files are required for `demux` mode.

For the current sciCUT&Tag test data, the index reads are arranged as:

- `I1`: `i7 ... j7`
- `I2`: `j5 ... i5`

while `sciCTextract` expects:

- `I1`: `j7 ... i7`
- `I2`: `i5 ... j5`

The pipeline therefore swaps the first and last 8 bases of each index read by default in `demux` mode. Disable this only if your run already matches the native `sciCTextract` layout:

```bash
--demux_swap_index_ends false
```

Primer annotation required columns:

- `Sample` or `ID`
- `i7_index_seq`
- `i5_index_seq`

The preprocessing step converts this to:

- `ID`
- `i7_index_seq`
- `i5_index_seq`
- `i7_index_id`
- `i5_index_id`

Tn5 annotation required columns:

- `Sample Name`
- `Tn5_s7`
- `Tn5_s7_seq`
- `Tn5_s5`
- `Tn5_s5_seq`

Barcode rewriting input:

- `--barcode_matrix`: CSV with columns `PAGE-1-s7,PAGE-1-s5,PAGE-2-s7,PAGE-2-s5,Well-ID`

## Sample Filtering (Optinal remove in the future?)


Filename-based sample filtering is configurable.

By default, no samples are excluded by filename.

Disable sample filtering completely:

```bash
--enable_sample_filter false
```

Provide your own comma-separated exclusion patterns:

```bash
--skip_patterns PosCtrl,NegCtrl,MyControl
```

## Runtime Options

### Conda

Use the included Conda environment:

```bash
nextflow run main.nf -profile conda \
  --input_mode paired_fastq \
  --input_dir /path/to/fastq \
  --barcode_matrix /path/to/barcode_matrix.csv \
  --ref /path/to/bowtie2/index_basename \
  --chrom_sizes /path/to/genome.chrom.sizes \
  --out_dir results
```

For SLURM:

```bash
nextflow run main.nf -profile slurm,conda \
  --input_mode paired_fastq \
  --input_dir /path/to/fastq \
  --barcode_matrix /path/to/barcode_matrix.csv \
  --ref /path/to/bowtie2/index_basename \
  --chrom_sizes /path/to/genome.chrom.sizes \
  --out_dir results
```

Example sciCT demultiplexing mode:

```bash
nextflow run main.nf -profile slurm,conda \
  --input_mode demux \
  --primer_annot /path/to/Primer_Annotation.csv \
  --tn5_annot /path/to/Tn5_Barcode_Annotation.csv \
  --fastq1 /path/to/sample_R1.fastq.gz \
  --fastq2 /path/to/sample_R2.fastq.gz \
  --umi1 /path/to/sample_I1.fastq.gz \
  --umi2 /path/to/sample_I2.fastq.gz \
  --demux_min_reads 10000 \
  --barcode_matrix /path/to/barcode_matrix.csv \
  --ref /path/to/bowtie2/index_basename \
  --chrom_sizes /path/to/genome.chrom.sizes \
  --out_dir results
```

### Singularity / Apptainer

Build the image from the included definition file:

```bash
containers/build_singularity_image.sh
```

On a local Linux machine where you have sudo privileges, the helper runs:

```bash
sudo singularity build containers/cuttag-preprocess.sif containers/cuttag-preprocess.def
```

If you use Apptainer instead:

```bash
sudo apptainer build containers/cuttag-preprocess.sif containers/cuttag-preprocess.def
```

If you cannot build locally, use a remote builder or ask your HPC admins to build it:

```bash
singularity build --remote containers/cuttag-preprocess.sif containers/cuttag-preprocess.def
```

Equivalent helper command:

```bash
containers/build_singularity_image.sh remote
```

Then run:

```bash
nextflow run main.nf -profile singularity \
  --input_dir /path/to/fastq \
  --barcode_matrix /path/to/barcode_matrix.csv \
  --ref /path/to/bowtie2/index_basename \
  --chrom_sizes /path/to/genome.chrom.sizes \
  --out_dir results
```

To use a different image location:

```bash
nextflow run main.nf -profile singularity \
  --singularity_image /path/to/container.sif \
  --input_dir /path/to/fastq \
  --barcode_matrix /path/to/barcode_matrix.csv \
  --ref /path/to/bowtie2/index_basename \
  --chrom_sizes /path/to/genome.chrom.sizes \
  --out_dir results
```

If your references or data are outside the working directory, bind those filesystem roots into the container. The default bind path is `./data` relative to the directory where you launch Nextflow.

```bash
nextflow run main.nf -profile singularity \
  --container_bind_paths ./data \
  --input_dir /path/to/fastq \
  --barcode_matrix /path/to/barcode_matrix.csv \
  --ref /varidata/research/projects/bbc/versioned_references/latest/data/hg38_gencode/indexes/bowtie2/hg38_gencode \
  --chrom_sizes /varidata/research/projects/bbc/versioned_references/2024-10-31_10.56.03_v17/data/hg38_gencode/sequence/hg38_gencode.fa.fai \
  --out_dir results
```

For multiple roots, provide a comma-separated list:

```bash
--container_bind_paths ./data,/scratch,/home
```

## Defaults

- sample-name filtering is available but no filename patterns are excluded unless `--skip_patterns` is provided
- demux mode keeps only samples with total reads greater than `--demux_min_reads`, default `10000`
- demux mode swaps the first and last 8 bases of `I1/I2` by default before calling `sciCTextract`
- barcode rewriting is always applied and requires `--barcode_matrix`
- adapter sequence defaults to `CTGTCTCTTATACACATCT`
- alignment resources default to `16 CPUs`, `256 GB`, and `18h`

## Notes

- `nextflow` must be installed on the host system.
- The Conda profile creates the software environment automatically from `envs/cuttag-preprocess.yml`.
- The Singularity and Apptainer profiles use `containers/cuttag-preprocess.sif` by default.
- Singularity and Apptainer bind `./data` from the launch directory by default; override with `--container_bind_paths` if needed.
- Sample filtering is controlled by `--enable_sample_filter` and `--skip_patterns`.
