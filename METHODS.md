# CUT&Tag / CoCnT preprocessing in Nextflow

## Workflow structure

The merged workflow supports two input modes:

1. `paired_fastq` mode: discover paired FASTQ files matching `*_R1.fq.gz` and `*_R2.fq.gz`.
2. `demux` mode: preprocess primer and Tn5 annotation tables, normalize raw sciCUT&Tag FASTQ headers, and demultiplex `R1/R2/I1/I2` input files with `sciCTextract`.
3. Optionally exclude samples by filename using configurable patterns.
4. Optionally rewrite FASTQ headers to replace `(s7, s5)` barcode pairs with `Well-ID`.
5. Trim adapters with Cutadapt.
6. Align paired reads with Bowtie2.
7. Convert SAM to BAM.
8. Convert BAM to a deduplicated BED-like fragment file and compress it with `bgzip`.
9. Generate a normalized BigWig track from the BED file.

## Files

- `main.nf`: Nextflow workflow definition.
- `nextflow.config`: defaults for resources and executors.
- `bin/modify_scict_header.sh`: generic header normalization used in `demux` mode.
- `bin/rewrite_fastq_barcodes`: wrapper that prefers the compiled C++ binary.
- `bin/rewrite_fastq_barcodes.py`: barcode rewrite script extracted from the notebook.
- `src/rewrite_fastq_barcodes.cpp`: fast C++ implementation for barcode rewriting.
- `tools/build_rewrite_fastq_barcodes.sh`: build script for the C++ binary.
- `envs/cuttag-preprocess.yml`: Conda environment for all required tools.
- `containers/cuttag-preprocess.sif`: default Singularity/Apptainer image path used by the config.
- `containers/cuttag-preprocess.def`: definition file used to build the Singularity/Apptainer image.

## Required software

The pipeline expects these commands to be available in your environment:

- `nextflow`

If you use the `conda` profile, the workflow will create the software environment automatically from `envs/cuttag-preprocess.yml`.
The environment pins Python 3.10 for compatibility with Cutadapt 4.3.

For faster barcode rewriting, build the compiled helper once:

```bash
./tools/build_rewrite_fastq_barcodes.sh
```

The wrapper used by Nextflow prefers the compiled binary and falls back to Python only if the binary is unavailable.

If you use `singularity` or `apptainer`, the default image path is `containers/cuttag-preprocess.sif`.

You can build that image from the repository with:

```bash
containers/build_singularity_image.sh
```

or manually on a local Linux machine with root or sudo privileges:

```bash
sudo singularity build containers/cuttag-preprocess.sif containers/cuttag-preprocess.def
```

That image should contain:

- `python3`
- `cutadapt`
- `bowtie2`
- `samtools`
- `bedtools`
- `bgzip`
- `bedGraphToBigWig`

## Required inputs

- `--input_mode`: `paired_fastq` or `demux`.
- `--input_dir`: directory containing paired FASTQ files.
- `--ref`: Bowtie2 index basename.
- `--chrom_sizes`: chromosome sizes file for BigWig generation.

For `demux` mode, the required inputs are:

- `--primer_annot`
- `--tn5_annot`
- `--fastq1`
- `--fastq2`
- `--umi1`
- `--umi2`
- `--demux_min_reads`
- `--demux_swap_index_ends`

The `--tn5_annot` input should be supplied as a CSV file.

Barcode rewriting is enabled by default and also requires:

- `--barcode_matrix`: CSV with columns
  `PAGE-1-s7,PAGE-1-s5,PAGE-2-s7,PAGE-2-s5,Well-ID`

For the current sciCUT&Tag data layout used in testing, the index reads are arranged as `I1 = i7 ... j7` and `I2 = j5 ... i5`, while `sciCTextract` expects `I1 = j7 ... i7` and `I2 = i5 ... j5`. The workflow therefore swaps the first and last 8 bases of `I1/I2` by default before demultiplexing. Disable this only for data already matching the native `sciCTextract` layout:

```bash
--demux_swap_index_ends false
```

## Sample filtering options

The filename filter is enabled by default and excludes sample names containing:
By default, no filename patterns are excluded.

Disable filename-based filtering entirely:

```bash
--enable_sample_filter false
```

Override the default exclusion patterns with a comma-separated list:

```bash
--skip_patterns PosCtrl,NegCtrl,MyControl
```

## Example runs

Run locally:

```bash
nextflow run main.nf \
  --input_mode paired_fastq \
  --input_dir /path/to/fastq \
  --barcode_matrix /path/to/barcode_matrix.csv \
  --ref /path/to/bowtie2/index_basename \
  --chrom_sizes /path/to/genome.chrom.sizes \
  --out_dir results
```

Run with Conda:

```bash
nextflow run main.nf -profile conda \
  --input_mode paired_fastq \
  --input_dir /path/to/fastq \
  --barcode_matrix /path/to/barcode_matrix.csv \
  --ref /path/to/bowtie2/index_basename \
  --chrom_sizes /path/to/genome.chrom.sizes \
  --out_dir results
```

Run on SLURM:

```bash
nextflow run main.nf -profile slurm \
  --input_mode paired_fastq \
  --input_dir /path/to/fastq \
  --barcode_matrix /path/to/barcode_matrix.csv \
  --ref /path/to/bowtie2/index_basename \
  --chrom_sizes /path/to/genome.chrom.sizes \
  --out_dir results
```

Run on SLURM with Conda:

```bash
nextflow run main.nf -profile slurm,conda \
  --input_mode paired_fastq \
  --input_dir /path/to/fastq \
  --barcode_matrix /path/to/barcode_matrix.csv \
  --ref /path/to/bowtie2/index_basename \
  --chrom_sizes /path/to/genome.chrom.sizes \
  --out_dir results
```

Run in `demux` mode:

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

Run with Singularity or Apptainer:

```bash
nextflow run main.nf -profile singularity \
  --input_dir /path/to/fastq \
  --barcode_matrix /path/to/barcode_matrix.csv \
  --ref /path/to/bowtie2/index_basename \
  --chrom_sizes /path/to/genome.chrom.sizes \
  --out_dir results
```

Override the default Singularity image location if needed:

```bash
nextflow run main.nf -profile singularity \
  --singularity_image /path/to/container.sif \
  --input_dir /path/to/fastq \
  --barcode_matrix /path/to/barcode_matrix.csv \
  --ref /path/to/bowtie2/index_basename \
  --chrom_sizes /path/to/genome.chrom.sizes \
  --out_dir results
```

## Notes

- The pipeline keeps the same adapter sequence and Bowtie2 arguments used in the existing shell script.
- The merged workflow can run either from already demultiplexed paired FASTQs or from raw sciCUT&Tag demultiplexing inputs.
- The `demux` mode uses a generic FASTQ header normalizer instead of the previous `ModifyHeader.sh` logic that depended on a specific instrument prefix.
- The `demux` handoff flattens the `sciCTextract` file list and matches `*_R1.fq.gz` and `*_R2.fq.gz` outputs explicitly before downstream processing.
- The barcode rewrite step is a required part of the workflow and prefers a compiled C++ implementation for speed, while preserving the original Python code as a fallback.
- The alignment step preserves the original high-memory setting (`256 GB`, `16 CPUs`, `18h`) but these can be changed in `nextflow.config`.
- Intermediate files are published into subdirectories under `--out_dir`.
- The Singularity and Apptainer profiles use `containers/cuttag-preprocess.sif` by default and can be overridden with `--singularity_image`.
- The Singularity and Apptainer profiles bind `./data` from the Nextflow launch directory by default. Override with `--container_bind_paths` if references or data are stored elsewhere.
- Sample-name filtering is configurable through `--enable_sample_filter` and `--skip_patterns`, but no samples are excluded unless patterns are provided explicitly.
