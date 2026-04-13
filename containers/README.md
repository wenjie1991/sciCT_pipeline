Place the Singularity or Apptainer image for this pipeline at:

`containers/cuttag-preprocess.sif`

The image should include:

- python3
- cutadapt
- bowtie2
- samtools
- bedtools
- bgzip
- bedGraphToBigWig

The container build uses `envs/cuttag-preprocess-container.yml`, a smaller runtime environment than the Nextflow Conda profile. The Conda profile environment keeps the compiler package needed to build the optional fast barcode rewrite binary on HPC.

Build the image from the included definition file on a local Linux machine where you have root or sudo privileges:

```bash
containers/build_singularity_image.sh
```

Equivalent manual command:

```bash
sudo singularity build containers/cuttag-preprocess.sif containers/cuttag-preprocess.def
```

If you use Apptainer instead:

```bash
sudo apptainer build containers/cuttag-preprocess.sif containers/cuttag-preprocess.def
```

If you cannot build locally, use your site's remote builder or ask the cluster admins to build the image:

```bash
singularity build --remote containers/cuttag-preprocess.sif containers/cuttag-preprocess.def
```

Equivalent helper command:

```bash
containers/build_singularity_image.sh remote
```

If you store the image elsewhere, override it with:

```bash
nextflow run main.nf -profile singularity --singularity_image /path/to/container.sif ...
```
