nextflow.enable.dsl = 2

params.input_mode = params.input_mode ?: 'paired_fastq'
params.input_dir = params.input_dir ?: null
params.primer_annot = params.primer_annot ?: null
params.tn5_annot = params.tn5_annot ?: null
params.fastq1 = params.fastq1 ?: null
params.fastq2 = params.fastq2 ?: null
params.umi1 = params.umi1 ?: null
params.umi2 = params.umi2 ?: null
params.demux_min_reads = params.demux_min_reads ?: 10000
// When true, swap the first and last 8 bases of I1/I2 during demux FASTQ
// normalization so the index reads match sciCTextract's expected layout.
params.demux_swap_index_ends = params.demux_swap_index_ends == null ? true : params.demux_swap_index_ends
params.out_dir = params.out_dir ?: 'results'
params.barcode_matrix = params.barcode_matrix ?: null
params.barcode_suffix = params.barcode_suffix ?: ''
// When false, skip the barcode-rewrite step because the input FASTQ headers
// already carry the rewritten Well-ID barcodes.
params.enable_barcode_rewrite = params.enable_barcode_rewrite == null ? true : params.enable_barcode_rewrite
params.enable_sample_filter = params.enable_sample_filter == null ? true : params.enable_sample_filter
params.skip_patterns = params.skip_patterns ?: []
params.adapter_seq = params.adapter_seq ?: 'CTGTCTCTTATACACATCT'
params.ref = params.ref ?: null
params.chrom_sizes = params.chrom_sizes ?: null
params.publish_mode = params.publish_mode ?: 'copy'

def validInputModes = ['paired_fastq', 'demux']
if (!validInputModes.contains(params.input_mode)) {
    error "Unsupported --input_mode '${params.input_mode}'. Valid modes: ${validInputModes.join(', ')}"
}

if (params.enable_barcode_rewrite && !params.barcode_matrix) {
    error "Missing required parameter for barcode rewriting: --barcode_matrix (or set --enable_barcode_rewrite false if the FASTQ headers are already rewritten)"
}

if (params.input_mode == 'paired_fastq' && !params.input_dir) {
    error "Missing required parameter for paired_fastq mode: --input_dir"
}

if (params.input_mode == 'demux') {
    ['primer_annot', 'tn5_annot', 'fastq1', 'fastq2', 'umi1', 'umi2'].each { key ->
        if (!params[key]) {
            error "Missing required parameter for demux mode: --${key}"
        }
    }
}

if (!params.ref) {
    error "Missing required parameter: --ref"
}

if (!params.chrom_sizes) {
    error "Missing required parameter: --chrom_sizes"
}

def normalizePatternList(patterns) {
    if (patterns == null) {
        return []
    }
    if (patterns instanceof CharSequence) {
        return patterns
            .toString()
            .split(',')
            .collect { it.trim() }
            .findAll { it }
    }
    return patterns.collect { it.toString().trim() }.findAll { it }
}

def skipPatternList = normalizePatternList(params.skip_patterns)
def skipRegex = skipPatternList ? skipPatternList.collect { java.util.regex.Pattern.quote(it) }.join('|') : null
def barcodeMatrixFile = params.enable_barcode_rewrite ? file(params.barcode_matrix, checkIfExists: true) : null

process PROCESS_PRIMER_ANNOT {
    tag "${primer_annot.simpleName}"
    cpus 1
    memory '4 GB'
    publishDir "${params.out_dir}/demux", mode: params.publish_mode, pattern: '*_primer.csv'

    input:
    path primer_annot

    output:
    path "primer_annotation.csv"

    script:
    """
    python ${projectDir}/bin/AddPrimerAnnotationCols.py \
      -i ${primer_annot} \
      -o primer_annotation.csv
    sed -i 's/\r\$//' primer_annotation.csv
    """
}

process PROCESS_TN5_ANNOT {
    tag "${tn5_annot.simpleName}"
    cpus 1
    memory '4 GB'
    publishDir "${params.out_dir}/demux", mode: params.publish_mode, pattern: '*_tn5.csv'

    input:
    path tn5_annot

    output:
    path "tn5_annotation.csv"

    script:
    """
    python ${projectDir}/bin/Check_Tn5_file.py \
      -i ${tn5_annot} \
      -o tn5_annotation.csv
    sed -i 's/\r\$//' tn5_annotation.csv
    """
}

process PROCESS_DEMUX_FASTQ {
    tag "${fastq.simpleName}"
    cpus 1
    memory '4 GB'
    publishDir "${params.out_dir}/demux/modified_fastq", mode: params.publish_mode, pattern: '*.fastq.gz'

    input:
    tuple val(idx), path(fastq)
    path helper_script

    output:
    tuple val(idx), path("${fastq.simpleName}.mod.fastq.gz")

    script:
    """
    bash ${helper_script} ${fastq} ${fastq.simpleName}.mod.fastq.gz ${params.demux_swap_index_ends}
    """
}

process RUN_DEMUX {
    tag "sciCTextract"
    cpus 4
    memory '32 GB'
    publishDir "${params.out_dir}/demux", mode: params.publish_mode, pattern: 'demux_out/*.fq.gz'

    input:
    path primer_annot2
    path tn5_annot2
    val processed_fastqs

    output:
    path "demux_out/*.fq.gz"

    script:
    """
    mkdir -p demux_out
    sciCTextract --forward-mode \
      --outdir demux_out \
      --Primer_Barcode ${primer_annot2} \
      --Tn5_Barcode ${tn5_annot2} \
      ${processed_fastqs.join(' ')}
    """
}

process COUNT_READS {
    tag "${sid}"
    cpus 1
    memory '4 GB'

    input:
    tuple val(sid), path(r1), path(r2)

    output:
    tuple val(sid), path(r1), path(r2), path("reads.txt")

    script:
    """
    n1=\$(zcat ${r1} | awk 'END{print NR/4}')
    n2=\$(zcat ${r2} | awk 'END{print NR/4}')
    echo \$(( \${n1:-0} + \${n2:-0} )) > reads.txt
    """
}

process WRITE_DEMUX_QC {
    tag "demux-qc"
    cpus 1
    memory '4 GB'
    publishDir "${params.out_dir}/demux", mode: params.publish_mode

    input:
    val qc_rows

    output:
    path "demux_qc.tsv"

    script:
    def lines = qc_rows.join('\n')
    """
    cat > demux_qc.tsv <<'EOF'
    status\tsample_id\treads
    ${lines}
    EOF
    """
}

process REWRITE_BARCODES {
    tag "${sample_id}:${read_label}"
    cpus 1
    memory '4 GB'
    publishDir "${params.out_dir}/rewritten_fastq", mode: params.publish_mode

    input:
    tuple val(sample_id), val(read_label), path(read_file)
    path barcode_matrix

    output:
    tuple val(sample_id), val(read_label), path("${sample_id}_${read_label}.rewritten.fq.gz")

    script:
    """
    rewrite_fastq_barcodes \
      --input ${read_file} \
      --output ${sample_id}_${read_label}.rewritten.fq.gz \
      --matrix ${barcode_matrix} \
      --suffix '${params.barcode_suffix}'
    """
}

process TRIM_ADAPTERS {
    tag "${sample_id}"
    cpus 8
    memory '32 GB'
    publishDir "${params.out_dir}/trimmed_fastq", mode: params.publish_mode, pattern: '*.fastq.gz'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.trimmed.fastq.gz"), path("${sample_id}_R2.trimmed.fastq.gz")

    script:
    """
    cutadapt \
      -j ${task.cpus} \
      -m 20 \
      -a ${params.adapter_seq} \
      -A ${params.adapter_seq} \
      -o ${sample_id}_R1.trimmed.fastq.gz \
      -p ${sample_id}_R2.trimmed.fastq.gz \
      ${r1} ${r2}
    """
}

process ALIGN_READS {
    tag "${sample_id}"
    cpus 16
    memory '256 GB'
    time '18h'
    publishDir "${params.out_dir}/sam", mode: params.publish_mode, pattern: '*.sam'
    publishDir "${params.out_dir}/bam", mode: params.publish_mode, pattern: '*.bam'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}.sam"), path("${sample_id}.bam")

    script:
    """
    bowtie2 \
      -p ${task.cpus} \
      --very-sensitive-local \
      --soft-clipped-unmapped-tlen \
      --no-mixed \
      --no-discordant \
      --dovetail \
      --phred33 \
      -I 10 \
      -X 1000 \
      -x ${params.ref} \
      -1 ${r1} \
      -2 ${r2} \
      > ${sample_id}.sam

    samtools view -S -b ${sample_id}.sam -o ${sample_id}.bam
    """
}

process BAM_TO_BED {
    tag "${sample_id}"
    cpus 4
    memory '32 GB'
    publishDir "${params.out_dir}/bed", mode: params.publish_mode, pattern: '*.bed.gz'

    input:
    tuple val(sample_id), path(sam_file), path(bam_file)

    output:
    tuple val(sample_id), path("${sample_id}.bed.gz")

    script:
    """
    bedtools bamtobed -bedpe -i ${bam_file} | \
      cut -f1,2,6,7 | \
      awk -F'\\t' '\$1 != "." || \$2 != -1 || \$3 != -1' | \
      awk 'BEGIN{OFS="\\t"} {barcode=\$4; sub(/\\/[12]\$/, "", barcode); ncolon=split(barcode, c, ":"); barcode=c[ncolon]; n=split(barcode, q, "_"); if (n >= 3) name=q[n-2] "_" q[n-1] "_" q[n]; else name=barcode; print \$1, \$2, \$3, name}' | \
      sort -k1,1 -k2,2n -k4,4 | \
      uniq -c | \
      awk '{print \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 "\\t" \$1}' | \
      sort -k1,1 -k2,2n -k3,3n \
      > ${sample_id}.bed

    bgzip -c ${sample_id}.bed > ${sample_id}.bed.gz
    """
}

process BED_TO_BIGWIG {
    tag "${sample_id}"
    cpus 4
    memory '32 GB'
    publishDir "${params.out_dir}/bigwig", mode: params.publish_mode, pattern: '*.bigwig'

    input:
    tuple val(sample_id), path(bed_gz)

    output:
    tuple val(sample_id), path("${sample_id}.bigwig")

    script:
    """
    coverage=\$(zcat ${bed_gz} | awk '{s+=(\$3-\$2)} END {print s}')
    if [ -z "\$coverage" ] || [ "\$coverage" = "0" ]; then
      echo "Coverage is zero for ${sample_id}, cannot scale BigWig." >&2
      exit 1
    fi

    scaling_factor=\$(awk -v c="\$coverage" 'BEGIN {print (1/c)*10^10}')

    zcat ${bed_gz} | \
      bedtools genomecov -bg -i stdin -g ${params.chrom_sizes} -scale "\$scaling_factor" \
      > ${sample_id}.bedgraph

    bedGraphToBigWig ${sample_id}.bedgraph ${params.chrom_sizes} ${sample_id}.bigwig
    """
}

workflow {
    if (params.input_mode == 'paired_fastq') {
        source_fastq_pairs = Channel
            .fromFilePairs("${params.input_dir}/*_{R1,R2}.fq.gz", flat: true, checkIfExists: true)
            .filter { row ->
                if (!params.enable_sample_filter || !skipRegex) {
                    return true
                }
                return !row[0].matches(".*(${skipRegex}).*")
            }
    } else {
        primer_annot_ch = Channel.value(file(params.primer_annot, checkIfExists: true))
        tn5_annot_ch = Channel.value(file(params.tn5_annot, checkIfExists: true))
        raw_demux_fastqs = Channel.of(
            tuple(0, file(params.fastq1, checkIfExists: true)),
            tuple(1, file(params.fastq2, checkIfExists: true)),
            tuple(2, file(params.umi1, checkIfExists: true)),
            tuple(3, file(params.umi2, checkIfExists: true))
        )
        demux_helper_script = Channel.value(file("${projectDir}/bin/modify_scict_header.sh", checkIfExists: true))

        primer_annot_out = PROCESS_PRIMER_ANNOT(primer_annot_ch)
        tn5_annot_out = PROCESS_TN5_ANNOT(tn5_annot_ch)

        // Keep R1/R2/I1/I2 in a stable order for sciCTextract positional args.
        ordered_fastqs = PROCESS_DEMUX_FASTQ(raw_demux_fastqs, demux_helper_script)
            .toList()
            .map { items -> items.sort { a, b -> a[0] <=> b[0] } }
            .map { items -> items.collect { it[1] } }

        demux_fastqs = RUN_DEMUX(primer_annot_out, tn5_annot_out, ordered_fastqs)
            .flatten()

        source_fastq_pairs = demux_fastqs
            .map { f ->
                def m = (f.name =~ /^(.+)_R([12])\.fq\.gz$/)
                assert m : "Unexpected demux filename: ${f.name}"
                tuple(m[0][1], m[0][2] as int, f)
            }
            .groupTuple(by: 0)
            .map { sid, mates, files ->
                def file_map = [:]
                mates.eachWithIndex { mate, idx -> file_map[mate as int] = files[idx] }
                assert file_map[1] && file_map[2] : "Missing mate for ${sid}"
                tuple(sid, file_map[1], file_map[2])
            }

        demux_read_counts = COUNT_READS(source_fastq_pairs)
            .map { sid, r1, r2, reads_file ->
                tuple(sid, r1, r2, reads_file.text.trim() as long)
            }

        demux_qc = demux_read_counts
            .map { sid, r1, r2, n ->
                def status = n >= (params.demux_min_reads as long) ? 'KEPT' : 'LOW_READS'
                "${status}\t${sid}\t${n}"
            }
            .collect()
            | WRITE_DEMUX_QC

        source_fastq_pairs = demux_read_counts
            .filter { sid, r1, r2, n -> n >= (params.demux_min_reads as long) }
            .map { sid, r1, r2, n -> tuple(sid, r1, r2) }
            .filter { row ->
                if (!params.enable_sample_filter || !skipRegex) {
                    return true
                }
                return !row[0].matches(".*(${skipRegex}).*")
            }
    }

    if (params.enable_barcode_rewrite) {
        rewrite_reads = source_fastq_pairs
            .flatMap { row ->
                def sample_id = row[0]
                def r1 = row[1]
                def r2 = row[2]
                [
                    tuple(sample_id, 'R1', r1),
                    tuple(sample_id, 'R2', r2)
                ]
            }

        rewritten_fastq = REWRITE_BARCODES(rewrite_reads, barcodeMatrixFile)
            .groupTuple(by: 0)
            .map { sample_id, read_labels, files ->
                def read_map = [:]
                read_labels.eachWithIndex { label, idx -> read_map[label] = files[idx] }
                tuple(sample_id, read_map['R1'], read_map['R2'])
            }
    } else {
        // FASTQ headers already carry rewritten Well-ID barcodes; pass through.
        rewritten_fastq = source_fastq_pairs
    }

    trimmed = rewritten_fastq | TRIM_ADAPTERS
    aligned = trimmed | ALIGN_READS
    bed = aligned | BAM_TO_BED
    bigwig = bed | BED_TO_BIGWIG

    bigwig.view { sample_id, bw -> "Completed ${sample_id}: ${bw}" }
}
