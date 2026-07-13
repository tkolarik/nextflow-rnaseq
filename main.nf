#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.samplesheet = params.samplesheet ?: 'data/samplesheet.csv'
params.transcriptome = params.transcriptome ?: 'https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/reference/prokaryotic/SL1344_sub.fasta'
params.outdir = params.outdir ?: 'results'
params.salmon_libtype = params.salmon_libtype ?: 'A'
params.salmon_min_assigned_frags = params.salmon_min_assigned_frags ?: 10

workflow {
    main:
    raw_samples = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            tuple(row.sample_id, row.condition ?: 'unknown', row.fastq_1, row.fastq_2)
        }

    url_samples = raw_samples
        .filter { sample_id, condition, fastq_1, fastq_2 -> fastq_1 ==~ /^https?:\/\/.*/ && fastq_2 ==~ /^https?:\/\/.*/ }

    local_samples = raw_samples
        .filter { sample_id, condition, fastq_1, fastq_2 -> !(fastq_1 ==~ /^https?:\/\/.*/) && !(fastq_2 ==~ /^https?:\/\/.*/) }
        .map { sample_id, condition, fastq_1, fastq_2 -> tuple(sample_id, condition, file(fastq_1, checkIfExists: true), file(fastq_2, checkIfExists: true)) }

    transcriptome_ch = Channel.value(params.transcriptome)

    FETCH_URL_FASTQ(url_samples)
    STAGE_LOCAL_FASTQ(local_samples)

    fastqs = FETCH_URL_FASTQ.out.fastqs.mix(STAGE_LOCAL_FASTQ.out.fastqs)

    FASTP(fastqs)
    SALMON_INDEX(transcriptome_ch)
    SALMON_QUANT(FASTP.out.trimmed, SALMON_INDEX.out.index)

    quant_files = SALMON_QUANT.out.quants
        .map { sample_id, condition, quant_file -> quant_file }
        .collect()

    sample_metadata = SALMON_QUANT.out.quants
        .map { sample_id, condition, quant_file -> "${sample_id}\t${condition}" }
        .collect()

    COUNT_MATRIX(quant_files, sample_metadata)
    MULTIQC(FASTP.out.json.mix(SALMON_QUANT.out.logs).collect())

}

process FETCH_URL_FASTQ {
    tag "$sample_id"
    label 'process_low'

    publishDir "${params.outdir}/fastq", mode: 'copy'

    input:
    tuple val(sample_id), val(condition), val(fastq_1), val(fastq_2)

    output:
    tuple val(sample_id), val(condition), path("${sample_id}_R1.fastq*"), path("${sample_id}_R2.fastq*"), emit: fastqs

    script:
    """
    curl -L --retry 3 --fail -o ${sample_id}_R1.fastq.gz '${fastq_1}'
    curl -L --retry 3 --fail -o ${sample_id}_R2.fastq.gz '${fastq_2}'
    """
}

process STAGE_LOCAL_FASTQ {
    tag "$sample_id"
    label 'process_low'

    publishDir "${params.outdir}/fastq", mode: 'copy'

    input:
    tuple val(sample_id), val(condition), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), val(condition), path("${sample_id}_R1.fastq*"), path("${sample_id}_R2.fastq*"), emit: fastqs

    script:
    """
    case '${reads_1}' in
      *.gz) cp ${reads_1} ${sample_id}_R1.fastq.gz ;;
      *) cp ${reads_1} ${sample_id}_R1.fastq ;;
    esac

    case '${reads_2}' in
      *.gz) cp ${reads_2} ${sample_id}_R2.fastq.gz ;;
      *) cp ${reads_2} ${sample_id}_R2.fastq ;;
    esac
    """
}

process FASTP {
    tag "$sample_id"
    label 'process_medium'

    publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
    tuple val(sample_id), val(condition), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), val(condition), path("${sample_id}.trimmed_R1.fastq.gz"), path("${sample_id}.trimmed_R2.fastq.gz"), emit: trimmed
    path "${sample_id}.fastp.html", emit: reports
    path "${sample_id}.fastp.json", emit: json

    script:
    """
    fastp \
      --in1 ${reads_1} \
      --in2 ${reads_2} \
      --out1 ${sample_id}.trimmed_R1.fastq.gz \
      --out2 ${sample_id}.trimmed_R2.fastq.gz \
      --html ${sample_id}.fastp.html \
      --json ${sample_id}.fastp.json \
      --thread ${task.cpus}
    """
}

process SALMON_INDEX {
    tag 'transcriptome'
    label 'process_medium'

    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    val transcriptome_url

    output:
    path 'salmon_index', emit: index

    script:
    def isUrl = transcriptome_url ==~ /^https?:\/\/.*/
    def isGzip = transcriptome_url.endsWith('.gz')
    def localTranscriptome = transcriptome_url.startsWith('/') ? transcriptome_url : "${projectDir}/${transcriptome_url}"
    def stageReference = isUrl && isGzip
        ? "curl -L --retry 3 --fail -o transcriptome.fa.gz '${transcriptome_url}'\ngzip -dc transcriptome.fa.gz > transcriptome.fa"
        : isUrl
            ? "curl -L --retry 3 --fail -o transcriptome.fa '${transcriptome_url}'"
            : "cp '${localTranscriptome}' transcriptome.fa"
    """
    ${stageReference}
    salmon index \
      --transcripts transcriptome.fa \
      --index salmon_index \
      --threads ${task.cpus}
    """
}

process SALMON_QUANT {
    tag "$sample_id"
    label 'process_medium'

    publishDir "${params.outdir}/salmon", mode: 'copy'

    input:
    tuple val(sample_id), val(condition), path(reads_1), path(reads_2)
    path index

    output:
    tuple val(sample_id), val(condition), path("${sample_id}.quant.sf"), emit: quants
    path "${sample_id}/logs", emit: logs

    script:
    """
    salmon quant \
      --index ${index} \
      --libType ${params.salmon_libtype} \
      --mates1 ${reads_1} \
      --mates2 ${reads_2} \
      --validateMappings \
      --minAssignedFrags ${params.salmon_min_assigned_frags} \
      --threads ${task.cpus} \
      --output ${sample_id}
    cp ${sample_id}/quant.sf ${sample_id}.quant.sf
    """
}

process COUNT_MATRIX {
    tag 'transcript-count-matrix'
    label 'process_low'

    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    path quant_files
    val metadata_rows

    output:
    path 'transcript_counts.tsv', emit: matrix
    path 'sample_metadata.tsv', emit: metadata

    script:
    def metadata = metadata_rows.join('\n')
    """
    mkdir quants
    cp *.quant.sf quants/
    printf 'sample_id\tcondition\n${metadata}\n' > sample_metadata.tsv
    make_count_matrix.py quants transcript_counts.tsv
    """
}

process MULTIQC {
    tag 'multiqc'
    label 'process_low'

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path qc_files

    output:
    path 'multiqc_report.html', emit: report

    script:
    """
    multiqc . --filename multiqc_report.html
    """
}
