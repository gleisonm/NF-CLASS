/*
 * pipeline input parameters
 */
params.reads = "$projectDir/data/ggal/gut_{1,2}.fq"
params.transcriptome_file = "$projectDir/data/ggal/transcriptome.fa"
params.multiqc = "$projectDir/multiqc"
params.outdir = "results"
log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    transcriptome: ${params.transcriptome_file}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent()

/*
 * define the `index` process that creates a binary index
 * given the transcriptome file
 */
process INDEX_SL {
    input:
    path transcriptome

    output:
    path 'salmon_index'

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i salmon_index
    """
}

process INDEX_BT2 {
    container 'quay.io/biocontainers/bowtie2:2.5.1--py38he00c5e5_2'
    input:
    path transcriptome

    output:
    path "bowtie2"

    script:
    """
    mkdir bowtie2
    bowtie2-build -f $transcriptome bowtie2/human_idx
    """
}

process BT2_ALIGN { 
    tag "Bowtie2 on $sample_id"
    publishDir params.outdir, mode:'copy'
    container 'quay.io/biocontainers/bowtie2:2.5.1--py38he00c5e5_2'

    input:
    path bowtie2_index
    tuple val(sample_id), path(reads)

    output:
    path "$sample_id"

    script:
    """
    bowtie2 -t -p $task.cpus -x $bowtie2_index/human_idx -1 ${reads[0]} -2 ${reads[1]} -S $sample_id 2> ${sample_id}.log
    """
}

process QUANTIFICATION {
    tag "Salmon on $sample_id"
    publishDir params.outdir, mode:'copy'

    input:
    path salmon_index
    tuple val(sample_id), path(reads)

    output:
    path "$sample_id"

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $salmon_index -1 ${reads[0]} -2 ${reads[1]} -o $sample_id
    """
}

process FASTQC {
    tag "FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

process MULTIQC {
    publishDir params.outdir, mode:'copy'

    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}


workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    index_bt2_ch = INDEX_BT2(params.transcriptome_file)
    algn = BT2_ALIGN(index_bt2_ch, read_pairs_ch)
    index_sl_ch = INDEX_SL(params.transcriptome_file)
    quant_ch = QUANTIFICATION(index_sl_ch, read_pairs_ch)
    fastqc_ch = FASTQC(read_pairs_ch)
    MULTIQC(quant_ch.mix(fastqc_ch).collect())
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
