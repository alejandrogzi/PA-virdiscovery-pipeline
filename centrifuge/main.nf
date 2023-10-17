#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// nextflow run main.nf --reads  --outdir 'results' --index index/p_compressed+h+v --db kraken_db


outdir = file(params.outdir)
outdir.mkdir()


process FASTQC {

    publishDir "${outdir}/fastqc/${sampleId}", mode: 'copy', overwrite: 'false'
    // tag "fastQC on $sampleId"
    // conda 'fastqc=0.11.9'
    executor 'local'
    cpus 4

    input:
        tuple val(sampleId), file(trim_paired)
    
    output:
        path("*.paired.trim*_fastqc.zip"), emit: fastqc_zip
        path("*.paired.trim*_fastqc.html")

    script:
    """
    mkdir -p ${outdir}/fastqc/${sampleId}

    fastqc ${trim_paired[0]} ${trim_paired[1]} \\
    --threads ${task.cpus} \\
    $params.fastqc_args
    """
}



process TRIMMOMATIC {

    publishDir "${outdir}/trimmomatic/${sampleId}", mode: 'copy', overwrite: 'false'
    // tag "Trimmomatic on $sampleId"
    executor 'local'
    cpus 6

    input:
        tuple val(sampleId), path(reads_ch)

    output:
        tuple val(sampleId), path("*.paired.trim*.fastq.gz"), emit: trim_paired 
        tuple val(sampleId), path("*.unpaired.trim_*.fastq.gz"), emit: trim_unpaired

    script:
    """
    mkdir -p ${outdir}/trimmomatic/${sampleId}

    trimmomatic \\
    PE \\
    -threads ${task.cpus} \\
    ${reads_ch[0]} ${reads_ch[1]} \\
    ${sampleId}.paired.trim_1.fastq.gz ${sampleId}.unpaired.trim_1.fastq.gz \\
    ${sampleId}.paired.trim_2.fastq.gz ${sampleId}.unpaired.trim_2.fastq.gz \\
    SLIDINGWINDOW:4:20 MINLEN:36
    """
}



process SPADES {
    publishDir "${outdir}/spades/${sampleId}", mode: 'copy', overwrite: 'false'
    executor 'local'
    cpus 6

    input:
        tuple val(sampleId), file(trim_paired)
        tuple val(sampleId), file(trim_unpaired)

    output:
        tuple val(sampleId), path("*scaffolds.fasta"), emit: contigs

    script:
    """
    mkdir -p ${outdir}/spades/${sampleId}

    spades.py \\
    -1 ${trim_paired[0]} \\
    -2 ${trim_paired[1]} \\
    -s ${trim_unpaired[0]} \\
    -o . \\
    -t ${task.cpus} \\
    --only-assembler
    """
}


process CENTRIFUGE {
    publishDir "${outdir}/centrifuge/${sampleId}", mode: 'copy', overwrite: 'false'
    executor 'local'
    cpus 6

    input:
        tuple val(sampleId), path(reads_ch)

    script:
    """
    mkdir -p ${outdir}/centrifuge/${sampleId}

    centrifuge \\
    -x ${params.index} \\
    -1 ${reads_ch[0]} \\
    -2 ${reads_ch[1]} \\
    -S ${outdir}/centrifuge/${sampleId}/${sampleId}.centrifuge.report.large.tsv \\
    --report-file ${outdir}/centrifuge/${sampleId}/${sampleId}.centrifuge.report.tsv \\
    --threads ${task.cpus}
    """
}


process KRAKEN2 {
    publishDir "${outdir}/kraken2/${sampleId}", mode: 'copy', overwrite: 'false'
    executor 'local'
    cpus 6

    input:
        tuple val(sampleId), path(reads_ch)

    output:
        tuple val(sampleId), path("*.kraken2.k2report"), emit: kraken2_report

    script:
    """
    mkdir -p ${outdir}/kraken2/${sampleId}

    kraken2 \\
    --db ${params.db} \\
    --paired \\
    --threads ${task.cpus} \\
    --report ${sampleId}.kraken2.k2report \\
    ${reads_ch[0]} ${reads_ch[1]}
    """
}


process BRACKEN {
    publishDir "${outdir}/bracken/${sampleId}", mode: 'copy', overwrite: 'false'
    executor 'local'
    cpus 6

    input:
        tuple val(sampleId), path(kraken2_report)

    output:
        tuple val(sampleId), path("*.bracken.report"), emit: bracken_report

    script:
    """
    mkdir -p ${outdir}/kraken2/${sampleId}

    bracken \\
    -d ${params.db} \\
    -i ${kraken2_report} \\
    -o ${sampleId}.bracken.report \\
    -r 100 \\
    -l G \\
    -t 10
    """
}


workflow {

    Channel.fromFilePairs("${params.reads}/*_R{1,2}_*.fastq.gz", type:'file').set{reads_ch}
    TRIMMOMATIC(reads_ch)
    // FASTQC(trim_paired)
    // SPADES(TRIMMOMATIC.out.trim_paired, TRIMMOMATIC.out.trim_unpaired)
    CENTRIFUGE(TRIMMOMATIC.out.trim_paired)
    KRAKEN2(TRIMMOMATIC.out.trim_paired)
    BRACKEN(KRAKEN2.out.kraken2_report)
}