nextflow.enable.dsl = 2

// Change is required! Specify your projectDir here
projectDir = "/scratch/rykalinav/rki_recency/Pipeline"


// Parameters for kraken2
params.krakendb = "/scratch/databases/kraken2_20230314/"

// Parameters for shiver
params.trimmomatic = "${projectDir}/Scripts/bin/trimmomatic-0.36.jar"
params.gal_primers = "${projectDir}/DataShiverInit/primers_GallEtAl2012.fasta"
params.illumina_adapters = "${projectDir}/DataShiverInit/adapters_Illumina.fasta"
params.alignment = "${projectDir}/DataShiverInit/HIV1_COM_2012_genome_DNA_NoGaplessCols.fasta"
params.config = "${projectDir}/Scripts/bin/config.sh"
params.remove_whitespace = "${projectDir}/Scripts/bin/tools/RemoveTrailingWhitespace.py"

log.info """
====================================================
                  TSI PIPELINE
====================================================
             Author: Vera Rykalina
       Affiliation: Robert Koch Institute 
        Acknowledgement: Tanay Golubchik
              Created: 17 July 2023
           Last Updated: 5 April 2024
====================================================
         """

// error codes
params.profile = null
if (params.profile) {
    exit 1, "--profile is WRONG use -profile" }

params.outdir = null
if (!params.outdir) {
  println "outdir: $params.outdir"
  error "Missing output directory!"
}


Set modes = ['paired', 'single']
if ( ! (params.mode in modes) ) {
    exit 1, "Unknown mode. Choose from " + modes
}

process FASTQC {
  conda "${projectDir}/Environments/fastqc.yml"
  publishDir "${params.outdir}/01_raw_fastqc", mode: "copy", overwrite: true
  debug true

  input:
    tuple val(id), path(reads)
  output:
    path "${id}*_fastqc.html", emit: html
    path "${id}*_fastqc.zip",  emit: zipped
  script:
    """
    fastqc ${reads}
    """

}

process FASTP {
  label "fastp"
  conda "${projectDir}/Environments/fastp.yml"
  publishDir "${params.outdir}/02_trimmed", mode: "copy", overwrite: true
  debug true

  input:
    tuple val(id), path(reads)

  output:
    tuple val(id), path("${id}.fastp.R{1,2}.fastq.gz"), emit: reads
    tuple val(id), path("${id}.fastp.json"),            emit: json
    tuple val(id), path("${id}.fastp.html"),            emit: html

 script:
    set_paired_reads = params.mode == 'single' ? '' : "--in2 ${reads[1]} --out2 ${id}.fastp.R2.fastq.gz --unpaired1 ${id}.SE.R1.fastq.gz --unpaired2 ${id}.SE.R2.fastq.gz"
    """
    fastp \
        --in1 ${reads[0]} \
        --out1 ${id}.fastp.R1.fastq.gz \
        ${set_paired_reads} \
        --adapter_fasta ${params.illumina_adapters} \
        --json ${id}.fastp.json \
        --html ${id}.fastp.html \
        --low_complexity_filter \
        --overrepresentation_analysis \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --thread ${task.cpus}

    """
}


// **************************************INPUT CHANNELS***************************************************
ch_ref_hxb2 = Channel.fromPath("${projectDir}/References/HXB2_refdata.csv", checkIfExists: true)

if (params.mode == 'paired') {
        ch_input_fastq = Channel
        .fromFilePairs( "${projectDir}/RawData/*_R{1,2}*.fastq.gz", checkIfExists: true )
} else { ch_input_fastq = Channel
        .fromPath( "${projectDir}/RawData/*.fastq.gz", checkIfExists: true )
        .map { file -> [file.simpleName, [file]]}
        .view()
}



workflow {
    ch_raw_reads = FASTQC ( ch_input_fastq )
    //ch_fastp = FASTP( ch_input_fastq)

}