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
if (params.profile) {
    exit 1, "--profile is WRONG use -profile" }

// warnings
def folder = new File(params.output)
if ( folder.exists() ) { 
    println ""
    println "\033[0;33mWARNING: Output folder already exists. Results might be overwritten! You can adjust the output folder via [--output]\033[0m"
}

Set modes = ['paired', 'single']
if ( ! (params.mode in modes) ) {
    exit 1, "Unknown mode. Choose from " + modes
}



process FASTP {
  label "fastp"
  conda "${projectDir}/Environments/fastp.yml"
  publishDir "${params.outdir}/01_trimmed", mode: "copy", overwrite: true


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
        --json ${id}.fastp.json \
        --html ${id}.fastp.html \
        --adapter_fasta ${params.illumina_adapters} \
        --low_complexity_filter \
        --overrepresentation_analysis \
        --thread $task.cpus \
        ${params.fastp_additional_parameters}
    """
}


// **************************************INPUT CHANNELS***************************************************
ch_ref_hxb2 = Channel.fromPath("${projectDir}/References/HXB2_refdata.csv", checkIfExists: true)


if (params.mode == 'paired') {
        ch_input_fastq = Channel
        .fromFilePairs( "${projectDir}/RawData/*_R{1,2}*.fastq.gz", checkIfExists: true ).view()
} else { ch_input_fastq = Channel
        .fromPath( "${projectDir}/RawData/*.fastq.gz", checkIfExists: true )
        .map { file -> [file.simpleName, [file]]}
}
