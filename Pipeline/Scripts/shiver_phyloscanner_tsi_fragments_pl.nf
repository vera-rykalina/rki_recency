
nextflow.enable.dsl = 2

// Change is required! Specify your projectDir here
projectDir = "/home/rykalinav/scratch/rki_recency/Pipeline"


// ********************** NEXTFLOW RUN **************************
// Activate nextflow environment
// cd to Pipeline folder and then type:
// nextflow run /Scripts/shiver_phyloscanner_tsi_fragments_pl.nf \
// -c Script/rki_profile.config \
// -profile rki_slurm,rki_mamba \
// --outdir ResultsFragments \
// -with-report HTML/report_$(date +%T).html
// # use this paramete to build a flowchart
// -with-dag pipeline_flowchart.png (pipeline_flowchart.mmd)
// # use this parameter for an empty test run
// -stub

params.outdir = null
if (!params.outdir) {
  println "outdir: $params.outdir"
  error "Missing output directory!"
}


log.info """
====================================================
            RECENCY PIPELINE for FRAGMENTS
====================================================
             Author: Vera Rykalina
       Affiliation: Robert Koch Institute 
              Created: 31 July 2023
           Last Updated: 31 July 2023
====================================================
         """



// **************************************INPUT CHANNELS***************************************************
ch_ref_hxb2 = Channel.fromPath("${projectDir}/References/HXB2_refdata.csv", checkIfExists: true)
ch_prrt_int_pairs = Channel.fromFilePairs("${projectDir}/PRRT_INT_FASTA/*_{INT,PRRT}*.fasta", checkIfExists: true)

workflow {
ch_prrt_int_pairs.view()
ch_consensus_from_prrt_int = MERGE_FRAGMENTS(ch_prrt_int_pairs)
}


process MERGE_FRAGMENTS {
  //conda "/home/beast2/anaconda3/envs/python3"
  conda "${projectDir}/Environments/python3.yml"
  publishDir "${projectDir}/${params.outdir}/01_merged_prrt_int", mode: "copy", overwrite: true

  input:
      tuple val(id), path(fastas) 
     
  output:
     path "${id}_prrt_int.fasta"
  
  script:
  """
    join_prrt_int.py \
       -p ${fastas[1]} \
       -i ${fastas[0]} \
       -o ${id}_prrt_int.fasta
  """  
}