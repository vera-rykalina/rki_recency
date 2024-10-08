nextflow.enable.dsl = 2

// Change is required! Specify your projectDir here
projectDir = "/scratch/rykalinav/rki_recency/Pipeline"

// ********************** NEXTFLOW RUN **************************
// Activate nextflow environment
// cd to Pipeline folder and then type:
// nextflow run /Scripts/shiver_phyloscanner_tsi_pipeline.nf \
// -c Script/rki_profile.config \
// -profile rki_slurm,rki_mamba \
// --outdir Results \
// -with-report HTML/report_$(date +%T).html
// # use this paramete to build a flowchart
// -with-dag pipeline_flowchart.png (pipeline_flowchart.mmd)
// # use this parameter for an empty test run
// -stub


// Parameters for kraken2
params.krakendb = "/scratch/databases/kraken2_20230314/"

// Parameters for shiver
params.trimmomatic = "${projectDir}/Scripts/bin/trimmomatic-0.36.jar"
params.gal_primers = "${projectDir}/DataShiverInit/primers_GallEtAl2012.fasta"
params.illumina_adapters = "${projectDir}/DataShiverInit/adapters_Illumina.fasta"
params.alignment = "${projectDir}/DataShiverInit/HIV1_COM_2012_genome_DNA_NoGaplessCols.fasta"
params.config = "${projectDir}/Scripts/bin/config.sh"
params.remove_whitespace = "${projectDir}/Scripts/bin/tools/RemoveTrailingWhitespace.py"

// Parameters for phyloscanner
params.raxmlargs = "raxmlHPC-SSE3 -m GTRCAT -p 1 --no-seq-check"
params.two_refs = "${projectDir}/References/2refs_HXB2_C.BW.fasta"
params.excision_coordinates = "${projectDir}/PhyloscannerExtra/DrugResistancePositionsInHXB2.txt"
params.windows_oneline = "${projectDir}/PhyloscannerExtra/windows250_VR_norms_oneline.txt"
params.hiv_distance_normalisation = "${projectDir}/PhyloscannerExtra/HIV_DistanceNormalisationOverGenome.csv"
params.k = 15

// Parameters for HIV-PhyloTSI
params.model = "${projectDir}/Scripts/bin/Model"

params.outdir = null
if (!params.outdir) {
  println "outdir: $params.outdir"
  error "Missing output directory!"
}


log.info """
====================================================
               PhyloTSI PIPELINE
====================================================
             Author: Vera Rykalina
       Affiliation: Robert Koch Institute 
        Acknowledgement: Tanay Golubchik
              Created: 17 July 2023
           Last Updated: 23 July 2023
====================================================
         """

// KRAKEN2
process KRAKEN2_CLASSIFY {

    label "kraken2"
    conda "${projectDir}/Environments/kraken2.yml"
    publishDir "${params.outdir}/00_kraken_report/${id}", pattern: "*.txt", mode: "copy", overwrite: true

    // SLURM cluster options
    // cpus 10,  memory "150 GB", time "4h"
    // clusterOptions "--job-name=classify_${sample}"
    //tag "${id}_kraken2_classify"

    input:
        tuple val(id), path(reads)
        val (db)

    output:
        tuple val(id), path("${id}.classified.R*.fastq"),     emit: fastq
        tuple val(id), path("${id}.kraken.out.txt"),          emit: kraken_output
        tuple val(id), path("${id}.kraken.report.txt"),       emit: kraken_report
  
 script:
        """
            kraken2 \
                --threads 10 \
                --db ${db} \
                --paired \
                --classified-out ${id}.classified.R#.fastq \
                --output ${id}.kraken.out.txt \
                --report ${id}.kraken.report.txt \
                ${reads[0]} ${reads[1]}
        """
}

// SHIVER PART (including IVA and KALLISTO)
process INITIALISATION {
  //conda "/home/beast2/anaconda3/envs/shiver"
  conda "${projectDir}/Environments/shiver.yml"
  publishDir "${projectDir}/${params.outdir}/01_init_dir", mode: "copy", overwrite: true

  input:
     
  output:
     path "InitDir", emit: InitDir
     path "InitDir/ExistingRefsUngapped.fasta", emit: ExistingRefsUngapped
     path "InitDir/IndividualRefs/*.fasta", emit: IndividualRefs
  script:
  
  """
  shiver_init.sh \
    InitDir \
    ${params.config} \
    ${params.alignment} \
    ${params.illumina_adapters} \
    ${params.gal_primers}
  """  
}

process FASTQ_ID_HEADER {
  //conda "/home/beast2/anaconda3/envs/shiver"
  conda "${projectDir}/Environments/shiver.yml"
  //publishDir "${params.outdir}/N_id_fastq", mode: "copy", overwrite: true
  debug true

  input:
    tuple val(id), path(reads)
  output:
    // RWS - removed white space (change fastq header here as well)
    tuple val("${id}"), path("${reads[0].getBaseName().split("_R")[0]}_RWS_1.fastq"), path("${reads[1].getBaseName().split("_R")[0]}_RWS_2.fastq")
  
  script:
   """
   zcat ${reads[0]} |\
      awk '{if (NR%4 == 1) {print \$1 "/" \$2} else print}' |\
      sed 's/:N:.*//' > ${reads[0].getBaseName().split("_R")[0]}_1.fastq
   rm ${reads[0]}

   python \
     ${params.remove_whitespace} \
     ${reads[0].getBaseName().split("_R")[0]}_1.fastq > ${reads[0].getBaseName().split("_R")[0]}_RWS_1.fastq

   zcat ${reads[1]} |\
      awk '{if (NR%4 == 1) {print \$1 "/" \$2} else print}' |\
      sed 's/:N:.*//' > ${reads[1].getBaseName().split("_R")[0]}_2.fastq
   rm ${reads[1]}
   
   python \
     ${params.remove_whitespace} \
     ${reads[0].getBaseName().split("_R")[0]}_2.fastq > ${reads[0].getBaseName().split("_R")[0]}_RWS_2.fastq
   """
}

process KALLISTO_INDEX {
  //conda "/home/beast2/anaconda3/envs/kalisto"
  conda "${projectDir}/Environments/kallisto.yml"
  publishDir "${projectDir}/${params.outdir}/02_kallisto_idx", mode: "copy", overwrite: true

  input:
     path fasta
  
  output:
     path "*.idx"
 
  script:
  """
  kallisto index --index ExistingRefsUngapped.idx ${fasta}
  """
}

process KALLISTO_QUANT {
  //conda "/home/beast2/anaconda3/envs/kalisto"
  conda "${projectDir}/Environments/kallisto.yml"
  publishDir "${projectDir}/${params.outdir}/03_kallisto_quant", mode: "copy", overwrite: true
  debug true

  input:
     tuple path(index), val(id), path(reads)

  output:
     tuple val(id), path("${id}/${id}_abundance.tsv")
   
  script:
   """
   kallisto quant \
    -i ${index} \
    -o ${id} \
    --plaintext ${reads[0]} ${reads[1]}

    mv ${id}/abundance.tsv ${id}/${id}_abundance.tsv
  """
}

process BEST_ALIGNMENT {
  publishDir "${projectDir}/${params.outdir}/04_best_ref", mode: "copy", overwrite: true
  debug true

  input:
     path alignments
     tuple val(id), path(abundancetsv)

  output:
     tuple val(id), path("${id}_bestRef.fasta")
   
  script:
   """
   BestRef=\$(sort -k5 -g ${abundancetsv} | tail -n1 | cut -f1)
   echo "Sample ID: " ${id} "\nBest Kallisto Reference: " \${BestRef} 
   mv \${BestRef}.fasta  ${id}_bestRef.fasta
   """
}

process IVA_CONTIGS {
  label "iva"
  //conda "/home/beast2/anaconda3/envs/iva"
  conda "${projectDir}/Environments/iva.yml"
  publishDir "${params.outdir}/05_iva_contigs", mode: "copy", overwrite: true
  errorStrategy "retry"

  input:
    tuple val(id), path(reads)

  output:
    tuple val("${id}"), path("${id}/${id}_contigs.fasta")

  script:
    """
    iva \
      -f ${reads[0]} \
      -r ${reads[1]} \
      --threads 16 \
      --pcr_primers ${params.gal_primers} \
      --adapters ${params.illumina_adapters} \
      --trimmomatic ${params.trimmomatic} \
      ${id} || mkdir -p ${id} && touch ${id}/contigs.fasta


    mv ${id}/contigs.fasta ${id}/${id}_contigs.fasta
    """
}

process ALIGN_CONTIGS {
  //conda "/home/beast2/anaconda3/envs/shiver"
  conda "${projectDir}/Environments/shiver.yml"
  publishDir "${params.outdir}/06_alignments/${id}", mode: "copy", overwrite: true
  //debug true
  
  input:
    path initdir
    tuple val(id), path(contigs)

  output:
    tuple val("${id}"), path("${id}_wRefs.fasta"), path("${id}.blast")

  script:
    if ( contigs.size() > 0 ) {
    """
    shiver_align_contigs.sh \
      ${initdir} \
      ${params.config} \
      ${contigs} \
      ${id}

    rm temp_*
    rm *_MergedHits.blast*
    mv ${id}_cut_wRefs.fasta ${id}_wRefs.fasta || mv ${id}_raw_wRefs.fasta ${id}_wRefs.fasta 
    """
  } else {
    """
     printf "There is no contig for sample with ID: ${id}"
     touch ${id}.blast
     touch ${id}_wRefs.fasta
    """
  }
}

process MAP {
  //conda "/home/beast2/anaconda3/envs/shiver"
  conda "${projectDir}/Environments/shiver.yml"
  publishDir "${params.outdir}/07_mapped/${id}", mode: "copy", overwrite: true
  //debug true

  input:
    path initdir
    tuple val(id), path(kallistoRef), path(read1), path(read2), path(contigs), path(shiverRef), path(blast)
  
  output:
    tuple val("${id}"), path("${id}*ref.fasta"), path("${id}*.bam"), path("${id}*.bam.bai"), path("${id}*WithHXB2.csv")
    
  script:
    if ( shiverRef.size() > 0 ) {
    """
    shiver_map_reads.sh \
        ${initdir} \
        ${params.config} \
        ${contigs} \
        ${id} \
        ${blast} \
        ${shiverRef} \
        ${read1} \
        ${read2}
    rm temp_* 
    rm *PreDedup.bam
    
    """ 
    } else {
     """
    touch ${id}.blast
    touch ${id}_contigs.fasta

    shiver_map_reads.sh \
        ${initdir} \
        ${params.config} \
        ${id}_contigs.fasta \
        ${id} \
        ${id}.blast\
        ${kallistoRef} \
        ${read1} \
        ${read2}
    rm temp_* 
    rm *PreDedup.bam
     """
  }
}

process MAF {
 //conda "/home/beast2/anaconda3/envs/python3"
 conda "${projectDir}/Environments/python3.yml"
 publishDir "${params.outdir}/08_maf", mode: "copy", overwrite: true
 //debug true

 input:
  tuple val(id), path(ref), path(bam), path(bai), path(basefreqs)
    
 output:
  path "${id}.csv"
    
 script:
  if (basefreqs instanceof List) {
  """
  produce_maf.py ${basefreqs[1]} ${id}.csv
  """ 
  } else {
  """
  produce_maf.py ${basefreqs} ${id}.csv
  """
   }
}

process JOIN_MAFS {
  //conda "/home/beast2/anaconda3/envs/python3"
  conda "${projectDir}/Environments/python3.yml"
  publishDir "${params.outdir}/09_joined_maf", mode: "copy", overwrite: true
  //debug true

  input:
    path mafcsv
    
  output:
    path "*.csv"
    
  script:
    """
    join_mafs.py ${mafcsv}
    """ 
  }


// PHYLOSCANNER PART (including IQTREE)

process BAM_REF_ID_CSV {
  publishDir "${params.outdir}/10_ref_bam_id", mode: "copy", overwrite: true
  //debug true

  input:
    tuple val(id), path(ref), path(bam), path(bai), path(basefreqs)
    
  output:
    path "*_bam_ref_id.csv"
  
  script:
    if (bam instanceof List) {
    """
    for bamfile in *_remap.bam; do
      echo ${id}_remap.bam,${id}_remap_ref.fasta,${id}
    done > ${id}_bam_ref_id.csv
    """ 
  } else {
     """
    for bamfile in *.bam; do
      echo ${id}.bam,${id}_ref.fasta,${id}  
    done > ${id}_bam_ref_id.csv
     """
  }
}

process MAKE_TREES {
 label "phyloscanner_make_trees"
 //conda "/home/beast2/anaconda3/envs/phyloscanner"
 conda "${projectDir}/Environments/phyloscanner.yml"
 publishDir "${params.outdir}/12_phylo_aligned_reads", mode: "copy", overwrite: true 
 //debug true

 input:
  path bam_ref_id_csv, name: "phyloscanner_input.csv"
  path bam_fasta_bai_files

 output:
  path "AlignedReads/*.fasta", emit: AlignedReads
  path "Consensuses/*.fasta", emit: Consensuses
  path "ReadNames/*.csv.gz", emit: ReadsNames
  path "*.csv", emit: WindowCoordinateCorrespondence

 script:
 // remove 9470,9720,9480,9730,9490,9740 from windows
 """
  phyloscanner_make_trees.py \
       ${bam_ref_id_csv} \
       --quality-trim-ends 25 \
       --alignment-of-other-refs ${params.two_refs} \
       --pairwise-align-to B.FR.83.HXB2_LAI_IIIB_BRU.K03455 \
       --excision-ref B.FR.83.HXB2_LAI_IIIB_BRU.K03455 \
       --excision-coords \$(cat ${params.excision_coordinates}) \
       --merge-paired-reads \
       --dont-check-duplicates \
       --no-trees \
       --read-names-only \
       --merging-threshold-a 0 \
       --min-read-count 1 \
       --windows \$(cat ${params.windows_oneline}) \
       --x-raxml "${params.raxmlargs}"
 """ 
}

process IQTREE {
  label "iqtree"
  conda "${projectDir}/Environments/iqtree.yml"
  publishDir "${params.outdir}/13_iqtree_trees", mode: "copy", overwrite: true
  //debug true

 input:
  path fasta

 output:
  path "*.treefile", emit: treefile
  path "*.log", emit: iqtreelog
  //path "*.iqtree", emit: iqtree
 script:
 """
  iqtree \
     -s ${fasta} \
     -pre IQTREE_bestTree.InWindow_${fasta.getSimpleName().split("Excised_")[1]} \
     -m GTR+F+R6 \
     -nt ${task.cpus} \
     --seed 0
 """ 
}

process TREE_ANALYSIS {
 label "phyloscanner_tree_analysis"
 publishDir "${params.outdir}/14_analysed_trees", mode: "copy", overwrite: true
 //debug true

 input:
  path treefile

 output:
   path "*patStats.csv", emit: patstat_csv
   path "*blacklistReport.csv", emit: blacklist_csv
   path "*patStats.pdf", emit: patstat_pdf
   //path "*.rda", emit: rda   
   //path "*.nex", emit: nex

 script:
 """
  phyloscanner_analyse_trees.R \
    --skipSummaryGraph \
    --overwrite \
    --outputRDA \
    --outputNexusTree \
    --verbose 1 \
    --windowThreshold 0.5 \
    --allowMultiTrans \
    --directionThreshold 0.33 \
    --readCountsMatterOnZeroLengthBranches \
    --blacklistReport \
    --parsimonyBlacklistK ${params.k} \
    --ratioBlacklistThreshold 0.005 \
    --rawBlacklistThreshold 3 \
    --multifurcationThreshold 1E-5 \
    --outgroupName B.FR.83.HXB2_LAI_IIIB_BRU.K03455 \
    --normRefFileName ${params.hiv_distance_normalisation} \
    --treeFileExtension .treefile IQTREE_bestTree.InWindow "k${params.k}" "s,${params.k}" 
 """ 
}

process PHYLO_TSI {
  conda "${projectDir}/Environments/phylo_tsi.yml"
  publishDir "${params.outdir}/15_phylo_tsi", mode: "copy", overwrite: true
  debug true

  input:
    path patstat
    path maf
    
  output:
    path "phylo_tsi.csv"
  
  script:
    """
    HIVPhyloTSI.py \
      -d ${params.model} \
      -p ${patstat} \
      -m ${maf} \
      -o phylo_tsi.csv \
      --amplicons True
    """ 
}

process PRETTIFY_AND_PLOT {
  conda "${projectDir}/Environments/python3.yml"
  publishDir "${params.outdir}/15_phylo_tsi", mode: "copy", overwrite: true
  debug true

  input:
    path phylo_tsi_csv

  output:
    path "phylo_tsi_prettified.csv"
    path "tsi_barplot.png"
  
  script:
    """
    prettify_plot_tsi.py ${phylo_tsi_csv} 
    """ 
}

process MAPPING_NOTES {
  //publishDir "${params.outdir}/15_phylo_tsi", mode: "copy", overwrite: true
  debug true

  input:
    tuple val(id), path(kallistoRef), path(read1), path(read2), path(contigs), path(shiverRef), path(blast)
    
  output:
    path "${id}_mapping_notes.csv"
  
  script:
    if (contigs.size() > 0) {
    """
    echo ${id},"Mapped with IVA contigs" > ${id}_mapping_notes.csv
    """ 
  } else {
     """
    bestref=\$(grep "^>" ${kallistoRef} | sed 's/>//g') 
    echo ${id},"Mapped with reference: \${bestref}" > ${id}_mapping_notes.csv
    """
  }
}
// **************************************INPUT CHANNELS***************************************************
ch_ref_hxb2 = Channel.fromPath("${projectDir}/References/HXB2_refdata.csv", checkIfExists: true)
ch_fastq_pairs = Channel.fromFilePairs("${projectDir}/RawData/*_R{1,2}*.fastq.gz", checkIfExists: true)

workflow {
  // ****************************************KRAKEN2******************************************************
  ch_kraken2_report = KRAKEN2_CLASSIFY(ch_fastq_pairs, params.krakendb)
  // *************************************SHIVER PART*****************************************************
  ch_initdir = INITIALISATION()
  ch_kallisto_index = KALLISTO_INDEX(ch_initdir.ExistingRefsUngapped)
  ch_kallisto_index_reads = ch_kallisto_index.combine(ch_fastq_pairs)
  ch_kallisto_quant = KALLISTO_QUANT(ch_kallisto_index_reads)
  ch_bestRef = BEST_ALIGNMENT(ch_initdir.IndividualRefs, ch_kallisto_quant)
  ch_fastq_id_header = FASTQ_ID_HEADER(ch_fastq_pairs)
  ch_iva_contigs = IVA_CONTIGS(ch_fastq_pairs)
  ch_wRef = ALIGN_CONTIGS(ch_initdir.InitDir, ch_iva_contigs)
  // Combine according to a key that is the first value of every first element, which is a list
  ch_map_args = ch_bestRef.combine(ch_fastq_id_header, by:0).combine(ch_iva_contigs, by:0).combine(ch_wRef, by:0)
  ch_map_out = MAP(ch_initdir.InitDir, ch_map_args)
  // ********************************************MAF*****************************************************
  ch_maf_out = MAF(ch_map_out)
  ch_hxb2_maf = ch_ref_hxb2.combine(ch_maf_out.collect())
  ch_joined_maf = JOIN_MAFS(ch_hxb2_maf)
  // ****************************************PHYLOSCANNER PART*******************************************
  ch_phyloscanner_csv = BAM_REF_ID_CSV(ch_map_out)
  // An easy way to concatinate bam_ref_id_csv files: use collectFile() operator
  ch_bam_ref_id_all = ch_phyloscanner_csv.collectFile(name: "phloscanner_input.csv", storeDir: "${projectDir}/${params.outdir}/11_bam_ref_id_all")
  // Exclude id and csv file/files from a channel (from an output tuple)
  ch_mapped_out_no_id = ch_map_out.map {id, fasta, bam, bai, csv -> [fasta, bam, bai]}
  ch_aligned_reads = MAKE_TREES(ch_bam_ref_id_all, ch_mapped_out_no_id.flatten().collect())
  ch_aligned_reads_positions_excised = ch_aligned_reads.AlignedReads.flatten().filter(~/.*PositionsExcised.*/)
  ch_iqtree = IQTREE(ch_aligned_reads_positions_excised)
  ch_analysed_trees = TREE_ANALYSIS(ch_iqtree.treefile.collect())
  ch_phylo_tsi = PHYLO_TSI(ch_analysed_trees.patstat_csv, ch_joined_maf)
  ch_prettified_tsi = PRETTIFY_AND_PLOT(ch_phylo_tsi)
   // Mapping notes
  ch_mapping_notes = MAPPING_NOTES(ch_map_args)
  ch_mapping_notes_all = ch_mapping_notes.collectFile(name: "mapping_report.csv", storeDir: "${projectDir}/${params.outdir}/15_phylo_tsi")
}

// ***********************************************Extras*************************************************
// Extra processes (can be used later)
process RENAME_FASTQ {
  conda "${projectDir}/Environments/python3.yml"
  publishDir "${params.outdir}/A_renamed_fastq", mode: "copy", overwrite: true
  input:
   path file
  output:
    path "*.fastq.gz"
  script:
   """
   regex-rename '\\d{6}_\\d{2}-\\d{5}_HIV(\\d{2}-\\d{5})_\\w{2,}_?\\d{2,}?_\\w{3}_\\w{4}_(R\\d)_\\d{3}.(\\w{5}.\\w{2})' '\\1_\\2.\\3' --rename
   """
}

//phyloscanner_input = PHYLOSCANNER_CSV(phyloscanner_csv.collect())
process PHYLOSCANNER_CSV {
  publishDir "${params.outdir}/B_phyloscanner_input", mode: "copy", overwrite: true
  input:
    path csv
  output:
    path "phyloscanner_input.csv"
  script:
    """
    cat *_bam_ref.csv > phyloscanner_input.csv
    """ 
}