nextflow.enable.dsl = 2

// Change is required! Specify your projectDir here
projectDir = "/home/rykalinav/scratch/rki_shiver/Pipeline"
//projectDir = "/home/beast2/rki_shiver/Pipeline"


// ********************** NEXTFLOW RUN **************************
// Activate nextflow environment
// cd to Pipeline folder and then type:
// nextflow run /Scripts/shiver_phyloscanner_tsi_pipeline.nf \
// -c Script/rki_profile.config \
// -profile rki_slurm,rki_mamba \
// --outdir Results \
// -with-report logs/test.$(date +%T).report.html
// # use this paramete to build a flowchart
// -with-dag flowchart.png
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
//params.extra_args = "-Q1 25 -A ${projectDir}/References/2refs_HXB2_C.BW.fasta -2 B.FR.83.HXB2_LAI_IIIB_BRU.K03455 -XR B.FR.83.HXB2_LAI_IIIB_BRU.K03455 -XC 823,824,825,892,893,894,907,908,909,1012,1013,1014,1156,1157,1158,1384,1385,1386,1444,1445,1446,1930,1931,1932,1957,1958,1959,2014,2015,2016,2023,2024,2025,2080,2081,2082,2134,2135,2136,2191,2192,2193,2280,2281,2282,2283,2284,2285,2298,2299,2300,2310,2311,2312,2316,2317,2318,2319,2320,2321,2322,2323,2324,2340,2341,2342,2346,2347,2348,2349,2350,2351,2352,2353,2354,2355,2356,2357,2358,2359,2360,2373,2374,2375,2379,2380,2381,2385,2386,2387,2388,2389,2390,2391,2392,2393,2394,2395,2396,2400,2401,2402,2409,2410,2411,2412,2413,2414,2415,2416,2417,2424,2425,2426,2430,2431,2432,2436,2437,2438,2439,2440,2441,2442,2443,2444,2457,2458,2459,2460,2461,2462,2463,2464,2465,2469,2470,2471,2472,2473,2474,2478,2479,2480,2481,2482,2483,2496,2497,2498,2499,2500,2501,2502,2503,2504,2505,2506,2507,2514,2515,2516,2517,2518,2519,2520,2521,2522,2526,2527,2528,2529,2530,2531,2535,2536,2537,2670,2671,2672,2679,2680,2681,2703,2704,2705,2709,2710,2711,2733,2734,2735,2742,2743,2744,2748,2749,2750,2751,2752,2753,2754,2755,2756,2757,2758,2759,2769,2770,2771,2772,2773,2774,2778,2779,2780,2811,2812,2813,2814,2815,2816,2817,2818,2819,2823,2824,2825,2841,2842,2843,2847,2848,2849,2850,2851,2852,2856,2857,2858,2865,2866,2867,2871,2872,2873,2892,2893,2894,2895,2896,2897,2901,2902,2903,2904,2905,2906,2952,2953,2954,2961,2962,2963,3000,3001,3002,3015,3016,3017,3018,3019,3020,3030,3031,3032,3042,3043,3044,3084,3085,3086,3090,3091,3092,3099,3100,3101,3111,3112,3113,3117,3118,3119,3135,3136,3137,3171,3172,3173,3177,3178,3179,3180,3181,3182,3189,3190,3191,3192,3193,3194,3204,3205,3206,3210,3211,3212,3222,3223,3224,3228,3229,3230,3237,3238,3239,3246,3247,3248,3249,3250,3251,3255,3256,3257,3261,3262,3263,3396,3397,3398,3501,3502,3503,3546,3547,3548,3705,3706,3707,4425,4426,4427,4449,4450,4451,4503,4504,4505,4518,4519,4520,4590,4591,4592,4641,4642,4643,4647,4648,4649,4656,4657,4658,4668,4669,4670,4671,4672,4673,4692,4693,4694,4722,4723,4724,4782,4783,4784,4974,4975,4976,5016,5017,5018,5067,5068,5069,7863,7864,7865,7866,7867,7868,7869,7870,7871,7872,7873,7874,7875,7876,7877,7881,7882,7883,7884,7885,7886"
params.windows_oneline = "${projectDir}/PhyloscannerExtra/windows250_VR_norms_oneline.txt"
params.hiv_distance_normalisation = "${projectDir}/PhyloscannerExtra/HIV_DistanceNormalisationOverGenome.csv"
params.k = 15

// Parameters for PhyloTSI
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
           Last Updated: 21 July 2023
====================================================
         """

// KRAKEN2
process KRAKEN2_CLASSIFY {

    label "kraken2"
    conda "${projectDir}/Environments/kraken2.yml"
    publishDir "${params.outdir}/00_kraken_report/${id}", pattern: "*.txt"

    // SLURM cluster options
    // cpus 10,  memory "150 GB", time "4h"
    // clusterOptions "--job-name=classify_${sample}"
    tag "${id}_kraken2_classify"

    input:
        tuple val(id), path(reads)

    output:
        tuple val(id), path("${id}.classified.R*.fastq"),     emit: fastq
        tuple val(id), path("${id}.kraken.out.txt"),          emit: kraken_output
        tuple val(id), path("${id}.kraken.report.txt"),       emit: kraken_report

 script:
        """
            kraken2 \
                --threads 10 \
                --db ${params.krakendb} \
                --paired \
                --classified-out ${sample}.classified.R#.fastq \
                --output ${sample}.kraken.out.txt \
                --report ${sample}.kraken.report.txt \
                ${reads[0]} ${reads[1]}
        """

    stub:
        """
            touch ${id}.classified.R_{1,2}.fastq ${id}.kraken.out.txt ${id}.kraken.report.txt
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
  publishDir "${projectDir}/${params.outdir}/04_bets_ref", mode: "copy", overwrite: true
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
      ${id}


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
    def ivacontig = contigs instanceof Path
    if ( ivacontig ) {
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
     printf "There is no conting for sample with ID: ${id}"
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
    def wRef = shiverRef instanceof Path
    if ( wRef ) {
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
        ${bestRef} \
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
  produce_maf.py ${basefreqs} ${id}.csv
  """ 
  } else {
  """
  produce_maf.py ${basefrecs} ${id}.csv
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
  debug true

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
 // --read-names-2 (Does not exists! - not used by me in the command; ask Tanya)
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
  publishDir "${params.outdir}/11_iqtree_trees", mode: "copy", overwrite: true
  //debug true

 input:
  path fasta

 output:
  path "*.treefile", emit: treefile
  path "*.iqtree", emit: iqtree
  path "*.log", emit: iqtreelog

 script:
 """
  iqtree \
     -s ${fasta} \
     -pre IQTREE_bestTree.InWindow_${fasta.getSimpleName().split("Excised_")[1]} \
     -m GTR+F+R6 \
     -nt 16
 """ 
}

process TREE_ANALYSIS {
 label "phyloscanner_tree_analysis"
 publishDir "${params.outdir}/12_analysed_trees", mode: "copy", overwrite: true
 //debug true

 input:
  path treefile

 output:
   path "*patStats.csv", emit: patstat_csv
   path "*blacklistReport.csv", emit: blacklist_csv
   path "*patStats.pdf", emit: patstat_pdf
   path "*.nex", emit: nex
   path "*.rda", emit: rda

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
    --readCountsMatterOnZeroLenghBranches \
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
  publishDir "${params.outdir}/13_phylo_tsi", mode: "copy", overwrite: true
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
      -o phylo_tsi.csv
    """ 
}

process PRETTIFY_AND_PLOT {
  conda "${projectDir}/Environments/python3.yml"
  publishDir "${params.outdir}/13_phylo_tsi", mode: "copy", overwrite: true
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

// **************************************INPUT CHANNELS***************************************************
ch_ref_hxb2 = Channel.fromPath("${projectDir}/References/HXB2_refdata.csv", checkIfExists: true)
ch_fastq_pairs = Channel.fromFilePairs("${projectDir}/RawData/*_R{1,2}*.fastq.gz", checkIfExists: true)

workflow {
  // ****************************************KRAKEN2******************************************************
  ch_kraken2_report = KRAKEN2_CLASSIFY(ch_fastq_pairs)
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
  //ch_analysided_trees = TREE_ANALYSIS(ch_iqtree.treefile.collect())
  //ch_phylo_tsi = PHYLO_TSI(ch_analysided_trees.patstat_csv, ch_joined_maf)
  //ch_prettified_tsi = PRETIFFY_AND_PLOT(ch_phylo_tsi)
}

// ***********************************************Extras*************************************************
// Extra process (can be used later)
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