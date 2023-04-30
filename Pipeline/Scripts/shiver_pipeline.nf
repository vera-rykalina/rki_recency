nextflow.enable.dsl = 2

projectDir = "/Users/vera/Projects/rki_shiver/Pipeline/"
params.trimmomatic = "${projectDir}/Scripts/bin/trimmomatic-0.36.jar"
params.gal_primers = "${projectDir}/DataShiverInit/primers_GallEtAl2012.fasta"
params.illumina_adapters = "${projectDir}/DataShiverInit/adapters_Illumina.fasta"
params.alignment = "${projectDir}/DataShiverInit/HIV1_COM_2012_genome_DNA_NoGaplessCols.fasta"
params.config = "${projectDir}/Scripts/bin/config.sh"


params.outdir = null
if (!params.outdir) {
  println "outdir: $params.outdir"
  error "Missing output directory!"
}

process RENAME_FASTQ {
  //conda "/home/beast2/anaconda3/envs/python3"
  publishDir "${params.outdir}/0_renamed_fastq", mode: "copy", overwrite: true

  input:
   path file

  output:
    path "*.fastq.gz"

  script:
   """
   regex-rename '\\d{6}_\\d{2}-\\d{5}_HIV(\\d{2}-\\d{5})_\\w{2,}_?\\d{2,}?_\\w{3}_\\w{4}_(R\\d)_\\d{3}.(\\w{5}.\\w{2})' '\\1_\\2.\\3' --rename
   """
}

process INITIALISATION {
  conda "/usr/local/Caskroom/miniconda/base/envs/shiver"
  //conda "/home/beast2/anaconda3/envs/shiver"
  publishDir "${projectDir}/${params.outdir}/1_init_dir", mode: "copy", overwrite: true

  input:
     
  output:
     path "InitDir"
  script:
  
  """
  shiver_init.sh InitDir ${params.config} ${params.alignment} ${params.illumina_adapters} ${params.gal_primers}
  """

}


process IVA_CONTIGS {
  //errorStrategy 'ignore'
  //conda "/home/beast2/anaconda3/envs/iva"
  conda "/usr/local/Caskroom/miniconda/base/envs/iva"
  publishDir "${params.outdir}/2_iva_contigs", mode: "copy", overwrite: true
  

  input:
    tuple val(id), path(reads)

  output:
    path "${id}"
    path "${id}/${id}_contigs.fasta", emit: fasta_contig

  script:
    """
    iva -f ${reads[0]} -r ${reads[1]} --pcr_primers ${params.gal_primers} --adapters ${params.illumina_adapters} --trimmomatic ${params.trimmomatic} ${id}
    mv ${id}/contigs.fasta ${id}/${id}_contigs.fasta
    """
}

process ALIGN_CONTIGS {
  errorStrategy "ignore"
  //conda "/home/beast2/anaconda3/envs/shiver"
  conda "/usr/local/Caskroom/miniconda/base/envs/shiver"
  publishDir "${params.outdir}/3_alignments", mode: "copy", overwrite: true
 
  
  input:
    path initdir
    path contigs

  output:
    path "*"
    //path "*.fasta", emit: fasta
    path "*cut*.fasta", emit: cut
    path "*raw*.fasta", emit: raw
    path "*.blast", emit: blast
    
  script:
    """
    shiver_align_contigs.sh ${initdir} ${params.config} ${contigs} ${contigs.getBaseName().split("_contigs")[0]}
    rm temp_*
    rm *_MergedHits.blast*
    """
}

process ID_BLAST {
  //conda "/home/beast2/anaconda3/envs/shiver"
  conda "/usr/local/Caskroom/miniconda/base/envs/shiver"
  publishDir "${params.outdir}/4_id_blast", mode: "copy", overwrite: true

  input:
    path blast
  
  output:
    tuple val("${blast.getBaseName()}"), path("${blast}")
    
    
  script:
    """
    echo "${blast.getBaseName()} ${blast}"
    """
}

process ID_CONTIGS {
  //conda "/home/beast2/anaconda3/envs/shiver"
  conda "/usr/local/Caskroom/miniconda/base/envs/shiver"
  publishDir "${params.outdir}/5_id_contigs", mode: "copy", overwrite: true

  input:
    path contigs
  
  output:
    tuple val("${contigs.getBaseName().split("_contigs")[0]}"), path("${contigs}")
    
    
  script:
    """
    echo "${contigs.getBaseName()} ${contigs}"
    """
}

process ID_REF {
  //conda "/home/beast2/anaconda3/envs/shiver"
  conda "/usr/local/Caskroom/miniconda/base/envs/shiver"
  publishDir "${params.outdir}/6_id_ref", mode: "copy", overwrite: true

  input:
    path refs
  
  output:
    tuple val("${refs.getBaseName().split("_cut")[0].split("_raw")[0]}"), path("${refs}")
    
  script:
   
    """
    echo "${refs}"
    """ 
   
}

process ID_FASTQ {
  conda "/usr/local/Caskroom/miniconda/base/envs/shiver"
  publishDir "${params.outdir}/7_id_fastq", mode: "copy", overwrite: true

  input:
    tuple val(id), path(fastq)
  output:
    tuple val("${id}"), path("${fastq[0]}"), path("${fastq[1]}")
  script:
   """
   echo "${id}"
   """
}

process MAP {
  //conda "/home/beast2/anaconda3/envs/shiver"
  conda "/usr/local/Caskroom/miniconda/base/envs/shiver"
  publishDir "${params.outdir}/8_mapped", mode: "copy", overwrite: true

  input:
    path initdir
    path args_list
    
  output:
    path "${args_list[1].getBaseName()}"
    
  script:
    if (args_list instanceof List) {
    """
    shiver_map_reads.sh ${initdir} ${params.config} ${args_list[0]} ${args_list[1].getBaseName()} ${args_list[1]} ${args_list[2]} ${args_list[3]} ${args_list[4]}
    """ 
    }
}


workflow {
  //fastq = channel.fromPath("${projectDir}/RawData/*.fastq.gz").collect()
  //renamed_fastq = RENAME_FASTQ(fastq.flatten())
  //id_fastq = channel.fromFilePairs("${projectDir}/${params.outdir}/1_renamed_fastq/*_R{1,2}.fastq.gz")
  fastq_pairs = channel.fromFilePairs("${projectDir}/RawData/*_R{1,2}*.fastq.gz")
  initdir = INITIALISATION()
  iva_contigs = IVA_CONTIGS(fastq_pairs)
  contigs_collected = iva_contigs.fasta_contig.collect()
  refs = ALIGN_CONTIGS(initdir, contigs_collected.flatten())
  id_blast = ID_BLAST(refs.blast)
  id_contigs = ID_CONTIGS(contigs_collected.flatten())
  //refs.cut.collect().view()
  id_refs = ID_REF(refs.cut.collect().flatten())
  id_fastq = ID_FASTQ(fastq_pairs)
  // Combine according to a key that is the first value of every first element, which is a list
  map_args = id_contigs.combine(id_blast, by:0).combine(id_refs, by:0).combine(id_fastq, by:0)
  // Get rid off id
  no_id_args = map_args.map { id, contigs, blast, ref_cut, read1, read2 -> [contigs, blast, ref_cut, read1, read2]}
  MAP(initdir,no_id_args)
  
}


