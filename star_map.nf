#!/usr/bin/env nextflow

/*
================================================================================
                           starmap/main
================================================================================
  Nextflow pipeline for RNA-seq mapping, STAR 2.7.5a
--------------------------------------------------------------------------------
 @Repository
 https://github.com/savytskanatalia
 @Author
 Savytska Natalia, 2020. Genome Biology of Neurodegenerative Diseases, DZNE TU
--------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run starmap.nf -c mappers_3_1.config --reads "/absolute/path/to/files/*.fastq.gz" -with-docker savytskanatalia/quantefication2 
    Mandatory arguments:
      --reads [str]                  Path to input data (must be surrounded with quotes). Default: "input/*.fastq.gz" 
      --outdir [str]                 Output directory path. Default: output/       
    References                        If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta [str]                  Path to fasta reference
      --gtf [str]                    Path to GTF annotation
      --star_index [str]             Path to STAR-Index reference

      --oufiltmm [int]               Numeric, parameter for STAR --outFilterMultimapNmax. Default: 1
      --winamm [int]                 Numeric, parameter for STAR --winAnchorMultimapNmax. Default: 1

      --singleEnd                    Flag is required, when single-end sequencing data is provided
      --staropt [str]                Additional STAR flags for processing. If any additional STAR parameters need to be included, they can be with this flag (e.g. '--alignIntronMin 20 ...').
    Options:
      --threads [int]                Number of threads to use. Default: 4
      --read_length [int]            Length of the reads. Default: 100
      --mism [int]                   outFilterMismatchNmax in STAR  alignment. Default: 999
      --misml [int]                  outFilterMismatchNoverLmax in STAR  alignment. Default: 0.1
    """.stripIndent()
}



params.singleEnd = false
Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .set { fastq_files }



// fastq_files=Channel.fromFilePairs(params.reads)



out=params.outdir
rfasta=params.fasta
rgtf=Channel.value(params.gtf)
read_length=Channel.value(params.readlen)
threads=Channel.value(params.threads)
mismatch_n=Channel.value(params.mism)
mismatch_nl=Channel.value(params.misml)
outdir=params.outdir
staropt=params.staropt
ofmm=Channel.value(params.oufiltmm)
wamm=Channel.value(params.winamm)

/*
================================================================================
                                  PREPROCESSING
================================================================================
*/


// CHECKING IF WE WANT TO MAP ONLY, QUANTIFY ONLY OR MAP AND QUANTIFY



/*
 * Build STAR index
 */

if (!params.star_index && (!params.fasta && !params.gtf)) exit 1, "Roses are red, Kaiju bloods blue, You forgot your STAR index, or give .fa for build! Either specify STAR-INDEX or Fasta and GTF!"
if (!params.star_index && (!params.fasta && !params.gtf)) exit 1, "Either specify STAR-INDEX or Fasta and GTF!"
if (!params.star_index){
    process build_star_index {
        echo true
        container='savytskanatalia/quantefication2'
        tag "${fasta}-${gtf}"
        label 'process_medium'
        publishDir params.outdir, mode: 'copy'
        input:
        path read from rfasta
        path gtf from rgtf
        val x from threads
        val y from read_length
        output:
        file("star-index") into star_index
        when: !(params.star_index)
        script:
        
        """
        source activate telocal
        echo TODAY WE ARE CANCELLING THE APOCALYPSE
        echo START BUILDING STAR INDEX....
        echo FASTA: $read
        echo GTF: $gtf
        mkdir star-index
        STAR \\
            --runMode genomeGenerate \\
            --runThreadN $x \\
            --sjdbGTFfile $gtf \\
            --sjdbOverhang ${y - 1} \\
            --genomeDir star-index/ \\
            --genomeFastaFiles $read \\
            --limitGenomeGenerateRAM 168632718037
            """
    }}

ch_star_index = params.star_index ? Channel.value(file(params.star_index)).ifEmpty{exit 1, "You thought there was index, but there was none. Sorry, hurts to be wrong. STAR index not found: ${params.star_index}" } : star_index

ch_star_index = ch_star_index.dump(tag:'ch_star_index')







process star_mapping {
    container='savytskanatalia/quantefication2'
     echo true
    publishDir "${params.outdir}/mapped", mode: 'copy'
    input:
    set sampleId, file(reads) from fastq_files
    path x from Channel.value(params.star_index)
    val y from threads
    val a from mismatch_n
    val b from mismatch_nl
    val s from staropt
    val d from ofmm
    val e from wamm
    output:
    file "*.bam"
    file "*"
    script:
    """
    source activate telocal
    pwd
    echo LETS GO FISHING
    echo HI THERE ${reads[0]} $sampleId

    STAR \\
        --outSAMstrandField intronMotif \\
        --outSAMunmapped Within \\
        --runThreadN $y \\
        --genomeDir $x \\
        --readFilesIn $reads  \\
        --readFilesCommand gunzip -c \\
        --outFilterMultimapNmax $d  \\
        --winAnchorMultimapNmax $e \\
        --outSAMtype BAM Unsorted \\
        --outFilterMismatchNmax $a \\
        --outFilterMismatchNoverLmax $b \\
        --outFileNamePrefix  ${sampleId} $s
        """
}


