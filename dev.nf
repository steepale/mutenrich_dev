#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/mutenrich
========================================================================================
 nf-core/mutenrich Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/mutenrich
----------------------------------------------------------------------------------------
*/

/*
    Steps in a typical somatic mutation enrichment analysis run:
    1) Perform MuSiC2 analysis
        a) Prepare files for analysis
            1) A Region Of Interest (ROI) file (example in roi.Rmd)
            2) A bamlist file (see example in bamlist.Rmd)
            3) Convert a VCF file to MAF format 
        b) Calculate per-base coverage across genome
        c) Calculate per-ROI coverage across genome
        d) Measure overall and per-ROI mutation rate
        e) Perform enrichment analysis for mutations in ROI's.
    2) Perform oncodriveCLUST analysis
    3) Perform MUFFINN analysis
    
*/

// Docker containers
container__roi = "steepale/20190819_roirscript:1.0"
container__bamlist = "steepale/20190819_roirscript:1.0"
container__maf = "steepale/vep_galgal5:release_92"

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/mutenrich \
    --ann genome.gtf \
    --bamlist bamlist.txt \
    --outdir = 'results'
    
    Mandatory arguments:
      --ann                         Path to the annotation file
      --bamlist                     Path to a tab seperated bamlist file

    Options:

    Other options:
      --outdir                      The output directory where the results will be saved
      --roi_script                  Path to the Region Of Interest R file to parse regions of interest from the annotation file
      -resume                       Resume the pipeline from prior runs
      -with-docker                  Launch docker containers
      --help                        Dsiplay this help message

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help || params.manifest == null){
    // Invoke the function above to display the help message
    helpMessage()
    // Exit out and do not run
    exit 0
}

/*
 * Defines some parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 *
 * params.ann = "data/Gallus_gallus.Gallus_gallus-5.0.87.gtf"
 * 
 */

// Hard coded parameters
params.workdir = '.'
params.roi_script = 'scripts/roi.R'
params.outdir = '.'

/* 
 * prints user convenience 
 */
println "M U T E N R I C H      D E V    "
println "================================="
println "Roi Script         : ${params.roi_script}"
println "annotation file    : ${params.ann}"
println ""
println ""
/*
 * get a file object for the given param string
 */
ann_file = file(params.ann)
roi_script = file(params.roi_script)

/*
 * STEP 1. Perform MuSiC2 Analysis
 * STEP 1.a. Prep input files
 * STEP 1.a.1. Create Region Of Interest (ROI) file
 */


// Create a 'sample_pairs' channel to emit tuples with 3 elements:
// the sample ID, the normal bam file, and the tumor bam file
bam_pairs = Channel.fromFilePairs(params.bamlist, flat: true)

/*
 * STEP 1.a. Create the MAF file from VCF files
 */
process maf {
    container "${container__maf}"
    label 'maf'
    publishDir "${params.workdir}/data", mode: 'copy'

    input:
    file 
    set sample_id, file(normal_bam), file(tumor_bam) from bam_pairs
    
    output:
    file 'roi.txt' into roi_file

    script:
    """
    Rscript ${r_roi} \
    --ann ${ann} \
    --output roi.txt \
    --roi exon
    """
}
















// A function for a header to the help message
def nfcoreHeader(){
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """    ${c_dim}----------------------------------------------------${c_reset}
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/mutenrich v${workflow.manifest.version}${c_reset}
    ${c_dim}----------------------------------------------------${c_reset}
    """.stripIndent()
}


