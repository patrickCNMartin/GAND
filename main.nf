//=============================================================================
// WORKFLOW
//=============================================================================
nextflow.enable.dsl=2

include { gand_bw } from './workflows/gand_bw.nf'
include { gand_scR } from './workflows/gand_scR.nf'
include { dwl_data } from './workflows/dwl_data.nf'

workflow {
    // data download is a seperate process
    // each data type can be downloaded if specified in config.
    // I didn't want to have this as part of a specific analyis process/workflow
    println "Check Data to be downloaded..."
    dwl_data()
    // Running ChIP analysis from BigWig Files
    if (params.bigwig.run) {
        println "Runing BigWig Analysis..."
        gand_bw()

    }
    if (params.scrna_r.run) {
        println "Runing Single Cell Analysis..."
        gand_scR()
    }
    println "Completed Run"
}