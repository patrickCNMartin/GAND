nextflow.enable.dsl=2
//=============================================================================
// PROCESSES
//=============================================================================
process seurat_integration {
    conda "${params.scrna_r.env}"
    publishDir params.tmp_scrna, 
        mode: 'copy',
        pattern: '*.rds' 

    input:
    path scrna_data        
    val manifest
    path ref_dir

    output:
    path "GAND_preprocessed.rds", emit: preprocessed
    path "GAND_seurat_integrated.rds", emit: integrated
    path "GAND_seurat_annotated.rds", emit: annotated

    // Note that we are going to use conda for now
    // we could also parse renv.lock files
    script:
    """
    seurat_scRNA.R --input_dir ${scrna_data} --manifest ${manifest} --ref_dir ${ref_dir}
    """
}


//=============================================================================
// MAIN WORKFLOW
//=============================================================================
include { RMARKDOWNNOTEBOOK } from  '../modules/nf-core/rmarkdownnotebook'
workflow gand_scR {
    // Create channels properly
    scrna_directory = Channel.fromPath(params.scrna_r.input)
    println "${params.scrna_r.input}"
    ref_directory = Channel.fromPath(params.scrna_r.ref)
    println "${params.scrna_r.ref}"
    manifest = params.scrna_r.manifest
    // Call the process
    integrated = seurat_integration(
        scrna_directory,
        manifest,
        ref_directory)
    
}