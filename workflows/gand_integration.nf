nextflow.enable.dsl=2
//=============================================================================
// PROCESSES
//=============================================================================
process seurat_integration {
    conda "${params.integration.env}"
    publishDir params.integration.tmp, 
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
    seurat_scRNA.R --input_dir ${scrna_data}\
        --manifest ${manifest} \
        --ref_dir ${ref_dir} \ 
    """
}

//=============================================================================
// MAIN WORKFLOW
//=============================================================================

workflow gand_integration {
    take:
    status
    main:
    gated_status = status.map { it }
    scrna_dir = gated_status.map { file(params.integration.input) }
    ref_dir   = gated_status.map { file(params.integration.ref) }
  
    integrated = seurat_integration(scrna_dir, params.integration.manifest, ref_dir)
        
    status_ch = integrated.annotated
                   .map { true }
                   .first()

    emit:
        preprocessed = integrated.preprocessed
        integrated   = integrated.integrated
        annotated    = integrated.annotated
        status       = status_ch
}