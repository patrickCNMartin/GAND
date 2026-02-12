nextflow.enable.dsl=2
//=============================================================================
// PROCESSES
//=============================================================================
process seurat_integration {
    publishDir params.integration.tmp, 
        mode: 'copy',
        pattern: '*.rds' 

    input:
    path scrna_data        
    val manifest
    path ref_dir
    val number_pcs
    val min_features
    val max_features
    val percent_mt
    val n_var_features
    val cluster_resolution
    val integration_tag
    val integration_method

    output:
    path "GAND_preprocessed.rds", emit: preprocessed
    path "GAND_seurat_integrated.rds", emit: integrated
    path "GAND_seurat_annotated.rds", emit: annotated


    // Note that we are going to use conda for now
    // we could also parse renv.lock files
    script:
    """
    seurat_scRNA.R --input_dir ${scrna_data} \
        --manifest ${manifest} \
        --ref_dir ${ref_dir} \
        --number_pcs ${number_pcs} \
        --min_features ${min_features} \
        --max_features ${max_features} \
        --percent_mt ${percent_mt} \
        --n_var_features ${n_var_features} \
        --cluster_resolution ${cluster_resolution} \
        --integration_tag ${integration_tag} \
        --integration_method ${integration_method}
    """
}

//=============================================================================
// MAIN WORKFLOW
//=============================================================================
workflow gand_integration {
    take:
    status

    main:
    // It's cleaner to use the process name .out rather than variable assignment
    seurat_integration(
        status.map { file(params.integration.input) },
        params.integration.manifest,
        status.map { file(params.integration.ref) },
        params.integration.number_pcs,
        params.integration.min_features,
        params.integration.max_features,
        params.integration.percent_mt,
        params.integration.n_var_features,
        params.integration.cluster_resolution,
        params.integration.integration_tag,
        params.integration.integration_method
    )
        
    status_ch = seurat_integration.out.annotated
                   .map { true }
                   .first()

    emit:
    preprocessed = seurat_integration.out.preprocessed
    integrated   = seurat_integration.out.integrated
    annotated    = seurat_integration.out.annotated
    status       = status_ch
}