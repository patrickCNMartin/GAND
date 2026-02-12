nextflow.enable.dsl=2
import groovy.json.JsonOutput
//=============================================================================
// PROCESSES
//=============================================================================
process cell_models {
    publishDir params.modelling.tmp, 
        mode: 'copy',
        pattern: '*.{csv}' 

    input:
    path annotated
    val gene_sets
    val mut_genes
    val score_type
    val min_cells

    output:
    path "GAND_seurat_annotated.csv", emit: annotated_df
    path "mutually_exclusive_genes.csv", emit: mut_genes
    path "*_geneset_list.csv", emit: gene_sets

    // Note that we are going to use conda for now
    // we could also parse renv.lock files? 
    script:
    def gene_sets_json = JsonOutput.toJson(gene_sets)
    def mut_genes_json = JsonOutput.toJson(mut_genes)
    def score_type_json = JsonOutput.toJson(score_type)
    """
    model_scRNA.R \
        --annotated ${annotated} \
        --gene_sets '${gene_sets_json}' \
        --mut_genes '${mut_genes_json}' \
        --score_type '${score_type_json}' \
        --min_cells ${min_cells}
    """
}

//=============================================================================
// MAIN WORKFLOW
//=============================================================================

workflow gand_modelling {
    take:
    status
    main:
    gated_status = status.map { it }
    annotated = gated_status.map { file(params.modelling.annotated) }
    gene_sets = params.modelling.gene_sets
    mut_genes = params.modelling.mut_genes
    score_type = params.modelling.score_type
    min_cells = params.modelling.min_cells
    // Run the modelling
    model_out = cell_models(
            annotated,
            gene_sets,
            mut_genes,
            score_type,
            min_cells)
    status_ch = model_out.annotated_df.map { true }.first()
    emit:
    annotated_df = model_out.annotated_df
    mut_genes    = model_out.mut_genes
    gene_sets    = model_out.gene_sets
    status       = status_ch
}