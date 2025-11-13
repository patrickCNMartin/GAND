nextflow.enable.dsl=2
//=============================================================================
// PROCESSES
//=============================================================================
process seurat_integration {
    conda "${params.scrna_r.env}"
    publishDir params.scrna_r.tmp, 
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

process build_report {
    conda "${params.scrna_r.env}"
    publishDir params.scrna_r.output, 
        mode: 'copy',
        overwrite: true

    input:
    path input_data
    path template
    
    output:
    path "GAND_seurat_analysis_report.pdf", emit: report
    path "figure_1_sample_umap_plots.pdf", emit: fig1_a
    path "figure_1_condition_umap_plots.pdf", emit: fig1_b
    path "figure_1_condition_agg_umap_plots.pdf", emit: fig1_c
    path "figure_1_integrated_umap_plots.pdf", emit: fig1_d
    
    script:
    """
    export PATH=\$HOME/bin:\$PATH
    Rscript -e '
      options(repos = c(CRAN = "https://cloud.r-project.org"))
      library(tinytex)
      library(rmarkdown)
      tinytex::install_tinytex(force = TRUE)
      rmarkdown::render("${template}",
                  params = list(input_dir = "${input_data}"),
                  output_file = "GAND_seurat_analysis_report.pdf",
                  run_pandoc = TRUE)
    '
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
    integrated_out = seurat_integration(
        scrna_directory,
        manifest,
        ref_directory)

    // Report channel - now uses the emitted annotated channel
    if (params.scrna_r.build_report == true) {
        template = Channel.fromPath(params.scrna_r.template)
        //pandoc = Channel.fromPath(params.scrna_r.pandoc)
        println "${params.scrna_r.tmp}"
        build_report(integrated_out.annotated, template)
    }
}