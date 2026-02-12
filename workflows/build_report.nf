nextflow.enable.dsl=2
//=============================================================================
// PROCESSES
//=============================================================================

process build_rmd_report {
    publishDir params.report.output, 
        mode: 'move',
        overwrite: true

    input:
    path template
    val data_map
    
    output:
    path "GAND_seurat_analysis_report.pdf", emit: report
    path "*.pdf", emit: all_plots
    
    script:
    def r_data_map = data_map.findAll { k, v -> v != null }
        .collect { k, v -> "${k} = \"${v}\"" }
        .join(',\n        ')
    """
    Rscript -e '
      library(tinytex)
      library(rmarkdown)  
      # Create a list of optional inputs
      data_map <- list(
        ${r_data_map}
      )
      rmarkdown::render("${template}",
                  params = list(
                    data_map = data_map
                  ),
                  output_file = "GAND_seurat_analysis_report.pdf",
                  run_pandoc = TRUE)
    '
    """
}
//=============================================================================
// WORKFLOW
//=============================================================================

workflow build_report {
    take:
    data_map
    status

    main:
    gated_status = status.map { it }

    template = gated_status.map { file(params.report.template) }

    data_files = Channel.from(data_map)
        .map { m -> m.findAll { k, v -> v != null } }

    build_rmd_report(template, data_files)
}
