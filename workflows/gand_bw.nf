nextflow.enable.dsl=2
//=============================================================================
// PROCESSES
//=============================================================================
process merge_by_tag {
    publishDir "${params.tmp_bw}", mode: 'copy', overwrite: true
    input:
    tuple val(tag), path(bw_files)
    val conda_loc
    output:
    path "${tag}_merged_hg38_uniqnorm_signal.bg", emit: bedgraph
    script:
    """
    if [ "${params.dev_mode}" = "true" ]; then
        echo 'Development setup'
        source ${conda_loc}
        conda activate gand_bw
    else
        echo 'Production setup'
        export DEBUG=0
    fi
    bigWigMerge ${bw_files.join(' ')} ${tag}_merged_hg38_uniqnorm_signal.bg
    """
}

process peak_from_bdg {
    publishDir "${params.tmp_bw}", mode: 'copy', overwrite: true
    input:
    path bdg_file
    val conda_loc
    val min_l
    val cutoff
    output:
    path "${bdg_file.baseName}.narrowPeak"
    script:
    """
    if [ "${params.dev_mode}" = "true" ]; then
        echo 'Development setup'
        source ${conda_loc}
        conda activate gand_bw
    else
        echo 'Production setup'
        export DEBUG=0
    fi
    macs3 bdgpeakcall -i ${bdg_file} -o ${bdg_file.baseName}.narrowPeak -l ${min_l} -c ${cutoff}
    """
}


//=============================================================================
// WORKFLOW
//=============================================================================

workflow gand_bw {

    // Logging
    println "Input Directory BW: ${params.input_bw}"
    println "Tmp Directory BW: ${params.tmp_bw}"
    println "Export to: ${params.output_bw}"
    println "BigWig File Tags: ${params.tags}"

    // Input
    tags = Channel.fromList(params.tags)
    conda_loc = params.conda_loc

    // Get all bigwig files from input directory
    bw_files = Channel
        .fromPath("${params.input_bw}", checkIfExists: true)
        .ifEmpty { error "No BigWig files found at ${params.input_bw}" }

    // Create channel of [tag, files_for_tag]
    tagged_files = tags
    .combine(bw_files)                    // Create [tag, file] pairs
    .filter { tag, file -> 
        file.toString().toLowerCase().contains(tag.toLowerCase()) 
    }
    .groupTuple()                         // Groups into [tag: [files]]
    .map { tag, files -> 
        println "Tag: ${tag}, Files: ${files}"
        tuple(tag, files)
    }
    .ifEmpty { error "No files matched any of the tags in ${params.tags}" }
    // Run merging
    merged = merge_by_tag(tagged_files, conda_loc)

    // Run peak calling
    peak_from_bdg(merged.bedgraph, conda_loc, params.min_peak_length, params.peak_cutoff)
}