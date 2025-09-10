nextflow.enable.dsl=2
//=============================================================================
// PROCESSES
//=============================================================================
process merge_by_tag {
    conda "${params.bw_envs}"
    publishDir "${params.tmp_bw}", mode: 'copy', overwrite: true
    input:
    tuple val(tag), path(bw_files)
    output:
    path "${tag}_merged_hg38_uniqnorm_signal.bg", emit: bedgraph
    script:
    """
    echo "Merging BigWig Files into BedGraph"
    bigWigMerge ${bw_files.join(' ')} ${tag}_merged_hg38_uniqnorm_signal.bg
    """
}

process peak_from_bdg {
    conda "${params.bw_envs}"
    publishDir "${params.tmp_bw}", mode: 'copy', overwrite: true
    input:
    path bdg_file
    val min_l
    val cutoff
    output:
    path "${bdg_file.baseName}.narrowPeak"
    script:
    """
    echo "Calling peaks from BedGraph"
    macs3 bdgpeakcall -i ${bdg_file} -o ${bdg_file.baseName}.narrowPeak -l ${min_l} -c ${cutoff}
    """
}


//=============================================================================
// WORKFLOW
//=============================================================================

workflow gand_bw {
    // Input
    tags = Channel.fromList(params.tags)
    

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
    merged = merge_by_tag(tagged_files)

    // Run peak calling
    peak_from_bdg(merged.bedgraph, params.min_peak_length, params.peak_cutoff)
}