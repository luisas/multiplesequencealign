process MERGE_LIBRARIES_MANUAL{

    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/foldseek:9.427df8a--pl5321hb365157_0':
    'biocontainers/foldseek:9.427df8a--pl5321hb365157_0' }"

    input:
    tuple val(meta), file(libraries)
    val(aggfunc)

    output:
    tuple val(id), file("${prefix}.library"), emit: library

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    merge_libraries.py ${libraries} ${aggfunc} ${prefix}.library 
    """
}