process PREP_FS_SEQS{
    tag "$meta.id"
    label 'process_low'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0' :
        'biocontainers/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0' }"

    input:
    tuple val(meta), val(seqs), path(mapping)
    
    output:
    tuple val (meta), path("${prefix}_fs"), emit: fs_dir
    
    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    fs_prep.py $mapping ${prefix}_fs
    """
}