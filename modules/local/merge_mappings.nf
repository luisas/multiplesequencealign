process MERGE_MAPPINGS {

  label 'process_low'
  tag "$meta.id"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0' :
        'biocontainers/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0' }"

  input:
  tuple val(meta), file(files)

  output:
  tuple val(meta), file("${prefix}.mapping"), emit: mapping

  script:
  def args = task.ext.args ?: ''
  prefix = task.ext.prefix ?: "${meta.id}"
  """
  for file in $files; do cat \$file >>  "${prefix}.mapping"; done
  """
}