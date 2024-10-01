process STRUCTURE_TO_3DI{
    tag "$meta.id"
    label 'process_low'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/foldseek:9.427df8a--pl5321hb365157_0':
    'biocontainers/foldseek:9.427df8a--pl5321hb365157_0' }"
    
    
    input:
    tuple val(meta), path (structures)
    
    output:
    tuple val(meta), path ("*3di.out"), emit: mapping

    
    script:
    """
    # Convert structures to 3di
    
    for structure in *.pdb; do
        st_id=\$(echo \$structure | cut -d'.' -f1)
        foldseek structureto3didescriptor \$structure \${st_id}_3di
        cut -f1,2,3 \${st_id}_3di > \${st_id}_3di_temp.out
        echo -n `cut -f1 \${st_id}_3di_temp.out | cut -f1 -d' '` > \${st_id}_3di.out
        echo -e -n ' \t ' >> \${st_id}_3di.out
        echo -n `cut -f2 \${st_id}_3di_temp.out` >> \${st_id}_3di.out
        echo -e -n ' \t ' >> \${st_id}_3di.out
        echo `cut -f3 \${st_id}_3di_temp.out` >> \${st_id}_3di.out
    done
    """
}