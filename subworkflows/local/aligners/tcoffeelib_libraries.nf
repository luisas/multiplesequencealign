include { TCOFFEE_LIBRARY as SEQUENCE_LIBRARY         } from '../../../modules/local/tcoffee/library/main.nf'
include { TCOFFEE_LIBRARY as SAP_LIBRARY              } from '../../../modules/local/tcoffee/library/main.nf'
include { TCOFFEE_LIBRARY as TMALIGN_LIBRARY          } from '../../../modules/local/tcoffee/library/main.nf'
include { TCOFFEE_MERGELIBS as MERGE_LIBRARIES        } from '../../../modules/local/tcoffee/mergelibs/main.nf'
include { TCOFFEE_LIBRARY as FOLDSEEK_LALIGN_LIBRARY  } from '../../../modules/local/tcoffee/library/main.nf'

workflow TCOFFEELIB_LIBRARIES {

    take: 
    ch_fastas       //channel: [ meta, /path/to/file.fasta ]
    ch_trees        //channel: [ meta, /path/to/file.tree ]
    ch_dependencies //channel: [ meta, /path/to/file.dependencies ]

    main:

    ch_libs     = Channel.empty()
    ch_versions = Channel.empty()

    // Extract the part of the command line that specifies the library
    ch_fastas.join(ch_trees).join(ch_dependencies).set{ ch_inputs }
    //ch_inputs.view()
    ch_inputs.map{
        meta, fasta, tree, template, dependencies ->
            [ extract_lib_parameters(meta["args_aligner"]).split().toList(), meta, fasta, tree, template, dependencies ]
    }.set{ ch_lib_params }

    // ch_lib_params.view()
    // TODO : make it so that it overlapping libraries are not computed twice
 
    // SEQUENCE LIBRARY 
    ch_lib_params.filter{ it[0].contains("proba_pair") }.multiMap{
                            metalib, meta, fastafile, tree, template, dependencies ->
                                fasta: [ meta, fastafile]
                                tree:  [ meta, tree ]
                        }.set{ ch_fastas_forlibs_sequence }
    
    SEQUENCE_LIBRARY(ch_fastas_forlibs_sequence.fasta,
                     ch_fastas_forlibs_sequence.tree,
                     [[:],[],[]],    
                     [[:],[]])
    ch_libs = ch_libs.mix(SEQUENCE_LIBRARY.out.lib)
    ch_versions = ch_versions.mix(SEQUENCE_LIBRARY.out.versions)

    // SAP LIBRARY
    ch_lib_params.filter{ it[0].contains("sap_pair") }.multiMap{
                            metalib, meta, fastafile, tree, template, dependencies ->
                                fasta: [ meta, fastafile]
                                tree:  [ meta, tree ]
                                dependencies: [ meta, template, dependencies ]
                        }.set{ ch_fastas_forlibs_sap }
    SAP_LIBRARY(ch_fastas_forlibs_sap.fasta,
                     ch_fastas_forlibs_sap.tree,
                     ch_fastas_forlibs_sap.dependencies,    
                     [[:],[]])
    ch_libs = ch_libs.mix(SAP_LIBRARY.out.lib)
    ch_versions = ch_versions.mix(SAP_LIBRARY.out.versions)

    // TMALIGN LIBRARY
    ch_lib_params.filter{ it[0].contains("TMalign_pair") }.multiMap{
                            metalib, meta, fastafile, tree, template, dependencies ->
                                fasta: [ meta, fastafile]
                                tree:  [ meta, tree ]
                                dependencies: [ meta, template, dependencies ]
                        }.set{ ch_fastas_forlibs_tmalign }
    TMALIGN_LIBRARY(ch_fastas_forlibs_tmalign.fasta,
                        ch_fastas_forlibs_tmalign.tree,
                        ch_fastas_forlibs_tmalign.dependencies,    
                        [[:],[]])
    ch_libs = ch_libs.mix(TMALIGN_LIBRARY.out.lib)
    ch_versions = ch_versions.mix(TMALIGN_LIBRARY.out.versions)

    // FOLDSEEK
    matrix_3di = Channel.fromPath("${params.matrix_3di}").map{ it -> [[:], it]}
    ch_lib_params.filter{ it[0].contains("FSlalign_pair") }.multiMap{
                            metalib, meta, fastafile, tree, template, dependencies ->
                                fasta: [ meta, fastafile]
                                tree:  [ meta, tree ]
                                dependencies: [ meta, template, dependencies ]
                        }.set{ ch_fastas_forlibs_fslalign }
    // FOLDSEEK_LALIGN_LIBRARY(ch_fastas_forlibs_fslalign.fasta,
    //                     ch_fastas_forlibs_fslalign.tree,
    //                     ch_fastas_forlibs_fslalign.dependencies,    
    //                     matrix_3di)
    // ch_libs = ch_libs.mix(FOLDSEEK_LALIGN_LIBRARY.out.lib)
    // ch_versions = ch_versions.mix(FOLDSEEK_LALIGN_LIBRARY.out.versions)

    // MERGE LIBRARIES
    MERGE_LIBRARIES(ch_libs.groupTuple())
    merged_libraries = MERGE_LIBRARIES.out.lib
    
    emit:
    ch_libs      = merged_libraries
}


def extract_lib_parameters(args){
    def matcher = args =~ /-lib\s(.*?)(?=\s-|$)/
    if(matcher.find()){
        return(matcher.group(1))
    }
}

