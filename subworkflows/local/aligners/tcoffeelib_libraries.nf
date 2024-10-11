include { TCOFFEE_LIBRARY as SEQUENCE_LIBRARY         } from '../../../modules/local/tcoffee/library/main.nf'
include { TCOFFEE_LIBRARY as SAP_LIBRARY              } from '../../../modules/local/tcoffee/library/main.nf'
include { TCOFFEE_LIBRARY as TMALIGN_LIBRARY          } from '../../../modules/local/tcoffee/library/main.nf'
include { TCOFFEE_MERGELIBS as MERGE_LIBRARIES        } from '../../../modules/local/tcoffee/mergelibs/main.nf'
include { TCOFFEE_LIBRARY as FOLDSEEK_LALIGN_LIBRARY  } from '../../../modules/local/tcoffee/library/main.nf'
include { TCOFFEE_LIBRARY as FOLDSEEK_FULL_LIBRARY    } from '../../../modules/local/tcoffee/library/main.nf'

include { STRUCTURE_TO_3DI as STRUCTURE_TO_3DI_LOCAL  } from '../../../modules/local/structure_to_3di.nf'
include { STRUCTURE_TO_3DI as STRUCTURE_TO_3DI_FULL   } from '../../../modules/local/structure_to_3di.nf'
include { MERGE_MAPPINGS as MERGE_MAPPINGS_FULL       } from '../../../modules/local/merge_mappings.nf'
include { MERGE_MAPPINGS as MERGE_MAPPINGS_LOCAL      } from '../../../modules/local/merge_mappings.nf'
include { PREP_FS_SEQS                                } from '../../../modules/local/prep_fs_seq.nf'
include { MERGE_LIBRARIES_MANUAL                      } from '../../../modules/local/merge_libraries_manual.nf'
include { ENCODE_FASTA                                } from '../../../modules/local/encode_fasta.nf'
include { CHANGE_LIBRARY_SEQUENCES                    } from '../../../modules/local/change_library_sequences.nf'

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
 
    // -----------------------------------------------------------------
    // SEQUENCE LIBRARY 
    // -----------------------------------------------------------------
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

    //-----------------------------------------------------------------
    // SAP LIBRARY
    //-----------------------------------------------------------------
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

    //-----------------------------------------------------------------
    // TMALIGN LIBRARY
    //-----------------------------------------------------------------
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

    //-----------------------------------------------------------------
    // FOLDSEEK LOCAL LALIGN
    //-----------------------------------------------------------------
    matrix_3di = Channel.fromPath("${params.matrix_3di}").map{ it -> [[:], it]}.collect()
    ch_lib_params.filter{ it[0].contains("FSlalign_pair") }.multiMap{
                            metalib, meta, fastafile, tree, template, dependencies ->
                                fasta: [ meta, fastafile]
                                tree:  [ meta, tree ]
                                dependencies: [ meta, template, dependencies ]
                        }.set{ ch_fastas_forlibs_fslalign }
    // Convert each structure to 3di
    STRUCTURE_TO_3DI_LOCAL(ch_fastas_forlibs_fslalign.dependencies.map{ meta, template, structures -> [ meta, structures ]})
    // Merge all the 3di mappings into one file (id, fasta sequence, 3di sequence)
    MERGE_MAPPINGS_LOCAL(STRUCTURE_TO_3DI_LOCAL.out.mapping)
    ch_fastas_forlibs_fslalign.fasta.combine(
        MERGE_MAPPINGS_LOCAL.out.mapping, by:0
    ).set{ mappings_3di }
    // Prepare the sequences for the library
    PREP_FS_SEQS(mappings_3di)

    // Prepare channels in the right order for the library computation 
    ch_fastas_forlibs_fslalign.fasta.combine(
        ch_fastas_forlibs_fslalign.tree, by : 0
    ).combine(
        ch_fastas_forlibs_fslalign.dependencies, by : 0
    ).map{
        meta, fasta, tree, template, dependencies -> 
            [ meta, fasta, tree, template ]
    }.combine(
        PREP_FS_SEQS.out.fs_dir, by : 0
    ).set{ ch_for_fslalign_lib }

    ch_for_fslalign_lib.multiMap{
        meta, fasta, tree, template, fs_dir ->
            fasta: [ meta, fasta ]
            tree:  [ meta, tree ]
            dependencies: [ meta, template, fs_dir ]
    }.set{ ch_fastas_forlibs_fslalign}

    // Compute the library
    FOLDSEEK_LALIGN_LIBRARY(ch_fastas_forlibs_fslalign.fasta,
                        ch_fastas_forlibs_fslalign.tree,
                        ch_fastas_forlibs_fslalign.dependencies,    
                        matrix_3di)
    ch_libs = ch_libs.mix(FOLDSEEK_LALIGN_LIBRARY.out.lib)
    ch_versions = ch_versions.mix(FOLDSEEK_LALIGN_LIBRARY.out.versions)


    //-----------------------------------------------------------------
    // FOLDSEEK LIBRARY FULL 
    //-----------------------------------------------------------------
    ch_lib_params.filter{ it[0].contains("FSfull_pair") }.multiMap{
                            metalib, meta, fastafile, tree, template, dependencies ->
                                fasta: [ meta, fastafile]
                                tree:  [ meta, tree ]
                                dependencies: [ meta, template, dependencies ]
                        }.set{ ch_fastas_forlibs_fsfull }
    
    // Convert each structure to 3di
    STRUCTURE_TO_3DI_FULL(ch_fastas_forlibs_fsfull.dependencies.map{ meta, template, structures -> [ meta, structures ]})
    // Merge all the 3di mappings into one file (id, fasta sequence, 3di sequence)
    MERGE_MAPPINGS_FULL(STRUCTURE_TO_3DI_FULL.out.mapping)
    ch_fastas_forlibs_fsfull.fasta.combine(
        MERGE_MAPPINGS_FULL.out.mapping, by:0
    ).set{ mappings_3di_full }

    ENCODE_FASTA(mappings_3di_full)
    ENCODE_FASTA.out.encoded_fasta.combine(ch_fastas_forlibs_fsfull.tree, by:0).multiMap{
        meta, fasta, mapping, encoded_fasta, tree ->
            mapping: [ meta, mapping ]
            tree: [ meta, tree ]
            encoded_fasta: [ meta, encoded_fasta ]
    }.set{ ch_fastas_forlibs_fsfull_tree }

    FOLDSEEK_FULL_LIBRARY(
        ch_fastas_forlibs_fsfull_tree.encoded_fasta,
        ch_fastas_forlibs_fsfull_tree.tree,
        [[:],[],[]],
        matrix_3di
    )

    CHANGE_LIBRARY_SEQUENCES(ch_fastas_forlibs_fsfull_tree.mapping.combine(FOLDSEEK_FULL_LIBRARY.out.lib.view(), by:0))
    ch_libs = ch_libs.mix(CHANGE_LIBRARY_SEQUENCES.out.library)

    // -----------------------------------------------------------------
    // MERGE LIBRARIES
    // -----------------------------------------------------------------
    MERGE_LIBRARIES_MANUAL(ch_libs.groupTuple(), "max")
    merged_libraries = MERGE_LIBRARIES_MANUAL.out.lib
    
    emit:
    ch_libs      = merged_libraries
}


def extract_lib_parameters(args){
    def matcher = args =~ /-lib\s(.*?)(?=\s-|$)/
    if(matcher.find()){
        return(matcher.group(1))
    }
}

