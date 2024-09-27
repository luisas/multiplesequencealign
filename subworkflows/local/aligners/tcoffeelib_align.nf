include { TCOFFEELIB_LIBRARIES } from './tcoffeelib_libraries.nf'
include { TCOFFEE_ALIGN as TCOFFEE_ALIGN_FORLIB } from '../../../modules/nf-core/tcoffee/align/main'

workflow TCOFFEELIB_ALIGN {

    take: 
    ch_fastas     //channel: [ meta, /path/to/file.fasta ]
    ch_trees      //channel: [ meta, /path/to/file.tree ]
    ch_templates  //channel: [ meta, /path/to/file.template, /path/to/file.accessory_informations ]
    compress      //boolean


    main: 

    ch_alignments = Channel.empty()
    ch_versions   = Channel.empty()

    // Extract the part of the command line that specifies the library
    ch_fastas.map{
        meta, fasta ->
            [ extract_lib_parameters(meta["args_aligner"]).split().toList(), meta, fasta  ]
    }.transpose().map{
        lib_params, meta, fastas ->
            [ meta, ["id_lib":lib_params], fastas ]
    }.set{ lib_params }

    // --------------------
    // Get the libraries
    // --------------------
    //TCOFFEELIB_LIBRARIES(lib_params)

    // // Prepare channel for alignment 
    // ch_fastas.join(TCOFFEELIB_LIBRARIES.out.ch_fasta_libs.map{
    //     meta, lib_file -> 
    //         newmeta = meta.clone()
    //         newmeta.remove("id_lib")
    //         [ newmeta, lib_file ]
    // }, by:0).multiMap{
    //     meta, fasta_file, lib_file ->
    //         fasta:   [ meta, fasta_file ]
    //         library: [ meta, lib_file ]
    // }.set{ ch_libs }

    // ch_libs.fasta.view()
    // // Align with the libraries
    // // TODO: ADD TREE
    // // TODO: ADD TEMPLATE
    // TCOFFEE_ALIGN_FORLIB(ch_libs.fasta,
    //                      [[:],[]],
    //                      [[:],[],[]],
    //                      ch_libs.library,
    //                      compress)

    // ch_alignments = ch_alignments.mix(TCOFFEE_ALIGN_FORLIB.out.alignment)
    // ch_versions   = ch_versions.mix(TCOFFEE_ALIGN_FORLIB.out.versions)

    emit: 
    alignment    = ch_alignments
    versions     = ch_versions

}

def extract_lib_parameters(args){
    def matcher = args =~ /-lib\s(.*?)(?=\s-|$)/
    if(matcher.find()){
        return(matcher.group(1))
    }
}