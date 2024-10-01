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

    // --------------------
    // Get the libraries
    // --------------------
    TCOFFEELIB_LIBRARIES(ch_fastas, 
                         ch_trees, 
                         ch_templates)

    // Prepare channel for alignment 
    ch_fastas.join(ch_trees).join(TCOFFEELIB_LIBRARIES.out.ch_libs).set{ch_fasta_libs}
    ch_fasta_libs.map{
                    meta, fasta, tree, lib_file ->
                        [ remove_lib_parameters(meta), fasta, tree, lib_file ]
                    }.multiMap{
                        meta, fasta_file, tree_file,  lib_file ->
                            fasta:   [ meta, fasta_file ]
                            tree:    [ meta, tree_file ]
                            library: [ meta, lib_file ]
                    }.set{ ch_libs }

    TCOFFEE_ALIGN_FORLIB(ch_libs.fasta,
                         ch_libs.tree,
                         [[:],[],[]],
                         ch_libs.library,
                         compress)

    ch_alignments = ch_alignments.mix(TCOFFEE_ALIGN_FORLIB.out.alignment)
    ch_versions   = ch_versions.mix(TCOFFEE_ALIGN_FORLIB.out.versions)

    // Re-add the library parameters to the metadata
    ch_alignments.map{
        meta, alignment -> 
            newmeta = meta.clone()
            newmeta["args_aligner"] = meta["args_aligner"] + " " + meta["args_aligner_lib"]
            newmeta.remove("args_aligner_lib")
            [ newmeta, alignment ]
    }.set{ ch_alignments }

    emit: 
    alignment    = ch_alignments
    versions     = ch_versions

}


def remove_lib_parameters(meta){
    // from the key in meta, remove the flag
    def args = meta["args_aligner"]
    def matcher = args =~ /(-lib\s.*?)(?=\s-|$)/
    def new_args
    def lib_params
    if(args == null){
        new_args =  args
    }else{
        if(matcher.find()){
            lib_params = matcher.group(1)
            new_args =  args.replaceAll(/-lib\s(.*?)(?=\s-|$)/, "")
        } else {
            new_args = args
        }
    }
    meta["args_aligner"] = new_args
    meta["args_aligner_lib"] = lib_params
    return(meta)
}
