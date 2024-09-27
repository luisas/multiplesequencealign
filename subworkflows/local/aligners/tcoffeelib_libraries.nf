include { TCOFFEE_LIBRARY as SEQUENCE_LIBRARY   } from '../../../modules/local/tcoffee/library/main.nf'

workflow TCOFFEELIB_LIBRARIES {

    take: 
    ch_fastas     //channel: [ meta, /path/to/file.fasta ]

    main:

    ch_libs     = Channel.empty()
    ch_versions = Channel.empty()


    // remove the library parameters from the command line
    ch_fastas.branch {
        sequence: it[1]["id_lib"] == "proba_pair"
        foldseek: it[1]["id_lib"] == "fs_pair"
    }.map{
        meta, fasta ->
            [ meta, fasta ]
    }.set { ch_fastas_forlibs }


    // SEQUENCE LIBRARY 
    SEQUENCE_LIBRARY(ch_fastas_forlibs.sequence,
                     [[:],[]],
                     [[:],[],[]],    
                     [[:],[]])
    ch_libs = ch_libs.mix(SEQUENCE_LIBRARY.out.lib)
    ch_versions = ch_versions.mix(SEQUENCE_LIBRARY.out.versions)

    

    emit:
    ch_fasta_libs      = ch_libs


}