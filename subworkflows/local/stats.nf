//
// Compute stats about the input sequences
//

include {   TCOFFEE_SEQREFORMAT_SIM       } from '../../modules/local/tcoffee_seqreformat_sim.nf'


workflow STATS {
    take:
    ch_seqs                //      channel: meta, /path/to/file.fasta
   

    main:

    ch_versions = Channel.empty()
    TCOFFEE_SEQREFORMAT_SIM(ch_seqs)
    ch_versions = ch_versions.mix(TCOFFEE_SEQREFORMAT_SIM.out.versions.first())


    emit:
    tcoffee_seqreformat_sim            = TCOFFEE_SEQREFORMAT_SIM.out.perc_sim                      // TODO
    tcoffee_seqreformat_simtot            = TCOFFEE_SEQREFORMAT_SIM.out.perc_sim_tot                      // TODO
    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}