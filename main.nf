#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { LIANA } from "./modules/liana_cell2cell_communication"

workflow {
    adata = Channel.fromPath(params.adata)
    //adata_liana= Channel.fromPath(params.adata_liana)
    outdir = Channel.fromPath(params.outdir)

    LIANA(adata)

}