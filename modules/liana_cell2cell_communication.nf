nextflow.enable.dsl=2

//out_dir = file(params.resDir)
mode = params.publish_dir_mode

process LIANA{
    publishDir "${out_dir}", mode: "$mode"
    label "LIANA"

    input:
    path(adata)


    output:
    path("*"), emit: liana_output
 


	script:
	"""
    001_liana.py \\
    --adata=${adata} \\
    --adata_out adata_processed_LIANA_by_sample.h5ad\\
	"""

}
