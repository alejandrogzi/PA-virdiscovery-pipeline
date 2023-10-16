process hmm_postprocessing {
      publishDir "${params.output}/${name}/${params.hmmerdir}/", mode: 'copy', pattern: "${set_name}_modified.tsv"
      label 'python3'

    input:
      tuple val(name), val(set_name),file(faa)
      file(viphogs)
      file(rvdb)
      file(pvogs)
      file(vogsdb)
    
    output:
      tuple val(name), val(set_name), file("${set_name}_modified.tsv"), file(faa)
    
    script:
    """
    make_hmm_table.py -t ${viphogs} ${rvdb} ${pvogs} ${vogsdb} -o ${set_name}_modified
    """
}

/*
input: File_hmmer_ViPhOG.tbl
output: File_hmmer_ViPhOG_modified.tbl
*/
