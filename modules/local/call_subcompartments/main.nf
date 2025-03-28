process CALL_SUBCOMPARTMENTS {
    container 'docker.io/ciuki97/sammy_subcompartments_env:latest'

    input:
    val tsv_content

    output:
    path "test_compartments_output.txt", emit: result
    path "compartments_input.tsv", emit: tsv_file

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo "${tsv_content.join('\n')}" > compartments_input.tsv
    template 'call_subcompartments.R' compartments_input.tsv
    """
}
