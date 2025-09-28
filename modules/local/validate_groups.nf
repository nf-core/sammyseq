process VALIDATE_GROUPS {
    tag "validate_groups"

    input:
    val comparison_results
    val compare_groups

    output:
    val "validated", emit: validation

    exec:
    def available = []
    for (int i = 0; i < comparison_results.size(); i += 2) {
        available.add(comparison_results[i].sample_group)
    }
    available = available.unique()
    
    def specified = compare_groups.split(/[,vs]/).collect { it.trim() }.unique() - ['']
    
    def missing = specified - available
    if (missing) {
        error "ERROR: Groups ${missing} not found! Available: ${available.join(', ')}"
    }
}