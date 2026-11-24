nextflow.enable.dsl=2

process MAKE_DB {
    tag "$name"
    container 'quay.io/biocontainers/diamond:2.1.11--h13889ed_1'
    cpus 2
    memory '8 GB'
    input:
    tuple val(name), path(reference)
    output:
    tuple val(name), path("${name}.dmnd"), emit: db
    script:
    """
    diamond makedb --in '${reference}' --db '${name}' --threads ${task.cpus}
    """
    stub:
    """
    touch '${name}.dmnd'
    """
}

process SCREEN {
    tag "${sample}:${db_name}"
    container 'quay.io/biocontainers/diamond:2.1.11--h13889ed_1'
    cpus 8
    memory '24 GB'
    maxForks 8
    input:
    tuple val(sample), path(proteins), val(db_name), path(db)
    output:
    tuple val(sample), val(db_name), path("${sample}.${db_name}.diamond.tsv"), emit: hits
    script:
    """
    diamond blastp --query '${proteins}' --db '${db}' --threads ${task.cpus} \\
      --out '${sample}.${db_name}.diamond.tsv' --evalue ${params.evalue} \\
      --id ${params.min_identity} --query-cover ${params.min_query_coverage} \\
      --subject-cover ${params.min_subject_coverage} --max-target-seqs 1 \\
      --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle
    """
    stub:
    """
    touch '${sample}.${db_name}.diamond.tsv'
    """
}

process ANNOTATE {
    tag "$sample"
    container 'python:3.11-slim'
    cpus 2
    memory '8 GB'
    input:
    tuple val(sample), path(proteins), path(ncyc), path(bacmet)
    path mapping
    output:
    tuple val(sample), path("${sample}.functional_hits.tsv"), emit: hits
    tuple val(sample), path("${sample}.summary.json"), emit: summary
    script:
    """
    annotate_hits.py --sample '${sample}' --proteins '${proteins}' \\
      --ncyc '${ncyc}' --bacmet '${bacmet}' --ncyc-map '${mapping}' \\
      --output '${sample}.functional_hits.tsv' --summary '${sample}.summary.json'
    """
    stub:
    """
    printf 'sample\\tcontig\\tgene_id\\n' > '${sample}.functional_hits.tsv'
    printf '{}\\n' > '${sample}.summary.json'
    """
}

workflow {
    main:
    if (!params.input || !params.ncyc_fasta || !params.ncyc_map || !params.bacmet_fasta)
        error 'Required: --input, --ncyc_fasta, --ncyc_map and --bacmet_fasta'
    if (params.min_identity < 0 || params.min_identity > 100 || params.min_query_coverage < 0 || params.min_query_coverage > 100 || params.min_subject_coverage < 0 || params.min_subject_coverage > 100)
        error 'Identity and coverage thresholds must be percentages between 0 and 100'
    samples = channel.fromPath(params.input, checkIfExists: true).splitCsv(header: true).map { row ->
        if (!row.sample || !row.proteins) error "Input requires sample and proteins columns: ${row}"
        if (!(row.sample ==~ /[A-Za-z0-9_.-]+/)) error "Unsafe sample identifier: ${row.sample}"
        tuple(row.sample.toString(), file(row.proteins.toString(), checkIfExists: true))
    }
    refs = channel.of(tuple('ncyc', file(params.ncyc_fasta, checkIfExists: true)), tuple('bacmet', file(params.bacmet_fasta, checkIfExists: true)))
    MAKE_DB(refs)
    SCREEN(samples.combine(MAKE_DB.out.db))
    ncyc_hits = SCREEN.out.hits.filter { sample, db, hit -> db == 'ncyc' }.map { sample, db, hit -> tuple(sample, hit) }
    bacmet_hits = SCREEN.out.hits.filter { sample, db, hit -> db == 'bacmet' }.map { sample, db, hit -> tuple(sample, hit) }
    annotation_inputs = samples.join(ncyc_hits, failOnDuplicate: true, failOnMismatch: true).join(bacmet_hits, failOnDuplicate: true, failOnMismatch: true)
    ANNOTATE(annotation_inputs, file(params.ncyc_map, checkIfExists: true))
    publish:
    reference_databases = MAKE_DB.out.db
    raw_alignments = SCREEN.out.hits
    annotated_hits = ANNOTATE.out.hits
    summaries = ANNOTATE.out.summary
}

output {
    reference_databases { path 'reference_databases' }
    raw_alignments { path 'diamond' }
    annotated_hits { path 'functional_hits' }
    summaries { path 'summaries' }
}
