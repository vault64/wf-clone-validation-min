#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/fastqingress'



process assembleCore {
    errorStrategy = {task.attempt <= 4 ? 'retry' : 'ignore'}
    maxRetries 4
    label "wfplasmid"
    cpus params.threads
    input:
        tuple val(sample_id), file(fastq), val(approx_size)
    output:
        tuple val(sample_id), path("*.reconciled.fasta"), optional: true, emit: assembly
        tuple val(sample_id), path("*.downsampled.fastq"), optional: true, emit: downsampled
        tuple val(sample_id), env(STATUS), emit: status
    script:
        name = sample_id
        cluster_dir = "trycycler/cluster_001"
        int target = params.assm_coverage * 3
        int min_dep = (params.assm_coverage / 3) * 2
        int min_len = 100
        int max_len = approx_size.toInteger() * 1.2
        int min_q = 7
        int exit_number = task.attempt <= 4 ? 1 : 0
        def fast = params.fast == true ? '-fast' : ''
        def cluster_option = (params.canu_useGrid == false) ? """\
        -useGrid=false \
        -obtovlThreads=$task.cpus \
        -utgovlThreads=$task.cpus \
        -corThreads=$task.cpus \
        -redThreads=$task.cpus \
        -batThreads=$task.cpus """ : ""
        def windows_params = System.properties['os.version'].toLowerCase().contains("wsl") ? """\
        -mhapPipe=false \
        -purgeOverlaps=false \
        -saveOverlaps=true """ : ""

    """

    ############################################################
    # Trimming
    ############################################################
    STATUS="Failed to trim reads"
    (seqkit subseq -j $task.cpus -r $params.trim_length:-$params.trim_length $fastq | \
        seqkit subseq -j $task.cpus -r 1:$max_len | \
        seqkit seq -j $task.cpus -m $min_len -Q $min_q -g > ${name}.trimmed.fastq) \
        && STATUS="Failed to downsample reads" &&

    ############################################################
    # Downsampling
    ############################################################


    (rasusa \
        --coverage $target \
        --genome-size $approx_size \
        --input ${name}.trimmed.fastq > ${name}.downsampled.fastq) \
        && STATUS="Failed to Subset reads" &&

    ############################################################
    # Subsetting
    ############################################################

    (trycycler subsample \
        --count 3 \
        --min_read_depth $min_dep \
        --reads ${name}.downsampled.fastq \
        --out_dir sets \
        --genome_size $approx_size) \
        && STATUS="Failed to assemble using Canu" &&

    ############################################################
    # Assembly
    ############################################################

    (for SUBSET in \$(ls sets/sample_*.fastq)
    do
        SUBSET_NAME=\$(basename -s .fastq \$SUBSET)
        canu \
            -p \$SUBSET_NAME \
            -d assm_\${SUBSET_NAME} \
            -maxThreads=$task.cpus \
            genomeSize=$approx_size \
            $fast \
            -nanopore \$SUBSET \
            $cluster_option \
            $windows_params
    done) && STATUS="Failed to trim Assembly" &&

    ############################################################
    # Trim assemblies
    ############################################################

    (for ASSEMBLY in \$(ls assm_*/*.contigs.fasta)
    do
        ASSEMBLY_NAME=\$(basename -s .fasta \$ASSEMBLY)
        trim.py \
            \$ASSEMBLY \
            -o \${ASSEMBLY_NAME}.trimmed.fasta
        deconcatenate.py \
            \${ASSEMBLY_NAME}.trimmed.fasta \
            -o \${ASSEMBLY_NAME}.deconcat.fasta
    done
    ls *.deconcat.fasta 1> /dev/null 2>&1) \
    && STATUS="Failed to reconcile assemblies" &&

    ############################################################
    # Reconciliation
    ############################################################

    (trycycler cluster \
        --assemblies *.deconcat.fasta \
        --reads ${name}.downsampled.fastq \
        --out_dir trycycler) &&
    (trycycler reconcile \
        --reads ${name}.downsampled.fastq \
        --cluster_dir $cluster_dir \
        --max_trim_seq_percent 20 \
        --max_add_seq_percent 10) &&
    (trycycler msa --cluster_dir $cluster_dir) &&
    (trycycler partition --reads ${name}.downsampled.fastq --cluster_dirs $cluster_dir) &&
    (trycycler consensus --cluster_dir $cluster_dir)

    ############################################################
    # Exit handling
    ############################################################

    if [ ! -f "${cluster_dir}/7_final_consensus.fasta" ]; then
        if ls ${cluster_dir}/1_contigs/*.fasta 1> /dev/null 2>&1; then
            STATUS="Completed but failed to reconcile"
            (seqkit sort ${cluster_dir}/1_contigs/*.fasta --by-length \
                | seqkit head -n 1 > ${name}.reconciled.fasta) \
                && echo "Trycycler failed, outputting un-reconciled assembly"
        elif [ "$exit_number" == "1" ]; then
            echo \$STATUS
            echo "Assembly failed, retrying process"
            exit 1
        elif [ "$exit_number" == "0" ]; then
            echo \$STATUS
            echo "Failed final attempt"
        fi
    else
        mv ${cluster_dir}/7_final_consensus.fasta ${name}.reconciled.fasta
        STATUS="Completed successfully"
    fi
    """
}


process medakaPolishAssembly {
    label "wfplasmid"
    cpus params.threads
    input:
        tuple val(sample_id), file(draft), file(fastq)
    output:
        tuple val(sample_id), path("*.final.fasta"), emit: polished
        tuple val(sample_id), env(STATUS), emit: status
    """
    STATUS="Failed to polish assembly with Medaka"
    medaka_consensus -i $fastq -d $draft -m r941_min_high_g360 -o . -t $task.cpus -f
    echo ">${sample_id}" >> ${sample_id}.final.fasta
    sed "2q;d" consensus.fasta >> ${sample_id}.final.fasta
    STATUS="Completed successfully"
    """
}



process report {
    label "wfplasmid"
    cpus 1
    input:
        path "per_barcode_stats/*"
        path "host_filter_stats/*"
        path lengths
    output:
        path "wf-clone-validation-*.html", emit: html
        path "sample_status.txt", emit: sample_stat
        path "inserts/*", optional: true, emit: inserts
    script:
        report_name = "wf-clone-validation-" + params.report_name + '.html'
    """
    report.py \
    --downsampled_stats downsampled_stats/* \
    --revision $workflow.revision \
    --commit $workflow.commitId \
    --status $final_status \
    --per_barcode_stats per_barcode_stats/* \
    --host_filter_stats host_filter_stats/* \
    --params params.json \
    --versions versions \
    --report_name $report_name \
    --plannotate_json $plannotate_json \
    --lengths $lengths \
    --inserts_json $inserts_json
    """
}

workflow pipeline {
    take:
        samples
        host_reference
        regions_bedfile
        database
        primers
        align_ref

    main:
        // Combine fastq from each of the sample directories into
        // a single per-sample fastq file
        named_samples = samples.map { it -> return tuple(it[1],it[0])}
        if(params.approx_size_sheet != null) {
            approx_size = Channel.fromPath(params.approx_size_sheet) \
            | splitCsv(header:true) \
            | map { row-> tuple(row.sample_id, row.approx_size) }
            final_samples = named_samples.join(approx_size)}
        else {
            final_samples = samples.map  { it -> return tuple(it[1],it[0], params.approx_size)}
        }
        sample_fastqs = combineFastq(final_samples)
        // Optionally filter the data, removing reads mapping to
        // the host or background genome
        if (host_reference.name != "NO_HOST_REF") {
            filtered = filterHostReads(
                    sample_fastqs.sample, host_reference, regions_bedfile)
            samples_filtered = filtered.unmapped
            updated_status = filtered.status
            filtered_stats = filtered.host_filter_stats.collect()
                             .ifEmpty(file("$projectDir/data/OPTIONAL_FILE"))
        }
        else {
            samples_filtered = sample_fastqs.sample
            updated_status = sample_fastqs.status
            filtered_stats = file("$projectDir/data/OPTIONAL_FILE")
        }
       
        // Core assembly and reconciliation
        assemblies = assembleCore(samples_filtered)
        
        named_drafts = assemblies.assembly.groupTuple()
        named_samples = assemblies.downsampled.groupTuple()
        named_drafts_samples = named_drafts.join(named_samples)


        // Polish draft assembly
        polished = medakaPolishAssembly(named_drafts_samples)
       
        // Concat statuses and keep the last of each
        final_status = sample_fastqs.status.concat(updated_status)
        .concat(assemblies.status).concat(polished.status).groupTuple()
        .map { it -> it[0].toString() + ',' + it[1][-1].toString() }
        final_status = final_status.collectFile(name: 'final_status.csv', newLine: true)
    
        downsampled_stats = downsampledStats(assemblies.downsampled)

        primer_beds = findPrimers(primers, polished.polished)
        software_versions = getVersions()
        workflow_params = getParams()

        annotation = runPlannotate(
            database, polished.polished.map { it -> it[1] }.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
            final_status)

        insert = inserts(primer_beds.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
            polished.polished.map { it -> it[1] }.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
            align_ref)
        report = report(
            downsampled_stats.collect().ifEmpty(file("$projectDir/data/OPTIONAL_FILE")),
            final_status,
            sample_fastqs.stats.collect(),
            filtered_stats,
            software_versions.collect(),
            workflow_params,
            annotation.report,
            insert.json,
            annotation.json)

        results = polished.polished.map { it -> it[1] }.concat(
            report.html,
            report.sample_stat,
            annotation.feature_table,
            insert.inserts,
            annotation.json,
            annotation.annotations,
            workflow_params)
    emit:
        results
        telemetry = workflow_params
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    label "wfplasmid"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*", saveAs: {
        f -> params.prefix ? "${params.prefix}-${f}" : "${f}" }
    input:
        file fname
    output:
        file fname
    """
    echo "Writing output files"
    """
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    if (params.disable_ping == false) {
        try { 
            Pinguscript.ping_post(workflow, "start", "none", params.out_dir, params)
        } catch(RuntimeException e1) {
        }
    }
    
    samples = fastq_ingress([
        "input":params.fastq,
        "sample":params.sample,
        "sample_sheet":params.sample_sheet,
        "sanitize": params.sanitize_fastq,
        "output":params.out_dir,
        "min_barcode":params.min_barcode,
        "max_barcode":params.max_barcode])
    host_reference = file(params.host_reference, type: "file")
    regions_bedfile = file(params.regions_bedfile, type: "file")
    primer_file = file("$projectDir/data/OPTIONAL_FILE")
    if (params.primers != null){
        primer_file = file(params.primers, type: "file")
    }
    align_ref = file("$projectDir/data/OPTIONAL_FILE")
    if (params.reference != null){
        align_ref = file(params.reference, type: "file")
    }
    database = file("$projectDir/data/OPTIONAL_FILE")
    if (params.db_directory != null){
         database = file(params.db_directory, type: "dir")

    }

    // Run pipeline
    results = pipeline(samples, host_reference, regions_bedfile, database, primer_file, align_ref)

    output(results[0])
   
}


if (params.disable_ping == false) {
    workflow.onComplete {
        try{
            Pinguscript.ping_post(workflow, "end", "none", params.out_dir, params)
        }catch(RuntimeException e1) {
        }
    }
    
    workflow.onError {
        try{
            Pinguscript.ping_post(workflow, "error", "$workflow.errorMessage", params.out_dir, params)
        }catch(RuntimeException e1) {
        }
    }
}
