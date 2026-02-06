nextflow.enable.dsl = 2

// --------------------------
// INPUT / OUTPUT
// --------------------------
params.vcf_dir     = "${projectDir}"                 // where your S*.hard-filtered.vcf live
params.vcf_pattern = "S*.hard-filtered.vcf"         // matches S3.hard-filtered.vcf etc.

params.outdir      = "${projectDir}/results"        // all outputs here


// --------------------------
// ANNOVAR CONFIG
// --------------------------
params.annovar_dir = "${System.getenv('HOME')}/annovar"
params.humandb     = "${params.annovar_dir}/humandb"
params.build       = "hg38"

// you already have these dbs in ~/annovar/humandb
params.annovar_protocol  = "refGene,ensGene,gnomad211_exome,gnomad30_genome,clinvar_202409"
params.annovar_operation = "g,g,f,f,f"


// --------------------------
// VEP CONFIG (OPTIONAL)
//  - you are currently running with: --run_vep false
// --------------------------
params.run_vep   = false                            // default off; override with --run_vep true
params.vep_exe   = "vep"                            // change if you install VEP
params.vep_cache = "/path/to/vep/cache"             // update later when VEP is set up


// --------------------------
// ZIMBABWE / AFRICAN LAYER
// --------------------------
params.run_zim_layer = true
params.zim_db        = "/home/nyasha01/projects/zi-mendelian/refs/zimbabwe_cohort_freq.tsv"


// --------------------------
// RESOURCES
// --------------------------
// Your machine only had ~7.6 GB free earlier, so don't request 8 GB.
params.max_cpus = 4
params.max_mem  = "7 GB"


/*
 * INPUT CHANNEL
 *  - picks up all VCFs that match the pattern in params.vcf_dir
 */
Channel
    .fromPath("${params.vcf_dir}/${params.vcf_pattern}")
    .ifEmpty { error "No VCFs found matching: ${params.vcf_dir}/${params.vcf_pattern}" }
    .map { vcf_file -> tuple(vcf_file.baseName, vcf_file) }
    .set { ch_vcf }

/*
 * OPTIONAL: VEP ANNOTATION
 *  - takes input VCF and produces VEP-annotated VCF
 */
process VEP_ANNOTATION {

    tag "${sample_id}"

    publishDir "${params.outdir}/vep", mode: 'copy', overwrite: true

    cpus params.max_cpus
    memory params.max_mem

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path("${sample_id}.vep.vcf.gz")

    when:
    params.run_vep

    script:
    """
    ${params.vep_exe} \
        --input_file ${vcf} \
        --output_file ${sample_id}.vep.vcf \
        --vcf \
        --compress_output bgzip \
        --assembly GRCh38 \
        --cache \
        --offline \
        --dir_cache ${params.vep_cache} \
        --symbol \
        --canonical \
        --af_gnomad \
        --everything

    tabix -p vcf ${sample_id}.vep.vcf.gz
    """
}

/*
 * If VEP is disabled, bypass it and just send original VCF to ANNOVAR.
 */
ch_vcf_for_annovar = params.run_vep ?
    VEP_ANNOTATION.out.map { sample_id, vep_vcf -> tuple(sample_id, vep_vcf) } :
    ch_vcf

/*
 * CONVERT VCF -> ANNOVAR INPUT (.avinput)
 */
process ANNOVAR_CONVERT {

    tag "${sample_id}"

    publishDir "${params.outdir}/annovar/input", mode: 'copy', overwrite: true

    cpus 1
    memory "2 GB"

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path("${sample_id}.avinput")

    script:
    """
    perl ${params.annovar_dir}/convert2annovar.pl \
        -format vcf4 \
        -allsample \
        -withfreq \
        -includeinfo \
        ${vcf} \
        > ${sample_id}.avinput
    """
}

process ANNOVAR_TABLE {

    tag "${sample_id}"

    publishDir "${params.outdir}/annovar/annotated",
               mode: 'copy',
               overwrite: true

    cpus params.max_cpus ?: 4
    memory '6 GB'

    input:
    tuple val(sample_id), path(avinput)

    output:
    // e.g. S4.hard-filtered.hg38_multianno.txt
    tuple val(sample_id), path("${sample_id}.${params.build}_multianno.txt")

    script:
    """
    perl ${params.annovar_dir}/table_annovar.pl \
        ${avinput} \
        ${params.humandb} \
        -buildver ${params.build} \
        -out ${sample_id} \
        -remove \
        -protocol ${params.annovar_protocol} \
        -operation ${params.annovar_operation} \
        -nastring .

    # table_annovar will produce ${sample_id}.${params.build}_multianno.txt
    """
}

process ZIM_PRIORITISATION {

    tag "${sample_id}"

    publishDir "${params.outdir}/prioritisation", mode: 'copy', overwrite: true

    cpus 1
    memory "2 GB"

    input:
    tuple val(sample_id), path(multianno_txt)

    output:
    tuple val(sample_id), path("${sample_id}.zim_prioritised.tsv")

    when:
    params.run_zim_layer

    script:
    """
    python3 ${projectDir}/bin/zim_prioritisation.py \
      --multianno ${multianno_txt} \
      --zim_db ${params.zim_db} \
      --out_prefix ${sample_id}
    """
}


workflow {

    // 1. Find VCFs
    ch_vcfs = Channel
        .fromPath("${params.vcf_dir}/${params.vcf_pattern}")
        .map { vcf ->
            // e.g. sample_id = "S4.hard-filtered"
            tuple(vcf.baseName, vcf)
        }

    // 2. VCF → AVINPUT
    ch_avinput   = ANNOVAR_CONVERT(ch_vcfs)

    // 3. AVINPUT → ANNOVAR multianno
    ch_multianno = ANNOVAR_TABLE(ch_avinput)

    // 4. Optional Zimbabwe layer
    ch_zim = ZIM_PRIORITISATION(ch_multianno)

    emit:
        annotated = ch_multianno
        zim_layer = ch_zim
}

       
