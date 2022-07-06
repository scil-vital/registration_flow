#!/usr/bin/env nextflow

if(params.help) {
    usage = file("$baseDir/USAGE")
    cpu_count = Runtime.runtime.availableProcessors()

    bindings = ["register_processes":"$params.register_processes",
                "cpu_count":"$cpu_count",
                "resampling": "$params.resampling"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)
    print template.toString()
    return
}

log.info "FINTA Multibundles Flow"
log.info "==============================================="
log.info ""
log.info "Start time: $workflow.start"
log.info ""

log.debug "[Command-line]"
log.debug "$workflow.commandLine"
log.debug ""

log.info "[Git Info]"
log.info "$workflow.repository - $workflow.revision [$workflow.commitId]"
log.info ""

log.info "Options"
log.info "======="
log.info ""
log.info "[Target]"
log.info "Target anat: $params.target_anat"
log.info "Resampling: $params.resampling"
log.info ""

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

log.info "Input: $params.input"
root = file(params.input)

Channel
    .fromFilePairs("$root/**/bundles/*.trk", size: -1) { it.parent.parent.name }
    .set{ tractogram } // [sid, tractogram.trk]

Channel
    .fromPath("$root/**/metrics/*.nii.gz")
    .map{[it.parent.parent.name, it]}
    .set{ reference } // [sid, t1.nii.gz]

target_anat = Channel.fromPath("$params.target_anat")

reference
    .combine(target_anat)
    .set{reference_target_anat} 

process Register_Anat {
    cpus params.register_processes
    memory '2 GB'

    input:
    set sid, file(reference), file(target_anat) from reference_target_anat

    output:
    // [sid, affine.mat, inverseWarp.nii.gz, fixed_t1.nii.gz]
    set sid, "${sid}__output0GenericAffine.mat", "${sid}__output1InverseWarp.nii.gz", "${target_anat}" into transformation_for_tractogram
    file "${sid}__outputWarped.nii.gz"
    file "${sid}__output1Warp.nii.gz"

    script:
    """
    export ANTS_RANDOM_SEED=1234
    antsRegistrationSyNQuick.sh -d 3 -f ${target_anat} -m ${reference} -o ${sid}__output -t s -n ${params.register_processes}
    """
}

// [sid, tractogram.trk, affine.mat, inverseWarp.nii.gz, outputWarped.nii.gz]
tractogram.join(transformation_for_tractogram).set{tractogram_registration} 

process Register_Streamlines {
    memory '2 GB'

    input:
    set sid, file(tractogram), file(affine), file(inverse_warp), file(target_anat) from tractogram_registration
    val resampling from params.resampling

    output:
    file "${sid}_registered_*"

    script:
    """
    for f in ${tractogram}
    do 
        if [[ ${resampling} -ge 3 ]]; then
           scil_resample_streamlines.py \$f \$f --nb_pts_per_streamline ${resampling} -f
        fi

        scil_apply_transform_to_tractogram.py \$f ${target_anat} \
        ${affine} ${sid}_registered_\$f \
        --inverse --in_deformation ${inverse_warp} -f -vv
    done
    """
}