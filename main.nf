#!/usr/bin/env nextflow

if(params.help) {
    usage = file("$baseDir/USAGE")
    cpu_count = Runtime.runtime.availableProcessors()

    bindings = ["register_processes":"$params.register_processes",
                "cpu_count":"$cpu_count",
                "resampling_streamlines": "$params.resampling_streamlines",
                "resampling_tractograms": "$params.resampling_tractograms",
                "registration_speed":"$params.registration_speed"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)
    print template.toString()
    return
}

log.info "Registration Flow"
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
log.info "Resampling (Streamlines): $params.resampling_streamlines"
log.info "Resampling (Tractograms): $params.resampling_tractograms"

log.info ""

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

log.info "Input: $params.input"
root = file(params.input)

if (params.registration_speed == 1){
    registration_script = Channel.value("antsRegistrationSyNQuick.sh")
}
else if (params.registration_speed == 0){
    registration_script = Channel.value("antsRegistrationSyN.sh")
}
else {
    error "Registration speed must be 0 or 1"
}

Channel
    .fromFilePairs("$root/**/bundles/*.trk", size: -1) { it.parent.parent.name }
    .set{ tractogram } // [sid, tractogram.trk]

Channel
    .fromFilePairs("$root/**/metrics/*.nii.gz", size: -1) { it.parent.parent.name }
    .into{ reference; reference_for_check } // [sid, ref.nii.gz]

num_ref = reference_for_check.map { [it[0], it[1].size()] }

target_anat = Channel.fromPath("$params.target_anat")

reference
    .combine(target_anat)
    .combine(registration_script)
    .join(num_ref)
    .set{reference_target_anat} 

process Register_Anat {
    cpus params.register_processes
    memory '2 GB'

    input:
    set sid, file(reference), file(target_anat), val(registration_script), val(number_of_ref) from reference_target_anat

    output:
    // [sid, affine.mat, inverseWarp.nii.gz, fixed_ref.nii.gz]
    set sid, "${sid}__output0GenericAffine.mat", "${sid}__output1InverseWarp.nii.gz", "${target_anat}" into transformation_for_tractogram
    file "${sid}__outputWarped.nii.gz"
    file "${sid}__outputInverseWarped.nii.gz"
    file "${sid}__output1Warp.nii.gz"

    when:
    number_of_ref == 1

    script:
    """
    export ANTS_RANDOM_SEED=1234
    ${registration_script} -d 3 -f ${target_anat} -m ${reference} -o ${sid}__output -t s -n ${params.register_processes}
    """
}

// [sid, tractogram.trk, affine.mat, inverseWarp.nii.gz, outputWarped.nii.gz]
tractogram.join(transformation_for_tractogram).set{tractogram_registration} 

process Register_Streamlines {
    memory '2 GB'

    input:
    set sid, file(tractogram), file(affine), file(inverse_warp), file(target_anat) from tractogram_registration
    val resampling_streamlines from params.resampling_streamlines
    val resampling_tractograms from params.resampling_tractograms

    output:
    file "*_registered.trk"

    script:
    """
    for f in ${tractogram}
    do 
        filename=\$(basename -- "\$f")
        root_name="\${filename%.*}"
        if [[ ${resampling_tractograms} -ge 1 ]]; then
           scil_resample_tractogram.py \$f ${resampling_tractograms} \$f -f -v
        fi

        scil_apply_transform_to_tractogram.py \$f ${target_anat} \
        ${affine} \${root_name}_registered.trk \
        --inverse --in_deformation ${inverse_warp} -f -vv

        if [[ ${resampling_streamlines} -ge 3 ]]; then
           scil_resample_streamlines.py \${root_name}_registered.trk \${root_name}_registered.trk --nb_pts_per_streamline ${resampling_streamlines} -f
        fi
    done
    """
}