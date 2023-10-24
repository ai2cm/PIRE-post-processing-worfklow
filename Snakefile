configfile: "config/config.yaml"

include: "rules/common.smk"
include: "rules/mppnccombine.smk"


def combine_and_coarsen_restart_inputs():
    inputs = []
    for run in RUNS:
        output = rules.coarsen_restarts_run_timestamp.output
	timestamps = COARSE_RESTART_TIMESTAMPS[run]
        restarts = expand(output, run=run, timestamp=timestamps)
        inputs.extend(restarts)
    return inputs


def combine_and_coarsen_diagnostics_inputs():
    inputs = []
    for run in RUNS:
        output = rules.combine_and_coarsen_tape.output
        tapes = TAPES[run]
        diagnostics = expand(output, run=run, tape=tapes)
        inputs.extend(diagnostics)
    return inputs


rule combine_and_coarsen_restarts:
    input:
        combine_and_coarsen_restart_inputs()


rule combine_and_coarsen_diagnostics:
    input:
        combine_and_coarsen_diagnostics_inputs()


rule transfer:
    input:
        rules.combine_and_coarsen_restarts.input,
	rules.combine_and_coarsen_diagnostics.input
    output:
        transfer_sentinel()
    params:
        source=transfer_source(),
	gcs_path=transfer_gcs_path()
    shell:
        """
        gsutil -m cp -n -r {params.source} {params.gcs_path}
        """

rule all:
    input:
        rules.transfer.output
