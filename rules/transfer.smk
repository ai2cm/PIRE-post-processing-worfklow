rule transfer_restarts_for_timestamp:
    input:
        rules.coarsen_restarts_run_timestamp.output
    output:
        touch(f"transfer-sentinels/{CONFIG}/restarts/{{run}}/{{timestamp}}.out")
    params:
        gcs_path=f"{GCS_BUCKET}/{{run}}/restarts/{{timestamp}}/"
    shell:
        """
        gsutil -m cp -n {input} {params.gcs_path}
        """


rule transfer_tape:
    input:
        sentinel=rules.combine_and_coarsen_tape.output,
        store=rules.initialize_combine_and_coarsen_tape.output.store
    output:
        touch(f"transfer-sentinels/{CONFIG}/zarrs/{{run}}/{{tape}}.out")
    params:
        gcs_path=f"{GCS_BUCKET}/{{run}}/diagnostics/"
    shell:
        """
        gsutil -m cp -n -r {input.store} {params.gcs_path}
        """
