rule upload_restart_files_for_timestamp:
    input:
        files=get_restart_files_for_timestamp
    output:
        touch("uploaded-restart-sentinels/{run}-{timestamp}.out")
    params:
        gcs_url=f"{GCS_BUCKET}/{{run}}/coarse-restarts/{{timestamp}}",
        partition="rdtn"
    script:
        """../scripts/upload_restarts.py"""
