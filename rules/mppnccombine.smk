rule mppncombine_restart_file_run_tape_tile_timestamp:
    input:
        restart_file_subtiles
    output:
        mppnccombined_restart_file()
    wildcard_constraints:
        # Ensure snakemake does not combine timestamp and tape wildcards.
        timestamp="\d+.\d+"
    shell:
        """
        module load intel/2021.1.2
        module load openmpi/intel-2021.1/4.1.2
        module load netcdf/intel-2021.1/hdf5-1.10.6/4.7.4
        module load hdf5/intel-2021.1/1.10.6
        mppnccombine {output} {input}
        """


rule coarsen_restarts_run_timestamp:
    input:
        fv_core=restarts_as_inputs("fv_core_coarse.res", from_mppnccombine=True),
	fv_tracer=restarts_as_inputs("fv_tracer_coarse.res", from_mppnccombine=True),
	fv_srf_wnd=restarts_as_inputs("fv_srf_wnd_coarse.res", from_mppnccombine=True),
	sfc_data=restarts_as_inputs("sfc_data_coarse", from_mppnccombine=True),
	coupler_res=restarts_as_inputs("coupler.res", from_mppnccombine=True)
    output:
        fv_core=restarts_as_outputs("fv_core.res"),
	fv_tracer=restarts_as_outputs("fv_tracer.res"),
	fv_srf_wnd=restarts_as_outputs("fv_srf_wnd.res"),
	sfc_data=restarts_as_outputs("sfc_data"),
	coupler_res=restarts_as_outputs("coupler.res")
    priority: 1000
    params:
        output_root=restart_output_root()
    script:
        "../scripts/coarsen_restarts.py"


rule mppncombine_diagnostics_file_run_tape_tile_segment:
    input:
        diagnostics_file_subtiles
    output:
        mppnccombined_diagnostics_file()
    wildcard_constraints:
        # Ensure snakemake does not combine segment and tape wildcards.
        segment="\d+.\d+"
    params:
        time="00:30:00"
    shell:
        """
        module load intel/2021.1.2
        module load openmpi/intel-2021.1/4.1.2
        module load netcdf/intel-2021.1/hdf5-1.10.6/4.7.4
        module load hdf5/intel-2021.1/1.10.6
        mppnccombine {output} {input}
        """


rule initialize_combine_and_coarsen_tape:
    input:
        files=combine_tape_input
    output:
        store=diagnostics_store()
    priority: 1000
    params:
        step="initialize",
	time="15:00:00"
    script:
        "../scripts/write_coarse_store_via_xpartition.py"


rule write_partition_combine_and_coarsen_tape:
    input:
        files=combine_tape_input,
	store=rules.initialize_combine_and_coarsen_tape.output.store
    output:
        xpartition_coarse_sentinel()
    params:
        ranks=XPARTITION_RANKS,
        rank="{rank}",
	step="write",
	time="15:00:00"
    wildcard_constraints:
        rank="\d+"
    script:
        "../scripts/write_coarse_store_via_xpartition.py"


rule combine_and_coarsen_tape:
    input:
        files=combine_tape_input,
        sentinels=expand(
            rules.write_partition_combine_and_coarsen_tape.output,
            rank=xpartition_rank_labels(XPARTITION_RANKS),
            run="{run}",
            tape="{tape}",
        )
    output:
        xpartition_coarse_sentinel(complete=True)
