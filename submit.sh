#!/bin/bash
set -e

# conda environment: 2023-09-18-X-SHiELD-snakemake

# The multi-year X-SHiELD simulations were run in two stages, and so the raw
# data exists in two places.  Before running the snakemake workflow it is
# therefore important to merge these datasets together through symbolic links.
python misc/symlink_data.py

# We must run some additional scripts to work around diagnostic data that needed
# to be restored from ppan:/archive for the plus-4K simulation.  This data was
# already run through mppnccombine, and therefore the mppnccombine steps for
# just these segments can be skipped.  These scripts provide a crude way for
# snakemake to recognize this.
python misc/create_dummy_uncombined_plus_4K_data.py
python misc/copy_plus_4K_diagnostics.py

# Run this workflow in stages, since the graph will be large.  Also group
# mppnccombine jobs within batch jobs since there will be a huge number of them.
n_batches=64

for batch in $(seq 1 ${n_batches})
do
    snakemake \
       --groups mppncombine_restart_file_run_tape_tile_timestamp=mppnccombine coarsen_restarts_run_timestamp=coarsen \
       --group-components mppnccombine=576 coarsen=5 \
       --batch combine_and_coarsen_restarts=${batch}/${n_batches} \
       --cluster "sbatch --time=00:30:00 --output=slurm-logs/restarts-%j.out" \
       --jobs 10 \
       combine_and_coarsen_restarts
    rm -rf .snakemake  # Clean up .snakemake directory to prevent the accumulation of files
done

n_batches=4
for batch in $(seq 1 ${n_batches})
do
    snakemake \
       --groups mppncombine_diagnostics_file_run_tape_tile_segment=mppnccombine \
       --group-components mppnccombine=100 \
       --batch combine_and_coarsen_diagnostics=${batch}/${n_batches} \
       --cluster "sbatch --time=15:00:00 --output=slurm-logs/diagnostics-%j.out" \
       --jobs 10 \
       combine_and_coarsen_diagnostics
    rm -rf .snakemake  # Clean up .snakemake directory to prevent the accumulation of files
done

# We will hold off on transferring until everything is done and looks good.
# snakemake -j1 transfer &> transfer-log.out
