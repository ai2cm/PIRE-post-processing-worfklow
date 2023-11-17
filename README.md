# X-SHiELD simulation post-processing on Stellar

This repository defines a `snakemake` workflow for post-processing data from
four two-plus-year C3072 resolution simulations completed on Princeton's Stellar
computer.

## Installing `snakemake` and other dependencies

`snakemake` requires many dependencies, so trying to build an environment with
plain `conda` does not always work.  [Per the `snakemake`
documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html),
it is therefore recommended to install `mamba`, which has a more advanced
dependency solver.  Once `mamba` is installed, one can create the
post-processing environment using:

```
$ mamba env create --file envs/environment.yaml
```

## Processing the data

The multi-year X-SHiELD simulations were run in two stages, and so the raw data
exists in two places.  Before running the `snakemake` workflow it is therefore
important to merge these datasets together through symbolic links.  This can be
done by calling the included script:

```
$ conda run --name 2023-09-18-X-SHiELD-snakemake misc/symlink_dataset.py
```

To process the data, activate the environment from a `screen` session (this will
take a while), and call the included top-level `submit.sh` bash script.  This
script handles partitioning the work into a sequence of batch jobs, grouping
tasks into single jobs where appropriate to prevent overwhelming the SLURM
scheduler with many small jobs.

```
$ screen
$ conda run --name 2023-09-18-X-SHiELD-snakemake submit.sh
```

## High-level overview of the workflow

This workflow is geared to produce data compatible with AI2's corrective machine
learning workflow.  At a high level it does the following:

- Combines the subtiles of the raw diagnostic and restart netCDF files output
  from the simulation into cohesive tiles using GFDL's `mppnccombine` tool,
  since even the coarse data was output with a 2x2 I/O layout to ease I/O
  overhead in the simulations.
- Concatenates the diagnostics datasets along the time and tile dimensions, and
  coarsens the partially coarsened C384 output to C48 resolution before dumping
  out to zarr.
- Coarsens each set of C384 restart files into its own subdirectory labeled with
  the timestamp of the form `"%Y%m%d.%H%M%S"` using AI2's default pressure-level
  coarse-graining strategy, and specialized coarsening of land surface fields.
  This arrangement and naming of the restart files is exactly what is required
  for performing an fv3net nudged run.

In total this workflow processes over 140 TB of restart files and 7.3 TB of
diagnostics (we have ignored the 3D diagnostics for now, though they too could
be processed by this workflow).  By way of coarsening to C48 resolution, this
data is reduced by a factor of 64 to a more manageable ~2 TB.
