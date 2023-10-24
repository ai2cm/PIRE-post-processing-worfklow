import os

import dask
import dask.diagnostics
import vcm
import xarray as xr
import xpartition

from fv3dataset import HistoryDataset

from metadata import X_DIM, Y_DIM
from metadata import get_segment_directories
from metadata import split_data_vars_on_coarsenability


def open_dataset(files, tape, target_chunk_size):
    directories = get_segment_directories(files)
    return HistoryDataset(
        tape,
        directories,
        target_chunk_size=target_chunk_size
    ).to_dask()


files = snakemake.input.files
tape = snakemake.wildcards["tape"]
target_chunk_size = snakemake.config["target_chunk_size"]
coarsening_factor = snakemake.config["coarsening_factor"]
grid_spec_stem = snakemake.config["grid_spec_stem"]
step = snakemake.params["step"]

grid_spec = vcm.open_tiles(grid_spec_stem)
fine = open_dataset(
    files,
    tape,
    target_chunk_size
)

coarsenable, uncoarsenable = split_data_vars_on_coarsenability(fine)
coarsened = vcm.cubedsphere.weighted_block_average(
    coarsenable,
    grid_spec.area_coarse,
    coarsening_factor,
    x_dim=X_DIM,
    y_dim=Y_DIM,
)
result = xr.merge([coarsened, uncoarsenable])

desired_chunks = {"time": "auto", "tile": -1}
with dask.config.set({"array.chunk-size": target_chunk_size}):
    chunks = {}
    for dim, size in desired_chunks.items():
        if dim in result.dims:
            chunks[dim] = size
    result = result.chunk(chunks)

if step == "initialize":
    result.partition.initialize_store(snakemake.output.store)
elif step == "write":
    ranks = int(snakemake.params["ranks"])
    rank = int(snakemake.params["rank"])
    partition_dims = ["time"]  # TODO: could parametrize this
    with dask.diagnostics.ProgressBar():
        result.partition.write(
            snakemake.input.store,
            ranks,
            partition_dims,
            rank
    )
else:
    raise ValueError("step parameter must be either 'initialize' or 'write'")
