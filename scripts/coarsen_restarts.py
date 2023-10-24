import shutil

import dask.diagnostics
import vcm
import vcm.cubedsphere
import xarray as xr

from pathlib import Path
from metadata import remove_coarse_suffix, tile_list_to_stem


coarsening_factor = snakemake.config["coarsening_factor"]
grid_spec_stem = snakemake.config["grid_spec_stem"]
output_root = Path(snakemake.params.output_root)
tiles = range(1, 7)
timestamp = snakemake.wildcards.timestamp

# Copy coupler.res file to the desired destination
source_coupler_res = snakemake.input.coupler_res
destination_coupler_res = snakemake.output.coupler_res
shutil.copy(source_coupler_res, destination_coupler_res)

# Coarsen the restart files and save them to the desired destination
fv_core_stem = tile_list_to_stem(snakemake.input.fv_core)
fv_tracer_stem = tile_list_to_stem(snakemake.input.fv_tracer)
fv_srf_wnd_stem = tile_list_to_stem(snakemake.input.fv_srf_wnd)
sfc_data_stem = tile_list_to_stem(snakemake.input.sfc_data)

restart_data = {
    "fv_core.res": vcm.open_tiles(fv_core_stem),
    "fv_tracer.res": vcm.open_tiles(fv_tracer_stem),
    "fv_srf_wnd.res": vcm.open_tiles(fv_srf_wnd_stem),
    "sfc_data": vcm.open_tiles(sfc_data_stem)
}
grid_spec = remove_coarse_suffix(vcm.open_tiles(grid_spec_stem))
coarsened = vcm.cubedsphere.coarsen_restarts_on_pressure(
    coarsening_factor,
    grid_spec,
    restart_data,
    coarsen_agrid_winds=True,
    extrapolate=False
)

for tape, ds in coarsened.items():
    datasets = {}

    with dask.diagnostics.ProgressBar():
        ds = ds.compute()
    
    for tile in tiles:
        path = output_root / f"{timestamp}.{tape}.tile{tile}.nc"
        datasets[path] = ds.isel(tile=tile - 1).drop("tile", errors="ignore")

    xr.save_mfdataset(datasets.values(), datasets.keys())
