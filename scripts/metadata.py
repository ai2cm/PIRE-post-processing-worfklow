import os
import typing

import xarray as xr


X_DIM = "grid_xt_coarse"
Y_DIM = "grid_yt_coarse"
HORIZONTAL_DIMS = [X_DIM, Y_DIM]


def remove_coarse_suffix(ds: xr.Dataset) -> xr.Dataset:
    """Remove '_coarse' suffix from all names in an Dataset"""
    rename = {}
    for var in ds.variables:
        rename[var] = var.replace("_coarse", "")
    return ds.rename(rename)


def validate_tile_list(tiles: typing.List[str]) -> None:
    """Validate that a list contains valid tile files"""
    assert len(tiles) == 6, "tiles must contain exactly 6 files"
    assert all(
        filename.endswith(f".tile{tile}.nc")
        for tile, filename in enumerate(tiles, start=1)
    ), "tiles must be an ordered list of files ending in '.tile[1-6].nc'"


def tile_list_to_stem(tiles: typing.List[str]) -> str:
    """Convert a list of tile files to vcm.open_tiles stem"""
    validate_tile_list(tiles)
    sample, *_ = tiles
    return sample.replace(".tile1.nc", "")


def get_segment_directories(files: typing.List[str]) -> typing.List[str]:
    """Return a sorted list of directories containing the files"""
    directories = {os.path.dirname(file) for file in files}
    return sorted(list(directories))


def is_coarsenable(da: xr.DataArray) -> bool:
    """Return whether a DataArray is coarsenable"""
    return all(dim in da.dims for dim in HORIZONTAL_DIMS)


def split_data_vars_on_coarsenability(ds: xr.Dataset) -> typing.Tuple[xr.Dataset, xr.Dataset]:
    """Return coarsenable and uncoarsenable components of a Dataset"""
    coarsenable = [var for var, da in ds.data_vars.items() if is_coarsenable(da)]
    uncoarsenable = [var for var in ds.data_vars if var not in coarsenable]
    return ds[coarsenable], ds[uncoarsenable]
