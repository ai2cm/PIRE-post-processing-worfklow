# #!/usr/bin/env python
# There was missing data in the plus-4K case.  Kai copied over data that
# was already combined for the missing segments.  This does not fit in
# neatly with our snakemake workflow, which expects uncombined raw data,
# so we'll need to manually symlink the pre-combined data into the work
# directory for those segments.
import glob
import logging
import os
import shutil

from pathlib import Path


TRANSFERRED = "transferred"
ROOT = Path("/scratch/cimes/GLOBALFV3/stellar_run")
SIMULATION = "20191020.00Z.C3072.L79x2_pire_PLUS_4K"
SOURCE_ROOT = ROOT / TRANSFERRED / SIMULATION / "history"
TARGET_ROOT = Path("/scratch/cimes/skclark/2023-09-18-PIRE-X-SHiELD-post-processing/work/plus-4K/diagnostics")
SEGMENTS = sorted(os.listdir(SOURCE_ROOT))
RERUN_SEGMENTS = ["2019102000"]  # Needed to be re-run (so uncombined)
FIRST_UNRESTORED_SEGMENT = "2020022700"
TAPES = [
    "pire_atmos_dyn_3h_coarse_inst",
    "pire_atmos_dyn_plev_coarse_3h",
    "pire_atmos_phys_3h_coarse",
    "pire_atmos_static_coarse"
]


def copy(segment):
    source = SOURCE_ROOT / segment
    destination = TARGET_ROOT / segment
    logging.info(f"Copying data from {source} to {destination}")
    os.makedirs(destination, exist_ok=True)
    for tape in TAPES:
        files = source.glob(f"{tape}.tile*")
        for file in files:
            shutil.copy(file, destination)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    restored_segments = []
    for segment in SEGMENTS:
        if segment < FIRST_UNRESTORED_SEGMENT and segment not in RERUN_SEGMENTS:
            restored_segments.append(segment)

    for segment in restored_segments:
        copy(segment)
