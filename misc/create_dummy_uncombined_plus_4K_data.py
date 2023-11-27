# #!/usr/bin/env python
import logging
import os

from pathlib import Path

ROOT = Path("/scratch/cimes/GLOBALFV3/stellar_run")
SIMULATION = "20191020.00Z.C3072.L79x2_pire_PLUS_4K"
TARGET_ROOT = Path("/scratch/gpfs/skclark/2023-09-15-X-SHiELD-symlinks")
RERUN_SEGMENTS = ["2019102000"]  # Needed to be re-run (so uncombined)
FIRST_UNRESTORED_SEGMENT = "2020022700"
TAPES = [
    "pire_atmos_dyn_3h_coarse_inst",
    "pire_atmos_dyn_plev_coarse_3h",
    "pire_atmos_phys_3h_coarse",
    "pire_atmos_static_coarse"
]
TILES = range(1, 7)
SUBTILES = [f"{subtile:04d}" for subtile in range(4)]


def create_dummy_files(segment, tape):
    root = TARGET_ROOT / SIMULATION / "history" / segment
    if not os.path.exists(root):
        os.makedirs(root)

    for tile in TILES:
        for subtile in SUBTILES:
            filename = root / f"{tape}.tile{tile}.nc.{subtile}"
            logging.info(f"Creating dummy file at {filename}")
            open(filename, "a").close()


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    history = ROOT / "transferred" / SIMULATION / "history"
    segments = sorted(os.listdir(history))

    restored_segments = []
    for segment in segments:
        if segment < FIRST_UNRESTORED_SEGMENT and segment not in RERUN_SEGMENTS:
            restored_segments.append(segment)

    for segment in restored_segments:
        for tape in TAPES:
            create_dummy_files(segment, tape)
