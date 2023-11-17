# #!/usr/bin/env python
# The simulations are split across multiple directories.  According to Kai,
# when the transferred and transferred_new simulations overlap, we should
# prioritize the transferred_new simulation over the transferred simulation.
import datetime
import logging
import os

from pathlib import Path

HISTORY = "history"
RESTART = "restart"
TRANSFERRED = "transferred"
TRANSFERRED_NEW = "transferred_new"
ROOT = Path("/scratch/cimes/GLOBALFV3/stellar_run")
SIMULATIONS = [
    "20191020.00Z.C3072.L79x2_pire",
    "20191020.00Z.C3072.L79x2_pire_PLUS_4K",
    "20191020.00Z.C3072.L79x2_pire_CO2_1270ppmv",
    "20191020.00Z.C3072.L79x2_pire_PLUS_4K_CO2_1270ppmv"
]
TARGET_ROOT = Path("/scratch/gpfs/skclark/2023-09-15-X-SHiELD-symlinks")
EXPECTED_TIMESTAMP_PATTERN = "%Y%m%d%H"


def is_timestamp_like(string):
    try:
        datetime.datetime.strptime(string, EXPECTED_TIMESTAMP_PATTERN)
    except ValueError:
        return False
    else:
        return True


def retain_timestamp_like(directories):
    return filter(is_timestamp_like, directories)


def merge_output_category(simulation, output_category):
    transferred = ROOT / TRANSFERRED / simulation
    transferred_new = ROOT / TRANSFERRED_NEW / simulation

    timestamps = {}
    timestamps[TRANSFERRED] = sorted(
        retain_timestamp_like(
            os.listdir(transferred / output_category)
        )
    )
    timestamps[TRANSFERRED_NEW] = sorted(
        retain_timestamp_like(
            os.listdir(transferred_new / output_category)
        )
    )

    links = {}
    for timestamp in timestamps[TRANSFERRED]:
        # Note the transferred simulations were potentially run with a
        # different segment length than the transferred_new simulation, so we need
        # to be careful to not rely on strict segment overlap to take care of the
        # precedence (if we do that we can end up with extra, redundant, segments
        # from the transferred simulation contaminating our symlinked data).
        if timestamp < min(timestamps[TRANSFERRED_NEW]):
            links[timestamp] = transferred / output_category / timestamp

    for timestamp in timestamps[TRANSFERRED_NEW]:
        links[timestamp] = transferred_new / output_category / timestamp

    return links
    
    
def merge(simulation):
    history = merge_output_category(simulation, HISTORY)
    restart = merge_output_category(simulation, RESTART)
    return {HISTORY: history, RESTART: restart}


def symlink(source, simulation, output_category, timestamp):
    target_root = TARGET_ROOT / simulation / output_category
    destination = target_root / timestamp
    if not os.path.exists(target_root):
        os.makedirs(target_root)
    logging.info(f"Creating a symlink from {source} to {destination}")
    os.symlink(source, destination)

    
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    merged_links = {}
    for simulation in SIMULATIONS:
        merged_links[simulation] = merge(simulation)

    for simulation, output_categories in merged_links.items():
        for output_category, links in output_categories.items():
            for timestamp, source in links.items():
                symlink(source, simulation, output_category, timestamp)
