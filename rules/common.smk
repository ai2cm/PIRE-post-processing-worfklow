import cftime
import numpy as np
import xarray as xr

from glob import glob
from pathlib import Path


CONFIG = config["config_name"]
GCS_BUCKET = config["gcs_bucket"]
OUTPUT_DIRECTORY = Path(config["output_directory"])
WORKING_DIRECTORY = Path(config["working_directory"])
XPARTITION_RANKS = config["xpartition_ranks"]


RUNS = config["runs"].keys()
SEGMENT_FORMAT = "%Y%m%d%H"
TIMESTAMP_FORMAT = "%Y%m%d.%H%M%S"
TILES = range(1, 7)


def get_run_directory(run):
    return Path(config["runs"][run]["root"])


def get_tapes(run):
    return config["runs"][run].get("tapes")


def xpartition_rank_labels(ranks):
    return [f"{rank:04d}" for rank in range(ranks)]


# cftime version 1.6.2 includes an implementation of strptime, which
# would be nice to use here...
def segment_label_to_datetime(label):
    year = int(label[:4])
    month = int(label[4:6])
    day = int(label[6:8])
    hour = int(label[8:10])
    return cftime.DatetimeJulian(year, month, day, hour)


def timestamp_label_to_datetime(label):
    year = int(label[:4])
    month = int(label[4:6])
    day = int(label[6:8])
    hour = int(label[9:11])
    minute = int(label[11:13])
    second = int(label[13:15])
    return cftime.DatetimeJulian(year, month, day, hour, minute, second)


def segment_label_to_timestamp_label(label):
    datetime = segment_label_to_datetime(label)
    return datetime.strftime(TIMESTAMP_FORMAT)


##### Metadata #####

SEGMENTS = {}
TAPES = {}
for run in RUNS:
    run_directory = get_run_directory(run)
    pattern = run_directory / "history" / "{segment}" / "{tape}.tile1.nc.0000"
    segments, tapes = glob_wildcards(pattern, followlinks=True)
    SEGMENTS[run] = sorted(list(set(segments)))
    specified_tapes = get_tapes(run)
    if specified_tapes is None:
        TAPES[run] = sorted(list(set(tapes)))  # Use globbed tapes
    else:
        TAPES[run] = specified_tapes


RESTART_SEGMENTS = {}
for run in RUNS:
    run_directory = get_run_directory(run)
    pattern = run_directory / "restart" / "{restart_segment}" / "fv_core.res.nc"
    segments, = glob_wildcards(pattern, followlinks=True)
    RESTART_SEGMENTS[run] = sorted(list(set(segments)))


COARSE_RESTART_TIMESTAMPS = {}
for run in RUNS:
    run_directory = get_run_directory(run)
    times = []
    for segment in RESTART_SEGMENTS[run]:
        pattern = run_directory / "restart" / segment / "{timestamp}.coupler.res"
        segment_times, = glob_wildcards(pattern, followlinks=True)
        segment_timestamp = segment_label_to_timestamp_label(segment)
        segment_times.append(segment_timestamp)
        times.extend(segment_times)
    COARSE_RESTART_TIMESTAMPS[run] = sorted(times)


RESTART_SEGMENTS_AS_TIMESTAMPS = {}
for run, segments in RESTART_SEGMENTS.items():
    RESTART_SEGMENTS_AS_TIMESTAMPS[run] = list(map(
        segment_label_to_timestamp_label, segments
    ))


##### Restarts #####

def get_restart_segment_directory(run, timestamp):
    run_directory = get_run_directory(run)

    segment_timestamps = RESTART_SEGMENTS_AS_TIMESTAMPS[run]
    segment_index = np.searchsorted(segment_timestamps, timestamp, side="left")
    segment_timestamp = segment_timestamps[segment_index]
    segment = RESTART_SEGMENTS[run][segment_index]

    return run_directory / "restart" / segment, segment_timestamp
    

def restart_file_subtiles(wildcards):
    run = wildcards.run
    tape = wildcards.tape
    tile = wildcards.tile
    timestamp = wildcards.timestamp

    root, segment_timestamp = get_restart_segment_directory(run, timestamp)

    if timestamp == segment_timestamp:
        pattern = root / f"{tape}.tile{tile}.nc*"
    else:
        pattern = root / f"{timestamp}.{tape}.tile{tile}.nc*"

    return sorted(glob(str(pattern)))


def get_coupler_res_for_timestamp(wildcards):
    run = wildcards.run
    timestamp = wildcards.timestamp

    root, segment_timestamp = get_restart_segment_directory(run, timestamp)
    if timestamp == segment_timestamp:
        return root / "coupler.res"
    else:
        return root / f"{timestamp}.coupler.res"


def combine_tape_input(wildcards):
    segments = SEGMENTS[wildcards.run]
    root = WORKING_DIRECTORY / "{{run}}" / "diagnostics" / "{segment}"
    filename = "{{tape}}.tile{tile}.nc"
    pattern = root / filename
    return expand(pattern, segment=segments, tile=TILES)


def mppnccombined_restart_file():
    root = WORKING_DIRECTORY / "{run}" / "restarts" / "{timestamp}"
    filename = "{timestamp}.{tape}.tile{tile}.nc"
    return temporary(root / filename)


def _tiled_restarts_as_inputs(tape, from_mppnccombine, wildcards):
    run = wildcards.run
    timestamp = wildcards.timestamp

    run_directory = get_run_directory(run)

    if from_mppnccombine:
        root = WORKING_DIRECTORY / run / "restarts" / timestamp
	filename = f"{timestamp}.{{tape}}.tile{{tile}}.nc"
    else:
        get_restart_segment_directory(run, timestamp)
	root, segment_timestamp = get_restart_segment_directory(run, timestamp)
	if timestamp == segment_timestamp:
            filename = "{tape}.tile{tile}.nc"
        else:
            filename = f"{timestamp}.{{tape}}.tile{{tile}}.nc"

    pattern = root / filename
    return expand(pattern, tape=tape, tile=TILES)


def restarts_as_inputs(tape, from_mppnccombine):
    def func(wildcards):
        if tape == "coupler.res":
            return get_coupler_res_for_timestamp(wildcards)
        else:
            return _tiled_restarts_as_inputs(tape, from_mppnccombine, wildcards)
    return func
    

def _tiled_restarts_as_outputs(tape):
    root = OUTPUT_DIRECTORY / "{{run}}" / "restarts" / "{{timestamp}}"
    filename = "{{timestamp}}.{tape}.tile{tile}.nc"
    pattern = root / filename
    return expand(pattern, tape=tape, tile=TILES)


def _coupler_res_as_output():
    root = OUTPUT_DIRECTORY / "{run}" / "restarts" / "{timestamp}"
    filename = "{timestamp}.coupler.res"
    return root / filename


def restarts_as_outputs(tape):
    if tape == "coupler.res":
        return _coupler_res_as_output()
    else:
        return _tiled_restarts_as_outputs(tape)


def restart_output_root():
    path = OUTPUT_DIRECTORY / "{run}" / "restarts" / "{timestamp}"
    return str(path)


##### Diagnostics #####

def diagnostics_file_subtiles(wildcards):
    run_directory = get_run_directory(wildcards.run)
    root = run_directory / "history" / wildcards.segment
    pattern = root / f"{wildcards.tape}.tile{wildcards.tile}.nc*"
    return sorted(glob(str(pattern)))


def mppnccombined_diagnostics_file():
    root = WORKING_DIRECTORY / "{run}" / "diagnostics" / "{segment}"
    filename = "{tape}.tile{tile}.nc"
    return root / filename


def diagnostics_store():
    root = OUTPUT_DIRECTORY / "{run}" / "diagnostics"
    store = "{tape}.zarr"
    return directory(root / store)


def xpartition_coarse_sentinel(complete=False):
    root = Path("xpartition-coarse-sentinels") / CONFIG / "{run}"
    if complete:
        filename = "{tape}.complete.out"
    else:
        filename = "{tape}.{rank}.out"
    return touch(root / filename)


##### Transfer #####

def transfer_sentinel():
    path = Path("transfer-sentinels") / f"{CONFIG}.out"
    return touch(path)


def transfer_source():
    path = OUTPUT_DIRECTORY / "*"
    return str(path)


def transfer_gcs_path():
    return f"{GCS_BUCKET}/"
