import multiprocessing
import os
import subprocess


def strip_coarse(file):
    return file.replace("_coarse", "")


files = snakemake.input.files
gcs_url = snakemake.params["gcs_url"]
timestamp = snakemake.wildcards["timestamp"]


calls = []
for file in files:
    call = ["gsutil", "cp"]
    root, filename = os.path.split(file)
    if timestamp in filename:
        new_filename = f"{strip_coarse(filename)}"
    else:
        new_filename = f"{timestamp}.{strip_coarse(filename)}"
    destination = os.path.join(gcs_url, new_filename)
    call.append(file)
    call.append(destination)
    calls.append(call)

with multiprocessing.Pool(6) as pool:
    pool.map(subprocess.call, calls)
