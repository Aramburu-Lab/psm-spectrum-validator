#!/usr/bin/env python3
# /// codebase-meta
# exempt: >-
#   Where a launcher run keeps its files: the run directory beside its results directory, and the
#   manifest naming its inputs. Split out of run_report so that export_mgf reads the same rule
#   without the report. Reads only the manifest it is given; writes nothing.
# ///
"""A launcher run's layout: <proj>/runs/<name> beside <proj>/results/<name>, and its manifest.

The manifest in the run directory names each input file with its experiment, bioreplicate and
data type. run_report reads it for the page's inputs, and export_mgf for where the spectra are.
"""

import os

MANIFEST = "fragpipe-files.fp-manifest"
MANIFEST_TYPE_COL = 3  # path, experiment, bioreplicate, data type


def _path_size(p):
    """Size of a file, or of a whole directory for Bruker .d inputs."""
    try:
        if os.path.isdir(p):
            return sum(os.path.getsize(os.path.join(dp, n)) for dp, _, ns in os.walk(p) for n in ns)
        return os.path.getsize(p)
    except OSError:
        return None


def load_manifest(path):
    """Read the manifest into dicts of path/name/size/experiment/type."""
    rows = []
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                rows.append(
                    {
                        "path": f[0],
                        "name": os.path.basename(f[0]),
                        "size": _path_size(f[0]),
                        "experiment": f[1] if len(f) > 1 else "",
                        "type": f[MANIFEST_TYPE_COL] if len(f) > MANIFEST_TYPE_COL else "",
                    }
                )
    except OSError:
        pass
    return rows


def input_dirs(files):
    """Directories the run's input files sit in, in manifest order.

    MSFragger writes *_uncalibrated.mzML beside a raw input, not into the results
    tree, and an mzML input sits there itself, so anything that reads spectra looks
    where the inputs came from. export_mgf uses this too, so the two can never
    disagree on where a scan lives.
    """
    dirs = []
    for f in files:
        d = os.path.dirname(f["path"])
        if d and d not in dirs:
            dirs.append(d)
    return dirs


def results_dir_for(run_dir, ru):
    """Where FragPipe wrote its tables for this run.

    The launcher records results_dir in the resource log; runs from before that
    are derived from its layout, <proj>/runs/<name> beside <proj>/results/<name>.
    Deriving it from the run directory is the point — reading a config file
    instead would resolve to whatever project that config currently names, which
    is not necessarily the project this run belongs to.
    """
    if ru.get("results_dir"):
        return ru["results_dir"]
    parent, name = os.path.split(os.path.abspath(run_dir).rstrip(os.sep))
    proj, runs = os.path.split(parent)
    cand = os.path.join(proj, "results", name)
    return cand if runs == "runs" and os.path.isdir(cand) else ""


def run_dir_for(work_dir):
    """The run directory behind a results directory: results_dir_for read backwards.

    Same layout, same reason: the manifest of THIS run is beside its results,
    whatever project a config file names today.
    """
    parent, name = os.path.split(os.path.abspath(work_dir).rstrip(os.sep))
    proj, results = os.path.split(parent)
    cand = os.path.join(proj, "runs", name)
    return cand if results == "results" and os.path.isdir(cand) else ""
