#!/usr/bin/env python3
# /// codebase-meta
# exempt: >-
#   psm.tsv reading for the spectrum section and export_mgf, split out of report_spectra when
#   ruff format pushed it past the 600-line ratchet (ADR-12). Reads only the files it is given.
# ///
"""FragPipe psm.tsv files and FASTA sequences, read for the PSMs a page or an export needs."""

import os
import re

from report_accession import fasta_accession, is_contaminant_row, is_decoy_row

MAX_PSMS_PER_PROTEIN = 8  # spectra you can page through for one protein
_RANK_FILE = re.compile(r"_rank(\d+)\.pep\.xml$")
_SPECTRUM = re.compile(r"^(?P<run>.+)\.(?P<scan>\d+)\.(?P<scan2>\d+)\.(?P<z>\d+)$")


def psm_rank(row):
    """A psm.tsv row's DDA+ rank, or None when its Spectrum File names none, as under DDA.

    The rank is in the file name and nowhere else in psm.tsv. FragPipe 24 converts Percolator's
    output per rank, to interact-<base>_rankN.pep.xml, for any file whose data type is not
    exactly DDA (CmdPercolator.java line 227, PercolatorOutputToPepXML.java lines 403 to 406),
    and Philosopher copies that path into Spectrum File. Every reader of a rank calls this.
    """
    found = _RANK_FILE.search(row.get("Spectrum File") or "")
    return int(found.group(1)) if found else None


def spectrum_parts(row):
    """(run, scan, charge text) of the row's Spectrum, <run>.<scan>.<scan>.<charge>, or None.

    (run, scan) is one spectrum whatever charge a PSM gave it: DDA+ can report one scan under
    two charges, and those are two Spectrum strings.
    """
    found = _SPECTRUM.match(row.get("Spectrum") or "")
    return (found["run"], int(found["scan"]), found["z"]) if found else None


def _psm_rows(path, ragged=None):
    """Yield psm.tsv rows as dicts.

    A non-blank line whose field count differs from the header, such as the last
    line of a file cut short, is skipped. Pass a list as `ragged` to receive its
    line number, so a caller that accounts for every row can count it.
    """
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            yield from psm_rows_in(fh, ragged)
    except OSError:
        return


def psm_rows_in(fh, ragged=None):
    """_psm_rows() of a psm.tsv already open as text, for a reader that hashes what it reads."""
    hdr = fh.readline().rstrip("\n").split("\t")
    for number, line in enumerate(fh, start=2):
        f = line.rstrip("\n").split("\t")
        if len(f) == len(hdr):
            yield dict(zip(hdr, f, strict=True))
        elif ragged is not None and line.strip():
            ragged.append(number)


def _num(row, key):
    try:
        return float(row.get(key, ""))
    except (TypeError, ValueError):
        return 0.0


def fasta_sequences(path, wanted):
    """Fetch sequences for a handful of accessions. Second pass over the FASTA.

    Needed for the coverage ribbon: peptide positions are found by locating the
    peptide in the protein sequence, which works regardless of whether the
    result tables happen to carry Start/End columns.

    Decoys are skipped, and that is not cosmetic. Philosopher writes them as
    `rev_sp|P12345|NAME`, which splits to the SAME accession as the target — so
    without this the decoy silently overwrites the real sequence and every
    coverage ribbon is drawn against a reversed protein.
    """
    out, acc, buf = {}, None, []
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    if acc in wanted:
                        out[acc] = "".join(buf)
                    hdr = line[1:].strip()
                    first = hdr.split()[0] if hdr.split() else hdr
                    acc = fasta_accession(first)
                    if first.startswith(("rev_", "decoy_")):
                        acc = None
                    buf = []
                elif acc in wanted:
                    buf.append(line.strip())
        if acc in wanted:
            out[acc] = "".join(buf)
    except OSError:
        return {}
    return out


def psm_tables(workdir):
    """Paths of <workdir>/<experiment>/psm.tsv, sorted by experiment name.

    <workdir>/psm.tsv is read only when no experiment has one, the fallback
    report_proteins.find_protein_tables applies to protein.tsv.
    """
    try:
        entries = sorted(os.listdir(workdir))
    except OSError:
        return []
    paths = (os.path.join(workdir, e, "psm.tsv") for e in entries)
    found = [p for p in paths if os.path.isfile(p)]
    root = os.path.join(workdir, "psm.tsv")
    return found or ([root] if os.path.isfile(root) else [])


def dda_plus_note(ranks):
    """The unit a DDA+ run's presence and spectra are counted in, or "" when no file names a rank.

    Args:
        ranks: {psm_rank: target PSMs} over the run's psm.tsv files, the key None for a row whose
            Spectrum File names no rank.

    Returns:
        One sentence of plain text, to be escaped by the page that shows it.
    """
    named = [r for r in ranks if r is not None]
    if not named:
        return ""
    higher = sum(n for r, n in ranks.items() if r is not None and r > 1)
    return (
        f"This is a DDA+ run. FragPipe pools ranks 1 to {max(named)} of each spectrum into one "
        "PSM list, up to msfragger.output_report_topN_dda_plus ranks (5 in this launcher's "
        f"base.workflow), and this section counts that pooled list: {higher:,} of "
        f"{sum(ranks.values()):,} PSMs are rank 2 or above, a peptide co-isolated with a "
        "better-scoring one on the same spectrum."
    )


def collect_psms(workdir, wanted, ranks=None):
    """One pass over every psm.tsv, keeping several PSMs per accession.

    A protein normally has many PSMs — several distinct peptides, each re-observed
    in every sample. Keeping only the best one throws away the evidence a reader
    most wants to leaf through, so up to MAX_PSMS_PER_PROTEIN are retained, one
    per distinct peptide (the highest-scoring scan for that peptide), ordered by
    score. Decoy rows are skipped: a decoy's PSM is no evidence for any protein. So are
    contaminant rows, whose Protein ID is the bare accession the report would class.

    Pass a Counter as `ranks` to receive the psm_rank of every target row read.

    -> {accession: {"psms": [row, ...], "peptides": set()}}
    """
    per_pep = {}  # (acc, peptide) -> best row
    peptides = {}  # acc -> set of peptides
    for psm in psm_tables(workdir):
        for row in _psm_rows(psm):
            if is_decoy_row(row) or is_contaminant_row(row):
                continue
            if ranks is not None:
                ranks[psm_rank(row)] += 1
            acc = fasta_accession(row.get("Protein ID", "").strip())
            if acc not in wanted:
                continue
            pep = row.get("Peptide", "").strip()
            if not pep:
                continue
            peptides.setdefault(acc, set()).add(pep)
            key = (acc, pep)
            prev = per_pep.get(key)
            if prev is None or _num(row, "Hyperscore") > _num(prev, "Hyperscore"):
                per_pep[key] = row

    found = {}
    for (acc, _pep), row in per_pep.items():
        found.setdefault(acc, []).append(row)
    return {
        acc: {
            "psms": sorted(rows, key=lambda r: -_num(r, "Hyperscore"))[:MAX_PSMS_PER_PROTEIN],
            "peptides": peptides.get(acc, set()),
        }
        for acc, rows in found.items()
    }
