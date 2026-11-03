#!/usr/bin/env python3
# /// codebase-meta
# exempt: >-
#   export_mgf's self-test, split out of it when it reached the 600-line ratchet (ADR-12). Run by
#   export_mgf --selftest, which hands it the module under test; reads nothing real and writes
#   only into a temporary directory.
# ///
"""export_mgf end to end on the planted run, then the DDA+ rank cases the audit found (ACC-2).

Scan 107 of the rank cases holds a canonical MDKNELVQK at rank 1 and a novel IEEDGEVR at rank 2,
the Spectrum File naming each rank as FragPipe 24 writes it under DDA+. Before the rank rule the
documented novel --match wrote the rank 2 alone as the spectrum's annotation, and with no --match
the two charges of one scan were written as two records.
"""

import contextlib
import io
import os
import re
import tempfile
import warnings
from collections import Counter

import export_mgf_fixture as fx

NOVEL_MATCH = r"label=novel\b"
_CANON = ("rs|GENE04300|GENE04300_ORGN", "label=canonical class=annotated_cds")
_NOVEL = ("rs|XORF01924|NOVEL_ORF", "label=novel class=uoORF")
_RANK_ROWS = {  # Spectrum, charge, peptide, mods, Observed M/Z, mass, Calculated M/Z, score, entry
    "r1": ("runA.00107.00107.2", "2", "MDKNELVQK", "N-term(42.0106), 1M(15.9949)", "581.7931",
           "1161.5700", "581.7923", "30.5000", _CANON, "interact-runA_rank1.pep.xml"),
    "r2": ("runA.00107.00107.2", "2", "IEEDGEVR", "", "473.7290", "945.4403", "473.7274",
           "9.0000", _NOVEL, "interact-runA_rank2.pep.xml"),
    "r2z3": ("runA.00107.00107.3", "3", "IEEDGEVR", "", "316.1550", "945.4403", "316.1540",
             "9.0000", _NOVEL, "interact-runA_rank2.pep.xml"),
}  # fmt: skip


def _rank_tree(tmp, name, keys):
    """The fixture run with its psm.tsv replaced by the rank case's rows."""
    work, data = fx.write_tree(os.path.join(tmp, name), "fragpipe-files.fp-manifest")
    lines = ["\t".join(fx.PSM_HEADER)]
    for key in keys:
        spec, z, pep, mods, omz, mass, cmz, score, entry, sfile = _RANK_ROWS[key]
        row = dict.fromkeys(fx.PSM_HEADER, "")
        row.update({"Spectrum": spec, "Spectrum File": sfile, "Peptide": pep, "Charge": z,
                    "Retention": "600.0", "Observed M/Z": omz, "Calculated Peptide Mass": mass,
                    "Calculated M/Z": cmz, "Hyperscore": score, "Probability": "0.9990",
                    "Assigned Modifications": mods, "Protein": entry[0],
                    "Protein ID": f"{entry[0]} {entry[1]} OS=Example organism"})  # fmt: skip
        lines.append("\t".join(row[c] for c in fx.PSM_HEADER))
    with open(os.path.join(work, "expA", "psm.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return work, data


def _records(path):
    """[{KEY: value}] of each BEGIN IONS block, the peak lines left out."""
    with open(path, encoding="utf-8") as fh:
        blocks = fh.read().split("BEGIN IONS\n")[1:]
    return [dict(ln.split("=", 1) for ln in b.splitlines() if "=" in ln) for b in blocks]


def check_rank_cases(em, tmp):
    """Both of the audit's cases, by default and under --all-ranks."""
    novel = re.compile(NOVEL_MATCH)
    out = os.path.join(tmp, "rank.mgf")
    work, data = _rank_tree(tmp, "same_charge", ("r1", "r2"))
    got, _ = em.export(work, [data], novel, out)
    # By default the rank 2 is a named drop, never the spectrum's annotation.
    assert (got[em.RANKED], got[em.NO_MATCH], got[em.WRITTEN]) == (1, 1, 0), got
    got, _ = em.export(work, [data], novel, out, em.Limits(all_ranks=True))
    assert got[em.WRITTEN] == 1 and got[em.RANKED] == 0, got
    # The scan's best is recorded before --match, so the canonical rank 1 still names it.
    (rec,) = _records(out)
    assert (rec["SEQ"], rec["RANK"]) == ("IEEDGEVR", "2"), rec
    assert (rec["SCAN_BEST_PEPTIDE"], rec["SCAN_BEST_HYPERSCORE"]) == ("MDKNELVQK", "30.5000")
    work, data = _rank_tree(tmp, "other_charge", ("r1", "r2z3"))
    got, _ = em.export(work, [data], novel, out, em.Limits(all_ranks=True))
    (rec,) = _records(out)
    assert (rec["TITLE"], rec["RANK"], rec["SCAN_BEST_PEPTIDE"]) == (
        "runA.00107.00107.3",
        "2",
        "MDKNELVQK",
    ), rec
    # Taken spectra are keyed on (run, scan): two charges of scan 107 are one spectrum.
    got, _ = em.export(work, [data], None, out, em.Limits(all_ranks=True))
    assert (got[em.WRITTEN], got[em.SHARED_SCAN]) == (1, 1), got
    assert [r["SEQ"] for r in _records(out)] == ["[+42.0106]-M[+15.9949]DKNELVQK"]
    got, _ = em.export(work, [data], None, out)
    assert (got[em.WRITTEN], got[em.RANKED], got[em.SHARED_SCAN]) == (1, 1, 0), got
    assert "RANK" not in _records(out)[0]


def check_planted(em, tmp):
    """The planted run: every count by default, then the same rows under --all-ranks."""
    # The "=" puts a workdir path in the header that pyteomics would read as a param.
    work, data = fx.write_tree(os.path.join(tmp, "batch=7"), em.MANIFEST)
    fx.check_index(os.path.join(data, "runA_uncalibrated.mzML"))
    dirs, _src = em.resolve_data_dirs(work, [])
    assert dirs == [data], dirs
    out = os.path.join(tmp, "x.mgf")
    expect = {em.READ: 14, em.DECOY: 1, em.RAGGED: 1, em.RANKED: 2, em.NO_MATCH: 1,
              em.NO_CALC: 1, em.MASS: 1, em.UNPARSED: 0, em.NO_PEPMASS: 1, em.NO_MZML: 1,
              em.MAPPED_ONLY: 2, em.NOT_BEST: 0, em.SHARED_SCAN: 0, em.BEYOND: 0, em.CAP: 0,
              em.UNREADABLE: 1, em.FALLBACK: 1, em.NO_RT: 1, em.WRITTEN: 4}  # fmt: skip
    for limits, ranked, shared, extra in (
        (em.Limits(), 2, 0, frozenset()),
        # Under --all-ranks the two rank 2 picks are on spectra a better peptide took.
        (em.Limits(all_ranks=True), 0, 2, em.RANK_KEYS),
    ):
        counts, runs = em.export(work, dirs, re.compile("XORF"), out, limits)
        want = {**expect, em.RANKED: ranked, em.SHARED_SCAN: shared}
        assert {k: counts[k] for k in want} == want, counts
        assert counts[em.READ] == sum(counts[k] for k in em.DROPS) + counts[em.WRITTEN], counts
        assert runs == {"runA": os.path.join(data, "runA_uncalibrated.mzML"), "runB": None}, runs
        fx.check_mgf_records(out, counts, work, extra)
    # VSTEVDAR's rank 2 is beyond the limit here, as VADLTNSFEK, dropped by it, takes no
    # spectrum; XORF01924's IEEDGEVR, beyond it too, counts as on scan 101, the first test.
    tight, _ = em.export(work, dirs, re.compile("XORF"), out, em.Limits(1, 1, True))
    got = tight[em.SHARED_SCAN], tight[em.BEYOND], tight[em.CAP], tight[em.WRITTEN]
    assert got == (1, 3, 2, 1), tight
    tight, _ = em.export(work, dirs, re.compile("XORF"), out, em.Limits(1, 1))
    got = tight[em.RANKED], tight[em.SHARED_SCAN], tight[em.BEYOND], tight[em.CAP]
    assert got == (2, 0, 2, 2), tight
    return work, dirs, out


def check_read_order(em):
    """Both halves of the order: reads grouped by file, records back in Hyperscore order.

    The fixture has one readable mzML, so the end-to-end export cannot tell the two apart.
    """
    keep = em.read_spectrum, em.mgf_record
    reads = []
    em.read_spectrum = lambda mzml, scan: reads.append(mzml) or ([scan], [1.0])
    em.mgf_record = lambda row, *_: row
    try:
        files = ["b", "a", "c", "a", "b", "c"]  # Hyperscore order, interleaving three files
        picks = [(k, (0,) * 5 + (k,), f, None) for k, f in enumerate(files)]
        got = em.read_records(picks, Counter(), lambda _pick: [])
    finally:
        em.read_spectrum, em.mgf_record = keep
    assert reads == sorted(files), reads
    assert got == list(range(len(files))), got


def check_spectra_dir(em, tmp):
    """--spectra-dir on a run with no launcher run directory writes what the manifest gives.

    Returns:
        The workdir of that run, which has no manifest.
    """
    records = []
    for name in ("launcher", "bare"):
        work, data = fx.write_tree(os.path.join(tmp, name), em.MANIFEST)
        run_dir = os.path.join(tmp, name, "proj", "runs", "fp")
        out = os.path.join(tmp, f"{name}.mgf")
        argv, source = ["--workdir", work, "--match", "XORF", "--out", out], run_dir
        if name == "bare":
            os.remove(os.path.join(run_dir, em.MANIFEST))
            os.rmdir(run_dir)
            argv, source = [*argv, "--spectra-dir", data], "--spectra-dir"
        with contextlib.redirect_stdout(io.StringIO()) as said:
            code = em.main(argv)
        assert code == 0 and f"  [{source}" in said.getvalue(), said.getvalue()
        records.append(_records(out))
    assert [len(r) for r in records] == [4, 4] and records[1] == records[0], records
    assert em._parser().parse_args(["--mzml-dir", data]).spectra_dir == [data]  # the old name
    return work


def check_mzml_names(em, tmp):
    """<run>.mzML is read where there is no <run>_uncalibrated.mzML, and loses to it where both are.

    A search of mzML files leaves <run>.mzML; a .raw run has MSFragger's uncalibrated copy, and a
    <run>.mzML beside that, here the same scans at twice the intensity, must not be read instead.
    """
    bodies = {}
    for name in ("uncalibrated", "plain", "both"):
        work, data = fx.write_tree(os.path.join(tmp, name), em.MANIFEST)
        uncal, plain = os.path.join(data, "runA_uncalibrated.mzML"), os.path.join(data, "runA.mzML")
        if name == "plain":
            os.rename(uncal, plain)
        elif name == "both":
            em.write_mzml(plain, [(s, mz, [2 * i for i in inten]) for s, mz, inten in fx.SCANS])
        out = os.path.join(tmp, f"{name}.mgf")
        _counts, runs = em.export(work, [data], re.compile("XORF"), out)
        read = plain if name == "plain" else uncal
        with open(out, encoding="utf-8") as fh:
            head, bodies[name] = fh.read().split("\n\n", 1)
        assert runs["runA"] == read and em._comment(f"run runA: {read}") in head, (name, runs)
    assert [b.count("BEGIN IONS") for b in bodies.values()] == [4, 4, 4], bodies
    assert bodies["plain"] == bodies["uncalibrated"] == bodies["both"], bodies


def check_refusals(em, tmp, work, bare):
    """Bad paths, and spectra that neither a manifest nor --spectra-dir finds, stop main.

    bare is the planted run with no launcher run directory, so with no manifest.
    """
    kept = os.path.join(tmp, "kept.mgf")
    with open(kept, "w", encoding="utf-8") as fh:
        fh.write("an earlier export\n")
    nope, empty = os.path.join(tmp, "nope"), os.path.join(tmp, "empty")
    os.mkdir(empty)
    for args, why in (
        (["--workdir", nope, "--spectra-dir", tmp], f"--workdir {nope} is not a directory"),
        (["--workdir", tmp, "--spectra-dir", tmp], "no psm.tsv in"),
        (["--workdir", work, "--out", os.path.join(tmp, "nodir", "x.mgf")], "does not exist"),
        (["--workdir", bare], "no launcher run directory"),
        (["--workdir", bare, "--spectra-dir", empty], "no *.mzML directly in"),
        (["--workdir", bare, "--spectra-dir", nope], f"--spectra-dir {nope} is not a directory"),
    ):
        err = io.StringIO()
        try:
            with contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
                em.main(["--out", kept, *args])  # a later --out in args wins
        except SystemExit as exc:
            assert exc.code == em.USAGE_ERROR and why in err.getvalue(), (why, err.getvalue())
        else:
            raise AssertionError(f"main accepted {args}")
    with open(kept, encoding="utf-8") as fh:
        assert fh.read() == "an earlier export\n"


def check_nothing_written(em, tmp, bare):
    """A spectra folder of other runs only: no MGF is written, and an earlier one stays as it was.

    The folder passes the check for spectra, so every row is dropped for a missing mzML. Written,
    the header-only file was refused by the validator after the export had said "wrote".
    """
    kept, elsewhere = os.path.join(tmp, "kept_empty.mgf"), os.path.join(tmp, "elsewhere")
    os.mkdir(elsewhere)
    em.write_mzml(os.path.join(elsewhere, "runZ_uncalibrated.mzML"), fx.SCANS)
    with open(kept, "w", encoding="utf-8") as fh:
        fh.write("an earlier export\n")
    err = io.StringIO()
    with contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
        code = em.main(["--workdir", bare, "--spectra-dir", elsewhere, "--out", kept])
    said = err.getvalue()
    assert code == 1 and f"{kept} was not written" in said and f" {em.NO_MZML}" in said, said
    with open(kept, encoding="utf-8") as fh:
        assert fh.read() == "an earlier export\n"


def run(em):
    """Every check, on the export_mgf module handed in."""
    em.require_mzml_reader()
    res, nt, ct = em.parse_mods("QMK", "1Q(-17.0265), 2M(15.9949), 2M(0.9840), C-term(0.9840)")
    assert em.proforma("QMK", res, nt, ct) == "Q[-17.0265]M[+16.9789]K-[+0.9840]"
    check_read_order(em)
    with tempfile.TemporaryDirectory() as tmp, warnings.catch_warnings():
        # Raised, not warned: read_spectrum then counts every scan unreadable.
        warnings.filterwarnings("error", fx.INDEX_FALLBACK)
        work, dirs, out = check_planted(em, tmp)
        bare = check_spectra_dir(em, tmp)
        check_refusals(em, tmp, work, bare)
        check_nothing_written(em, tmp, bare)
        check_mzml_names(em, tmp)
        os.replace(os.path.join(work, "expA", "psm.tsv"), os.path.join(work, "psm.tsv"))
        top, _ = em.export(work, dirs, re.compile("XORF"), out)
        assert (top[em.READ], top[em.WRITTEN]) == (14, 4), top
        check_rank_cases(em, tmp)
    print(
        "export_mgf self-test ok: 4 of 14 PSM rows written; dropped 1 decoy, 1 ragged, 2 at DDA+ "
        "rank 2, 1 by --match, 1 with no mass to check, 1 by the mass check, 1 with no precursor "
        "m/z, 1 by a missing mzML, 1 by an unreadable scan; under --all-ranks the 2 rank 2 picks "
        "are on spectra taken and fill no protein's place; per-protein and cap limits hold; scan "
        "107's rank 2 is dropped by default and written with RANK and its scan's best peptide "
        "under --all-ranks, and its two charges are one spectrum; a top-level psm.tsv is read; "
        "--spectra-dir with no manifest writes the manifest's records; a run's <run>.mzML is "
        "read where it has no <run>_uncalibrated.mzML, which wins where it has both; bad paths, "
        "and spectra neither finds, stop before --out is opened; a spectra folder of other runs "
        "writes nothing and leaves --out as it was; scans are read file by file and written in "
        "Hyperscore order"
    )
