#!/usr/bin/env python3
# /// codebase-meta
# exempt: >-
#   Self-test fixture for export_mgf, split out of it when that file crossed the 600-line ratchet
#   (ADR-12): a psm.tsv, a manifest and an indexed mzML written into a temporary directory.
#   It also holds the checks on the MGF read back. Imported by export_mgf's self-test only;
#   reads nothing real and writes only under the path it is given.
# ///
"""A finished FragPipe run in miniature, for export_mgf's self-test.

Every offset in the mzML is counted in bytes, and check_index holds the file to
it, because pyteomics hides a wrong or missing index: it re-parses the file
instead, and a fixture with a broken index would pass every other check.
"""

import os
import re
import warnings

from psm_spectrum import AA, H2O, SELFTEST_MASS_TOL, mzml_reader, write_mzml

INDEX_FALLBACK = "Could not extract the embedded offset index"

# The self-test's psm.tsv carries the header of a real run's psm.tsv. VADLTNSFEK and YSPGAVEK
# are made up: their masses are psm_spectrum.peptide_mass's, written to 4 decimals as FragPipe
# writes them, and each Observed M/Z keeps its row's old step from the Calculated M/Z. Every
# other mass is FragPipe's own from two real runs, so the mass check is held to numbers this
# module did not compute. MDKNELVQK is n-acetylated in one of them (1145.5751); the oxidised
# form adds 15.9949. The accessions are made up: XORF entries stand for novel ORFs, GENE ones
# for others.
PSM_HEADER = (
    "Spectrum",
    "Spectrum File",
    "Peptide",
    "Modified Peptide",
    "Extended Peptide",
    "Prev AA",
    "Next AA",
    "Peptide Length",
    "Charge",
    "Retention",
    "Observed Mass",
    "Calibrated Observed Mass",
    "Observed M/Z",
    "Calibrated Observed M/Z",
    "Calculated Peptide Mass",
    "Calculated M/Z",
    "Delta Mass",
    "SpectralSim",
    "RTScore",
    "Expectation",
    "Hyperscore",
    "Nextscore",
    "Probability",
    "Qvalue",
    "Number of Enzymatic Termini",
    "Number of Missed Cleavages",
    "Protein Start",
    "Protein End",
    "Intensity",
    "Assigned Modifications",
    "Observed Modifications",
    "Purity",
    "Is Decoy",
    "Is Contaminant",
    "Is Unique",
    "Protein",
    "Protein ID",
    "Entry Name",
    "Gene",
    "Protein Description",
    "Mapped Genes",
    "Mapped Proteins",
    "Parent Scan Number",
    "Apex Retention Time",
    "Apex Scan Number",
    "Retention Time Start",
    "Retention Time End",
    "Retention Time FWHM",
    "Traced Scans",
)
_NOVEL = "rs|XORF01762|NOVEL_ORF"
_CANON = "rs|GENE04300|GENE04300_ORGN"
_FIELDS = (
    "Spectrum",
    "Peptide",
    "Modified Peptide",
    "Charge",
    "Retention",
    "Observed M/Z",
    "Calculated Peptide Mass",
    "Calculated M/Z",
    "Hyperscore",
    "Probability",
    "Assigned Modifications",
    "Protein",
    "Mapped Proteins",
    "Spectrum File",
)
_A, _A2, _B = "interact-runA.pep.xml", "interact-runA_rank2.pep.xml", "interact-runB.pep.xml"
PSMS = (
    # Written: N-term acetyl and M1 oxidation, as one SEQ.
    (
        "runA.00101.00101.2",
        "MDKNELVQK",
        "n[43]M[147]DKNELVQK",
        "2",
        "1033.2210",
        "581.7931",
        "1161.5700",
        "581.7923",
        "30.5000",
        "0.9990",
        "N-term(42.0106), 1M(15.9949)",
        _NOVEL,
        "",
        _A,
    ),
    # No mass to check and an oxidation parse_mods cannot read: unguarded, it would take
    # the slot above with the oxidation missing from SEQ.
    (
        "runA.00108.00108.2",
        "MDKNELVQK",
        "",
        "2",
        "1040.1000",
        "581.7931",
        "",
        "581.7923",
        "35.0000",
        "0.9990",
        "N-term(42.0106), M1(15.9949)",
        _NOVEL,
        "",
        _A,
    ),
    # Written with no RTINSECONDS, although a better-scoring runB PSM of the same peptide
    # exists: runB has no mzML.
    (
        "runA.00102.00102.2",
        "VADLTNSFEK",
        "",
        "2",
        "",
        "562.2854",
        "1122.5557",
        "562.2851",
        "16.4710",
        "0.9998",
        "",
        _NOVEL,
        "",
        _A,
    ),
    (
        "runB.00201.00201.2",
        "VADLTNSFEK",
        "",
        "2",
        "3511.3474",
        "562.7902",
        "1122.5557",
        "562.2851",
        "18.6980",
        "0.9999",
        "",
        _NOVEL,
        "",
        _B,
    ),
    # Rank 2 of the spectrum written above for VADLTNSFEK, as DDA+ writes it. One TITLE holds
    # one record, and with the protein's limit applied before the spectrum test, this pick
    # would fill XORF01762's last place.
    (
        "runA.00102.00102.2",
        "VSTEVDAR",
        "",
        "2",
        "",
        "438.7252",
        "875.4348",
        "438.7247",
        "14.6190",
        "0.9989",
        "",
        _NOVEL,
        "",
        _A2,
    ),
    # XORF01762's fourth peptide: written only because the rank-2 pick above takes no place.
    (
        "runA.00107.00107.2",
        "IEEDGEVR",
        "",
        "2",
        "915.7432",
        "473.7290",
        "945.4403",
        "473.7274",
        "12.7390",
        "0.8499",
        "",
        _NOVEL,
        "",
        _A,
    ),
    # A decoy of a novel entry: it matches XORF, and is no --match hit.
    (
        "runA.00111.00111.2",
        "LIPGSHER",
        "",
        "2",
        "873.0659",
        "454.7516",
        "907.4876",
        "454.7511",
        "11.6110",
        "0.9992",
        "",
        "rev_rs|XORF01762|NOVEL_ORF",
        "",
        _A,
    ),
    # Wrong Calculated Peptide Mass, and the best Hyperscore of its peptide: it must not
    # take the slot of the correct PSM below.
    (
        "runA.00103.00103.2",
        "GDDGAHSINEK",
        "",
        "2",
        "641.3631",
        "571.7577",
        "1157.4948",
        "571.7572",
        "40.0000",
        "1.0000",
        "",
        _CANON,
        "rs|XORF01293|SHARES_STOP_ORF",
        _A,
    ),
    # Kept only through Mapped Proteins, with no Observed M/Z, so PEPMASS falls back.
    (
        "runA.00104.00104.2",
        "GDDGAHSINEK",
        "",
        "2",
        "639.7367",
        "",
        "1141.4999",
        "571.7572",
        "23.9530",
        "1.0000",
        "",
        _CANON,
        "rs|XORF01293|SHARES_STOP_ORF",
        _A,
    ),
    # Chosen, but scan 105 is not in the mzML.
    (
        "runA.00105.00105.2",
        "YSPGAVEK",
        "",
        "2",
        "670.5567",
        "425.7189",
        "849.4232",
        "425.7189",
        "9.5090",
        "0.9914",
        "",
        "rs|XORF01924|NOVEL_ORF",
        "",
        _A,
    ),
    # No precursor m/z at all: dropped before it can take the slot above.
    (
        "runA.00110.00110.2",
        "YSPGAVEK",
        "",
        "2",
        "690.0000",
        "",
        "849.4232",
        "",
        "12.0000",
        "0.9914",
        "",
        "rs|XORF01924|NOVEL_ORF",
        "",
        _A,
    ),
    # Rank 2 on scan 101, which MDKNELVQK takes, of XORF01924, whose one place under
    # --per-protein 1 YSPGAVEK fills first. It is beyond that limit and on a spectrum taken, and
    # counts as the spectrum taken, the test select() makes first.
    (
        "runA.00101.00101.2",
        "IEEDGEVR",
        "",
        "2",
        "915.7432",
        "473.7290",
        "945.4403",
        "473.7274",
        "9.0000",
        "0.8499",
        "",
        "rs|XORF01924|NOVEL_ORF",
        "",
        _A2,
    ),
    # No XORF anywhere.
    (
        "runA.00106.00106.2",
        "SNPSHNPK",
        "n[43]SNPSHNPK",
        "2",
        "572.5961",
        "461.7231",
        "921.4305",
        "461.7225",
        "20.3350",
        "0.9999",
        "N-term(42.0106)",
        "rs|GENE04720|GENE04720_ORGN",
        "",
        _A,
    ),
)
SCANS = (
    (101, [401.28707, 175.11895, 288.20301], [50.0, 1000.0, 250.5]),
    (102, [147.1128, 262.1397], [10.0, 20.0]),
    (104, [120.0, 130.5, 140.25, 150.125], [1.0, 2.0, 3.0, 4.0]),
    (107, [175.11895, 274.18737, 403.22996], [300.0, 120.0, 60.0]),
)


def write_tree(tmp, manifest):
    """<tmp>/proj/{runs,results}/fp, the manifest naming <tmp>/data, and the mzML there."""
    data, run_dir = os.path.join(tmp, "data"), os.path.join(tmp, "proj", "runs", "fp")
    work = os.path.join(tmp, "proj", "results", "fp")
    for d in (data, run_dir, os.path.join(work, "expA")):
        os.makedirs(d)
    with open(os.path.join(run_dir, manifest), "w", encoding="utf-8") as fh:
        fh.writelines(f"{data}/{r}.raw\texpA\t\tDDA\n" for r in ("runA", "runB"))
    write_mzml(os.path.join(data, "runA_uncalibrated.mzML"), SCANS)
    with open(os.path.join(work, "expA", "psm.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(PSM_HEADER) + "\n")
        for psm in PSMS:
            row = dict(zip(_FIELDS, psm, strict=True))
            row["Protein ID"] = row["Protein"] + " label=x class=y OS=Example organism"
            fh.write("\t".join(row.get(c, "") for c in PSM_HEADER) + "\n")
        fh.write("\nrunA.00109.00109.2\tMDKNEL")  # a blank line, then a file cut mid-write
    return work, data


def check_index(path):
    """Each embedded offset lands on its own spectrum, and pyteomics reads that index.

    Neither shows in what export reads back: get_by_id quietly re-parses from the top on
    a bad offset, and PreIndexedMzML only warns, with INDEX_FALLBACK, when it has no index.
    """
    from pyteomics import mzml  # noqa: PLC0415  optional, as in psm_spectrum

    with open(path, "rb") as fh:
        raw = fh.read()
    with warnings.catch_warnings():
        warnings.filterwarnings("error", INDEX_FALLBACK)
        with mzml_reader(path, mzml) as reader:
            index = dict(reader.index["spectrum"])
    assert len(index) == len(SCANS), index
    for i, (scan, _mz, _inten) in enumerate(SCANS):
        sid = f"controllerType=0 controllerNumber=1 scan={scan}"
        assert raw.startswith(f'<spectrum index="{i}" id="{sid}"'.encode(), index[sid]), sid


_PROFORMA = (
    r"(?:\[([+-]\d+\.\d{4})\]-)?((?:[A-Z](?:\[[+-]\d+\.\d{4}\])?)+)"
    r"(?:-\[([+-]\d+\.\d{4})\])?"
)


def _mass_from_proforma(seq):
    m = re.fullmatch(_PROFORMA, seq)
    assert m, seq
    nterm, body, cterm = m.groups()
    mass = H2O + float(nterm or 0) + float(cterm or 0)
    for aa, delta in re.findall(r"([A-Z])(?:\[([+-]\d+\.\d{4})\])?", body):
        mass += AA[aa] + float(delta or 0)
    return mass


def check_mgf_records(out, counts, work, extra=frozenset()):
    """The MGF read back with pyteomics, field by field; extra names the --all-ranks keys."""
    from pyteomics import mgf  # noqa: PLC0415  deferred like psm_spectrum's mzml import

    want = {  # title -> (SEQ, PEPMASS, PROTEIN, peaks, first peak)
        "runA.00101.00101.2": (
            "[+42.0106]-M[+15.9949]DKNELVQK",
            581.7931,
            "XORF01762",
            3,
            (175.11895, 1000.0),
        ),
        "runA.00104.00104.2": ("GDDGAHSINEK", 571.7572, "GENE04300", 4, (120.0, 1.0)),
        "runA.00102.00102.2": ("VADLTNSFEK", 562.2854, "XORF01762", 2, (147.1128, 10.0)),
        "runA.00107.00107.2": ("IEEDGEVR", 473.7290, "XORF01762", 3, (175.11895, 300.0)),
    }
    # Any header line pyteomics took for a param would show up here as an extra key.
    keys = {
        "title",
        "pepmass",
        "charge",
        "seq",
        "scans",
        "rtinseconds",
        "protein",
        "hyperscore",
        "probability",
    } | extra
    odd = {
        "runA.00102.00102.2": keys - {"rtinseconds"},
        "runA.00104.00104.2": keys | {"pepmass_source", "matched_in"},
    }
    # Keyed by peptide as well: a DDA+ spectrum carries a PSM per rank.
    calc = {(p[0], p[1]): float(p[6]) for p in PSMS if p[0] in want}
    spectra = list(mgf.read(out))
    assert [s["params"]["title"] for s in spectra] == list(want), spectra
    for s in spectra:
        prm = s["params"]
        assert set(prm) == odd.get(prm["title"], keys), prm
        seq, pepmass, prot, n_peaks, first = want[prm["title"]]
        assert prm["seq"] == seq, prm["seq"]
        mass = calc[prm["title"], re.sub("[^A-Z]", "", seq)]
        assert abs(_mass_from_proforma(seq) - mass) < SELFTEST_MASS_TOL, seq
        assert list(prm["charge"]) == [2] and prm["pepmass"][0] == pepmass, prm
        assert prm["protein"] == prot and prm["scans"] == prm["title"].split(".")[1].lstrip("0")
        assert len(s["m/z array"]) == n_peaks, s
        assert (s["m/z array"][0], s["intensity array"][0]) == first, s
        assert list(s["m/z array"]) == sorted(s["m/z array"])
    with open(out, encoding="utf-8") as fh:
        text = fh.read()
    assert not re.search(r"\[(?![+-])", text), "unsigned bracket"
    assert (spectra[1]["params"]["pepmass_source"], spectra[1]["params"]["matched_in"]) == (
        "Calculated M/Z",
        "Mapped Proteins only",
    ), spectra[1]
    assert f"# export_mgf.py from {work.replace('=', '%3D')}\n" in text, text
    assert "\n# --match XORF\n# --per-protein 3, --max-spectra 500, " in text, text
    assert all(s["params"].get("rank", 1) in (1, "1") for s in spectra), spectra
    assert counts["PROBABILITY from psm.tsv column Probability"] == len(want)
