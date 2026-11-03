#!/usr/bin/env python3
# /// codebase-meta
# exempt: >-
#   One PSM's peptide and its scan, for report_spectra and export_mgf: residue masses, the
#   modification parse and mass check, and the run's mzML found, read and, for a reader probe or
#   a fixture, written. Split out of report_spectra so that export_mgf needs none of the
#   report's drawing. Reads only the mzML it is given and the PSI-MS vocabulary psims ships,
#   never the network; writes only the path it is given.
# ///
"""One PSM's peptide, checked against FragPipe's own mass, and its scan in the run's mzML.

report_spectra draws these PSMs and export_mgf writes them as MGF. Both read the peptide, the
mass check and the scan here, so the two never disagree about the same PSM.

Correctness decisions, each of which was a real trap:

* Peptides come from `Peptide` + `Assigned Modifications`, never `Modified
  Peptide`. The latter's brackets are absolute residue+mod in some MSFragger
  workflows and plain deltas in others; `Assigned Modifications` is always a
  delta, written `%.4f`.
* The mzML spectrum id format changed in MSFragger 4.4, so it is read from the
  file's own index rather than assumed.
* Every annotation is checked against the `Calculated Peptide Mass` column
  before it is used. A peptide whose reconstructed mass disagrees is skipped,
  because a plausible-looking wrong annotation is worse than none.
"""

import base64
import functools
import hashlib
import os
import re
import struct
import sys

from report_psms import _num, spectrum_parts

# Monoisotopic residue masses. Hard-coded rather than imported so this module
# keeps working if pyteomics is ever absent — only mzML reading truly needs it.
AA = {
    "G": 57.02146,
    "A": 71.03711,
    "S": 87.03203,
    "P": 97.05276,
    "V": 99.06841,
    "T": 101.04768,
    "C": 103.00919,
    "L": 113.08406,
    "I": 113.08406,
    "N": 114.04293,
    "D": 115.02694,
    "Q": 128.05858,
    "K": 128.09496,
    "E": 129.04259,
    "M": 131.04049,
    "H": 137.05891,
    "F": 147.06841,
    "R": 156.10111,
    "Y": 163.06333,
    "W": 186.07931,
}
H2O = 18.0105646837
PROTON = 1.00727646677
MASS_CHECK_TOL = 0.01  # Da; reconstructed vs FragPipe's own Calculated Peptide Mass
SELFTEST_MASS_TOL = 0.001
MZML_SUFFIX = "_uncalibrated.mzML"  # what MSFragger writes beside each Thermo .raw input
MZML_NAMES = (MZML_SUFFIX, ".mzML")  # looked for in this order in each folder, after <run>
PSI_MS_URI = "http://purl.obolibrary.org/obo/ms/psi-ms.obo"  # psims maps it to its own copy


def parse_mods(peptide, assigned):
    """Return (residue masses, n-term delta, c-term delta).

    `assigned` is FragPipe's Assigned Modifications: '1M(15.9949), N-term(42.0106)'
    — position, residue, delta. Terminal entries carry no position.
    """
    res = [AA.get(a) for a in peptide]
    if any(r is None for r in res):
        return None, 0.0, 0.0  # non-standard residue (B/X/Z): skip
    nterm = cterm = 0.0
    for m in re.finditer(r"(?:(\d+)([A-Z])|(N-term|C-term))\(([-+]?\d*\.?\d+)\)", assigned or ""):
        pos, aa_, term, delta = m.groups()
        delta = float(delta)
        if term == "N-term":
            nterm += delta
        elif term == "C-term":
            cterm += delta
        else:
            i = int(pos) - 1
            if not (0 <= i < len(peptide)) or peptide[i] != aa_:
                return None, 0.0, 0.0  # parse disagrees with the sequence
            res[i] += delta
    return res, nterm, cterm


def peptide_mass(res, nterm, cterm):
    """Neutral monoisotopic mass of the modified peptide."""
    return sum(res) + H2O + nterm + cterm


def _validated_psm(row):
    """Parse and sanity-check one PSM. -> (peptide, res, nt, ct, run, scan) or None.

    The mass check is the important part: FragPipe already computed the peptide
    mass, so any disagreement means our modification handling is wrong and the
    annotation would be confidently incorrect.
    """
    pep = row.get("Peptide", "").strip()
    if not pep:
        return None
    res, nt, ct = parse_mods(pep, row.get("Assigned Modifications", ""))
    if res is None:
        return None
    calc = _num(row, "Calculated Peptide Mass")
    if calc and abs(calc - peptide_mass(res, nt, ct)) > MASS_CHECK_TOL:
        return None
    spectrum = spectrum_parts(row)
    return None if spectrum is None else (pep, res, nt, ct, *spectrum)


def _find_mzml(data_dirs, run):
    """The run's spectra: in each folder in turn, <run>_uncalibrated.mzML, else <run>.mzML.

    MSFragger writes <run>_uncalibrated.mzML beside a Thermo .raw input, not into results. A
    search of mzML files, from msconvert for example, has its input <run>.mzML and no such copy.
    The order keeps a .raw run on its uncalibrated file where a <run>.mzML sits beside it too.
    A Bruker .d folder is not read.
    """
    for d in data_dirs:
        for suffix in MZML_NAMES:
            cand = os.path.join(d, run + suffix)
            if os.path.isfile(cand):
                return cand
    return None


def _spectrum_key(index_keys, scan):
    """Build the mzML lookup key for `scan` in whichever id format the file uses.

    MSFragger 4.4 changed Thermo ids from 'scan=N' to the controllerType form,
    so this is read from the file's own index rather than assumed.
    """
    sample = next(iter(index_keys), "")
    if sample.startswith("controllerType"):
        return f"controllerType=0 controllerNumber=1 scan={scan}"
    return f"scan={scan}"


@functools.lru_cache(maxsize=1)
def offline_vocabulary():
    """The PSI-MS vocabulary from the copy psims ships, built once, or None where that cannot be.

    pyteomics 5 builds the vocabulary for every mzML reader, and psims by default downloads it
    from purl.obolibrary.org each time, using its own copy only when that fails: 11 downloads in
    one export_mgf self-test. Private data read on a laptop, and a compute node with no internet,
    must wait on no such call. Where this pyteomics or psims lacks the API, the readers are built
    as before, and one line on stderr says so.
    """
    try:
        from psims.controlled_vocabulary.controlled_vocabulary import OBOCache  # noqa: PLC0415
        from pyteomics import mzml, xml  # noqa: PLC0415

        if not issubclass(mzml.PreIndexedMzML, xml.CVParamParser):
            raise TypeError("this pyteomics' PreIndexedMzML takes no cv")
        return OBOCache(enabled=False, use_remote=False).load(PSI_MS_URI)
    except Exception as exc:  # noqa: BLE001  any older or newer API keeps the old call
        print(
            f"psm_spectrum: no offline PSI-MS vocabulary ({exc!r}); pyteomics loads its own, "
            "from the network where it can",
            file=sys.stderr,
        )
        return None


def mzml_reader(mzml_path, pymzml):
    """A pyteomics PreIndexedMzML over mzml_path, given offline_vocabulary() where there is one."""
    cv = offline_vocabulary()
    if cv is None:
        return pymzml.PreIndexedMzML(mzml_path)
    return pymzml.PreIndexedMzML(mzml_path, cv=cv)


@functools.lru_cache(maxsize=4)
def _open_mzml(mzml_path, pymzml):
    """An open reader per file, so its index is parsed once rather than once per scan.

    Opening reads the file's whole offset index. Reopened per scan, a 5,422-scan export of one
    real run had not finished after 10 minutes; with this cache it took 6 s. It holds 4 files,
    so a caller reading more must go one file at a time, as export_mgf.read_records does.
    A file that fails to open raises, and lru_cache does not keep a raised call, so it is retried.
    """
    return mzml_reader(mzml_path, pymzml)


def read_spectrum(mzml_path, scan):
    """Return (mz, intensity) arrays for one scan, or (None, None)."""
    try:
        # Deferred: pyteomics is optional, and this module must import cleanly
        # without it so the rest of the report still builds.
        from pyteomics import mzml as pymzml  # noqa: PLC0415
    except ImportError:
        return None, None
    try:
        reader = _open_mzml(mzml_path, pymzml)
        spec = reader.get_by_id(_spectrum_key(reader.index["spectrum"].keys(), scan))
        return spec["m/z array"], spec["intensity array"]
    except Exception:  # noqa: BLE001
        # Deliberately broad. pyteomics parses with lxml, whose XMLSyntaxError
        # descends from SyntaxError — outside every exception class this used to
        # name — so a zero-length or head-truncated *_uncalibrated.mzML (left by
        # a preempted job, or by a second run truncating the file as it reopens
        # it for writing) escaped all the way out of build() and cost a
        # multi-hour run its entire report. One unreadable scan is one missing
        # slide; it is never worth the page.
        return None, None


def _b64(values):
    return base64.b64encode(struct.pack(f"<{len(values)}d", *values)).decode()


def _spectrum(i, sid, mz, inten):
    arrays = "".join(
        f'<binaryDataArray encodedLength="{len(_b64(v))}">'
        '<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>'
        '<cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>'
        f'<cvParam cvRef="MS" accession="{acc}" name="{name}" value=""/>'
        f"<binary>{_b64(v)}</binary></binaryDataArray>"
        for acc, name, v in (
            ("MS:1000514", "m/z array", mz),
            ("MS:1000515", "intensity array", inten),
        )
    )
    return (
        f'<spectrum index="{i}" id="{sid}" defaultArrayLength="{len(mz)}">'
        '<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2"/>'
        f'<binaryDataArrayList count="2">{arrays}</binaryDataArrayList></spectrum>\n'
    )


def write_mzml(path, scans):
    """Indexed mzML with MSFragger 4.4 ids; every offset is counted in bytes, not characters.

    Args:
        path: The file to write.
        scans: [(scan number, m/z values, intensities)] of MS2 scans, in file order.
    """
    body = (
        '<?xml version="1.0" encoding="utf-8"?>\n'
        '<indexedmzML xmlns="http://psi.hupo.org/ms/mzml">\n'
        '<mzML xmlns="http://psi.hupo.org/ms/mzml" id="runA" version="1.1.0">\n'
        f'<run id="runA"><spectrumList count="{len(scans)}">\n'
    ).encode()
    offsets = []
    for i, (scan, mz, inten) in enumerate(scans):
        sid = f"controllerType=0 controllerNumber=1 scan={scan}"
        offsets.append((sid, len(body)))
        body += _spectrum(i, sid, mz, inten).encode()
    body += b"</spectrumList></run>\n</mzML>\n"
    index_at = len(body)
    body += (
        '<indexList count="1"><index name="spectrum">\n'
        + "".join(f'<offset idRef="{s}">{o}</offset>\n' for s, o in offsets)
        + f"</index></indexList>\n<indexListOffset>{index_at}</indexListOffset>\n"
        "<fileChecksum>"
    ).encode()
    body += hashlib.sha1(body).hexdigest().encode() + b"</fileChecksum>\n</indexedmzML>\n"
    with open(path, "wb") as fh:
        fh.write(body)
