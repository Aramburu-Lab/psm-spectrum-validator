#!/usr/bin/env python3
# /// codebase-meta
# exempt: >-
#   Export utility, not a pipeline step: point it at a finished FragPipe results directory and it
#   writes the chosen PSMs' spectra as annotated MGF for the spectrum-validation browser tool.
#   Reads psm.tsv, each run's mzML and, without --spectra-dir, the launcher's run manifest;
#   writes only its --out.
# ///
r"""Chosen PSMs of a FragPipe run as annotated MGF, small enough to pull to a laptop.

    python3 scripts/python/export_mgf.py --workdir <proj>/results/<run> \
        --match 'label=novel\b' --out <run>_novel.mgf

A run the launcher did not make has no manifest to say where its spectra are, so name
their folder with --spectra-dir. Each run is read from <run>_uncalibrated.mzML, else
from <run>.mzML, as psm_spectrum._find_mzml says; a Bruker .d folder is not read.

Match a label, not an accession prefix that Mapped Proteins also carries: on one
real run such a prefix matched 2,541 rows, 2,064 of them only through Mapped
Proteins, so labelled by an annotated Protein ID. Such a record carries
MATCHED_IN=Mapped Proteins only.

Nothing that decides whether an annotation is right is written here. The
peptide, its modifications and the mass check are psm_spectrum's
`_validated_psm`, the scan is its `read_spectrum`, where the mzML lives is
run_layout's rule unless --spectra-dir says, and the protein is
report_accession's `fasta_accession`. A second copy of any of them is how this
file and the report's spectra section would come to disagree about the same
PSM. The self-test and its fixture are imported only under --selftest, so an
export needs neither.

Filter order matters. Decoy rows go first, so a decoy of a matching entry is
never counted a --match hit. The mass check and the mzML lookup run BEFORE the
best-PSM-per-peptide pick, so a PSM that can never be written cannot take the
slot of one that can. Every drop is printed and stamped in the MGF header, and
rows read always equal the drops plus the spectra written.

A DDA+ psm.tsv holds up to five ranks of one scan, the rank named only by the
row's Spectrum File (report_psms.psm_rank). Only rank 1 is exported by default:
a rank 2 to 5 PSM is a co-isolated peptide a better-scoring one outranks on the
same spectrum, and written unmarked it would hand the validator a chimeric
spectrum annotated with the weaker peptide. Ranks 2 and up are a named drop.
--all-ranks exports them, and marks every record with RANK= and with the best
non-decoy peptide of its scan, taken over every row before --match, so a rank 2
whose rank 1 is canonical still says so. A scan is (run, scan), whatever charge
a PSM gave it, so a pick on a scan a better-scoring peptide took is dropped
before its protein's limit and never fills a protein's place. A scan taken that
way and an unreadable scan are the only unwritable picks found after the pick;
they are counted, not replaced.

SEQ is ProForma 2 with SIGNED deltas, never FragPipe's `Modified Peptide`,
whose brackets hold absolute masses (`n[43]M[147]...`). The browser tool
rejects an unsigned bracket for exactly that reason.
"""

import argparse
import os
import re
import sys
import tempfile
from collections import Counter, namedtuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from psm_spectrum import (  # noqa: E402
    AA,
    MASS_CHECK_TOL,
    MZML_SUFFIX,
    _find_mzml,
    _validated_psm,
    mzml_reader,
    parse_mods,
    peptide_mass,
    read_spectrum,
    write_mzml,
)
from report_accession import fasta_accession, is_decoy_row  # noqa: E402
from report_psms import _num, _psm_rows, psm_rank, psm_tables, spectrum_parts  # noqa: E402
from run_layout import MANIFEST, input_dirs, load_manifest, run_dir_for  # noqa: E402

PER_PROTEIN = 3
MAX_SPECTRA = 500
MATCH_PRIMARY = ("Protein ID", "Protein")
MATCH_MAPPED = "Mapped Proteins"
USAGE_ERROR = 2  # argparse's exit status for p.error
# On real runs with PeptideProphet off and Percolator on, psm.tsv has only
# "Probability", so it is not PeptideProphet's number. The column actually used is
# stamped in the MGF header so the two are never confused.
PROB_COLS = ("PeptideProphet Probability", "Probability")

READ = "psm.tsv rows read"
DECOY = "rows of decoy PSMs"
RAGGED = "rows whose field count differs from the header"
RANKED = "rows at DDA+ rank 2 or above, exported only with --all-ranks"
NO_MATCH = "rows not matching --match"
NO_CALC = "rows with no Calculated Peptide Mass, so no mass check"
MASS = "rows failing the mass check"
UNPARSED = "rows with an unparseable peptide, modification or Spectrum id"
NO_PEPMASS = "rows with neither Observed nor Calculated M/Z"
NO_MZML = "rows whose run has no <run>_uncalibrated.mzML or <run>.mzML"
MAPPED_ONLY = "rows matching --match only through Mapped Proteins, labelled by Protein ID"
NOT_BEST = "rows not the best-Hyperscore PSM of their peptide"
SHARED_SCAN = "picks whose scan is taken by a better-scoring peptide"
BEYOND = "peptides beyond --per-protein"
CAP = "spectra dropped by --max-spectra"
UNREADABLE = "scans unreadable in their mzML"
FALLBACK = "PEPMASS from Calculated M/Z, no Observed M/Z"
NO_RT = "spectra with no RTINSECONDS, Retention blank or not a number"
WRITTEN = "spectra written"
PRE_PICK = (DECOY, RAGGED, RANKED, NO_MATCH, NO_CALC, MASS, UNPARSED, NO_PEPMASS, NO_MZML)
POST_PICK = (SHARED_SCAN, BEYOND, CAP, UNREADABLE)
DROPS = (*PRE_PICK, NOT_BEST, *POST_PICK)
ORDER = (READ, *PRE_PICK, MAPPED_ONLY, NOT_BEST, *POST_PICK, FALLBACK, NO_RT, WRITTEN)
RANK_KEYS = frozenset({"rank", "scan_best_peptide", "scan_best_hyperscore", "scan_best_rank"})
Limits = namedtuple(
    "Limits", "per_protein max_spectra all_ranks", defaults=(PER_PROTEIN, MAX_SPECTRA, False)
)
DEFAULTS = Limits()


def proforma(pep, res, nterm, cterm):
    """ProForma 2 with signed 4-decimal deltas, read off parse_mods' output.

    Args:
        pep: Bare peptide sequence.
        res: Per-residue masses with their modifications added, from parse_mods.
        nterm: N-terminal delta.
        cterm: C-terminal delta.

    Returns:
        The ProForma string. A delta is a residue's mass minus the bare residue,
        so two modifications on one residue arrive already summed into one bracket.
    """
    out = [f"[{nterm:+.4f}]-"] if round(nterm, 4) else []
    for aa, mass in zip(pep, res, strict=True):
        delta = mass - AA[aa]
        out.append(f"{aa}[{delta:+.4f}]" if round(delta, 4) else aa)
    if round(cterm, 4):
        out.append(f"-[{cterm:+.4f}]")
    return "".join(out)


def _why_invalid(row):
    """Which of _validated_psm's refusals this was, so the mass check is counted alone."""
    pep = row.get("Peptide", "").strip()
    res, nt, ct = parse_mods(pep, row.get("Assigned Modifications", ""))
    calc = _num(row, "Calculated Peptide Mass")
    if pep and res is not None and calc and abs(calc - peptide_mass(res, nt, ct)) > MASS_CHECK_TOL:
        return MASS
    return UNPARSED


def _match(row, pattern):
    """'primary', 'mapped' when only Mapped Proteins matched, or '' for no match."""
    if pattern is None or any(pattern.search(row.get(c, "")) for c in MATCH_PRIMARY):
        return "primary"
    return "mapped" if pattern.search(row.get(MATCH_MAPPED, "")) else ""


def _number(row, key):
    """The cell as written when float() reads it, else ''; MGF readers convert with float()."""
    text = row.get(key, "").strip()
    try:
        float(text)
    except ValueError:
        return ""
    return text


def _writable(row, pattern, data_dirs, mzml_of, counts):
    """(parsed PSM, mzML path, match) for a writable row, else None with its reason counted."""
    hit = _match(row, pattern)
    if not hit:
        counts[NO_MATCH] += 1
        return None
    # Counted before any other drop: it describes what --match selects, not what survives.
    counts[MAPPED_ONLY] += hit == "mapped"
    # _validated_psm passes a row with no mass to check, and parse_mods skips a token it
    # cannot read, so such a row would be written with a modification missing from SEQ.
    if not _num(row, "Calculated Peptide Mass"):
        counts[NO_CALC] += 1
        return None
    parsed = _validated_psm(row)
    if parsed is None:
        counts[_why_invalid(row)] += 1
        return None
    if not (_number(row, "Observed M/Z") or _number(row, "Calculated M/Z")):
        counts[NO_PEPMASS] += 1
        return None
    run = parsed[4]
    if run not in mzml_of:
        mzml_of[run] = _find_mzml(data_dirs, run)
    if not mzml_of[run]:
        counts[NO_MZML] += 1
        return None
    return parsed, mzml_of[run], hit


def _scan_best(row, scans):
    """Keep the row in scans if it is its scan's best so far, the tie to the lower peptide."""
    spectrum = spectrum_parts(row)
    if spectrum is None:
        return
    key, score = spectrum[:2], _num(row, "Hyperscore")
    here = (score, row.get("Peptide", "").strip(), psm_rank(row) or 1, row.get("Hyperscore", ""))
    prev = scans.get(key)
    if prev is None or score > prev[0] or (score == prev[0] and here[1] < prev[1]):
        scans[key] = here


def best_per_peptide(workdir, pattern, data_dirs, counts, all_ranks=False):
    """Each peptide's best writable PSM, from every <experiment>/psm.tsv under workdir.

    Args:
        workdir: FragPipe results directory.
        pattern: Compiled --match regex, or None to keep every protein.
        data_dirs: Directories searched for each run's mzML, by psm_spectrum._find_mzml.
        counts: Counter that receives every reason a row was not kept.
        all_ranks: Keep DDA+ ranks 2 and up instead of counting them a drop.

    Returns:
        ({(accession, peptide): (row, parsed, mzml, match)}, {run: its mzML, or None},
        {(run, scan): (Hyperscore, peptide, rank, Hyperscore text)} of each scan's best
        non-decoy row, read before --match and the rank rule).
    """
    best, mzml_of, scans = {}, {}, {}
    for path in psm_tables(workdir):
        ragged = []
        for row in _psm_rows(path, ragged):
            counts[READ] += 1
            if is_decoy_row(row):
                counts[DECOY] += 1
                continue
            _scan_best(row, scans)
            if not all_ranks and (psm_rank(row) or 1) > 1:
                counts[RANKED] += 1
                continue
            ok = _writable(row, pattern, data_dirs, mzml_of, counts)
            if ok is None:
                continue
            key = (fasta_accession(row.get("Protein ID", "").strip()), ok[0][0])
            prev = best.get(key)
            if prev is None or _num(row, "Hyperscore") > _num(prev[0], "Hyperscore"):
                best[key] = (row, *ok)
        counts[READ] += len(ragged)
        counts[RAGGED] += len(ragged)
    return best, mzml_of, scans


def _rank(pick):
    """Hyperscore, highest first, the Spectrum id breaking ties so reruns write one order."""
    return -_num(pick[0], "Hyperscore"), pick[0].get("Spectrum", "")


def select(best, per_protein, max_spectra, counts):
    """Each scan once, then the per_protein best peptides of each protein, then max_spectra.

    One walk in Hyperscore order tests each pick for a scan, (run, scan), already taken before its
    protein's limit, so a pick that cannot be written never fills a protein's last place, and
    a pick its protein's limit drops never takes a spectrum.

    Args:
        best: Output of best_per_peptide.
        per_protein: Distinct peptides kept per protein.
        max_spectra: Hard cap on spectra written.
        counts: Counter that receives what each limit dropped.

    Returns:
        [(row, parsed, mzml, match)] in descending Hyperscore.
    """
    taken, kept, chosen = set(), Counter(), []
    # (accession, peptide) breaks a tie, so which pick takes a spectrum never follows file order.
    for (acc, _pep), pick in sorted(best.items(), key=lambda item: (*_rank(item[1]), item[0])):
        scan = pick[1][4:6]
        if scan in taken:
            counts[SHARED_SCAN] += 1
        elif kept[acc] >= per_protein:
            counts[BEYOND] += 1
        else:
            taken.add(scan)
            kept[acc] += 1
            chosen.append(pick)
    counts[CAP] += max(0, len(chosen) - max_spectra)
    return chosen[:max_spectra]


def read_records(picks, counts, keys_of):
    """Each pick's MGF record in the picks' own order, the scans read one mzML file at a time.

    select() returns Hyperscore order, which interleaves a run's files, and psm_spectrum's
    reader cache holds 4. Read in that order, one real run's controls, 6,841 scans in 14
    files, cost 3,451 index parses; read file by file, they cost 14.

    Args:
        picks: select()'s [(row, parsed, mzml, match)].
        counts: Counter that receives unreadable scans and mgf_record's counts.
        keys_of: pick -> the record's own KEY=value lines, written after HYPERSCORE.

    Returns:
        The readable picks' records, in the order of picks.
    """
    records = {}
    for i in sorted(range(len(picks)), key=lambda k: picks[k][2]):
        row, parsed, mzml, _match = picks[i]
        mz, inten = read_spectrum(mzml, parsed[5])
        if mz is None:
            counts[UNREADABLE] += 1
            continue
        records[i] = mgf_record(row, parsed, mz, inten, counts, keys_of(picks[i]))
    return [records[i] for i in sorted(records)]


def _plain(value):
    """An intensity without exponent or trailing zeros."""
    return f"{float(value):.4f}".rstrip("0").rstrip(".")


def record_keys(pick, scans, all_ranks):
    """MATCHED_IN for a Mapped-Proteins-only match; under --all-ranks, the rank and scan's best.

    Args:
        pick: select()'s (row, parsed, mzml, match).
        scans: best_per_peptide's {(run, scan): best non-decoy row's fields}.
        all_ranks: Whether --all-ranks was given.

    Returns:
        The KEY=value lines, keys inside the block because a parser drops every line outside it.
    """
    row, parsed, _mzml, match = pick
    keys = ["MATCHED_IN=Mapped Proteins only"] if match == "mapped" else []
    if all_ranks:
        _score, pep, rank, score_text = scans[parsed[4:6]]
        keys += [f"RANK={psm_rank(row) or 1}", f"SCAN_BEST_PEPTIDE={pep}"]
        keys += [f"SCAN_BEST_HYPERSCORE={score_text.strip()}", f"SCAN_BEST_RANK={rank}"]
    return keys


def mgf_record(row, parsed, mz, inten, counts, keys):
    """One BEGIN IONS ... END IONS block in the browser tool's contract.

    Args:
        row: The psm.tsv row.
        parsed: Its _validated_psm tuple.
        mz: m/z array of the scan.
        inten: Intensity array of the scan.
        counts: Counter that records PEPMASS fallbacks, missing RT and the probability column.
        keys: record_keys' lines for this pick.

    Returns:
        The record as text. A theoretical PEPMASS is marked by a key inside the block, because
        a parser drops every line outside it.
    """
    pep, res, nt, ct, _run, scan, z_str = parsed
    lines = ["BEGIN IONS", f"TITLE={row['Spectrum']}"]
    pepmass = _number(row, "Observed M/Z")
    if pepmass:
        lines.append(f"PEPMASS={pepmass}")
    else:
        counts[FALLBACK] += 1
        lines += [f"PEPMASS={_number(row, 'Calculated M/Z')}", "PEPMASS_SOURCE=Calculated M/Z"]
    lines += [
        f"CHARGE={int(_num(row, 'Charge') or z_str)}+",
        f"SEQ={proforma(pep, res, nt, ct)}",
        f"SCANS={scan}",
    ]
    # An empty RTINSECONDS= makes pyteomics refuse the whole file, not one record.
    rt = _number(row, "Retention")
    if rt:
        lines.append(f"RTINSECONDS={rt}")
    else:
        counts[NO_RT] += 1
    lines += [
        f"PROTEIN={fasta_accession(row.get('Protein ID', '').strip())}",
        f"HYPERSCORE={row.get('Hyperscore', '').strip()}",
    ]
    lines += keys
    col = next((c for c in PROB_COLS if row.get(c, "").strip()), None)
    if col:
        lines.append(f"PROBABILITY={row[col].strip()}")
        counts[f"PROBABILITY from psm.tsv column {col}"] += 1
    lines += [f"{m:.5f} {_plain(i)}" for m, i in sorted(zip(mz, inten, strict=True))]
    lines.append("END IONS")
    return "\n".join(lines) + "\n"


def export(workdir, data_dirs, pattern, out, limits=DEFAULTS):
    """Write the MGF.

    Args:
        workdir: FragPipe results directory holding <experiment>/psm.tsv.
        data_dirs: Directories searched for each run's mzML, by psm_spectrum._find_mzml.
        pattern: Compiled --match regex, or None.
        out: Path of the .mgf to write.
        limits: --per-protein, --max-spectra and --all-ranks.

    Returns:
        (Counter of what was read, dropped and written, {run: the mzML it was read from, or
        None}), each run also named with its file in the MGF's header.
    """
    counts = Counter()
    best, runs, scans = best_per_peptide(workdir, pattern, data_dirs, counts, limits.all_ranks)
    counts[NOT_BEST] = counts[READ] - len(best) - sum(counts[k] for k in PRE_PICK)
    picks = select(best, limits.per_protein, limits.max_spectra, counts)
    records = read_records(picks, counts, lambda p: record_keys(p, scans, limits.all_ranks))
    counts[WRITTEN] = len(records)
    ranks = "--all-ranks, every record marked RANK" if limits.all_ranks else "DDA+ rank 1 only"
    head = [
        f"export_mgf.py from {workdir}",
        f"--match {pattern.pattern if pattern else '(none, every protein kept)'}",
        f"--per-protein {limits.per_protein}, --max-spectra {limits.max_spectra}, {ranks}",
        f"mzML from {', '.join(data_dirs)}",
    ]
    head += [f"run {run}: {path or 'no mzML found'}" for run, path in sorted(runs.items())]
    head += [f"{counts[k]:,} {k}" for k in ORDER]
    head += [f"{n:,} {k}" for k, n in sorted(counts.items()) if k not in ORDER]
    head = [_comment(h) for h in head]
    head.insert(1, "# Here an equals sign is written %3D and a percent sign %25")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(head) + "\n\n" + "\n".join(records))
    return counts, runs


def _comment(text):
    """A header line with no "=": pyteomics makes a line with one "=" a param of every spectrum."""
    return "# " + text.replace("%", "%25").replace("=", "%3D")


def resolve_data_dirs(workdir, given):
    """--spectra-dir as given, else run_layout's rule: the input dirs in this run's manifest.

    Args:
        workdir: FragPipe results directory, <proj>/results/<name>.
        given: Directories from --spectra-dir, possibly empty.

    Returns:
        (directories, where they came from, or why there are none).
    """
    if given:
        return given, "--spectra-dir"
    run_dir = run_dir_for(workdir)
    if not run_dir:
        return [], "no launcher run directory, <proj>/runs/<name>, beside this results directory"
    manifest = os.path.join(run_dir, MANIFEST)
    dirs = input_dirs(load_manifest(manifest))
    return dirs, manifest if dirs else f"no input file named in {manifest}"


def _holds_spectra(folder):
    """Whether folder directly holds an .mzML file, as each name _find_mzml opens ends."""
    try:
        names = os.listdir(folder)
    except OSError:
        return False
    return any(n.endswith(".mzML") and os.path.isfile(os.path.join(folder, n)) for n in names)


def _spectra_dirs(p, args):
    """The folders each run's mzML is looked for in, or a stop that says what to pass instead."""
    data_dirs, source = resolve_data_dirs(args.workdir, args.spectra_dir)
    print(f"[export-mgf] results : {args.workdir}")
    print(f"[export-mgf] mzML in : {', '.join(data_dirs) or '(none)'}  [{source}]")
    hint = f"pass --spectra-dir DIR, the folder holding each <run>{MZML_SUFFIX} or <run>.mzML"
    if not data_dirs:
        p.error(f"{source}; {hint}")
    if not any(_holds_spectra(d) for d in data_dirs):
        p.error(f"no *.mzML directly in {', '.join(data_dirs)}; {hint}")
    return data_dirs


def _check_paths(p, args):
    """Stop on a path that cannot work before any reading, and before --out is overwritten."""
    if not os.path.isdir(args.workdir):
        p.error(f"--workdir {args.workdir} is not a directory")
    if not psm_tables(args.workdir):
        p.error(f"no psm.tsv in {args.workdir} or in any directory directly under it")
    out_dir = os.path.dirname(os.path.abspath(args.out))
    if not os.path.isdir(out_dir):
        p.error(f"the directory of --out, {out_dir}, does not exist")
    for folder in args.spectra_dir:
        if not os.path.isdir(folder):
            p.error(f"--spectra-dir {folder} is not a directory")


def require_mzml_reader():
    """Stop unless this Python's pyteomics can open an indexed mzML.

    read_spectrum turns every failure into "no spectrum", which is right for one bad
    scan in the report. Here a missing pyteomics, or pyteomics 5 without psims, would
    be counted as every scan unreadable and blamed on the mzML. The probe is a
    one-scan file psm_spectrum writes, opened the way every scan is, so a run reads
    no fixture and the network is no more asked here than for a scan.
    """
    try:
        from pyteomics import auxiliary, mzml  # noqa: PLC0415  optional, as in psm_spectrum
    except ImportError as exc:
        sys.exit(f"export_mgf.py: {sys.executable} cannot import pyteomics.mzml: {exc}")
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "probe.mzML")
        write_mzml(path, [(1, [100.0], [1.0])])
        try:
            with mzml_reader(path, mzml):
                pass
        except auxiliary.PyteomicsError as exc:
            sys.exit(f"export_mgf.py: {sys.executable} cannot read mzML: {exc}")


def _parser():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--workdir", help="FragPipe results dir holding <experiment>/psm.tsv")
    p.add_argument(
        "--spectra-dir",
        "--mzml-dir",
        action="append",
        default=[],
        metavar="DIR",
        help=f"folder holding each run's <run>{MZML_SUFFIX}, which MSFragger writes beside a "
        "Thermo .raw file, else its <run>.mzML, as a search of mzML files leaves; a Bruker .d "
        "folder is not read. <run> is the run in psm.tsv's Spectrum column. Repeatable. Given, "
        "no manifest is read. Default: the input folders the launcher's run manifest names, "
        f"{MANIFEST} in <proj>/runs/<name> beside <proj>/results/<name>. --mzml-dir is its old "
        "name",
    )
    p.add_argument(
        "--match",
        help=r"regex on Protein ID, Protein and Mapped Proteins, e.g. "
        r"'label=novel\b'. Match a label, not an accession prefix, which can also "
        "match annotated proteins through Mapped Proteins",
    )
    p.add_argument(
        "--per-protein",
        type=int,
        default=PER_PROTEIN,
        help=f"best-Hyperscore distinct peptides per protein (default {PER_PROTEIN})",
    )
    p.add_argument(
        "--max-spectra",
        type=int,
        default=MAX_SPECTRA,
        help=f"hard cap on spectra written (default {MAX_SPECTRA})",
    )
    p.add_argument(
        "--all-ranks",
        action="store_true",
        help="export DDA+ ranks 2 to 5 too, each record marked RANK= and with its scan's best "
        "non-decoy peptide; by default only rank 1 is exported",
    )
    p.add_argument("--out", help="the .mgf to write")
    p.add_argument("--selftest", action="store_true", help="run the self-test and exit")
    return p


def main(argv=None):
    """Parse the command line, export, and print every count.

    Args:
        argv: Arguments, default sys.argv[1:].

    Returns:
        0 when at least one spectrum was written, 1 when none was.
    """
    p = _parser()
    args = p.parse_args(argv)
    if args.selftest:
        import export_mgf_selftest  # noqa: PLC0415  test-only, and it imports this module

        export_mgf_selftest.run(sys.modules[__name__])
        return 0
    if not (args.workdir and args.out) or min(args.per_protein, args.max_spectra) < 1:
        p.error("--workdir and --out are required; --per-protein and --max-spectra must be >= 1")
    try:
        pattern = re.compile(args.match) if args.match else None
    except re.error as exc:
        p.error(f"--match is not a valid regex: {exc}")
    _check_paths(p, args)
    require_mzml_reader()
    data_dirs = _spectra_dirs(p, args)
    limits = Limits(args.per_protein, args.max_spectra, args.all_ranks)
    counts, runs = export(args.workdir, data_dirs, pattern, args.out, limits)
    for k in ORDER:
        print(f"[export-mgf] {counts[k]:>9,}  {k}")
    for k, n in sorted(counts.items()):
        if k not in ORDER:
            print(f"[export-mgf] {n:>9,}  {k}")
    for run, path in sorted(runs.items()):
        print(f"[export-mgf] run {run}: {path or 'no mzML found'}")
    print(f"[export-mgf] wrote {args.out}")
    return 0 if counts[WRITTEN] else 1


if __name__ == "__main__":
    sys.exit(main())
