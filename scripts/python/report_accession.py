#!/usr/bin/env python3
# /// codebase-meta
# exempt: >-
#   How to read a protein identifier: the accession, gene and description behind whatever the
#   search tool wrote. Pure parsing, no I/O. Imported by run_report, report_fasta, report_groups,
#   report_quant and report_spectra; writes nothing.
# ///
"""One rule for reading a protein identifier, because two rules do not join.

The FASTA and the result tables have to agree on what an entry is called. Every
time they have not, the failure was silent: the report rendered, the tables
filled, and the entries the search existed to find were the ones that matched
nothing. It has happened three times, each with a different database:

    sp|P62805|H4_HUMAN            the accession is field 2
    alt|pid=1302903|sc=ATG|...    a custom id, kept WHOLE — splitting it
                                  unmatched 146 microproteins in a real run
    rs|GENE02080|HSP12_ORGN       a tag Philosopher does not know, so it wrote
      label=housekeeping_sanity   the ENTIRE header into the id column and
      pathway=chaperones ...      repeated it in gene and description

So this module is the only place that answers the question, and the self-test
carries one case from each of those databases.
"""

MIN_ACC_PARTS = 2  # sp|P12345|NAME -> the accession is field 2
DB_TAGS = ("sp", "tr", "rs")
DECOY_PREFIXES = ("rev_", "decoy_")
# Philosopher's database --contamprefix writes it before each cRAP header (dat.go lines 344 to 346).
# Its own psm.tsv and protein.tsv test the same prefix for Is Contaminant (psm.go line 160).
CONTAM_PREFIX = "contam_"


def fasta_accession(token):
    """Accession from a FASTA header, the way Philosopher assigns it.

    `sp|P62805|H4_HUMAN` -> P62805, but a custom identifier is kept WHOLE.
    Splitting on "|" unconditionally turned every alt-ORF entry —
    `alt|pid=1302903|sc=ATG|...` — into `pid=1302903`, which matches nothing in
    any result table, so 146 microproteins in a real run were silently
    unmatched: no species, no length, absent from the short-protein table. The
    entries the search exists to find were the ones the parser dropped.

    The whitespace split is not defensive tidying. For a database tag it does
    not recognise, Philosopher writes the ENTIRE header line into the protein
    table's id column — description and all — so the id arriving here is
    `rs|GENE02080|HSP12_ORGN label=housekeeping_sanity class=...`. Without
    this, that never matches the FASTA side, and a housekeeping chaperone
    rendered as a 127-character accession classed as a novel ORF.
    """
    first = token.split()[0] if token.split() else token
    parts = first.split("|")
    return parts[1] if parts[0] in DB_TAGS and len(parts) >= MIN_ACC_PARTS + 1 else first


def split_header(raw):
    """-> (accession, gene, description) from an id that may be a whole header.

    When Philosopher dumps the header into the id column, FragPipe repeats it in
    the gene and description columns too, so all three have to be recovered from
    the one string. The gene comes from the entry name — `HSP12_ORGN` -> HSP12 —
    which is only ever used as a fallback, because a database whose own gene
    column is populated gives a better answer than trimming a species suffix.
    """
    first = raw.split()[0] if raw.split() else raw
    parts = first.split("|")
    entry = parts[2] if len(parts) >= MIN_ACC_PARTS + 1 else ""
    # `alt|pid=34598|sc=ATG|...` has no entry name — field 3 is another key=value
    # of the identifier, and reading it as a gene labelled alt-ORF rows "sc=ATG".
    gene = "" if "=" in entry else (entry.rsplit("_", 1)[0] if "_" in entry else entry)
    rest = raw[len(first) :].strip()
    # `NOVEL_ORF` is a class, not a gene: reading it as one labelled every
    # novel ORF in the table "NOVEL". A stem that merely repeats the entry's own
    # label carries no per-entry information, so the accession is the better name.
    if gene and f"label={gene.lower()}" in rest.lower():
        gene = ""
    return fasta_accession(first), gene, rest


def accession_keys(token):
    """Every accession Philosopher might write for this entry. -> (primary, {keys}).

    A database with a tag this parser does not know — `rs|GENE00110|CCT8_ORGN`
    from one custom build — leaves a choice that cannot be made by
    reading the FASTA: Philosopher may key the protein table on the whole token
    or on the middle field, and if the report guesses the other one then nothing
    joins and every length, species and class silently goes missing.

    So it does not guess. Indexing under both costs one dict entry and makes the
    join correct either way.
    """
    keys = {token}
    parts = token.split("|")
    if len(parts) >= MIN_ACC_PARTS + 1 and parts[1]:
        keys.add(parts[1])
    return fasta_accession(token), keys


def is_decoy(token):
    """True for a decoy entry, wherever the builder put the tag.

    Philosopher writes `rev_sp|P12345|NAME`, prefixing the whole header. One
    custom build writes `rs|rev_XORF00064|DECOY_REV`, prefixing the
    accession field instead. Checking only the first token counted 6,239 decoys
    as targets and reported a database twice its real size.
    """
    return any(part.startswith(DECOY_PREFIXES) for part in token.split("|"))


def is_decoy_row(row):
    """True for a psm.tsv or protein-table row of a decoy, the one test every reader applies.

    FragPipe writes the entry as the FASTA names it under `Protein`, beside `Protein ID`, and
    Philosopher tags both for a decoy. report_demo writes `Protein ID` alone, where a test on
    `Protein` alone raised KeyError. A cell csv.DictReader left as None reads as empty.
    """
    return is_decoy(row.get("Protein") or row.get("Protein ID") or "")


def is_contaminant(token):
    """True for an entry Philosopher's --contamprefix marked, and for its decoy, rev_contam_.

    The one test of a contaminant. Without the prefix a cRAP entry is an ordinary sp| header,
    which is how hemoglobin was classed Globin in every short-protein table of one real run: a
    database meant for these readers is built with --contamprefix (decision D-F24).
    """
    first = token.split()[0] if token.split() else token
    for tag in DECOY_PREFIXES:
        if first.startswith(tag):
            first = first[len(tag) :]
            break
    return first.startswith(CONTAM_PREFIX)


def is_contaminant_row(row):
    """is_contaminant for a psm.tsv or protein-table row, read from Protein first.

    Philosopher reads Protein ID from the header past the prefix (getIDUniprot, db.go line 319),
    so a contaminant's Protein ID is the bare P68871 and only Protein carries the tag.
    Protein ID is read when Protein is absent or empty, as report_demo writes rows.
    """
    return is_contaminant(row.get("Protein") or row.get("Protein ID") or "")


def display_name(acc, gene):
    """A label short enough to read. Alt-ORF entries have no gene, only an id.

    `alt|pid=34598|orfs=5|pi=NP_001289194.1{OpenProt};XP_011533042.1{OpenProt}`
    is 140 characters of provenance around one identifier. Truncating it to fit
    keeps the prefix every such entry shares and cuts the part that tells them
    apart, so three entries render as the same name.
    """
    if gene:
        return gene
    for part in acc.split("|"):
        if part.startswith("pid="):
            return part
    return acc


def _selftest_headers():
    """Accession, gene and description out of whatever the tool wrote."""
    # Philosopher wrote the WHOLE header into the id column for the rs| tag, so
    # a housekeeping chaperone classed as a novel ORF with a 127-char accession.
    hk = (
        "rs|GENE02080|HSP12_ORGN label=housekeeping_sanity class=annotated_cds "
        "pathway=chaperones OS=Example organism OX=0 PE=1 SV=1"
    )
    assert fasta_accession(hk) == "GENE02080", fasta_accession(hk)
    acc, gene, desc = split_header(hk)
    assert (acc, gene) == ("GENE02080", "HSP12"), (acc, gene)
    assert desc.startswith("label=housekeeping_sanity") and "|" not in desc

    # A class name is not a gene: reading NOVEL_ORF as one labelled every
    # novel ORF in the table "NOVEL".
    assert split_header("rs|XORF01762|NOVEL_ORF label=novel aa=72")[1] == ""
    # A custom identifier has no entry name to take a gene from.
    assert split_header("alt|pid=34598|sc=ATG|orfs=2")[1] == ""
    assert display_name("alt|pid=34598|sc=ATG|orfs=2", "") == "pid=34598"

    # The rules that already existed must not move.
    assert fasta_accession("sp|P05109|S10A8_HUMAN") == "P05109"
    assert fasta_accession("alt|pid=7|sc=ATG|x") == "alt|pid=7|sc=ATG|x"
    assert fasta_accession("contam_sp|P00761|TRYP_PIG").startswith("contam")
    assert is_contaminant("contam_sp|P00761|TRYP_PIG") and is_contaminant("rev_contam_sp|P1|X")
    assert not is_contaminant("sp|P00761|TRYP_PIG") and not is_contaminant("rs|XORF00001|X")
    assert is_contaminant_row({"Protein": "contam_sp|P68871|HBB_HUMAN", "Protein ID": "P68871"})
    assert not is_contaminant_row({"Protein": "sp|P68871|HBB_HUMAN", "Protein ID": "P68871"})
    assert is_decoy("rev_sp|P05109|X") and is_decoy("rs|rev_XORF00064|DECOY_REV")
    assert not is_decoy("rs|GENE02080|HSP12_ORGN")
    # Protein first, Protein ID when Protein is absent or empty, as report_demo writes it.
    assert is_decoy_row({"Protein": "rev_sp|P05109|X", "Protein ID": "P05109"})
    assert not is_decoy_row({"Protein": "sp|P05109|X", "Protein ID": "rev_P05109"})
    assert is_decoy_row({"Protein ID": "rev_P05109"}) and is_decoy_row(
        {"Protein": "", "Protein ID": "rev_P05109"}
    )
    assert not is_decoy_row({"Protein ID": "P05109"}) and not is_decoy_row({"Protein": None})
    assert accession_keys("rs|GENE02080|HSP12_ORGN")[1] == {
        "rs|GENE02080|HSP12_ORGN",
        "GENE02080",
    }
    print("report_render header self-test ok")


if __name__ == "__main__":
    _selftest_headers()
