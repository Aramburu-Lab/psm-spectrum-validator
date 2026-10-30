# PSM Spectrum Validator

A single, self-contained web page for judging peptide-spectrum matches (PSMs) by eye. It walks one annotated MS/MS
spectrum at a time through a short flowchart of questions: coverage by b and y ions, unannotated ions, the most
prominent ions, extra peaks, gaps and peaks beyond the end. It compares each of your answers with a rule-based
reference model.

## Using it

Download `spectrum_validator.html` and open it in a browser. Nothing is installed and nothing leaves your computer.

- **Load MGF** reads centroided spectra, one record per spectrum, with the peptide on a `SEQ=` line in ProForma with
  signed mass deltas, such as `[+42.0106]-LVNELTEFAK` or `PEPTM[+15.9949]K`.
- **Load controls** reads PSMs you trust, each with a `HYPERSCORE=` line, so every spectrum is set beside the controls
  nearest it in score.
- **Information** in the page shows the format and offers a small example file.
- **Export results to csv** saves every answer, mark and verdict.

## Feedback

Use **Help make this tool better** at the bottom of the page, or **Have a different opinion?** under any reply of the
reference model. Each saves a text note with the tool's version and settings. Post it as an
[issue here](https://github.com/Aramburu-Lab/psm-spectrum-validator/issues), or email it to rytis.stakaitis@umu.se.

## License

MIT, see [LICENSE](LICENSE).
