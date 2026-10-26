# 🧬 DNA Viewer

A small Streamlit app that turns FASTA files and raw DNA sequences into something easier to read than endless lines of A, T, G and C.

## Features

- **Double-helix view:** a scrolling B-DNA helix with base pairs, hydrogen bonds, sugar-phosphate backbone, minor and major grooves, and turn markers. A reading head scans along the sequence (or click, scroll, or use the arrow keys), with a live sequence window and helix facts for the selected pair.
- **Input:** upload a FASTA file (multi-record supported), paste a sequence, or try the built-in example. RNA input (U) is converted to T automatically.
- **Overview:** length, GC content, GC skew, approximate dsDNA weight and base composition charts.
- **Coloured sequence view:** colour-coded bases with position numbers, forward or reverse-complement strand, optional double-stranded view, and FASTA download.
- **GC content:** sliding-window GC % and GC skew plots with adjustable window and step, plus cumulative skew.
- **ORF finder:** ATG-to-stop open reading frames on both strands and all six frames, shown on a map and a table, with protein FASTA download.
- **Six-frame translation:** with start (M) and stop (*) residues highlighted.
- **Motif search:** IUPAC codes supported (e.g. `TATAWAW`), overlapping matches, both strands, with matches outlined in context.

## Run locally

```bash
git clone https://github.com/Tanurgvns/dna-viewer.git
cd dna-viewer
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
streamlit run app.py
```

Then open http://localhost:8501.

## Run the tests

```bash
pip install pytest
pytest
```

## Project structure

```
dna-viewer/
├── app.py              # Streamlit interface
├── helix_view.html     # Double-helix viewer (HTML/JavaScript)
├── helix_view.py       # Puts the sequence into the helix viewer
├── dna_utils.py        # Sequence analysis (parsing, ORFs, GC, translation, motifs)
├── tests/              # Unit tests for dna_utils
├── requirements.txt
└── .streamlit/config.toml
```

## Notes

- Uses the standard genetic code (NCBI table 1).
- ORFs are reported from the first ATG to the next in-frame stop codon; nested ORFs are not listed separately.
- The coloured view shows up to 6,000 bases at a time; use the slider to move through longer sequences.

## Roadmap

- Circular plasmid map
- Restriction enzyme presets
- Codon usage table
- Alternative genetic codes (e.g. mitochondrial)

## License

MIT
