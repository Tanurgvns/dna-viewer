"""Sequence analysis helpers for DNA Viewer.

Kept separate from the Streamlit UI so the logic can be unit-tested.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

VALID_BASES = set("ACGTN")

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

_BASES = "TCAG"
_AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {
    a + b + c: _AMINO[i]
    for i, (a, b, c) in enumerate((a, b, c) for a in _BASES for b in _BASES for c in _BASES)
}
STOP_CODONS = {"TAA", "TAG", "TGA"}

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}


@dataclass
class Record:
    header: str
    sequence: str


@dataclass
class ORF:
    strand: str        # "+" or "-"
    frame: int         # 1, 2 or 3
    start: int         # 1-based, forward-strand coordinate
    end: int           # 1-based, inclusive
    length_nt: int
    protein: str

    @property
    def length_aa(self) -> int:
        return len(self.protein)


def clean_sequence(raw: str) -> tuple[str, int]:
    """Uppercase, strip whitespace/digits, convert U to T.

    Returns the cleaned sequence and the number of invalid characters
    that were replaced with N.
    """
    seq = re.sub(r"[\s\d]", "", raw).upper().replace("U", "T")
    invalid = sum(1 for b in seq if b not in VALID_BASES)
    seq = "".join(b if b in VALID_BASES else "N" for b in seq)
    return seq, invalid


def parse_fasta(text: str) -> list[Record]:
    """Parse FASTA text. Plain sequence with no header becomes one record."""
    text = text.strip()
    if not text:
        return []
    if not text.startswith(">"):
        seq, _ = clean_sequence(text)
        return [Record("Pasted sequence", seq)] if seq else []

    records: list[Record] = []
    header, chunks = None, []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if header is not None:
                records.append(Record(header, clean_sequence("".join(chunks))[0]))
            header, chunks = line[1:].strip() or "Unnamed", []
        elif line and not line.startswith(";"):
            chunks.append(line)
    if header is not None:
        records.append(Record(header, clean_sequence("".join(chunks))[0]))
    return [r for r in records if r.sequence]


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)


def base_counts(seq: str) -> dict[str, int]:
    return {b: seq.count(b) for b in "ACGTN"}


def gc_content(seq: str) -> float:
    """GC percentage, ignoring N."""
    acgt = sum(seq.count(b) for b in "ACGT")
    if acgt == 0:
        return 0.0
    return 100 * (seq.count("G") + seq.count("C")) / acgt


def gc_skew(seq: str) -> float:
    g, c = seq.count("G"), seq.count("C")
    return 0.0 if g + c == 0 else (g - c) / (g + c)


def sliding_window(seq: str, window: int, step: int):
    """Yield (centre position, GC %, GC skew) for each window."""
    window = max(1, min(window, len(seq)))
    step = max(1, step)
    for i in range(0, len(seq) - window + 1, step):
        chunk = seq[i:i + window]
        yield i + window // 2 + 1, gc_content(chunk), gc_skew(chunk)


def translate(seq: str, frame: int = 0, stop_at_stop: bool = False) -> str:
    """Translate from a 0-based frame offset. Unknown codons become X."""
    protein = []
    for i in range(frame, len(seq) - 2, 3):
        aa = CODON_TABLE.get(seq[i:i + 3], "X")
        if aa == "*" and stop_at_stop:
            break
        protein.append(aa)
    return "".join(protein)


def six_frame_translation(seq: str) -> dict[str, str]:
    rc = reverse_complement(seq)
    frames = {}
    for f in range(3):
        frames[f"+{f + 1}"] = translate(seq, f)
    for f in range(3):
        frames[f"-{f + 1}"] = translate(rc, f)
    return frames


def find_orfs(seq: str, min_aa: int = 30) -> list[ORF]:
    """Find ATG-to-stop ORFs on both strands (non-nested, per frame)."""
    length = len(seq)
    orfs: list[ORF] = []
    for strand, s in (("+", seq), ("-", reverse_complement(seq))):
        for frame in range(3):
            i = frame
            while i <= len(s) - 3:
                if s[i:i + 3] == "ATG":
                    j = i
                    while j <= len(s) - 3 and s[j:j + 3] not in STOP_CODONS:
                        j += 3
                    if j > len(s) - 3:      # no stop codon before the end
                        break
                    protein = translate(s[i:j], 0)
                    if len(protein) >= min_aa:
                        end_excl = j + 3
                        if strand == "+":
                            start, end = i + 1, end_excl
                        else:
                            start, end = length - end_excl + 1, length - i
                        orfs.append(ORF(strand, frame + 1, start, end,
                                        end_excl - i, protein))
                    i = j + 3
                else:
                    i += 3
    return sorted(orfs, key=lambda o: o.length_nt, reverse=True)


def motif_to_regex(motif: str) -> str:
    motif = motif.upper().replace("U", "T")
    try:
        return "".join(IUPAC[ch] for ch in motif)
    except KeyError as exc:
        raise ValueError(f"'{exc.args[0]}' is not a DNA or IUPAC code") from exc


def find_motif(seq: str, motif: str, both_strands: bool = True) -> list[dict]:
    """Find overlapping motif matches. Positions are 1-based forward-strand."""
    pattern = re.compile(f"(?=({motif_to_regex(motif)}))")
    hits = [
        {"strand": "+", "start": m.start() + 1, "end": m.start() + len(motif),
         "match": m.group(1)}
        for m in pattern.finditer(seq)
    ]
    if both_strands:
        rc = reverse_complement(seq)
        n = len(seq)
        for m in pattern.finditer(rc):
            start = n - (m.start() + len(motif)) + 1
            hits.append({"strand": "-", "start": start,
                         "end": start + len(motif) - 1, "match": m.group(1)})
    return sorted(hits, key=lambda h: (h["start"], h["strand"]))


def molecular_weight_dsdna(seq: str) -> float:
    """Approximate double-stranded DNA weight in daltons (~617.96 Da/bp + 36)."""
    return len(seq) * 617.96 + 36.04
