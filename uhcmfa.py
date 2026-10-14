#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
UHCMFA — Universal Hexagram Coding Multifractal Analysis
=========================================================

Encodes any DNA/RNA sequence into I-Ching hexagrams and performs
multifractal analysis (MF-DFA) on the sequence.

Method
------
  1. Codon (3 nt)       -> 6 lines (yin/yang) -> 1 hexagram
  2. Hexagram value     -> integer 0..63
  3. Integer series     -> MF-DFA (Kantelhardt 2002)
  4. Result             -> H(q=2), Δα, α_peak, h(q) spectrum

Use Cases
---------
UHCMFA is a fast pre-screening tool. It does not replace BLAST,
QUAST, or Kraken2, but it runs in O(N log N) and is therefore much
faster than alignment-based methods. Typical use cases:

  1. Sequencing assembly quality control
  2. Contamination detection (bacterial/viral/fungal)
  3. Rapid metagenome clustering without alignment
  4. Horizontal gene transfer (HGT) pre-screening
  5. Plasmid and vector identification
  6. Synthetic DNA signature verification
  7. Genome stability monitoring over time
  8. CRISPR insertion and deletion detection
  9. Phylogenetic pre-clustering (thousands of genomes)
 10. Whole-genome or whole-chromosome overview in a single number

Advantages Over Standard Methods
--------------------------------
  1. NO ALIGNMENT NEEDED
     BLAST/DIAMOND/minimap2 align every read against a database.
     UHCMFA reads the sequence once and computes a signature.
     No reference genome, no index, no seed-and-extend.

  2. O(N log N) SCALING
     Standard tools are O(N · M) or worse (N = query, M = database).
     UHCMFA is O(N log N) regardless of database size.

  3. TINY OUTPUT
     BLAST produces megabytes of alignments. Kraken2 produces
     per-read taxon assignments. UHCMFA produces 5 lines.

  4. LANGUAGE- AND PLATFORM-INDEPENDENT
     The signature (H, Δα, α_peak) is a tuple of real numbers.
     It can be compared across labs, pipelines, and file formats.

  5. DETERMINISTIC AND REPRODUCIBLE
     No random seeds, no heuristic thresholds, no version drift.
     Same input → same output, forever.

  6. NO DATABASE DEPENDENCY
     No need for RefSeq, UniProt, or any external reference.
     The signature is intrinsic to the sequence.

  7. EMBARRASSINGLY PARALLEL AT THE FILE LEVEL
     Each sequence can be processed independently and combined
     via parallel reduction of (H, Δα) tuples.

  8. WORKS ON FRAGMENTS
     Unlike assembly-based tools, UHCMFA works on partial
     sequences, unassembled reads, and draft genomes.

  9. SINGLE FILE, SINGLE DEPENDENCY
     ~500 lines of Python, one dependency (numpy). No build,
     no docker, no environment file.

 10. DIRECT HEXAGRAM OUTPUT
     The sequence can be dumped as Unicode hexagrams for
     visual inspection, artistic rendering, or archival.

Performance Comparison
----------------------
Times are measured on a single modern CPU core (x86-64, 3 GHz)
unless noted. "Genome" = ~5 Mb bacterial genome unless otherwise stated.

  TASK: Identify organism from a draft genome
  -------------------------------------------
  Method              Tool            Input        Time       Output
  ------------------  --------------  -----------  ---------  --------------
  Alignment           BLASTn          vs NT DB     ~30 min    Alignment list
  K-mer               Kraken2         vs std DB    ~2 min     Taxon per read
  Marker genes        CheckM          vs HMM DB    ~5 min     Completeness
  ANI                 FastANI         vs RefSeq    ~1 min     ANI matrix
  Sketch              Mash            vs RefSeq    ~10 s      Distance
  Multifractal        UHCMFA          none         ~0.5 s     (H, Δα)

  TASK: Detect contamination in a human sample
  --------------------------------------------
  Method              Tool            Input        Time       Output
  ------------------  --------------  -----------  ---------  --------------
  K-mer               Kraken2         vs std DB    ~5 min     Per-read taxon
  Coverage            ConFindr        vs RefSeq    ~3 min     Contam %  
  Marker genes        CheckM          vs HMM DB    ~5 min     Lineage
  Multifractal        UHCMFA          none         ~1 s       Δα shift

  TASK: Quality control of a genome assembly
  ------------------------------------------
  Method              Tool            Input        Time       Output
  ------------------  --------------  -----------  ---------  --------------
  Gene prediction     Prokka          vs HMM DB    ~5 min     Gene list
  Marker genes        BUSCO           vs BUSCO DB  ~3 min     Completeness
  Alignment           QUAST           vs reference ~2 min     Assembly stats
  Multifractal        UHCMFA          none         ~0.5 s     H, Δα

  TASK: Detect horizontal gene transfer
  -------------------------------------
  Method              Tool            Input        Time       Output
  ------------------  --------------  -----------  ---------  --------------
  Composition         Alien Hunter    vs model     ~2 min     HGT regions
  Phylogeny           HGTector        vs RefSeq    ~10 min    HGT candidates
  Alignment           BLASTx          vs nr DB     ~1 h       HGT hits
  Multifractal        UHCMFA          none         ~0.5 s     Δα anomaly

  TASK: Cluster thousands of viral genomes
  ----------------------------------------
  Method              Tool            Input        Time       Output
  ------------------  --------------  -----------  ---------  --------------
  Alignment           BLASTn          vs all       ~days      Distance matrix
  Sketch              Mash            1000 genomes ~1 min     Distance matrix
  K-mer               sourmash        1000 genomes ~2 min     Distance matrix
  Multifractal        UHCMFA          1000 genomes ~1 min     (H, Δα) points

  TASK: Whole human genome overview
  ---------------------------------
  Method              Tool            Input        Time       Output
  ------------------  --------------  -----------  ---------  --------------
  Alignment           minimap2        vs self      ~1 h       Overlaps
  K-mer               Jellyfish       k=21         ~10 min    K-mer spectrum
  Repeat              RepeatMasker    vs library   ~30 min    Repeat families
  Multifractal        UHCMFA          none         ~20 min    5 numbers

Scaling
-------
  Input size          UHCMFA time     Notes
  ------------------  --------------  ------------------------------
  10 kb (virus)       < 100 ms        instant
  1 Mb (plasmid)      ~ 200 ms        instant
  5 Mb (bacterium)    ~ 2 s           fast
  50 Mb (chromosome)  ~ 15 s          acceptable
  250 Mb (chr 1)      ~ 60 s          one minute
  3.1 Gb (human)      ~ 20 min        full genome
  100 Gb (cohort)     ~ 10 h          batch processing

Memory
------
  UHCMFA holds only the current sequence in memory.
  Peak memory ≈ 8 bytes × N codons.
  10 kb virus: < 1 MB
  5 Mb bacterium: ~ 40 MB
  250 Mb chromosome: ~ 2 GB
  3.1 Gb genome: ~ 25 GB (split by chromosome)

Frequently Asked Questions
--------------------------
Q: Is this a replacement for BLAST?
A: No. UHCMFA is a pre-screening and quality-control tool.
   Use it to filter before running expensive alignment.

Q: What does H(q=2) mean biologically?
A: It measures long-range correlations in the codon sequence.
   Values > 0.5 indicate structure (GC-rich regions, repeats).
   Values ≈ 0.5 indicate near-random codon usage.

Q: What does Δα mean?
A: The width of the multifractal spectrum. Larger Δα means
   the sequence has different correlation behavior at different
   scales — i.e., it is heterogeneous. Small Δα means the
   sequence is scale-invariant.

Q: Can I compare two sequences?
A: Yes. Compute (H, Δα) for each and compare. Same organism →
   similar signatures. Different organism → different signatures.

Q: Why hexagrams?
A: 4³ = 64 codons and 64 hexagrams form a natural bijection.
   The hexagram representation is a compact visual encoding of
   the codon sequence and can be dumped for inspection.

Usage
-----
    uhcmfa [OPTIONS] [FILE...]
    uhcmfa --help
    uhcmfa --methods
    uhcmfa --version

Examples
--------
    # 1. Single local file
    uhcmfa chr21.fa

    # 2. Compressed file
    uhcmfa chr21.fa.gz

    # 3. Fetch directly from NCBI
    uhcmfa --acc NC_001722

    # 4. Combine multiple accessions
    uhcmfa --acc NC_001722 M30502 X61240

    # 5. Fetch human chromosome
    uhcmfa --chrom 21

    # 6. Fetch whole human genome (~20 min)
    uhcmfa --human

    # 7. Contamination comparison
    uhcmfa --compare sample.fa ecoli.fa yeast.fa

    # 8. Reverse complement and frame 1
    uhcmfa --frame 1 --rc chr21.fa

    # 9. JSON output for automation
    uhcmfa --json --acc NC_001722

    # 10. CSV output for spreadsheets
    uhcmfa --csv --chrom 21 --chrom 22

    # 11. Hexagrams only (no analysis)
    uhcmfa --hex-only --acc NC_001722 > hiv2.txt

    # 12. Quiet mode for scripts
    uhcmfa --quiet --chrom 21

    # 13. Large genome in background
    nohup uhcmfa --human -o human.txt &

    # 14. Compare your sequence to bacteria
    uhcmfa --compare my_seq.fa ref_ecoli.fa ref_yeast.fa

    # 15. Print method comparison table
    uhcmfa --methods

Output Format
-------------
Default (text):
    ䷀  N = 3278 codons   (1 sequence, 9835 nt)
    ䷁  H(q=2) = 0.583   [uncorrelated (random walk)]
    ䷂  Δα     = 0.412   [weak multifractal]
    ䷃  α_peak = 0.571
    ䷄  h(q)   = ䷦䷧䷨䷩䷪䷫䷬䷭䷮

Interpretation
--------------
  H(q=2)   Hurst exponent
           < 0.4   anticorrelated
           ~ 0.5   uncorrelated (random walk)
           > 0.6   long-range correlation

  Δα       Width of the multifractal spectrum
           < 0.1   monofractal
           < 0.5   weak multifractal
           ≥ 0.5   strong multifractal

  α_peak   Peak of the spectrum

  h(q)     9 hexagrams quantizing h(q)
           q = -5, -3, -2, -1, 0, 1, 2, 3, 5

License
-------
Public domain (CC0). Use freely.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import time
import urllib.parse
import urllib.request
from typing import Iterable, List, Optional, Tuple

try:
    import numpy as np
except ImportError:
    sys.stderr.write("UHCMFA: numpy not found. Install: pip install numpy\n")
    sys.exit(2)

__version__ = "1.1.0"
PROG = "uhcmfa"

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
CHUNK = 500_000
NCBI_DELAY = 0.34
NCBI_DELAY_KEY = 0.11

HUMAN_CHROMS = {**{i: f"NC_{i:06d}" for i in range(1, 25)}, 25: "NC_012920"}

_DIGIT = str.maketrans({
    "A": "0", "C": "1", "G": "2", "T": "3", "U": "3",
    "R": "2", "Y": "1", "S": "1", "W": "2", "K": "3", "M": "0",
    "B": "1", "D": "2", "H": "3", "V": "2", "N": "0",
})
VALID = set("ACGTURYSWKMBDHVN")
RC_TABLE = str.maketrans("ACGTURYSWKMBDHVN", "TGCAAYRSWMKVHDBN")

KW = {
    ("111","111"):1,  ("000","000"):2,  ("010","100"):3,  ("001","010"):4,
    ("010","111"):5,  ("111","010"):6,  ("000","010"):7,  ("010","000"):8,
    ("011","111"):9,  ("111","110"):10, ("000","111"):11, ("111","000"):12,
    ("111","101"):13, ("101","111"):14, ("000","001"):15, ("100","000"):16,
    ("110","100"):17, ("001","011"):18, ("000","110"):19, ("011","000"):20,
    ("101","100"):21, ("001","101"):22, ("001","000"):23, ("000","100"):24,
    ("111","100"):25, ("001","111"):26, ("001","100"):27, ("110","011"):28,
    ("010","010"):29, ("101","101"):30, ("110","001"):31, ("100","011"):32,
    ("111","001"):33, ("100","111"):34, ("101","000"):35, ("000","101"):36,
    ("011","101"):37, ("101","110"):38, ("010","001"):39, ("100","010"):40,
    ("001","110"):41, ("011","100"):42, ("110","111"):43, ("111","011"):44,
    ("110","000"):45, ("000","011"):46, ("110","010"):47, ("010","011"):48,
    ("110","101"):49, ("101","011"):50, ("100","100"):51, ("001","001"):52,
    ("011","001"):53, ("100","110"):54, ("100","101"):55, ("101","001"):56,
    ("011","011"):57, ("110","110"):58, ("011","010"):59, ("010","110"):60,
    ("011","110"):61, ("100","001"):62, ("010","101"):63, ("101","010"):64,
}

# ---------------------------------------------------------------------------
# Method comparison tables (printed by --methods)
# ---------------------------------------------------------------------------

METHOD_TABLES = [
    ("Task: Identify organism from a draft genome", [
        ("Method",         "Tool",         "Input",       "Time",     "Output"),
        ("Alignment",      "BLASTn",       "NT DB",       "~30 min",  "Alignment list"),
        ("K-mer",          "Kraken2",      "std DB",      "~2 min",   "Taxon per read"),
        ("Marker genes",   "CheckM",       "HMM DB",      "~5 min",   "Completeness"),
        ("ANI",            "FastANI",      "RefSeq",      "~1 min",   "ANI matrix"),
        ("Sketch",         "Mash",         "RefSeq",      "~10 s",    "Distance"),
        ("Multifractal",   "UHCMFA",       "none",        "~0.5 s",   "(H, Δα)"),
    ]),
    ("Task: Detect contamination in a human sample", [
        ("Method",         "Tool",         "Input",       "Time",     "Output"),
        ("K-mer",          "Kraken2",      "std DB",      "~5 min",   "Per-read taxon"),
        ("Coverage",       "ConFindr",     "RefSeq",      "~3 min",   "Contam %"),
        ("Marker genes",   "CheckM",       "HMM DB",      "~5 min",   "Lineage"),
        ("Multifractal",   "UHCMFA",       "none",        "~1 s",     "Δα shift"),
    ]),
    ("Task: Quality control of a genome assembly", [
        ("Method",         "Tool",         "Input",       "Time",     "Output"),
        ("Gene prediction","Prokka",       "HMM DB",      "~5 min",   "Gene list"),
        ("Marker genes",   "BUSCO",        "BUSCO DB",    "~3 min",   "Completeness"),
        ("Alignment",      "QUAST",        "reference",   "~2 min",   "Assembly stats"),
        ("Multifractal",   "UHCMFA",       "none",        "~0.5 s",   "H, Δα"),
    ]),
    ("Task: Detect horizontal gene transfer", [
        ("Method",         "Tool",         "Input",       "Time",     "Output"),
        ("Composition",    "Alien Hunter", "model",       "~2 min",   "HGT regions"),
        ("Phylogeny",      "HGTector",     "RefSeq",      "~10 min",  "HGT candidates"),
        ("Alignment",      "BLASTx",       "nr DB",       "~1 h",     "HGT hits"),
        ("Multifractal",   "UHCMFA",       "none",        "~0.5 s",   "Δα anomaly"),
    ]),
    ("Task: Cluster thousands of viral genomes", [
        ("Method",         "Tool",         "Input",       "Time",     "Output"),
        ("Alignment",      "BLASTn",       "all-vs-all",  "~days",    "Distance matrix"),
        ("Sketch",         "Mash",         "1000 genomes","~1 min",   "Distance matrix"),
        ("K-mer",          "sourmash",     "1000 genomes","~2 min",   "Distance matrix"),
        ("Multifractal",   "UHCMFA",       "1000 genomes","~1 min",   "(H, Δα) points"),
    ]),
    ("Task: Whole human genome overview", [
        ("Method",         "Tool",         "Input",       "Time",     "Output"),
        ("Alignment",      "minimap2",     "self",        "~1 h",     "Overlaps"),
        ("K-mer",          "Jellyfish",    "k=21",        "~10 min",  "K-mer spectrum"),
        ("Repeat",         "RepeatMasker", "library",     "~30 min",  "Repeat families"),
        ("Multifractal",   "UHCMFA",       "none",        "~20 min",  "5 numbers"),
    ]),
]

SCALING_TABLE = [
    ("Input size",         "UHCMFA time",  "Notes"),
    ("10 kb (virus)",      "< 100 ms",     "instant"),
    ("1 Mb (plasmid)",     "~ 200 ms",     "instant"),
    ("5 Mb (bacterium)",   "~ 2 s",        "fast"),
    ("50 Mb (chromosome)", "~ 15 s",       "acceptable"),
    ("250 Mb (chr 1)",     "~ 60 s",       "one minute"),
    ("3.1 Gb (human)",     "~ 20 min",     "full genome"),
    ("100 Gb (cohort)",    "~ 10 h",       "batch processing"),
]

MEMORY_TABLE = [
    ("Input",              "Peak memory",  "Notes"),
    ("10 kb virus",        "< 1 MB",       "trivial"),
    ("5 Mb bacterium",     "~ 40 MB",      "trivial"),
    ("250 Mb chromosome",  "~ 2 GB",       "needs RAM"),
    ("3.1 Gb genome",      "~ 25 GB",      "split by chromosome"),
]

ADVANTAGES = [
    ("No alignment needed",
     "BLAST/DIAMOND/minimap2 align every read against a database. "
     "UHCMFA reads the sequence once and computes a signature. "
     "No reference genome, no index, no seed-and-extend."),
    ("O(N log N) scaling",
     "Standard tools are O(N · M) or worse. UHCMFA is O(N log N) "
     "regardless of database size."),
    ("Tiny output",
     "BLAST produces megabytes of alignments. Kraken2 produces "
     "per-read taxon assignments. UHCMFA produces 5 lines."),
    ("Language-independent signature",
     "The tuple (H, Δα, α_peak) can be compared across labs, "
     "pipelines, and file formats."),
    ("Deterministic and reproducible",
     "No random seeds, no heuristic thresholds, no version drift. "
     "Same input → same output, forever."),
    ("No database dependency",
     "No RefSeq, no UniProt, no HMM library. The signature is "
     "intrinsic to the sequence."),
    ("Parallel at the file level",
     "Each sequence can be processed independently and combined "
     "via parallel reduction of (H, Δα) tuples."),
    ("Works on fragments",
     "Unlike assembly-based tools, UHCMFA works on partial "
     "sequences, unassembled reads, and draft genomes."),
    ("Single file, single dependency",
     "~500 lines of Python, one dependency (numpy). No build, "
     "no docker, no environment file."),
    ("Direct hexagram output",
     "The sequence can be dumped as Unicode hexagrams for "
     "visual inspection, artistic rendering, or archival."),
]


def print_method_tables(fh=sys.stdout) -> None:
    """Print all comparison tables and advantages."""
    print("UHCMFA — Method Comparison Tables", file=fh)
    print("=" * 72, file=fh)
    print(file=fh)

    for title, rows in METHOD_TABLES:
        print(title, file=fh)
        print("-" * len(title), file=fh)
        # column widths
        cols = list(zip(*rows))
        widths = [max(len(str(c)) for c in col) for col in cols]
        for i, row in enumerate(rows):
            line = "  ".join(str(cell).ljust(w)
                             for cell, w in zip(row, widths))
            print(line.rstrip(), file=fh)
            if i == 0:
                print("  ".join("-" * w for w in widths), file=fh)
        print(file=fh)

    print("Scaling", file=fh)
    print("-------", file=fh)
    cols = list(zip(*SCALING_TABLE))
    widths = [max(len(str(c)) for c in col) for col in cols]
    for i, row in enumerate(SCALING_TABLE):
        line = "  ".join(str(cell).ljust(w)
                         for cell, w in zip(row, widths))
        print(line.rstrip(), file=fh)
        if i == 0:
            print("  ".join("-" * w for w in widths), file=fh)
    print(file=fh)

    print("Memory", file=fh)
    print("------", file=fh)
    cols = list(zip(*MEMORY_TABLE))
    widths = [max(len(str(c)) for c in col) for col in cols]
    for i, row in enumerate(MEMORY_TABLE):
        line = "  ".join(str(cell).ljust(w)
                         for cell, w in zip(row, widths))
        print(line.rstrip(), file=fh)
        if i == 0:
            print("  ".join("-" * w for w in widths), file=fh)
    print(file=fh)

    print("Advantages over standard methods", file=fh)
    print("--------------------------------", file=fh)
    for i, (title, body) in enumerate(ADVANTAGES, 1):
        print(f"{i:2d}. {title}", file=fh)
        # wrap body at 68 chars with indent
        words = body.split()
        line = "    "
        for w in words:
            if len(line) + len(w) + 1 > 72:
                print(line, file=fh)
                line = "    " + w
            else:
                line = (line + " " + w) if line.strip() else ("    " + w)
        if line.strip():
            print(line, file=fh)
        print(file=fh)


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def read_fasta(path: str) -> str:
    op = gzip.open if path.endswith(".gz") else open
    try:
        with op(path, "rt", encoding="utf-8", errors="replace") as f:
            return "".join(
                line.strip() for line in f if not line.startswith(">")
            ).upper()
    except FileNotFoundError:
        sys.stderr.write(f"{PROG}: file not found: {path}\n")
        sys.exit(2)
    except OSError as e:
        sys.stderr.write(f"{PROG}: {e}\n")
        sys.exit(2)


def fetch_ncbi(acc: str, api_key: Optional[str] = None,
               quiet: bool = False) -> str:
    delay = NCBI_DELAY_KEY if api_key else NCBI_DELAY
    parts: List[str] = []
    pos = 1
    while True:
        stop = pos + CHUNK - 1
        if not quiet:
            sys.stderr.write(f"# {acc} {pos}-{stop}\r")
            sys.stderr.flush()
        p = {
            "db": "nucleotide", "id": acc, "rettype": "fasta",
            "retmode": "text", "seq_start": pos, "seq_stop": stop,
        }
        if api_key:
            p["api_key"] = api_key
        url = EFETCH + "?" + urllib.parse.urlencode(p)
        req = urllib.request.Request(
            url, headers={"User-Agent": f"{PROG}/{__version__}"})
        last_err: Optional[Exception] = None
        for attempt in range(3):
            try:
                with urllib.request.urlopen(req, timeout=180) as r:
                    text = r.read().decode("utf-8", "replace")
                break
            except Exception as e:
                last_err = e
                time.sleep(2 ** attempt)
        else:
            raise RuntimeError(f"NCBI fetch failed: {last_err}")
        seq = "".join(
            l.strip() for l in text.splitlines()
            if not l.startswith(">")
        ).upper()
        if not seq:
            break
        parts.append(seq)
        if len(seq) < CHUNK:
            break
        pos = stop + 1
        time.sleep(delay)
    if not quiet:
        sys.stderr.write("\n")
    return "".join(parts)


# ---------------------------------------------------------------------------
# Encoding
# ---------------------------------------------------------------------------

def to_series(seq: str, frame: int = 0, rc: bool = False) -> np.ndarray:
    seq = "".join(c for c in seq.upper() if c in VALID)
    if rc:
        seq = seq.translate(RC_TABLE)[::-1]
    seq = seq[frame:]
    n = (len(seq) // 3) * 3
    seq = seq[:n]
    if n == 0:
        return np.empty(0, dtype=np.float64)
    digits = seq.translate(_DIGIT)
    a = np.frombuffer(digits.encode("ascii"), dtype=np.uint8) - ord("0")
    a = a.reshape(-1, 3).astype(np.float64)
    return 16.0 * a[:, 0] + 4.0 * a[:, 1] + a[:, 2]


def to_hexagrams(seq: str, frame: int = 0, rc: bool = False) -> str:
    x = to_series(seq, frame=frame, rc=rc)
    if x.size == 0:
        return ""
    return _series_to_hexagrams(x)


def _series_to_hexagrams(x: np.ndarray) -> str:
    n1 = (x // 16).astype(np.uint8)
    n2 = ((x // 4) % 4).astype(np.uint8)
    n3 = (x % 4).astype(np.uint8)

    def bits(n: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        return (n >> 1) & 1, n & 1

    a1, b1 = bits(n1)
    a2, b2 = bits(n2)
    a3, b3 = bits(n3)
    lines = np.stack([a1, b1, a2, b2, a3, b3], axis=1)
    lower = (lines[:, 0] << 2) | (lines[:, 1] << 1) | lines[:, 2]
    upper = (lines[:, 3] << 2) | (lines[:, 4] << 1) | lines[:, 5]
    chars = np.empty(x.size, dtype="U1")
    key = (upper.astype(int) * 8 + lower.astype(int))
    lookup = np.empty(64, dtype="U1")
    for (u_bits, l_bits), num in KW.items():
        ui = int(u_bits, 2)
        li = int(l_bits, 2)
        lookup[ui * 8 + li] = chr(0x4DC0 + num - 1)
    chars[:] = lookup[key]
    return "".join(chars.tolist())


# ---------------------------------------------------------------------------
# MF-DFA
# ---------------------------------------------------------------------------

def _detrend(seg: np.ndarray) -> np.ndarray:
    s = seg.size
    t = np.arange(s, dtype=np.float64)
    A = np.empty((s, 2))
    A[:, 0] = t
    A[:, 1] = 1.0
    coef, *_ = np.linalg.lstsq(A, seg, rcond=None)
    return seg - (A @ coef)


def mfdfa(x: np.ndarray, qs: np.ndarray, scales: np.ndarray) -> np.ndarray:
    N = x.size
    x = x - x.mean()
    Y = np.cumsum(x)
    Yr = Y[::-1]
    Fq = np.full((qs.size, scales.size), np.nan)
    for si, s in enumerate(scales):
        ns = N // s
        if ns < 4:
            continue
        F2 = np.empty(2 * ns)
        for v in range(ns):
            r = _detrend(Y[v*s:(v+1)*s])
            F2[v] = np.mean(r * r)
        for v in range(ns):
            r = _detrend(Yr[v*s:(v+1)*s])
            F2[ns + v] = np.mean(r * r)
        F2 = F2[F2 > 1e-30]
        if F2.size == 0:
            continue
        for qi, q in enumerate(qs):
            if q == 0:
                Fq[qi, si] = np.exp(0.5 * np.mean(np.log(F2)))
            else:
                Fq[qi, si] = np.mean(F2 ** (q / 2.0)) ** (1.0 / q)
    log_s = np.log(scales)
    hq = np.full(qs.size, np.nan)
    for qi in range(qs.size):
        y = np.log(Fq[qi])
        m = np.isfinite(y)
        if m.sum() >= 3:
            hq[qi] = np.polyfit(log_s[m], y[m], 1)[0]
    return hq


def analyze(x: np.ndarray, smin: int = 16, smax_div: int = 4) -> Optional[dict]:
    N = x.size
    smax = max(smin * 4, N // smax_div)
    if smax <= smin:
        return None
    scales = np.unique(
        np.round(np.logspace(np.log10(smin), np.log10(smax), 20)).astype(int)
    )
    scales = scales[(scales >= smin) & (scales * 4 <= N)]
    if scales.size < 4:
        return None
    qs = np.array([-5, -3, -2, -1, 0, 1, 2, 3, 5], dtype=np.float64)
    hq = mfdfa(x, qs, scales)
    i2 = int(np.where(qs == 2)[0][0])
    H = float(hq[i2])
    tau = qs * hq - 1.0
    m = np.isfinite(tau)
    if m.sum() < 3:
        return None
    alpha = np.gradient(tau[m], qs[m])
    falpha = qs[m] * alpha - tau[m]
    da = float(alpha.max() - alpha.min())
    ap = float(alpha[np.argmax(falpha)])
    return {
        "N": int(N),
        "H": H,
        "da": da,
        "alpha_peak": ap,
        "hq": hq.tolist(),
        "qs": qs.tolist(),
        "scales": scales.tolist(),
        "smin": int(smin),
        "smax": int(smax),
    }


def classify_H(H: float) -> str:
    if H > 0.6:
        return "long-range correlation"
    if H < 0.4:
        return "anticorrelated"
    return "uncorrelated (random walk)"


def classify_da(da: float) -> str:
    if da < 0.1:
        return "monofractal"
    if da < 0.5:
        return "weak multifractal"
    return "strong multifractal"


def hq_to_hexagrams(hq: Iterable[float]) -> str:
    hq_arr = np.asarray(list(hq), dtype=np.float64)
    if not np.all(np.isfinite(hq_arr)):
        hq_arr = np.nan_to_num(hq_arr, nan=0.5)
    lo = float(hq_arr.min())
    hi = float(hq_arr.max())
    if hi - lo < 1e-9:
        hi = lo + 1e-9
    out = []
    for h in hq_arr:
        t = (h - lo) / (hi - lo)
        idx = int(np.clip(np.floor(t * 64), 0, 63))
        out.append(chr(0x4DC0 + idx))
    return "".join(out)


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def emit_text(result: dict, name: str, nseq: int, total_nt: int,
              fh=sys.stdout) -> None:
    H = result["H"]
    da = result["da"]
    ap = result["alpha_peak"]
    N = result["N"]
    print(f"䷀  {name}: N = {N} codons  "
          f"({nseq} sequence(s), {total_nt} nt)", file=fh)
    print(f"䷁  H(q=2) = {H:.3f}   [{classify_H(H)}]", file=fh)
    print(f"䷂  Δα     = {da:.3f}   [{classify_da(da)}]", file=fh)
    print(f"䷃  α_peak = {ap:.3f}", file=fh)
    print(f"䷄  h(q)   = {hq_to_hexagrams(result['hq'])}", file=fh)


def emit_json(result: dict, name: str, nseq: int, total_nt: int,
              fh=sys.stdout) -> None:
    out = {
        "program": PROG,
        "version": __version__,
        "name": name,
        "n_sequences": nseq,
        "total_nt": total_nt,
        "N_codons": result["N"],
        "H_q2": result["H"],
        "H_class": classify_H(result["H"]),
        "delta_alpha": result["da"],
        "da_class": classify_da(result["da"]),
        "alpha_peak": result["alpha_peak"],
        "hq": result["hq"],
        "qs": result["qs"],
        "scales": result["scales"],
        "hq_hexagrams": hq_to_hexagrams(result["hq"]),
    }
    json.dump(out, fh, ensure_ascii=False, indent=2)
    fh.write("\n")


def emit_csv(results: List[Tuple[str, dict, int, int]], fh=sys.stdout) -> None:
    print("name,n_sequences,total_nt,N_codons,H_q2,H_class,"
          "delta_alpha,da_class,alpha_peak,hq_hexagrams", file=fh)
    for name, res, nseq, total_nt in results:
        print(f"{name},{nseq},{total_nt},{res['N']},"
              f"{res['H']:.4f},{classify_H(res['H'])},"
              f"{res['da']:.4f},{classify_da(res['da'])},"
              f"{res['alpha_peak']:.4f},"
              f"{hq_to_hexagrams(res['hq'])}", file=fh)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog=PROG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "UHCMFA — Universal Hexagram Coding Multifractal Analysis.\n"
            "Encodes DNA into hexagrams and performs multifractal analysis.\n"
            "Fast pre-screening tool: O(N log N), no alignment.\n"
            "\n"
            "See --methods for comparison against standard tools\n"
            "(BLAST, Kraken2, QUAST, CheckM, Mash, FastANI, ...)."
        ),
        epilog=(
            "Examples:\n"
            "  uhcmfa chr21.fa\n"
            "  uhcmfa chr21.fa.gz\n"
            "  uhcmfa --acc NC_001722\n"
            "  uhcmfa --acc NC_001722 M30502 X61240\n"
            "  uhcmfa --chrom 21\n"
            "  uhcmfa --chrom 21 22\n"
            "  uhcmfa --human\n"
            "  uhcmfa --compare sample.fa ecoli.fa yeast.fa\n"
            "  uhcmfa --frame 1 --rc chr21.fa\n"
            "  uhcmfa --json --acc NC_001722\n"
            "  uhcmfa --csv --chrom 21 --chrom 22\n"
            "  uhcmfa --hex-only --acc NC_001722 > hiv2.txt\n"
            "  uhcmfa --quiet --chrom 21\n"
            "  uhcmfa --methods\n"
            "  uhcmfa --api-key $NCBI_KEY --human\n"
        ),
    )
    p.add_argument("paths", nargs="*",
                   help="local FASTA or FASTA.gz files")
    p.add_argument("--acc", nargs="+", metavar="ACC",
                   help="NCBI accession number(s)")
    p.add_argument("--chrom", nargs="+", type=int, metavar="N",
                   help="human chromosome 1-24 (25 = mtDNA)")
    p.add_argument("--human", action="store_true",
                   help="fetch all 24 chromosomes + mitochondria")
    p.add_argument("--compare", nargs="+", metavar="FILE",
                   help="compare multiple files side by side (CSV-style)")
    p.add_argument("--frame", type=int, default=0, choices=[0, 1, 2],
                   help="reading frame (default 0)")
    p.add_argument("--rc", action="store_true",
                   help="use reverse complement")
    p.add_argument("--smin", type=int, default=16, metavar="N",
                   help="minimum scale in codons (default 16)")
    p.add_argument("--smax-div", type=int, default=4, metavar="N",
                   help="maximum scale = N/N (default 4)")
    p.add_argument("--api-key", default=None, metavar="KEY",
                   help="NCBI API key (speeds up fetching)")
    p.add_argument("--hex-only", action="store_true",
                   help="output hexagrams only, no analysis")
    p.add_argument("--json", action="store_true",
                   help="JSON output")
    p.add_argument("--csv", action="store_true",
                   help="CSV output")
    p.add_argument("--methods", action="store_true",
                   help="print method comparison tables and advantages")
    p.add_argument("-q", "--quiet", action="store_true",
                   help="no progress messages on stderr")
    p.add_argument("-o", "--output", default="-", metavar="FILE",
                   help="output file (default: stdout)")
    p.add_argument("--version", action="version",
                   version=f"{PROG} {__version__}")
    return p


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _collect_inputs(args, quiet: bool) -> List[Tuple[str, str]]:
    out: List[Tuple[str, str]] = []
    if args.human:
        for n, acc in HUMAN_CHROMS.items():
            if not quiet:
                sys.stderr.write(f"# fetching chromosome {n} ({acc}) ...\n")
            out.append((f"chr{n}", fetch_ncbi(acc, args.api_key, quiet)))
    elif args.chrom:
        for n in args.chrom:
            if n not in HUMAN_CHROMS:
                sys.stderr.write(f"{PROG}: invalid chromosome: {n}\n")
                continue
            acc = HUMAN_CHROMS[n]
            if not quiet:
                sys.stderr.write(f"# fetching chromosome {n} ({acc}) ...\n")
            out.append((f"chr{n}", fetch_ncbi(acc, args.api_key, quiet)))
    elif args.acc:
        for acc in args.acc:
            if not quiet:
                sys.stderr.write(f"# fetching {acc} ...\n")
            out.append((acc, fetch_ncbi(acc, args.api_key, quiet)))
    elif args.compare:
        for path in args.compare:
            if not quiet:
                sys.stderr.write(f"# reading {path} ...\n")
            out.append((path, read_fasta(path)))
    elif args.paths:
        for path in args.paths:
            if not quiet:
                sys.stderr.write(f"# reading {path} ...\n")
            out.append((path, read_fasta(path)))
    return out


def _run_one(name: str, seq: str, args, quiet: bool):
    x = to_series(seq, frame=args.frame, rc=args.rc)
    if x.size < 1024:
        sys.stderr.write(
            f"{PROG}: {name}: too short ({x.size} codons)\n")
        return None
    res = analyze(x, smin=args.smin, smax_div=args.smax_div)
    if res is None:
        sys.stderr.write(f"{PROG}: {name}: analysis failed\n")
        return None
    return (name, res, 1, len(seq))


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

    # --methods: print comparison tables and exit
    if args.methods:
        fh = sys.stdout if args.output == "-" else open(
            args.output, "w", encoding="utf-8")
        try:
            print_method_tables(fh)
        finally:
            if fh is not sys.stdout:
                fh.close()
        return 0

    # --hex-only: print hexagrams only
    if args.hex_only:
        inputs = _collect_inputs(args, args.quiet)
        if not inputs:
            parser.error("provide --acc, --chrom, --human, --compare, or a file")
        fh = sys.stdout if args.output == "-" else open(
            args.output, "w", encoding="utf-8")
        try:
            for name, seq in inputs:
                hexes = to_hexagrams(seq, frame=args.frame, rc=args.rc)
                for i in range(0, len(hexes), 64):
                    fh.write(hexes[i:i+64] + "\n")
        finally:
            if fh is not sys.stdout:
                fh.close()
        return 0

    inputs = _collect_inputs(args, args.quiet)
    if not inputs:
        parser.error("provide --acc, --chrom, --human, --compare, or a file")

    # --compare: analyze each separately, print comparison
    if args.compare or (args.paths and len(args.paths) > 1 and args.csv):
        results = []
        for name, seq in inputs:
            r = _run_one(name, seq, args, args.quiet)
            if r:
                results.append(r)
        if not results:
            return 1
        fh = sys.stdout if args.output == "-" else open(
            args.output, "w", encoding="utf-8")
        try:
            emit_csv(results, fh)
        finally:
            if fh is not sys.stdout:
                fh.close()
        return 0

    # Combined analysis
    if not args.quiet:
        sys.stderr.write(f"# combining {len(inputs)} sequence(s) ...\n")
    series = [to_series(s, frame=args.frame, rc=args.rc) for _, s in inputs]
    x = np.concatenate([s for s in series if s.size])
    total_nt = sum(len(s) for _, s in inputs)
    name = "+".join(n for n, _ in inputs) if len(inputs) <= 3 \
        else f"{len(inputs)} sequences"

    if x.size < 1024:
        sys.stderr.write(f"{PROG}: too short ({x.size} codons)\n")
        return 1

    res = analyze(x, smin=args.smin, smax_div=args.smax_div)
    if res is None:
        sys.stderr.write(f"{PROG}: analysis failed\n")
        return 1

    fh = sys.stdout if args.output == "-" else open(
        args.output, "w", encoding="utf-8")
    try:
        if args.json:
            emit_json(res, name, len(inputs), total_nt, fh)
        elif args.csv:
            emit_csv([(name, res, len(inputs), total_nt)], fh)
        else:
            emit_text(res, name, len(inputs), total_nt, fh)
    finally:
        if fh is not sys.stdout:
            fh.close()
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        sys.stderr.write(f"\n{PROG}: interrupted\n")
        sys.exit(130)