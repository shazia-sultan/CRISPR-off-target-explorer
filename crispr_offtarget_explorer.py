# crispr_offtarget_explorer.py
# Paste-only CRISPR (SpCas9) off-target explorer.
# What this does:
# - User pastes a genome sequence (A/C/G/T only, no FASTA headers etc.)
# - User provides a 20-nt CRISPR guide sequence.
# - Program scans the genome for possible target sites:
#       * On the + strand: looks for [20nt protospacer] + [NGG PAM].
#       * On the - strand: looks for [CCN PAM upstream] + [20nt protospacer].
# - Reports perfect matches (on-targets, 0 mismatches).
# - Reports off-targets up to a mismatch threshold (default ≤3).
# - Provides a "minesweeper"-style grid where each cell shows the number
#   of off-targets in that chunk of the genome.
# - Interactive inspector lets the user explore hits.

from dataclasses import dataclass
from typing import List, Tuple
import math
import re
import sys


PAM_PATTERN = "NGG"   # SpCas9 PAM (forward strand)
MM_THRESHOLD = 3      # report off-targets with ≤3 mismatches
CHUNK_BP    = 200     # genome is divided into chunks of 200 bp for grid display
DEBUG_SCAN  = False   # set True to print candidate windows during scan


@dataclass
class Hit:
    """
    A single CRISPR site hit (either on-target or off-target).
    """
    pos0: int            # 0-based index of protospacer start on + reference
    strand: str          # '+' (forward) or '-' (reverse strand)
    mm: int              # number of mismatches compared to guide
    seed_mm: int         # mismatches in the "seed" region (first 12 bases)
    prot: str            # 20-nt protospacer sequence (in guide orientation)
    pam_triplet: str     # PAM triplet adjacent to protospacer (e.g., TGG/CGG)
    pam_obs_plus: str    # observed PAM triplet on the + reference (useful for - strand hits)



def revcomp(s: str) -> str:
    """Return reverse complement of DNA sequence."""
    return s.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]

def pam_match(pattern: str, triplet: str) -> bool:
    """
    Check if a PAM matches the pattern (e.g., 'NGG').
    'N' is treated as a wildcard.
    """
    return len(pattern) == len(triplet) and all(p == "N" or p == b for p, b in zip(pattern, triplet))

def hamming(a: str, b: str) -> int:
    """Compute Hamming distance (number of mismatches) between two equal-length strings."""
    return sum(x != y for x, y in zip(a, b))

def seed_mismatches(prot: str, guide: str, seed_len: int = 12) -> int:
    """Count mismatches in the seed region (positions 1–12 of the protospacer)."""
    return hamming(prot[:seed_len], guide[:seed_len])


def read_guide() -> str:
    """Ask user for a 20-nt guide RNA sequence."""
    g = input("Enter 20-nt guide (A/C/G/T): ").strip().upper()
    if len(g) != 20 or re.search(r"[^ACGT]", g):
        print("ValueError: guide must be exactly 20 nt of A/C/G/T.")
        sys.exit(1)
    return g

def read_genome_from_paste() -> str:
    """Ask user to paste genome sequence, stop on blank line."""
    print("Paste genome sequence (A/C/G/T). Press ENTER on a blank line to finish:")
    lines: List[str] = []
    while True:
        try:
            line = input()
        except EOFError:
            break
        if line.strip() == "":
            if lines:
                break
            continue
        lines.append(line)

    raw = "".join(lines)
    seq = re.sub(r"[^ACGTacgt]", "", raw).upper()

    if not seq:
        print("Error: no A/C/G/T bases detected.")
        sys.exit(1)

    compact_raw = re.sub(r"\s+", "", raw)
    if len(seq) != len(compact_raw):
        print("[info] Removed non-ACGT characters (whitespace/others).")
    return seq



def scan(genome: str, guide: str, pam: str = PAM_PATTERN, mm_thr: int = MM_THRESHOLD) -> Tuple[List[Hit], List[Hit]]:
    """
    Scan both strands for potential CRISPR targets.
    Returns (on_targets, off_targets).
    """
    g, L, k, p = genome, len(genome), len(guide), len(pam)
    on:  List[Hit] = []
    off: List[Hit] = []

    # -------- Forward strand (protospacer followed by NGG PAM) --------
    for i in range(0, L - k - p + 1):
        prot = g[i:i+k]                  # 20-nt protospacer
        pam_plus = g[i+k:i+k+p]          # next 3 bases (should be NGG)
        if DEBUG_SCAN and (guide in prot or hamming(prot, guide) <= mm_thr):
            print(f"[dbg +] i={i} prot={prot} pam={pam_plus} mm={hamming(prot, guide)}")

        if pam_match(pam, pam_plus):
            mm = hamming(prot, guide)
            hit = Hit(
                pos0=i, strand="+", mm=mm, seed_mm=seed_mismatches(prot, guide),
                prot=prot, pam_triplet=pam_plus, pam_obs_plus=pam_plus
            )
            # Add hit to on-targets (0 mismatches) or off-targets (≤mm_thr mismatches).
            (on if mm == 0 else off).append(hit) if mm <= mm_thr else None

    # -------- Reverse strand (CCN PAM upstream on +) --------
    rev_pam = revcomp(pam)  # NGG → CCN
    for end in range(k + p, L + 1):
        prot_plus = g[end-k:end]       # 20 bases on + (but actually protospacer on - strand)
        pam_up    = g[end-k-p:end-k]   # 3 bases upstream on + (should be CCN)
        prot_rc   = revcomp(prot_plus) # protospacer in guide orientation

        if DEBUG_SCAN and (guide in prot_rc or hamming(prot_rc, guide) <= mm_thr):
            print(f"[dbg -] end={end} prot_plus={prot_plus} pam_up={pam_up} mm={hamming(prot_rc, guide)}")

        if pam_match(rev_pam, pam_up):
            mm = hamming(prot_rc, guide)
            hit = Hit(
                pos0=end-k, strand="-", mm=mm, seed_mm=seed_mismatches(prot_rc, guide),
                prot=prot_rc, pam_triplet=revcomp(pam_up), pam_obs_plus=pam_up
            )
            (on if mm == 0 else off).append(hit) if mm <= mm_thr else None

    return on, off



def build_grid(genome_len: int, off_targets: List[Hit], chunk_bp: int = CHUNK_BP):
    """
    Divide genome into chunks and count off-targets per chunk.
    Returns (rows, cols, cells, chunks).
    """
    chunks = math.ceil(genome_len / chunk_bp)
    cols   = math.ceil(math.sqrt(chunks))
    rows   = math.ceil(chunks / cols)
    cells  = [[{"count": 0, "idxs": []} for _ in range(cols)] for _ in range(rows)]

    def chunk_id(pos0: int) -> int:
        return min(pos0 // chunk_bp, chunks - 1)

    for idx, site in enumerate(off_targets):
        ci = chunk_id(site.pos0)
        r, c = divmod(ci, cols)
        cells[r][c]["count"] += 1
        cells[r][c]["idxs"].append(idx)

    return rows, cols, cells, chunks

def print_grid(rows: int, cols: int, cells, width: int = 3) -> None:
    """Print grid header and counts (· means zero)."""
    print("   " + " ".join(f"{c:>{width}d}" for c in range(cols)))
    for r in range(rows):
        row = [f"{r:>2d} "]
        for c in range(cols):
            cnt = cells[r][c]["count"]
            row.append(f"{(str(cnt) if cnt else '·'):>{width}s}")
        print(" ".join(row))


def main() -> None:
    print("=== CRISPR Off-Target Explorer  ===")
    guide  = read_guide()
    genome = read_genome_from_paste()

    print("\nScanning…")
    on, off = scan(genome, guide, pam=PAM_PATTERN, mm_thr=MM_THRESHOLD)

   
    print("\n--- Summary ---")
    print(f"Genome length: {len(genome)} bp")
    print(f"Guide: {guide} | PAM: {PAM_PATTERN} | mismatch ≤ {MM_THRESHOLD}")
    print(f"On-target (0 mm): {len(on)}")
    print(f"Off-targets (1..{MM_THRESHOLD} mm): {len(off)}")

    rows, cols, cells, _ = build_grid(len(genome), off, CHUNK_BP)
    print(f"Grid: {rows} × {cols} cells | {CHUNK_BP} bp per cell\n")
    print("Each cell shows # of OFF-TARGETS in that genome chunk (· means zero).")
    print_grid(rows, cols, cells)

   
    print("\nCommands:")
    print("  i <row> <col>  : list off-targets in a cell")
    print("  on             : list on-targets (0 mm)")
    print("  q              : quit")

    while True:
        parts = input("> ").strip().split()
        if not parts:
            continue
        cmd = parts[0].lower()

        if cmd == "q":
            break

        elif cmd == "on":
            # Print all perfect matches
            if not on:
                print("No on-targets found.")
                continue
            for s in sorted(on, key=lambda d: (d.pos0, d.strand))[:200]:
                print(f" - pos={s.pos0+1} strand={s.strand} mm={s.mm} seed_mm={s.seed_mm} "
                      f"prot={s.prot} PAM={s.pam_triplet}"
                      + (f" (obs:+ {s.pam_obs_plus})" if s.strand == "-" else ""))

        elif cmd == "i" and len(parts) == 3:
            # Inspect a single grid cell
            try:
                r, c = int(parts[1]), int(parts[2])
            except ValueError:
                print("Use integers: i <row> <col>")
                continue
            if not (0 <= r < rows and 0 <= c < cols):
                print("Out of bounds.")
                continue

            idxs = cells[r][c]["idxs"]
            if not idxs:
                print("No off-targets in this chunk.")
                continue

            chunk_id = r * cols + c
            start = chunk_id * CHUNK_BP
            end   = min((chunk_id + 1) * CHUNK_BP, len(genome))
            print(f"Chunk [{r},{c}] spans {start+1}-{end} (len={end-start})")

            for s in sorted((off[i] for i in idxs), key=lambda d: (d.mm, d.pos0))[:50]:
                print(f" - pos={s.pos0+1} strand={s.strand} mm={s.mm} seed_mm={s.seed_mm} "
                      f"prot={s.prot} PAM={s.pam_triplet}"
                      + (f" (obs:+ {s.pam_obs_plus})" if s.strand == "-" else ""))

        else:
            print("Commands: i <row> <col>  |  on  |  q")

if __name__ == "__main__":
    main()
