#!/usr/bin/env python3
"""Rewrite reference-list compounds to the library's modification notation.

Context matches a reference compound to a library entry as a string after
stripping flanking residues, so a Skyline-style ``C[+57.0]`` never equals the
library's ``C[+57.0214635]`` and such peptides can never be reference hits.
Each compound is matched to a library entry by its stripped sequence and its
modification masses rounded to 0.1 Da at the same positions, and replaced by
that entry's PeptideModSeq. Compounds without a match are kept unchanged. The
file layout (delimiter, columns, column order) is preserved, since EncyclopeDIA
reads mass lists by position.

    normalize_reference_mods.py --reference in.tsv --library lib.dlib --out out.tsv
"""

import argparse
import csv
import re
import sqlite3
import sys

MOD = re.compile(r"\[([+-]?\d+(?:\.\d+)?)\]")


def key(modseq: str):
    """(stripped sequence, ((position, mass rounded to 0.1), ...)); position 0
    holds an n-terminal modification, residue i holds position i (1-based)."""
    seq, mods, pos, i = [], [], 0, 0
    while i < len(modseq):
        c = modseq[i]
        if c == "[":
            m = MOD.match(modseq, i)
            if not m:
                raise ValueError(f"bad modification in {modseq!r}")
            mods.append((pos, round(float(m.group(1)), 1)))
            i = m.end()
        elif c.isalpha():
            seq.append(c)
            pos += 1
            i += 1
        else:
            i += 1
    return "".join(seq), tuple(sorted(mods))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--library", required=True, help=".dlib or .elib")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    con = sqlite3.connect(f"file:{a.library}?mode=ro", uri=True)
    lib = {}
    for (modseq,) in con.execute("SELECT DISTINCT PeptideModSeq FROM entries"):
        lib.setdefault(key(modseq), modseq)
    con.close()

    lines = open(a.reference, newline="").read().splitlines()
    if not lines:
        sys.exit(f"normalize_reference_mods: {a.reference} is empty")
    delim = "\t" if lines[0].count("\t") >= lines[0].count(",") else ","
    rows = list(csv.reader(lines, delimiter=delim))
    changed = unmatched = 0
    for row in rows[1:]:
        if not row or not row[0].strip():
            continue
        new = lib.get(key(row[0].strip()))
        if new is None:
            unmatched += 1
        elif new != row[0]:
            row[0] = new
            changed += 1
    with open(a.out, "w", newline="") as f:
        csv.writer(f, delimiter=delim, lineterminator="\n").writerows(rows)
    print(f"normalize_reference_mods: {len(rows) - 1} rows, {changed} rewritten to the "
          f"library notation, {unmatched} without a library entry (kept as is)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
