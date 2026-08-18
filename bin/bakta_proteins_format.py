#!/usr/bin/env python3
"""
Build a combined --proteins file for Bakta from:
  - a type-strain GenBank file (.gbk/.gb/.gbff) -- gene/product/translation pulled
    directly from each CDS feature's qualifiers (reliable, no header guessing), and
  - any number of extra protein FASTA files (your own curated/additional proteins).

Output header schema (Bakta short form):
    >id gene~~~product~~~dbxrefs

Requires biopython for GenBank input:  pip install biopython

Usage:
    python bakta_proteins_format.py type_strain.gbk additional_proteins.fasta -o combined.proteins.fasta

Input type is auto-detected (GenBank vs FASTA) by content, regardless of extension.
Duplicate IDs are dropped (first occurrence wins); duplicates are reported to stderr.
"""
import argparse
import re
import sys

NCBI_ORG_RE = re.compile(r"\s*\[[^\[\]]+\]\s*$")
UNIPROT_HEADER_RE = re.compile(r"^(sp|tr)\|([^|]+)\|(\S+)\s+(.*)$")
GN_TAG_RE = re.compile(r"\bGN=(\S+)")
OS_TAG_RE = re.compile(r"\s+OS=.*$")


def is_genbank(path):
    with open(path) as fh:
        for line in fh:
            if line.strip():
                return line.startswith("LOCUS")
    return False


def wrap_seq(seq, width=70):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def parse_fasta_header(raw_header):
    if "~~~" in raw_header:
        parts = raw_header.split(None, 1)
        seq_id = parts[0]
        rest = parts[1] if len(parts) > 1 else "~~~~~"
        return seq_id, rest

    m = UNIPROT_HEADER_RE.match(raw_header)
    if m:
        _, acc, _, desc = m.groups()
        gene_match = GN_TAG_RE.search(desc)
        gene = gene_match.group(1) if gene_match else ""
        product = OS_TAG_RE.sub("", desc).strip()
        return acc, f"{gene}~~~{product}~~~"

    parts = raw_header.split(None, 1)
    seq_id = parts[0]
    desc = parts[1].strip() if len(parts) > 1 else ""
    desc = NCBI_ORG_RE.sub("", desc).strip()
    if not desc:
        desc = "hypothetical protein"
    return seq_id, f"~~~{desc}~~~"


def entries_from_fasta(path):
    header = None
    seq_lines = []

    def build():
        seq_id, formatted_desc = parse_fasta_header(header)
        return seq_id, formatted_desc, "".join(seq_lines)

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield build()
                header = line[1:]
                seq_lines = []
            elif line:
                seq_lines.append(line)
        if header is not None:
            yield build()


def entries_from_genbank(path):
    from Bio import SeqIO

    for record in SeqIO.parse(path, "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue
            q = feature.qualifiers
            translation = q.get("translation", [None])[0]
            if not translation:
                print(f"[skip] CDS with no translation (likely pseudogene) in {path}: "
                      f"{q.get('locus_tag', ['?'])[0]}", file=sys.stderr)
                continue
            seq_id = q.get("protein_id", [None])[0] or q.get("locus_tag", [None])[0]
            if not seq_id:
                continue
            gene = q.get("gene", [""])[0]
            product = q.get("product", ["hypothetical protein"])[0]
            dbxrefs = ",".join(q.get("db_xref", []))
            yield seq_id, f"{gene}~~~{product}~~~{dbxrefs}", translation


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="Input file(s): GenBank (.gbk/.gb/.gbff) and/or FASTA")
    ap.add_argument("-o", "--output", required=True, help="Combined output FASTA for --proteins")
    args = ap.parse_args()

    seen_ids = set()
    n_written = 0
    n_dupes = 0

    with open(args.output, "w") as out:
        for path in args.inputs:
            gen = entries_from_genbank(path) if is_genbank(path) else entries_from_fasta(path)
            for seq_id, formatted_desc, seq in gen:
                if seq_id in seen_ids:
                    n_dupes += 1
                    print(f"[skip] duplicate id '{seq_id}' from {path}", file=sys.stderr)
                    continue
                seen_ids.add(seq_id)
                out.write(f">{seq_id} {formatted_desc}\n{wrap_seq(seq)}\n")
                n_written += 1

    print(f"Wrote {n_written} sequences to {args.output} ({n_dupes} duplicate ids skipped)", file=sys.stderr)


if __name__ == "__main__":
    main()
