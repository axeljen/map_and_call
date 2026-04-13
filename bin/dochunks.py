#!/usr/bin/env python3

import sys

def usage() -> None:
    print("Usage: dochunks.py <fai_file> <chunk_size_mb> <optional:scaffold_list>", file=sys.stderr)


# fasta index and chunk size in Mb as arguments
if len(sys.argv) < 3 or len(sys.argv) > 4:
    usage()
    sys.exit(1)

fai_file = sys.argv[1]
chunk_size_mb = int(sys.argv[2])
chunk_size_bp = chunk_size_mb * 1_000_000

scaffold_sizes = []
scaffold_list = None
if len(sys.argv) == 4:
    scaffold_list = [scaffold.strip() for scaffold in open(sys.argv[3], 'r').readlines()]

with open(fai_file, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        scaffold = parts[0]
        size = int(parts[1])
        if scaffold_list is None or scaffold in scaffold_list:
            scaffold_sizes.append((scaffold, size))

chunks = []

current_chunk = []
current_size = 0


def flush_chunk() -> None:
    if current_chunk:
        chunks.append(list(current_chunk))


for scaffold, size in scaffold_sizes:
    start = 1
    while start <= size:
        remaining = size - start + 1
        space_left = chunk_size_bp - current_size

        if space_left <= 0:
            flush_chunk()
            current_chunk = []
            current_size = 0
            space_left = chunk_size_bp

        take = remaining if remaining < space_left else space_left
        end = start + take - 1

        current_chunk.append(f"{scaffold}:{start}-{end}")
        current_size += take
        start = end + 1

        if current_size == chunk_size_bp:
            flush_chunk()
            current_chunk = []
            current_size = 0

if current_chunk:
    flush_chunk()

for chunk in chunks:
    print(",".join(chunk))

