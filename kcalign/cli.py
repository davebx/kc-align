#!/usr/bin/env python
import argparse
from kcalign import genome_mode, gene_mode, mixed_mode

def main():
    long_description = 'Takes as input a single FASTA sequence of interest, the start and end points of the gene of interest in that sequence, and a multiFASTA of sequences to align to, returns a multiple alignment of the gene of interest aligned to the homologous gene in each other input sequence by performing multiple alignment on the protein sequences and then restoring the original nucleotide sequence while maintaining any gaps that may have been inserted.'
    parser = argparse.ArgumentParser(description='Align a sequence against multiple reference.', epilog=long_description)
    parser.add_argument('--reference', '-r', dest='reference', action='store', required=True, help='Reference to align against')
    parser.add_argument('--reads', '-R', dest='reads', action='store', required=True, help='Reads to align')
    parser.add_argument('--start', '-s', dest='start', type=str, action='store', help='Start position, required in genome mode')
    parser.add_argument('--end', '-e', dest='end', type=str, action='store', help='End position, required in genome mode')
    parser.add_argument('--mode', '-m', dest='mode', action='store', choices=['genome', 'gene', 'mixed'], required=True, help='Alignment mode')
    args = parser.parse_args()

    if args.mode == 'genome':
        genome_mode(args.reference, args.reads, args.start, args.end)
    elif args.mode == 'gene':
        gene_mode(args.reference, args.reads)
    else:
        mixed_mode(args.reference, args.reads)

if __name__ == '__main__':
    exit(main())
