#!/usr/bin/env python
# Takes as input a single FASTA sequence of interest, the start and end points
# of the gene of interest in that sequence, and a multiFASTA of sequences to
# align to, returns a multiple alignment of the gene of interest aligned to the
# homologous gene in each other input sequence by performing multiple alignment
# on the protein sequences and then restoring the original nucleotide sequence
# while maintaining any gaps that may have been inserted. Allows for selection
# analysis of small genomes (viruses and maybe some bacteria).
# Can also be run in 'gene' mode where the reference input is a single
# in-frame gene and does not require [start] and [end] parameters.
# USAGE:
# kc-align --mode genome --reference [reference FASTA] --reads [reads FASTA] --start [start] --end [end]
# kc-align --mode gene --reference [reference FASTA] --reads [reads FASTA]
# kc-align --mode mixed --reference [reference FASTA] --reads [reads FASTA]

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import subprocess
import pdb


# Given a result from the aligning with Kalign, trims residues from longer 2nd
# sequence so it matches the 1st
def trim(align, trans):
    count = 0
    for i in align[0]:
        if i == 'M':
            break
        else:
            count += 1
    start = count
    count = 0
    for i in reversed(align[0]):
        if i == '*':
            break
        else:
            count += 1
    end = count
    while align[1][start] != 'M':
        start -= 1
    gapped = align[1][start:len(align[1])-end]
    no_gaps = ''
    for i in gapped:
        if i != '-':
            no_gaps += i
    return no_gaps


# Reinserts '*' at end of gapped alignment to make finding the end of the
# aligned portion more accurate
def reinsert_star(seq, gapped_seq):
    length = len(seq)
    count = 0
    for i in range(len(gapped_seq)):
        if gapped_seq[i] == seq[count] and count == 0:
            start = i
        if count == length-1:
            break
        elif gapped_seq[i] == seq[count]:
            count += 1
    return gapped_seq[:start]+gapped_seq[start:i+1]+'*'+gapped_seq[i+1:]


# Given a single gene's protein sequence and its entire DNA genome, finds the
# DNA seqeunce that corresponds with the given gene
def extract_DNA_seq(seq, trans):
    Pseq = seq.translate()
    for i in range(len(Pseq)):
        if str(Pseq[i:i+len(trans)]) == trans:
            break
    return seq[i*3:(i+len(trans))*3+3]


# Calculates Levenshtein distance between two strings
def distance(s1, s2):
    if len(s1) > len(s2):
        s1, s2 = s2, s1
    distances = range(len(s1)+1)
    for i2, c2 in enumerate(s2):
        distances_ = [i2+1]
        for i1, c1 in enumerate(s1):
            if c1 == c2:
                distances_.append(distances[i1])
            else:
                distances_.append(1+min((distances[i1], distances[i1+1],
                	              distances_[-1])))
        distances = distances_
    return distances[-1]


# Checks if the sliced out peptide sequence is complete by comparing to the
# original translated DNA sequence and extending the peptide sequence until
# it reaches where a stop codon would be in the original sequence
def extend_to_stop(prot_seq, frame, DNA_seq):
    trans = DNA_seq[frame:].translate()
    length = len(prot_seq)
    for i in range(len(trans)):
        if trans[i:i+length] == prot_seq:
            count = 0
            while trans[i+length+count] != '*':
                prot_seq += trans[i+length+count]
                count += 1
    return prot_seq


# Given a protein sequence and genome sequence finds the protein sequence of
# the homologous protein in the genome (currently checks the 3 coding strand
# reading frames)
# Calculates the Levenshtein distance of each candidate peptide to determine
# the correct reading frame
# If the Levenshtein distance divided by the length of the peptide of interest
# is less than 0.5, that seqeunce is thrown out (removes poorly aligning
# sequences)
def find_homologs(seq1, seq2):
    trans = []
    distances = []
    for i in range(0, 3):
        trans0 = seq2[i:].translate(stop_symbol='')
        stops = seq2[i:].translate()
        to_align = [SeqRecord(seq1, id='Seq1'), SeqRecord(trans0, id='Seq2')]
        SeqIO.write(to_align, 'tmp.fasta', 'fasta')
        subprocess.call(['kalign', '-i', 'tmp.fasta', '-o', 'out.fasta', '-d',
                        'pair'])
        alignments = []
        for record in SeqIO.parse('out.fasta', 'fasta'):
            alignments.append(record.seq)
        alignments[0] = reinsert_star(seq1, alignments[0])
        trans.append(Seq(trim(alignments, stops)))
        distances.append(distance(str(seq1), str(trans[-1])))
    subprocess.call(['rm', 'tmp.fasta', 'out.fasta'])
    if distances.index(min(distances))/len(seq1) <= 0.5:
        return (trans[distances.index(min(distances))],
                distances.index(min(distances)))
    else:
        return 1


# Use pairwise alignment with Kalign to determine the frame that the homolog
# for the gene of interest is located
def create_lists(reads, seq, og_seqs):
    seqs = []
    names = []
    ids = []
    err = []
    for record in SeqIO.parse(reads, 'fasta'):
        print(record.id)
        result = find_homologs(seq[:-1], record.seq)
        if result != 1:
            ids.append(record.id)
            names.append(record.description)
            prot_seq = extend_to_stop(result[0], result[1], record.seq)
            seqs.append(prot_seq)
            og_seqs[record.id] = extract_DNA_seq(record.seq[result[1]:],
                                                 seqs[-1])
        else:
            err.append(record.id)
    return seqs, names, ids, og_seqs, err


# Create multiFASTA file with protein sequences of the gene of interest and the
# in-frame sequences of the genomes to be multiple aligned and then align with
# Kalign
def combine_align(records, ids, names, seqs):
    for i, n, s in zip(ids, names, seqs):
        records.append(SeqRecord(s, id=i, description=n))
    SeqIO.write(records, 'pre_align.fasta', 'fasta')
    subprocess.call(['kalign', '-in', 'pre_align.fasta', '-out',
                     'aligned.fasta'])
    subprocess.call(['rm', 'pre_align.fasta'])


# Restore original codon sequence while maintaining gaps and write to file
def restore_codons(og_seqs):
    records = []
    for record in SeqIO.parse('aligned.fasta', 'fasta'):
        prot_seq = record.seq
        prot_seq += '*'
        nuc_seq = og_seqs[record.id]
        new_seq = ''
        counter = 0
        for i in range(len(prot_seq)):
            if prot_seq[i] != '-':
                new_seq += nuc_seq[counter:counter+3]
                counter += 3
            else:
                new_seq += '---'
        records.append(SeqRecord(new_seq, id=record.id,
                                 description=record.description))
    subprocess.call(['rm', 'aligned.fasta'])
    SeqIO.write(records, 'codon_aligned.fasta', 'fasta')
    SeqIO.write(records, 'codon_aligned.clustal', 'clustal')


# For when inputs are both whole genomes
def genome_mode(reference, reads, start, end):
    start = int(start)-1
    end = int(end)
    # Find protein sequence of gene of interest and extract the original DNA
    # sequence (only for genome mode)
    og_seqs = {}
    for record in SeqIO.parse(reference, 'fasta'):
        idd = record.id
        name = record.description
        seq = record.seq[start % 3:].translate()[start//3:end//3]
        og_seqs[record.id] = record.seq[start:end]
    seqs, names, ids, og_seqs, err = create_lists(reads, seq, og_seqs)
    records = [SeqRecord(seq[:-1], id=idd, description=name)]
    combine_align(records, ids, names, seqs)
    restore_codons(og_seqs)
    if len(err) > 0:
        print('The following sequences did not align sufficiently well to the\
               chosen gene and were thrown out before multiple alignment:\n')
        for e in err:
            print(e+'\n')


# For when inputs are in-frame genes
def gene_mode(reference, reads):
    ids = []
    names = []
    seqs = []
    og_seqs = {}
    for record in SeqIO.parse(reference, 'fasta'):
        ids.append(record.id)
        names.append(record.description)
        seqs.append(record.seq.translate())
        og_seqs[record.id] = record.seq
    for record in SeqIO.parse(reads, 'fasta'):
        ids.append(record.id)
        names.append(record.description)
        seqs.append(record.seq.translate())
        og_seqs[record.id] = record.seq
    records = []
    combine_align(records, ids, names, seqs)
    restore_codons(og_seqs)


# For when reference input is an in-frame gene but the reads are whole genomes
def mixed_mode(reference, reads):
    og_seqs = {}
    for record in SeqIO.parse(reference, 'fasta'):
        idd = record.id
        name = record.description
        seq = record.seq.translate()
        og_seqs[record.id] = record.seq
    seqs, names, ids, og_seqs, err = create_lists(reads, seq, og_seqs)
    records = [SeqRecord(seq[:-1], id=idd, description=name)]
    combine_align(records, ids, names, seqs)
    restore_codons(og_seqs)
    if len(err) > 0:
        print('The following sequences did not align sufficiently well to the\
               chosen gene and were thrown out before multiple alignment:\n')
        for e in err:
            print(e+'\n')
