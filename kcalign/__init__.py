#!/usr/bin/env python

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import subprocess
import warnings
warnings.filterwarnings('ignore')


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
# the correct reading frame (lowest distance is correct)
# To detect frameshifts, the Levenshtein distances of each reading frame with
# the reference are compard. If the distance of the frame with the lowest 
# distance divided by either of the other reading frame's distances is 
# between 0.5 and 1.5 then a frameshift is likely. If there is no frameshift
# the distance between the correct reading frame and the reference should
# be considerably lower than the distance between the other reading frames
# and the reference.
def find_homologs(seq1, seq2):
    trans = []
    distances = []
    for i in range(0, 3):
        trans0 = seq2[i:].translate(stop_symbol='')
        stops = seq2[i:].translate()
        to_align = [SeqRecord(seq1, id='Seq1'), SeqRecord(trans0, id='Seq2')]
        SeqIO.write(to_align, 'tmp.fasta', 'fasta')
        subprocess.call(['kalign', '-i', 'tmp.fasta', '-o', 'out.fasta', '-d',
                        'pair'], stderr=subprocess.DEVNULL)
        alignments = []
        for record in SeqIO.parse('out.fasta', 'fasta'):
            alignments.append(record.seq)
        alignments[0] = reinsert_star(seq1, alignments[0])
        trans.append(Seq(trim(alignments, stops)))
        distances.append(distance(str(seq1), str(trans[-1])))
    subprocess.call(['rm', 'tmp.fasta', 'out.fasta'])
    inds = [0, 1, 2]
    del inds[distances.index(min(distances))]
    for i in inds:
        if min(distances) == 0:
            break
        if (distances[i]/min(distances) < 1.5 
             and distances[i]/min(distances) > 0.5):
            return 1
    return (trans[distances.index(min(distances))],
            distances.index(min(distances)))


# Determines if the homologous sequence found by find_homologs contains a
# frameshift mutation by first calculating a sliding window of edit distances
# in order to find the approximate location of the indel (done by finding the
# region with the largest increase in edit distance between windows). An 
# extra nucleotide is then inserted in that region, the sequence is converted
# to amino acids, realigned with the reference and the edit distance between 
# the alignments calculated. This is repeated one more time and if the edit 
# distance of the alignments after inserting on or two extra nucleotides is
# less than the distance without adding nucleotides then there is likely to
# be a frameshift mutation and the function will return a value of 1. If no 
# frameshift is detected, it returns 0. 
def detect_frameshift(ref, read):
    records = [SeqRecord(ref, id='reference'),
               SeqRecord(read.translate(), id='read')]
    SeqIO.write(records, 'tmp.fasta', 'fasta')
    subprocess.call(['kalign', '-in', 'tmp.fasta', '-out', 'out.fasta', '-d',
                     'pair'], stderr=subprocess.DEVNULL)
    seqs = {}
    for record in SeqIO.parse('out.fasta', 'fasta'):
        seqs[record.id] = str(record.seq)
    distances = []
    for i in range(0, len(seqs['reference'])-50, 20):
        distances.append(distance(seqs['reference'][i:i+50],
                                  seqs['read'][i:i+50]))
    diffs = []
    for i in range(1, len(distances)):
        diffs.append(distances[i]-distances[i-1])
    ind = diffs.index(max(diffs))*60
    new_distances = [distance(seqs['reference'], seqs['read'])]
    records = [SeqRecord(ref, id='reference'),
               SeqRecord((read[:ind]+'A'+read[ind:]).translate())]
    SeqIO.write(records, 'tmp.fasta', 'fasta')
    subprocess.call(['kalign', '-in', 'tmp.fasta', '-out', 'out.fasta', '-d',
                     'pair'], stderr=subprocess.DEVNULL)
    for record in SeqIO.parse('out.fasta', 'fasta'):
        seqs[record.id] = str(record.seq)
    new_distances.append(distance(seqs['reference'], seqs['read']))
    records = [SeqRecord(ref, id='reference'),
               SeqRecord((read[:ind]+'AA'+read[ind:]).translate())]
    SeqIO.write(records, 'tmp.fasta', 'fasta')
    subprocess.call(['kalign', '-in', 'tmp.fasta', '-out', 'out.fasta', '-d',
                     'pair'], stderr=subprocess.DEVNULL)
    for record in SeqIO.parse('out.fasta', 'fasta'):
        seqs[record.id] = str(record.seq)
    subprocess.call(['rm', 'out.fasta', 'tmp.fasta'])
    new_distances.append(distance(seqs['reference'], seqs['read']))
    if (new_distances[1] < new_distances[0] or 
        new_distances[2] < new_distances[0]):
        return 1
    else:
        return 0


# Use pairwise alignment with Kalign to determine the frame that the homolog
# for the gene of interest is located. Will also look for any potential
# frameshift mutations in the gene of interest and throw out the sequence
# if one is detected.
def create_lists(reads, seq, og_seqs):
    seqs = []
    names = []
    ids = []
    err = []
    for record in SeqIO.parse(reads, 'fasta'):
        print(record.id)
        result = find_homologs(seq[:-1], record.seq)
        if result == 1:
            err.append(record.id)
        else:
            ids.append(record.id)
            names.append(record.description)
            prot_seq = extend_to_stop(result[0], result[1], record.seq)
            seqs.append(prot_seq)
            og_seqs[record.id] = extract_DNA_seq(record.seq[result[1]:], 
                                                 seqs[-1])
            shift = detect_frameshift(seq, og_seqs[record.id])
            if shift == 1:
                seqs = seqs[:-1]
                names = names[:-1]
                del og_seqs[record.id]
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
    SeqIO.write(records, 'kc-align.fasta', 'fasta')
    SeqIO.write(records, 'kc-align.clustal', 'clustal')


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
        print('The following sequences were suspected of containing ' \
              'frameshifts and so were thrown out before multiple alignment:')
        for e in err:
            print(e+'\n')


# For when inputs are in-frame genes
def gene_mode(reference, reads):
    ids = []
    names = []
    seqs = []
    og_seqs = {}
    err = []
    for record in SeqIO.parse(reference, 'fasta'):
        ids.append(record.id)
        names.append(record.description)
        seqs.append(record.seq.translate())
        og_seqs[record.id] = record.seq
    for record in SeqIO.parse(reads, 'fasta'):
        if '*' in record.seq.translate()[:-1]:
            err.append(record.id)
        else:
            ids.append(record.id)
            names.append(record.description)
            seqs.append(record.seq.translate())
            og_seqs[record.id] = record.seq
    records = []
    combine_align(records, ids, names, seqs)
    restore_codons(og_seqs)
    if len(err) > 0:
        print('The following sequences were suspected of containing ' \
              'frameshifts and so were thrown out before multiple alignment:')
        for e in err:
            print(e+'\n')


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
        print('The following sequences were suspected of containing ' \
              'frameshifts and so were thrown out before multiple alignment:')
        for e in err:
            print(e+'\n')
