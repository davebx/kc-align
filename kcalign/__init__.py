#!/usr/bin/env python

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import pty
import shlex
import sys
import subprocess
import warnings
warnings.filterwarnings('ignore')
import concurrent.futures


# Given a result from the aligning with Kalign, trims residues from
# longer 2nd sequence so it matches the 1st
def trim(align, trans):
    count = 0
    for i in align[0]:
        if i != '-':
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
    gapped = align[1][start:len(align[1])-end]
    no_gaps = ''
    for i in gapped:
        if i != '-':
            no_gaps += i
    return no_gaps


# Kalign3 might take a few seconds on larger datasets. Throw in a wait()
# just to be safe.
def invoke_kalign(input_file, output_file):
    if not os.path.exists(input_file):
        print('Input file missing')
        exit(1)
    # If STDIN is not a tty, kalign3 assumes that STDIN is an input file
    # and fails to detect its type.
    # Pass in a fake pty to work around this behavior, which is appropriate
    # for command-line usage but doesn't work in a headless environment.
    fakepty, alsofakepty = pty.openpty()
    command = shlex.split(f'kalign -i {input_file} -o {output_file}')
    kaligner = subprocess.Popen(command, stderr=subprocess.PIPE,
                                stdout=subprocess.PIPE, stdin=fakepty)
    kaligner.wait()
    stdout, stderr = kaligner.communicate()
    os.close(fakepty)
    os.close(alsofakepty)
    # If Kalign fails try again with MAFFT (Kalign sometimes seg faults on some data)
    if kaligner.returncode != 0:
        with open(output_file, 'w') as out_fh:
            command = shlex.split('mafft %s' % input_file)
            mafft = subprocess.Popen(command, stdout=out_fh)
            mafft.wait()
            _, stderr = mafft.communicate()
        if mafft.returncode != 0:
            sys.stderr.write('Error running mafft')
            sys.stderr.write(stderr.decode('utf-8'))
            return 1
        else:
            return 0
    else:
        return 0


# Reinserts '*' at end of gapped alignment to make finding the end of
# the aligned portion more accurate
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


# Given a single gene's protein sequence and its entire DNA genome,
# finds the DNA seqeunce that corresponds with the given gene.
def extract_DNA_seq(seq, trans, tab):
    Pseq = seq.translate(table=tab)
    check = 0
    for i in range(len(Pseq)):
        if str(Pseq[i:i+len(trans)]) == trans:
            check = 1
            break
    if check == 1:
        return seq[i*3:(i+len(trans))*3]
    else:
        return 1


# Same as extract_DNA_seq() but for when the gene of interests is split
# between two different reading frames.
def join_extract_DNA_seq(seq, homologs, tab):
    Dseqs = []
    for homolog in homologs:
        shift = homolog[1]
        Pseq = seq[homolog[1]:].translate(table=tab)
        check = 0
        for i in range(len(Pseq)):
            if str(Pseq[i:i+len(homolog[0])]) == homolog[0]:
                check = 1
                break
        if check == 1:
            Dseqs.append(seq[i*3+shift:(i+len(homolog[0]))*3+shift])
        else:
            return 1
    return Dseqs[0]+Dseqs[1]


# Calculates Levenshtein distance between two strings.
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


# Given a protein sequence and genome sequence finds the protein
# sequence of the homologous protein in the genome (currently checks
# the 3 coding strand reading frames).
# Calculates the Levenshtein distance of each candidate peptide to
# determine the correct reading frame (lowest distance is correct).
# To detect frameshifts, the Levenshtein distances of each reading frame
# with the reference are compard. If the distance of the frame with the
# lowest distance divided by either of the other reading frame's
# distances is between 0.9 and 1.1 then a frameshift is likely. If there
# is no frameshift the distance between the correct reading frame and
# the reference should be considerably lower than the distance between
# the other reading frames and the reference.
def para_find_homologs(seq1, seq2, tab):
    trans = []
    distances = []
    frames2 = []
    with concurrent.futures.ProcessPoolExecutor() as executor:
        frames = [0, 1, 2]
        results = [executor.submit(test_frames, seq1,seq2,f,tab) for f in frames]
        for res in concurrent.futures.as_completed(results):
            if res.result() == 1:
                return 1
            else:
                trans.append(res.result()[0])
                distances.append(res.result()[1])
                frames2.append(res.result()[2])
    subprocess.call('rm tmpfilexyz*.fasta outfilexyz*.fasta', shell=True)
    del frames[distances.index(min(distances))]
    for i in frames:
        if min(distances) == 0:
            break
        if (distances[i]/min(distances) < 1.1 and
           distances[i]/min(distances) > 0.9):
            return 1
    return (trans[distances.index(min(distances))],
            frames2[distances.index(min(distances))])


# Same as find_homologs() but for when the gene of interest is split
# between two different reading frames.
def para_join_find_homologs(seqs, seq2, tab):
    final = []
    for seq in seqs:
        trans = []
        distances = []
        frames2 = []
        with concurrent.futures.ProcessPoolExecutor() as executor:
            frames = [0, 1, 2]
            results = [executor.submit(test_frames, seq,seq2,f,tab) for f in frames]
            for res in concurrent.futures.as_completed(results):
                if res.result() == 1:
                    return 1
                else:
                    trans.append(res.result()[0])
                    distances.append(res.result()[1])
                    frames2.append(res.result()[2])
        subprocess.call('rm tmpfilexyz*.fasta outfilexyz*.fasta', shell=True)
        del frames[distances.index(min(distances))]
        for i in frames:
            if min(distances) == 0:
                break
            if (distances[i]/min(distances) < 1.1 and
               distances[i]/min(distances) > 0.9):
                return 1
        final.append((trans[distances.index(min(distances))],
                      frames2[distances.index(min(distances))]))
    return final


# Function that finds the correct reading frame for homolog. Originally
# a part of find_homologs but was split to enable parallelization.
def test_frames(seq1, seq2, frame, tab):
    trans0 = seq2[frame:].translate(stop_symbol='', table=tab)
    stops = seq2[frame:].translate(table=tab)
    to_align = [SeqRecord(seq1, id='Seq1'), SeqRecord(trans0, id='Seq2')]
    SeqIO.write(to_align, 'tmpfilexyz'+str(frame)+'.fasta', 'fasta')
    kresult = invoke_kalign('tmpfilexyz'+str(frame)+'.fasta', 'outfilexyz'+str(frame)+'.fasta')
    if kresult == 1:
        return 1
    alignments = []
    for record in SeqIO.parse('outfilexyz'+str(frame)+'.fasta', 'fasta'):
        alignments.append(record.seq)
    alignments[0] = reinsert_star(seq1, alignments[0])
    trans = Seq(trim(alignments, stops))
    dist = distance(str(seq1), str(trans))
    return (trans, dist, frame)


# Original non parallel versions
def find_homologs(seq1, seq2, tab):
    trans = []
    distances = []
    for i in range(0, 3):
        trans0 = seq2[i:].translate(stop_symbol='', table=tab)
        stops = seq2[i:].translate(table=tab)
        to_align = [SeqRecord(seq1, id='Seq1'), SeqRecord(trans0, id='Seq2')]
        SeqIO.write(to_align, 'tmpfilexyz.fasta', 'fasta')
        kresult = invoke_kalign('tmpfilexyz.fasta', 'outfilexyz.fasta')
        if kresult == 1:
            return 1
        alignments = []
        for record in SeqIO.parse('outfilexyz.fasta', 'fasta'):
            alignments.append(record.seq)
        alignments[0] = reinsert_star(seq1, alignments[0])
        trans.append(Seq(trim(alignments, stops)))
        distances.append(distance(str(seq1), str(trans[-1])))
    subprocess.call(['rm', 'tmpfilexyz.fasta', 'outfilexyz.fasta'])
    inds = [0, 1, 2]
    del inds[distances.index(min(distances))]
    for i in inds:
        if min(distances) == 0:
            break
        if (distances[i]/min(distances) < 1.1 and
           distances[i]/min(distances) > 0.9):
            return 1
    return (trans[distances.index(min(distances))],
            distances.index(min(distances)))


def join_find_homologs(seqs, seq2, tab):
    final = []
    for seq in seqs:
        trans = []
        distances = []
        for i in range(0, 3):
            trans0 = seq2[i:].translate(stop_symbol='', table=tab)
            stops = seq2[i:].translate(table=tab)
            to_align = [SeqRecord(seq, id='Seq1'),
                        SeqRecord(trans0, id='Seq2')]
            SeqIO.write(to_align, 'tmpfilexyz.fasta', 'fasta')
            kresult = invoke_kalign('tmpfilexyz.fasta', 'outfilexyz.fasta')
            if kresult == 1:
                return 1
            alignments = []
            for record in SeqIO.parse('outfilexyz.fasta', 'fasta'):
                alignments.append(record.seq)
            alignments[0] = reinsert_star(seq, alignments[0])
            trans.append(Seq(trim(alignments, stops)))
            distances.append(distance(str(seq), str(trans[-1])))
        subprocess.call(['rm', 'tmpfilexyz.fasta', 'outfilexyz.fasta'])
        inds = [0, 1, 2]
        del inds[distances.index(min(distances))]
        for i in inds:
            if min(distances) == 0:
                break
            if (distances[i]/min(distances) < 1.1 and
               distances[i]/min(distances) > 0.9):
                return 1
        final.append((trans[distances.index(min(distances))],
                      distances.index(min(distances))))
    return final


# Determines if the homologous sequence found by find_homologs contains
# a frameshift mutation by first calculating a sliding window of edit
# distances in order to find the approximate location of the indel (done
# by finding the region with the largest increase in edit distance
# between windows). An extra nucleotide is then inserted in that region,
# the sequence is converted to amino acids, realigned with the reference
# and the edit distance between the alignments calculated. This is
# repeated one more time and if the edit distance of the alignments
# after inserting on or two extra nucleotides is less than the distance
# without adding nucleotides then there is likely to be a frameshift
# mutation and the function will return a value of 1. If no frameshift
# is detected, it returns 0. Not effective for frameshifts that occur
# at the end of a sequence.
def detect_frameshift(ref, read):
    records = [SeqRecord(ref, id='reference'),
               SeqRecord(read.translate(table=tab), id='read')]
    SeqIO.write(records, 'tmpfilexyz.fasta', 'fasta')
    invoke_kalign('tmpfilexyz.fasta', 'outfilexyz.fasta')
    seqs = {}
    for record in SeqIO.parse('outfilexyz.fasta', 'fasta'):
        seqs[record.id] = str(record.seq)
    distances = []
    for i in range(0, len(seqs['reference'])-30, 10):
        distances.append(distance(seqs['reference'][i:i+30],
                                  seqs['read'][i:i+30]))
    diffs = []
    for i in range(1, len(distances)):
        diffs.append(distances[i]-distances[i-1])
    ind = diffs.index(max(diffs))*30
    new_distances = [distance(seqs['reference'], seqs['read'])]
    records = [SeqRecord(ref, id='reference'),
               SeqRecord((read[:ind]+'A'+read[ind:]).translate(table=tab))]
    SeqIO.write(records, 'tmpfilexyz.fasta', 'fasta')
    invoke_kalign('tmpfilexyz.fasta', 'outfilexyz.fasta')
    for record in SeqIO.parse('outfilexyz.fasta', 'fasta'):
        seqs[record.id] = str(record.seq)
    new_distances.append(distance(seqs['reference'], seqs['read']))
    records = [SeqRecord(ref, id='reference'),
               SeqRecord((read[:ind]+'AA'+read[ind:]).translate(table=tab))]
    SeqIO.write(records, 'tmpfilexyz.fasta', 'fasta')
    invoke_kalign('tmpfilexyz.fasta', 'outfilexyz.fasta')
    for record in SeqIO.parse('outfilexyz.fasta', 'fasta'):
        seqs[record.id] = str(record.seq)
    subprocess.call(['rm', 'outfilexyz.fasta', 'tmpfilexyz.fasta'])
    new_distances.append(distance(seqs['reference'], seqs['read']))
    if (new_distances[1] <= new_distances[0] or
       new_distances[2] <= new_distances[0]):
        return 1
    else:
        return 0


# Calculates N cotent of sequences. Returns 0 if N content is less than 5%
# and 1 if it is above 5%
def check_n(seq):
    n = 0
    for i in seq:
        if i == 'N':
            n += 1
    if n/len(seq) < 0.05:
        return 0
    else:
        return 1


# Use pairwise alignment with Kalign to determine the frame that the
# homolog for the gene of interest is located. Will also look for any
# potential frameshift mutations in the gene of interest and throw out
# the sequence if one is detected.
def create_lists(reads, seq, og_seqs, join, para, tab):
    seqs = []
    names = []
    ids = []
    err = []
    if join == 0 and seq[-1] == '*':
        seq = seq[:-1]
    for record in SeqIO.parse(reads, 'fasta'):
        if para:
            if join == 0:
                result = para_find_homologs(seq, record.seq, tab)
            elif join == 1:
                result = para_join_find_homologs(seq, record.seq, tab)
        else:
            if join == 0:
                result = find_homologs(seq, record.seq, tab)
            elif join == 1:
                result = join_find_homologs(seq, record.seq, tab)
        if result == 1:
            err.append(record.id)
        else:
            ids.append(record.id)
            names.append(record.description)
            if join == 0:
                seqs.append(result[0])
                og_seqs[record.id] = extract_DNA_seq(record.seq[result[1]:],
                                                     seqs[-1], tab)
            elif join == 1:
                seqs.append(result[0][0]+result[1][0])
                og_seqs[record.id] = join_extract_DNA_seq(record.seq, result, tab)
            if og_seqs[record.id] == 1:
                shift = 1
            else:
                if len(og_seqs[record.id]) > 200 and join == 0:
                    try:
                        shift = detect_frameshift(seq, og_seqs[record.id], tab)
                    except:
                        shift = 1
                else:
                    shift = 0
            if shift == 0:
                n_check = check_n(og_seqs[ids[-1]])
            else:
                n_check = 0
            if shift == 1 or n_check == 1:
                ids = ids[:-1]
                names = names[:-1]
                seqs = seqs[:-1]
                del og_seqs[record.id]
                err.append(record.id)
    return seqs, names, ids, og_seqs, err


# Create multiFASTA file with protein sequences of the gene of interest
# and the in-frame sequences of the genomes to be multiple aligned and
# then align with Kalign.
def combine_align(records, ids, names, seqs):
    for i, n, s in zip(ids, names, seqs):
        records.append(SeqRecord(s, id=i, description=n))
    SeqIO.write(records, 'pre_align.fasta', 'fasta')
    invoke_kalign('pre_align.fasta', 'protein_align.fasta')
    subprocess.call(['rm', 'pre_align.fasta'])


# Restore original codon sequence while maintaining gaps and write to
# file.
def restore_codons(og_seqs, names):
    records = []
    for record in SeqIO.parse('protein_align.fasta', 'fasta'):
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
                                 description=names[record.id]))
    subprocess.call(['rm', 'protein_align.fasta'])
    SeqIO.write(records, 'codon_align.fasta', 'fasta')
    SeqIO.write(records, 'codon_align.clustal', 'clustal')


# Checks sequences for duplicates and compresses them into a single
# representaitve sequence if found
def compressor(seqs, names, ids, og_seqs):
    all_sames = [[]]
    for index, idd1 in enumerate(ids):
        for s in all_sames:
            if idd1 in s:
                break
            else:
                pass
        if idd1 in s:
            pass
        else:
            same = [idd1]
            for idd2 in ids[index+1:]:
                if og_seqs[idd1] == og_seqs[idd2]:
                    same.append(idd2)
            all_sames.append(same)
    all_sames = all_sames[1:]
    new_seqs = []
    new_names = []
    new_ids = []
    ref_id = str(set(og_seqs.keys())-set(ids))[2:-2]
    new_og_seqs = {}
    new_og_seqs[ref_id] = og_seqs[ref_id]
    count = 0
    for same in all_sames:
        if len(same) == 1:
            new_ids.append(same[0])
            new_names.append(names[ids.index(same[0])])
            new_seqs.append(seqs[ids.index(same[0])])
            new_og_seqs[new_ids[-1]] = og_seqs[new_ids[-1]]
        else:
            new_ids.append('MultiSeq'+str(count)+'_'+str(len(same)))
            count += 1
            name = ''
            for s in same:
                name += s
                name += ','
            new_names.append(name[:-1])
            new_seqs.append(seqs[ids.index(same[0])])
            new_og_seqs[new_ids[-1]] = og_seqs[same[0]]
    return new_seqs, new_names, new_ids, new_og_seqs


# Check input FASTA is occupied and doesn't use invalid characters
def check_input(reference, reads):
    records = []
    invalids = ['E', 'F', 'I', 'J', 'L', 'O', 'P', 'Q', 'X', 'Z']
    try:
        for record in SeqIO.parse(reference, 'fasta'):
            if any(x in invalids for x in record.seq):
                print(f'Input Error: Invalid characters detected in reference sequence')
                exit()
            records.append(records)
        if len(records) != 1:
            print('User Error: reference input should contain exactly one sequence')
            exit()
        records = []
        for record in SeqIO.parse(reads, 'fasta'):
            if any(x in invalids for x in record.seq):
                print(f'Input Error: Invalid characters detected in sequence ID: {record.id}')
                exit()
            records.append(record)
        if len(records) == 0:
            print('User Error: sequence FASTA file is empty')
            exit()
    except:
        print('User Error: improperly formated FASTA')
        exit()


# Check that translation table chosen is valid
def check_tab(tab):
    if tab == None:
        tab = 1
    else:
        if int(tab) in [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33]:
            tab = int(tab)
        else:
            print('User Error: Chosen translation table number is invalid. See: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for valid options')
            exit()
    return tab


# For when inputs are whole genomes
def genome_mode(reference, reads, start, end, compress, para, tab):
    check_input(reference, reads)
    tab = check_tab(tab)
    # Check that start and end coordinates are present and valid
    if ',' not in str(start) and ',' not in str(end):
        try:
            start = int(start)-1
            end = int(end)
        except:
            print('User Error: invalid start/end coordinate(s)')
            exit()
        join = 0
        if start >= end:
            print('User Error: start coordinate(s) must be less than the end coordinate(s)')
            exit()
    else:
        try:
            start = (int(start.split(',')[0])-1, int(start.split(',')[1])-1)
            end = (int(end.split(',')[0]), int(end.split(',')[1]))
        except:
            print('User Error: invalid start/end coordinate(s)')
            exit()
        join = 1
        if start[0] >= end[0] or start[1] >= end[1]:
            print('User Error: start coordinate(s) must be less than the end coordinate(s)')
            exit()

    # Find protein sequence of gene of interest and extract the original DNA
    # sequence (only for genome mode)
    og_seqs = {}
    for record in SeqIO.parse(reference, 'fasta'):
        idd = record.id
        name = record.description
        if isinstance(start, tuple):
            seq = (record.seq[start[0] % 3:].translate(table=tab)[start[0]//3:end[0]//3],
                   record.seq[start[1] % 3:].translate(table=tab)[start[1]//3:end[1]//3])
            og_seqs[record.id] = record.seq[start[0]:end[0]]+record.seq[start[1]:end[1]]
        else:
            seq = record.seq[start % 3:].translate(table=tab)[start//3:end//3]
            og_seqs[record.id] = record.seq[start:end]
    seqs, names, ids, og_seqs, err = create_lists(reads, seq, og_seqs, join, para, tab)
    if compress:
        seqs, names, ids, og_seqs = compressor(seqs, names, ids, og_seqs)
    if join == 1:
        seq = seq[0]+seq[1]
    elif join == 0 and seq[-1] == '*':
        seq = seq[:-1]
        if og_seqs[idd][-3:] in ['TAG', 'TAA', 'TGA']:
            og_seqs[idd] = og_seqs[idd][:-3]
    records = [SeqRecord(seq, id=idd, description=name)]
    if len(seqs) == 0:
        print('No homologous sequences were found')
        exit()
    combine_align(records, ids, names, seqs)
    names = dict(zip(ids, names))
    names[idd] = name
    restore_codons(og_seqs, names)
    if len(err) > 0:
        print('The following '+str(len(err))+' sequence(s) were suspected of '
              'containing frameshifts, an early stop codon, or contained too '
              'many Ns and so were thrown out before multiple alignment:')
        for e in err:
            print(e)
        print('\n')


# For when inputs are in-frame genes
def gene_mode(reference, reads, compress, tab):
    check_input(reference, reads)
    tab = check_tab(tab)
    ids = []
    names = []
    seqs = []
    og_seqs = {}
    err = []
    for record in SeqIO.parse(reference, 'fasta'):
        ids.append(record.id)
        names.append(record.description)
        seqs.append(record.seq.translate(table=tab))
        og_seqs[record.id] = record.seq
    for record in SeqIO.parse(reads, 'fasta'):
        if '*' in record.seq.translate(table=tab)[:-1]:
            err.append(record.id)
        else:
            ids.append(record.id)
            names.append(record.description)
            seqs.append(record.seq.translate(table=tab))
            og_seqs[record.id] = record.seq
    records = [SeqRecord(seqs[0], id=ids[0], description=names[0])]
    if compress:
        seqs, names, ids, og_seqs = compressor(seqs[1:], names[1:], ids[1:],
                                               og_seqs)
    combine_align(records, ids, names, seqs)
    new_names = dict(zip(ids, names))
    new_names[records[0].id] = records[0].description
    restore_codons(og_seqs, new_names)
    if len(err) > 0:
        print('The following '+str(len(err))+' sequence(s) were suspected of '
              'containing frameshifts, an early stop codon, or contained too '
              'many Ns and so were thrown out before multiple alignment:')
        for e in err:
            print(e)
        print('\n')


# For when reference input is an in-frame gene but the reads are whole genomes
def mixed_mode(reference, reads, compress, para, tab):
    check_input(reference, reads)
    tab = check_tab(tab)
    join = 0
    og_seqs = {}
    for record in SeqIO.parse(reference, 'fasta'):
        idd = record.id
        name = record.description
        seq = record.seq.translate(table=tab)
        og_seqs[record.id] = record.seq
    seqs, names, ids, og_seqs, err = create_lists(reads, seq, og_seqs, join, para, tab)
    if compress:
        seqs, names, ids, og_seqs = compressor(seqs, names, ids, og_seqs)
    if join == 1:
        seq = seq[0]+seq[1]
    elif join == 0 and seq[-1] == '*':
        seq = seq[:-1]
        og_seqs[idd] = og_seqs[idd][:-3]
    records = [SeqRecord(seq, id=idd, description=name)]
    if len(seqs) == 0:
        print('No homologous sequences were found')
        exit()
    combine_align(records, ids, names, seqs)
    names = dict(zip(ids, names))
    names[idd] = name
    restore_codons(og_seqs, names)
    if len(err) > 0:
        print('The following '+str(len(err))+' sequence(s) were suspected of '
              'containing frameshifts, an early stop codon, or contained too '
              'many Ns and so were thrown out before multiple alignment:')
        for e in err:
            print(e)
        print('\n')
