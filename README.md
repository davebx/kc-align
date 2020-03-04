# Kc-Align

Kc-Align is a codon-aware multiple aligner that uses Kalgin2 to produce in-frame gapped codon alignments for selection analysis of small genomes (mostly viral and some smaller bacterial genomes). Takes nucleotide seqeunces as inputs, converts them to their in-frame amino acid sequences, performs multiple alignment with Kalign, and then converts the alignments back to their original codon sequence while preserving the gaps. Produces three outputs: the gapped nucleotide alignments in FASTA format and in CLUSTAL format and the amino acid level alignment.

Kc-Align will also attempt to detect any frameshift mutations in the input reads. If a frameshift is detected, that sequence will not be included in the multiple alignment and its ID will be printed to stdout.

# Modes

Kc-Align can be run in three different modes, depending on your input data. 

In "genome" mode, the "reference" and "reads" input parameters are all full genome FASTA files. This mode also requires the 1-based start and end position numbers corresponding to the gene you are interested in aligning from the reference input.

If both the "reference" and "reads" inputs are already in-frame genes, the "gene" mode should be used. This mode does not require start and end position parameters as the reference is already in-frame.

For the case when your "reference" is an in-frame gene while the "reads" are whole genomes, the "mixed" mode can be used. Like gene mode, this mode does not require the start and end point position parameters.

## USAGE:

#### kc-align --mode genome --reference [reference FASTA] --reads [reads FASTA] --start [start] --end [end]

Ex: kc-align --mode genome --reference ref.fasta --reads reads.fasta --start 3512 --end 7831

#### kc-align --mode gene --reference [reference FASTA] --reads [reads FASTA]

Ex: kc-align --mode gene --reference ref.fasta --reads reads.fasta

#### kc-align --mode mixed --reference [reference FASTA] --reads [reads FASTA]

Ex: kc-align --mode mixed --reference ref.fasta --reads reads.fasta

For genes that are split into two parts and later joined together:

#### kc-align --mode genome --reference [reference FASTA] --reads [reads FASTA] --start [start1,start2] --end [end1,end2]

Ex: kc-align --mode genome --reference ref.fasta --reads reads.fasta --start 3512,3511 --end 7831,3721
