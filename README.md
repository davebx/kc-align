# Kc-Align

Kc-Align is a fast and accurate tool for performing codon-aware multiple sequence alignments (MSA). It makes use of the very fast multiple alignment program Kalign3 to ensure maximum speed. Kc-Align is a extremely extremely versatile tool, capable of taking a variety of inputs and achieving the same result. Every other aligner requires the sequence inputs to be the coding sequences of the gene/open reading frame (ORF) to be aligned, requiring curation from the whole-genome sequences and also preventing use of assemblies that may not be properly annotated. Kc-Align solves this problem by using pairwise alignments to extract the sequence from each whole genome that is homologous to the sequence of a high quality and well annotated reference sequence. This feautre may also be bypassed for those who already have curated data and simply desire a quick and accurate codon-aware multiple aligner (see Modes below).


## Obtaining Kc-Align

Kc-Align is availbe through PyPI (`pip install kcalign`) and through Bioconda (`conda install kcalign`). Alternatively, a GUI interface for Kc-Align is installed and available for use on Galaxy (http://usegalaxy.org).

## Using Kc-Align

### USAGE:

`kc-align --mode genome --reference [reference sequence] --sequences [other seqs to align] --start [start coordinate] --end [end coordinate]`

#### (or)

`kc-align --mode [gene | mixed] --reference [reference sequence] --sequences [other seqs to align]`

### Arguments:

```
--mode/-m         Alignment mode (genome, gene, or mixed) (required)

--reference/-r    Reference sequences to align against (required)

--sequences/-S    Other sequences to align (required)

--start/-s        Start position (required in genome mode)

--end/-e          End position (required in genome mode)

--compress/-c     Compress identical sequences

--parallel/-p     Enable parallelization of homology search
```

### Modes

Kc-Align can be run in three different modes, depending on your input data. These modes are: genome, gene, and mixed.

#### Genome

If the reference sequences and all other sequences to be aligned are full genomes, this mode should be used. This mode requires the start and end coordinates of the gene/ORF to be aligned with regard to the reference sequence. Kc-Align will excise the appropriate subsequence from the referencce and then use pairwise alignments to find the corresponding homologous sequences from each of the other input genomes. It will then perform the MSA using the extracted sequences.

If the gene/ORF you wish to align exists in multiple distinct segments in the genome (ribosomal frameshift during transcription), Kc-Align can find each homologous segment separately for each input sequence and then concatenate them before performing the final MSA. This requires the user to enter the start and end coordinates of each segment as a comma-separated list following their respective arguments.

##### Example

`kc-align -m genome -r reference.fasta -S sequences.fasta -s 3532,3892 -e 3894,5326`

In the above example, the two segments being aligned have the coordinates 3532-3894 and 3892-5326.

#### Gene

If the input sequences have already been trimmed to the coding sequences of the gene/ORF of interest, then this mode can be used and Kc-Align will not perform the homology search as it does in genome mode. It instead simply immediately performs the MSA, making this mode much faster than the others.

#### Mixed

For the case when your reference is a coding sequence while all other sequences are whole genomes. Like gene mode, this mode does not require the start and end point position parameters but like genome mode it will perform homology searching in order to extract the sequences homologous to the reference from the other input sequences.

### Outputs

Kc-Align will output two files: a FASTA format alignment and a CLUSTAL format alignment.

### Compress Identical Sequences

If the `--compress/-c` parameter is specified, Kc-Align will compress identical sequences into a single sequence. In the FASTA output, compressed sequences will have an ID of the form MultiSeq[incremental index]_[number of sequences that were compressed] (ex: MultiSeq3_321, third compression with 321 sequences having that same sequence) while the description field is a comma-separated list of every ID that was compressed into that single sequence. The reference sequence will not be compressed.

### Parallelization

If the `--parallel/-p` parameter is used in genome or mixed mode, the calculations for the homology search will be split between 3 cores (if possible), decreasing runtimes by up to 35%.
