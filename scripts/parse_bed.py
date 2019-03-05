from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict




def load_bed(bedfile):
    with open(bedfile) as f:
        for line in f:
            rchr, rstart, rend, qchr, qstart, qend, _, _, _ = line.split()
            ref_sequences[rchr].append((int(rstart), int(rend)))
            query_sequences[qchr].append((int(qstart), int(qend)))



def coords_to_cnes(fasta_file, sequences, out_file):
    """
    Function steps:
    1. parse ref/query .fa file and turn into dict
    2. loop through each `chr` in ref_sequences
    3. loop though each cne in `chr` + convert to ACGT alphabet
    4. write out
    """
    # source for basic structure:
    # https://stackoverflow.com/questions/30503543

    # misc related question
    # https://stackoverflow.com/questions/39249121

    ref_records = SeqIO.to_dict(SeqIO.parse(open(fasta_file), 'fasta'))

    short_seq_records = []
    for name in sequences:
        for (start, stop) in sequences[name]:
            long_seq_record = ref_records[name]
            long_seq = long_seq_record.seq
            alphabet = long_seq.alphabet
            short_seq = str(long_seq)[start:stop] #[start-1:stop]
            short_seq_record = SeqRecord(Seq(short_seq, alphabet), id=name, description='')
            short_seq_records.append(short_seq_record)

    # write to file
    with open(out_file, 'w') as f:
        SeqIO.write(short_seq_records, f, 'fasta')


def bed_row_to_list(index):
    pass

if __name__ == "__main__":
    ref_sequences = defaultdict(list)
    query_sequences = defaultdict(list)

    # handle these from argparsing
    bedfile = '../output/outfile.bed'
    ref_file = '../input/hg38.fa'
    query_file = '../input/galGal4.fa'

    # populate ref_sequences and query_sequences
    print("Converting .bed file is Python dict...")
    load_bed(bedfile)

    print("Generating .fasta file of CNEs...")
    coords_to_cnes(ref_file, ref_sequences, 'ref_cnes.fasta')
