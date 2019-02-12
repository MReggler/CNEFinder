from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

# source for basic structure:
# https://stackoverflow.com/questions/30503543

# misc related question
# https://stackoverflow.com/questions/39249121


ref_sequences = defaultdict(list)
query_sequences = defaultdict(list)

with open('../output/outfile.bed') as f:
    for line in f:
        rchr, rstart, rend, qchr, qstart, qend, rl, ql, sim = line.split()
        ref_sequences[rchr].append((int(rstart), int(rend)))
        query_sequences[qchr].append((int(qstart), int(qend)))

# parse fasta file and turn into dictionary
ref_records = SeqIO.to_dict(SeqIO.parse(open('../input/hg38.fa'), 'fasta'))

short_seq_records = []
for name in ref_sequences:
    for (start, stop) in ref_sequences[name]:
        long_seq_record = ref_records[name]
        long_seq = long_seq_record.seq
        alphabet = long_seq.alphabet
        short_seq = str(long_seq)[start-1:stop]
        short_seq_record = SeqRecord(Seq(short_seq, alphabet), id=name, description='')
        short_seq_records.append(short_seq_record)

# write to file
with open('output.fasta', 'w') as f:
    SeqIO.write(short_seq_records, f, 'fasta')
