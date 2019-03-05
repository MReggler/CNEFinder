from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import os
def running_local(filename, flag=True):
    if flag:
        return ('..' + filename).strip()
    else:
        return (filename).strip()


def parse_env_file(env_file, in_docker=False):
    var, value = "", ""
    ref_file, query_file = "", ""

    with open(env_file) as f:
        for line in f:
            if line.strip():
                try:
                    var, value = line.split("=")
                except ValueError as e:
                    pass

            if var in ["REF_GENOME_FILE"]:
                ref_file = value
            elif var in ["QUERY_GENOME_FILE"]:
                query_file = value

    if not in_docker:
        return (running_local(ref_file), running_local(query_file))
    else:
        return (ref_file, query_file)


def load_bed(bed_file):
    ref_sequences = defaultdict(list)
    query_sequences = defaultdict(list)

    with open(bed_file) as f:
        for line in f:
            rchr, rstart, rend, qchr, qstart, qend, _, _, _ = line.split()
            ref_sequences[rchr].append((int(rstart), int(rend)))
            query_sequences[qchr].append((int(qstart), int(qend)))

    return (ref_sequences, query_sequences)


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
    bed_file = '/output/outfile.bed'
    env_file = '/input/metadata.txt' # make sure this is here
    r_cnes = '/output/ref_cnes.fa'
    q_cnes = '/output/query_cnes.fa'
    json_table = '/output/table.json'

    in_docker = False
    if os.environ.get('APP_ENV') == 'docker':
        in_docker = True
    else:
        print("parse_bed.py is not running in a container.")

    # if not on docker, handle fact that `input` and `output`
    # are not mounted at root loc: /
    files = [bed_file, env_file, r_cnes, q_cnes, json_table]
    bed_file, env_file, r_cnes, q_cnes, json_table \
        = [running_local(f, not in_docker) for f in files]

    # get two fasta filenames from CNEfinder's .env
    ref_file, query_file = parse_env_file(env_file, in_docker)

    # populate ref_sequences and query_sequences dicts
    print("Converting .bed file into Python dicts...")
    ref_sequences, query_sequences = load_bed(bed_file)

    # create output .fa files (takes ages so check when local testing)
    if not os.path.isfile(ref_file):
        print("Generating .fasta file of CNEs for {}".format(ref_file))

    print("Generating .fasta file of CNEs...")
    coords_to_cnes(ref_file, ref_sequences, 'ref_cnes.fasta')
