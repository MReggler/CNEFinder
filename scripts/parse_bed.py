from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import os
import re
import json

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


def field_map(dictseq, name, func):
    for d in dictseq:
        d[name] = func(d[name])
        yield d


def objectify_bed(bed_file):
    """
    Converts .bed `table` into a list of dicts.
    """
    colnames = ('ref_chromosome', 'ref_start_coord', 'ref_end_coord',
                'query_chromosome', 'query_start_coord','query_end_coord',
                'ref_cne_length', 'query_cne_length', 'similarity')

    logpats = r'(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+' \
              r'(\S+)\s+(\S+)\s+(\S+)\s+(\S+)'
    logpat = re.compile(logpats)

    with open(bed_file) as f:
        rows   = (line.split() for line in f)
        groups = (logpat.match(line) for line in f)
        tuples = (g.groups() for g in groups if g)
        dicts  = (dict(zip(colnames, t)) for t in tuples)
        dicts  = field_map(dicts, 'ref_start_coord', int)
        dicts  = field_map(dicts, 'ref_end_coord', int)
        dicts  = field_map(dicts, 'query_start_coord', int)
        dicts  = field_map(dicts, 'query_end_coord', int)
        dicts  = field_map(dicts, 'ref_cne_length', int)
        dicts  = field_map(dicts, 'query_cne_length', int)

        return list(dicts)


def get_fasta_sequences(fasta_file):
    fasta_records = SeqIO.parse(open(fasta_file), 'fasta')
    return [str(rec.seq).upper() for rec in fasta_records]


def bed_to_json(bed_file, r_cnes, q_cnes, out_file):
    row_dicts = objectify_bed(bed_file)
    ref_seqs = get_fasta_sequences(r_cnes)
    query_seqs = get_fasta_sequences(q_cnes)

    for idx, row_dict in enumerate(list(row_dicts)):
        row_dict['ref_sequence'] = ref_seqs[idx]
        row_dict['query_sequence'] = query_seqs[idx]

    with open(out_file, 'w') as f:
        json.dump(row_dicts, f)


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
        coords_to_cnes(ref_file, ref_sequences, r_cnes)
    if not os.path.isfile(ref_file):
        print("Generating .fasta file of CNEs for {}".format(query_file))
        coords_to_cnes(query_file, ref_sequences, q_cnes)

    # convert bed file to json object with addition of actual sequences.
    print("Generating json file of .bed + CNE sequences")
    bed_to_json(bed_file, r_cnes, q_cnes, json_table)
