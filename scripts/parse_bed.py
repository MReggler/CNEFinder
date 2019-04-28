from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import os
import re
import json

def running_local(filename, flag=True):
    """Modifies a filename if script is not running in a container.

    Args:
        filename: the filename string to modify.
        flag: a boolean flag determining whether the modification
            should happen.

    Returns:
        A string containing the (adjusted) filename.
    """
    if flag:
        return ('.' + filename).strip()
    else:
        return (filename).strip()


def parse_env_file(env_file, in_docker=False):
    """Parses a files containing environment variables.

    Retrieves the environment variables that correspond to the
    reference and query genome file names.

    Args:
        env_file: a string holding the env file's name.
        in_docker: a boolean flag that holds whether the
            function is executing within a container.

    Returns:
        A tuple containing the reference and query genomes filenames
        as strings.
    """
    var, value = "", ""
    ref_file, query_file = "", ""

    with open(env_file) as f:
        for line in f:
            if line.strip():
                try:
                    var, value = line.split("=")
                except ValueError:
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
    """Converts a .bed file in to a pair of default dictionaries.

    Extracts the critical information - start and end indexs of CNEs -
    from the input .bed file and places them in a defaultdict where
    the key is the chromosome number.

    Args:
        bed_file: a string holding the .bed file's name.

    Returns:
        Two dict mapping chromosome name keys to a list of tuples
        containing the start and end index of CNEs described by the
        .bed file. For example:

        {'chr2': [(1224, 1452), (3456, 4134)], ...}
    """
    ref_sequences = defaultdict(list)
    query_sequences = defaultdict(list)

    with open(bed_file) as f:
        for line in f:
            rchr, rstart, rend, qchr, qstart, qend, _, _, _ = line.split()
            ref_sequences[rchr].append((int(rstart), int(rend)))
            query_sequences[qchr].append((int(qstart), int(qend)))

    return (ref_sequences, query_sequences)


def coords_to_cnes(fasta_file, sequences, out_file):
    """Converts CNE start and end coordinates into CNE FASTA sequences.

    Iterates through a defaultdict of CNE start and end index positions
    keyed to chromosome numbers, and uses a parsed FASTA file of the
    associated genome to map coordinates to ACTG sequences. Writes out the
    result to a file.

    Args:
        fasta_file: a string holding a genome-wide FASTA file's name.
        seqences: a defaultdict output of the load_bed() function.
        out_file: a string holding the name of the output file.
    """
    # Function steps:
    # 1. parse ref/query .fa file and turn into dict
    # 2. loop through each `chr` in ref_sequences
    # 3. loop though each cne in `chr` + convert to ACGT alphabet
    # 4. write out

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

    with open(out_file, 'w') as f:
        SeqIO.write(short_seq_records, f, 'fasta')


def field_map(dictseq, name, func):
    """Maps a function to an item value of a dict generator.

    Args:
        dictseq: the dict generator to whose item value the function
            will be mapped.
        name: the key of the item to be mapped.
        func: the function to be mapped.

    Yields:
        the mapped and unmapped items of the dict generator.
    """
    for d in dictseq:
        d[name] = func(d[name])
        yield d


def objectify_bed(bed_file):
    """Converts a .bed table into a list of dicts.

    Args:
        bed_file: a string holding the name of the .bed file

    Returns:
        a list of dicts holding the same information as the bed
        file, but with variable types instead of all strings.
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
    """Extracts all the sequences from a FASTA file.

    Args:
        fasta_file: a string holding the name of the FASTA file

    Returns:
        a list of upper-case strings holding the sequences present
        in the FASTA input file.
    """
    fasta_records = SeqIO.parse(open(fasta_file), 'fasta')
    return [str(rec.seq).upper() for rec in fasta_records]


def bed_to_json(bed_file, r_cnes, q_cnes, out_file):
    """Converts a .bed table into into a JSON object.

    Uses objectify_bed() to convert a .bed file into a list of dicts,
    then extracts the FASTA sequences from the reference and query
    CNE files to add additional information to the output file.

    Args:
        bed_file: a string holding the name of the .bed file.
        r_cnes: a string holding the name of reference CNEs FASTA file.
        q_cnes: a string holding the name of query CNEs FASTA file.
        out_file: a string holding the name the JSON output file.

    """
    row_dicts = objectify_bed(bed_file)
    ref_seqs = get_fasta_sequences(r_cnes)
    query_seqs = get_fasta_sequences(q_cnes)

    for idx, row_dict in enumerate(list(row_dicts)):
        row_dict['ref_sequence'] = ref_seqs[idx]
        row_dict['query_sequence'] = query_seqs[idx]

    with open(out_file, 'w') as f:
        json.dump(row_dicts, f)


def main(bed_file, env_file, r_cnes, q_cnes, json_table):
    """Handler function to parsing CNEFinder output bed file.

    Args:
        bed_file: the CNEFinder output file (exists).
        env_file: the run configuration .env file (exists).
        r_cnes: the reference cnes fasta file (to be made).
        q_cnes: the query cnes fasta file (to be made).
        json_table: the json table containing all data (to
            be made).
    """
    in_docker = False
    if os.environ.get('APP_ENV') == 'docker':
        in_docker = True
        print("parse_bed.py is running in a container.")
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
    print("Generating .fasta file of CNEs for {}".format(ref_file))
    coords_to_cnes(ref_file, ref_sequences, r_cnes)

    print("Generating .fasta file of CNEs for {}".format(query_file))
    coords_to_cnes(query_file, ref_sequences, q_cnes)

    # convert bed file to json object with addition of actual sequences.
    print("Generating json file of .bed + CNE sequences")
    bed_to_json(bed_file, r_cnes, q_cnes, json_table)

if __name__ == "__main__":
    bed_file = '/output/outfile.bed'
    env_file = '/input/metadata.txt' # make sure this is here
    r_cnes = '/output/ref_cnes.fa'
    q_cnes = '/output/query_cnes.fa'
    json_table = '/output/table.json'

    main(bed_file, env_file, r_cnes, q_cnes, json_table)