FROM python:3

# Add main python script
COPY scripts/parse_bed.py /
WORKDIR /

# Needed for SeqIO to parse .fasta files
RUN pip install biopython

CMD [ "python", "./parse_bed.py" ]
