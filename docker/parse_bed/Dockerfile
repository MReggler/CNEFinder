FROM python:3

# Add main python script
ADD scripts/parse_bed.py /

# Needed for SeqIO to parse .fasta files
RUN pip install biopython

CMD [ "python", "./my_script.py" ]
