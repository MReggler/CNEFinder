FROM bioconductor/release_core2

# Set up a basic python 3 install
RUN ["apt-get", "update"]
RUN ["apt-get", "upgrade", "-y", "--fix-missing"]
RUN apt-get update && apt-get -y install python3.5 python3-pip

ADD scripts/pre_process.py /

# Needed for the pre_process.py script to run
RUN pip install rpy2

CMD [ "python", "./pre_process.py" ]
