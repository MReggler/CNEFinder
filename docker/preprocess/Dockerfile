FROM bioconductor/release_core2

# Set up a basic python 3 install
RUN ["apt-get", "update"]
RUN ["apt-get", "upgrade", "-y", "--fix-missing"]
RUN apt-get update && apt-get -y --no-install-recommends install locate build-essential \
    libssl-dev libffi-dev libbz2-dev liblzma-dev libpcre3-dev \
    python3.5 python3-pip python3-setuptools python3-dev python3-wheel libreadline-dev

COPY ./scripts/pre_process.py /
COPY ./scripts/pre_process.sh /
COPY ./docker/preprocess/requirements.txt /
WORKDIR /

# Needed for the pre_process.py script to run
#RUN pip3 install pandas tzlocal
#RUN pip3 install rpy2==2.9.5
# install python requirements
RUN pip3 install -r requirements.txt

ENTRYPOINT ["./pre_process.sh"]
