FROM gcc:latest

# copy over all cnefinder source files
COPY ./source /usr/src/cnefinder

# temp copy of exon and genes files for galGal4 + hg38
# (these will eventually come from volume mount)
COPY ./Experiments /usr/src/cnefinder/Experiments

WORKDIR /usr/src/cnefinder

# Install prereqs and then compile CNEFinder
RUN bash -c ./pre-install.sh
RUN make -f Makefile

# run CNEFinder for search from index positions
# note: if a,b,c,d == 0 then searches within entire chromosomes

ENTRYPOINT ["./run.sh"]
