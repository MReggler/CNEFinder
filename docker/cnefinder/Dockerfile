FROM gcc:latest

# copy over all cnefinder source files
COPY ./source /usr/src/cnefinder

WORKDIR /usr/src/cnefinder

# Install prereqs and then compile CNEFinder
RUN bash -c ./pre-install.sh
RUN make -f Makefile

# run CNEFinder for search from index positions
# note: if a,b,c,d == 0 then searches within entire chromosomes

ENTRYPOINT ["./run.sh"]
