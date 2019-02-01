# GCC support can be specified at major, minor, or micro version
# (e.g. 8, 8.2 or 8.2.0).
# See https://hub.docker.com/r/library/gcc/ for all supported GCC
# tags from Docker Hub.
# See https://docs.docker.com/samples/library/gcc/ for more on how to use this image
FROM gcc:latest

# All Args from docker-compose.yml
ARG ref_genome_file
ARG query_genome_file
ARG exons_ref_file
ARG exons_query_file
ARG min_seq_length
ARG sim_threshold
ARG output_file
ARG ref_chrom
ARG query_chrom
ARG ref_start
ARG ref_end
ARG query_start
ARG query_end

# copy files to dir in the image
# and set that as the working location
COPY . /usr/src/myapp
WORKDIR /usr/src/myapp

# Install prereqs and then compile CNEFinder
RUN bash -c ./pre-install.sh
RUN make -f Makefile

# run CNEFinder for search from index positions
# note: if a,b,c,d == 0 then searches within entire chromosomes

CMD ./cnef -r $REF_GENOME_FILE -q $QUERY_GENOME_FILE \
    -e $EXONS_REF_FILE -f $EXONS_QUERY_FILE \
    -y $REF_CHROM -z $QUERY_CHROM \
    -a $REF_START -b $REF_END -c $QUERY_START -d $QUERY_END \
    -t $SIM_THRESHOLD -l $MIN_SEQ_LENGTH \
    -o $OUTPUT_FILE