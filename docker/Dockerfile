# GCC support can be specified at major, minor, or micro version
# (e.g. 8, 8.2 or 8.2.0).
# See https://hub.docker.com/r/library/gcc/ for all supported GCC
# tags from Docker Hub.
# See https://docs.docker.com/samples/library/gcc/ for more on how to use this image
FROM gcc:latest

# copy files to dir in the image
# and set that as the working location
COPY . /usr/src/myapp
WORKDIR /usr/src/myapp

# Install prereqs and then compile CNEFinder
RUN bash -c ./pre-install.sh
RUN make -f Makefile

# run CNEFinder for search from index positions
# note: if a,b,c,d == 0 then searches within entire chromosomes

ENTRYPOINT ["./run.sh"]
