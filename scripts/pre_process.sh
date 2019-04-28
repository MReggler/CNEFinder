#!/bin/sh

# LOC=$(locate libR.so)

# echo "libR.so found at:"
# echo ${LOC}

echo $(ls /usr/local/lib/R/lib)

# R_LIB_PATH=$(echo ${LOC} | sed -e "s,/[^/]*$,,")

# echo "r library found at:"
# echo ${R_LIB_PATH}

export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/usr/local/lib/R/lib"

# run the python script
python3 ./pre_process.py
