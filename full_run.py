import subprocess
import sys
import shutil

#---------------------------------------------------------
# PREREQS
# --------------------------------------------------------
# 1. Ensure necessary run configuration settings in .env
# 2. Download reference and query .fa files to ./input/
# 3. All docker containers have been built

# subprocess.call(['./build.sh'])



#---------------------------------------------------------
# STEP 1
# --------------------------------------------------------
# Copy the environment variables file into the ./input/
# and ./output volume mounts for later reference
shutil.copy2('./.env','./input/metadata.txt')
shutil.copy2('./.env','./ouput/metadata.txt')


#---------------------------------------------------------
# STEP 2
# --------------------------------------------------------
# This will populate the ./input/ volume mount with the
# gene and exon files needed by CNEFinder (it produces sets
# both to cover all options)
subprocess.call(['./run_preprocess.sh'])


#---------------------------------------------------------
# STEP 3
# --------------------------------------------------------
# This will run CNEFinder in a container and produce a
# .bed file of identified CNEs in the ./output/ volume mount
subprocess.call(['./run_cnefinder.sh'])


#---------------------------------------------------------
# STEP 4
# --------------------------------------------------------
# This will run parse_bed.py script to produce .fa files
# of the reference and query CNEs and will wrap all this
# information up into a single JSON file
subprocess.call(['./run_parse_bed.sh'])

# or run the python script locally (if biopython is installed)
#subprocess.Popen([sys.executable, "./scripts/parse_bed.py"])


#---------------------------------------------------------
# STEP 5
# --------------------------------------------------------
# Clean-up:
#   - Delete the contents of ./input to save spave
shutil.rmtree('./input') 
