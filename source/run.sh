#!/bin/sh

chmod +x ./cnef

# Example gene-based run:
# ./cnef -r hg38.fa  -q galGal4.fa -e hg38_exons -f galGal4_exons -g hg38_genes -n zeb2 -j galGal4_genes -m zeb2 -t 0.95 -l 150 -o output

# Example index-position-based run:
# ./cnef  -r hg38.fa  -q galGal4.fa  -e hg38_exons  -f galGal4_exons -y 2 -z 7  -a 144121063 -b 146282147 -c 31812269 -d 33916777 -t 0.95 -l 150 -o output

# Assumes `perfect` input. i.e., env vars exists and are valid for use
# FUTURE WORK: handle invalid REF_GENE and/or QUERY_GENE values
if [ -z "${REF_START}" ]; then
    echo "Running CNEFinder from gene names..."

    ./cnef -r ${REF_GENOME_FILE} -q ${QUERY_GENOME_FILE} \
        -e ${EXONS_REF_FILE} -f ${EXONS_QUERY_FILE} \
        -g ${GENES_REF_FILE} -n ${REF_GENE} \
        -j ${GENES_QUERY_FILE} -m ${QUERY_GENE} \
        -t ${SIM_THRESHOLD} -l ${MIN_SEQ_LENGTH} \
        -o ${OUTPUT_PATH}
else
    echo "Running CNEFinder from index positions..."

    ./cnef -r ${REF_GENOME_FILE} -q ${QUERY_GENOME_FILE} \
        -e ${EXONS_REF_FILE} -f ${EXONS_QUERY_FILE} \
        -y ${REF_CHROM} -z ${QUERY_CHROM} \
        -a ${REF_START} -b ${REF_END} -c ${QUERY_START} -d ${QUERY_END} \
        -t ${SIM_THRESHOLD} -l ${MIN_SEQ_LENGTH} \
        -o ${OUTPUT_PATH}
fi
