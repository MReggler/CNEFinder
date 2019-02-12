#!/bin/sh

chmod +x ./cnef

./cnef -r ${REF_GENOME_FILE} -q ${QUERY_GENOME_FILE} \
    -e ${EXONS_REF_FILE} -f ${EXONS_QUERY_FILE} \
    -y ${REF_CHROM} -z ${QUERY_CHROM} \
    -a ${REF_START} -b ${REF_END} -c ${QUERY_START} -d ${QUERY_END} \
    -t ${SIM_THRESHOLD} -l ${MIN_SEQ_LENGTH} \
    -o ${OUTPUT_PATH}
