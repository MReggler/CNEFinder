from rpy2.robjects.vectors import StrVector
from rpy2.robjects import pandas2ri
from rpy2.robjects import r as R
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
import os
import pandas as pd
import rpy2.robjects as ro

pandas2ri.activate()

def in_container():
    if os.environ.get('APP_ENV') == 'docker':
        return True
    return False


def parse_env_file(env_file, in_docker=False):
    var, value = "", ""
    ref_info, query_info = {}, {}

    with open(env_file) as f:
        for line in f:
            if line.strip():
                try:
                    var, value = line.split("=")
                except ValueError:
                    pass

            if var in ["REF_MART", "REF_HOST", "REF_DATASET"]:
                ref_info[var] = value.rstrip()
            elif var in ["QUERY_MART", "QUERY_HOST", "QUERY_DATASET"]:
                query_info[var] = value.rstrip()

    return (ref_info, query_info)


def use_ensembl(biomart, dataset, host, verbose=False):
    return R.useEnsembl(biomart, dataset, host, verbose=verbose)


def get_exons(mart):
    exons = R.getBM(attributes = StrVector(("chromosome_name",
                "exon_chrom_start", "exon_chrom_end", "strand")),
                mart=mart)

    print(exons.head())

    exons_ranges = R.GRanges(
        seqnames=exons.rx2('chromosome_name'),
        ranges=R.IRanges(start=exons.rx2('exon_chrom_start'),
                         end=exons.rx2('exon_chrom_end')),
        strand='+' if exons.rx2('strand') == '1L' else '-')

    # This was hell to find
    # https://stackoverflow.com/questions/38806898/
    set_method = R("`seqlevelsStyle<-`")
    exons_ranges = set_method(exons_ranges, "UCSC")

    as_data_frame = R("function(x) as.data.frame(x)")
    exons_ranges_df = as_data_frame(exons_ranges)

    print(exons_ranges_df.head())

    #exons_df = pandas2ri.ri2py(exons_ranges_df)

    return exons_ranges_df


def search_datasets(mart, pattern):
    return  R.searchDatasets(mart=mart, pattern=pattern)


def get_genes(mart):
    genes = R.getBM(
        attributes = StrVector(("external_gene_name", "chromosome_name",
            "start_position", "end_position")),
        mart=mart)

    # with localconverter(ro.default_converter + pandas2ri.converter):
    #     genes_df = ro.conversion.ri2py(genes)

    #genes_df = pandas2ri.ri2py_dataframe(genes)
    genes_df = pandas2ri.ri2py(genes)
    genes_df.rename(columns={'external_gene_name': 'gene_name'}, inplace=True)

    return genes_df


def tab_delim_file_pd(pd_df, filename, on_docker=False):
    if on_docker:
        filename = '/output/' + filename
    pd_df.to_csv(filename, sep='\t', encoding='utf-8', index=False)


def tab_delim_file_rpy2(rpy_df, filename, on_docker=False):
    if on_docker:
        filename = '/output/' + filename
    rpy_df.to_csvfile(filename, sep='\t', row_names=False)


if __name__ == "__main__":
    local_test = False

    base = importr('base')
    dollar = base.__dict__['$']

    # if testing locally get library path in which
    # biomaRt has been installed from env variable
    if local_test:
        base._libPaths(os.getenv("LIB_PATH"))

    # load+attach biomaRt package at `top` scope
    R.library("biomaRt")
    R.library("GenomicRanges")

    #---------------------------------------------------------------------------------------------
    # LOCAL TEST
    #---------------------------------------------------------------------------------------------
    if local_test:
        mart_test = "ENSEMBL_MART_ENSEMBL"
        host = "dec2015.archive.ensembl.org" # "jul2016.archive.ensembl.org" "ensembl.org"
        dataset = "ggallus_gene_ensembl" # "hsapiens_gene_ensembl"

        mart_obj = use_ensembl(mart_test, dataset, host, verbose=False)
        dataset_info = pandas2ri.ri2py(search_datasets(mart_obj, dataset))
        dataset_version = dataset_info.at[0, 'version']
        genes_filename = "{}_genes".format(dataset_version)
        exons_filename = "{}_exons".format(dataset_version)

        # get all genes in all chromosomes
        gene_ranges = get_genes(mart_obj)
        tab_delim_file_pd(gene_ranges, genes_filename, not local_test)

        # get all exonic coordinates in all chromosomes
        exons_ranges = get_exons(mart_obj)
        tab_delim_file_rpy2(exons_ranges, exons_filename, not local_test)

    #---------------------------------------------------------------------------------------------
    # ON DOCKER
    #---------------------------------------------------------------------------------------------
    else:
        ref_info, query_info = parse_env_file('/input/metadata.txt')

        marts = [ref_info.get('REF_MART'), query_info.get('QUERY_MART')]
        hosts = [ref_info.get('REF_HOST'), query_info.get('QUERY_HOST')]
        datasets = [ref_info.get('REF_DATASET'), query_info.get('QUERY_DATASET')]

        # weirdly pointless loop to handle both the ref + query data in same script
        for i in [0, 1]:
            # Set up the connection to the mart
            print("Connecting to mart {} at host {} for dataset {}".format(marts[i], hosts[i], datasets[i]))
            mart_obj = use_ensembl(marts[i], datasets[i], hosts[i], verbose=False)

            if i == 0:
                genes_filename = "ref_genes"
                exons_filename = "ref_exons"
            else:
                genes_filename = "query_genes"
                exons_filename = "query_exons"

            # get all genes in all chromosomes
            gene_ranges = get_genes(mart_obj)
            tab_delim_file_pd(gene_ranges, genes_filename, not local_test)

            # get all exonic coordinates in all chromosomes
            exons_ranges = get_exons(mart_obj)
            tab_delim_file_rpy2(exons_ranges, exons_filename, not local_test)
