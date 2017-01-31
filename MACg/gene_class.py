# -*- coding: utf-8
# pylint: disable=line-too-long

import csv
import numpy as np
import sys
import os
import anvio.utils as utils

def get_data_from_txt_file(file_path):
    """ Reads the coverage data from TAB delimited file """
    samples = utils.get_columns_of_TAB_delim_file(file_path)
    data = utils.get_TAB_delimited_file_as_dictionary(file_path, column_mapping=[int] + [float] * len(samples))
    return data, samples


def apply_func_to_genes_in_sample(data,samples, func, list_of_genes=None):
    """ Apply the give function on the list of genes in each sample. The function is expected to accept a list """
    if list_of_genes is None:
        list_of_genes = data.keys()
    d = dict(zip(samples, [next(map(func, [[data[gene_id][sample_id] for gene_id in list_of_genes]])) for
                       sample_id in samples]))
    return d


def get_mean_coverage_in_samples(data,samples,list_of_genes=None):
    """ Returns a dictionary with of the average coverage value of the list of genes per sample. if no list of genes is
    supplied then the average is calculated over all genes """
    mean_coverage_in_samples = apply_func_to_genes_in_sample(data, samples, np.mean, list_of_genes)
    return mean_coverage_in_samples


def get_std_in_samples(data,samples,list_of_genes=None):
    """ Returns a dictionary with of the standard deviation of the coverage values of the list of genes per sample.
    if no list of genes is supplied then the average is calculated over all genes """
    std_in_samples = apply_func_to_genes_in_sample(data, samples, np.std, list_of_genes)
    return std_in_samples


def get_detection_of_genes(data, samples, mean_coverage_in_samples, std_in_samples, gamma):
    """ Returns a dictionary (of dictionaries), where for each gene_id, and each sample_id the detection of the gene
    is determined. The criteria for detection is having coverage that is greater than 0 and also that is not more
    than 3 (assuming gamma=3 for example) standard deviations below the mean coverage in the sample.
    Notice that the mean coverage isn't the mean of all genes in the sample necesarilly. In fact it would be the mean of
    only the taxon-specific genes."""
    detection_of_genes = {}
    for gene_id in data:
        detection_of_genes[gene_id] = {}
        for sample in samples:
            detection_of_genes[gene_id][sample] = data[gene_id][sample] > max(0,mean_coverage_in_samples[sample] -
                                                                             gamma*std_in_samples[sample])
            if data[gene_id][sample] > 0 and data[gene_id][sample] < mean_coverage_in_samples[sample] - \
                    gamma*std_in_samples[sample]:
                print('gene %s, in sample %s has non-zero coverage %s, and it has been marked as not detected due to '
                      'the detection criteria' % (gene_id, sample, data[gene_id][sample]))
    return detection_of_genes


def get_detection_of_genome_in_samples(detection_of_genes, samples, alpha):
    detection_of_genome_in_samples = {}
    for sample_id in samples:
        number_of_detected_genes_in_sample = len([gene_id for gene_id in detection_of_genes if detection_of_genes[
            gene_id][sample_id]])
        print(number_of_detected_genes_in_sample)
        detection_of_genome_in_samples[sample_id] = number_of_detected_genes_in_sample > alpha * len(
            detection_of_genes)
        print(alpha * len(detection_of_genes))
        print(detection_of_genome_in_samples[sample_id])
    return detection_of_genome_in_samples

def get_gene_classes(data, samples, alpha, beta, gamma):
    """ returning the classification per gene along with detection in samples (i.e. for each sample, whether the
    genome has been detected in the sample or not """
    taxon_specific_genes = None
    converged = False
    loss = None
    while not converged:
        # mean of coverage of all TS genes in each sample
        mean_coverage_in_samples = get_mean_coverage_in_samples(data,samples,taxon_specific_genes)
        std_in_samples = get_std_in_samples(data, samples)
        detection_of_genes = get_detection_of_genes(data, samples, mean_coverage_in_samples, std_in_samples, gamma)
        detection_of_genome_in_samples = get_detection_of_genome_in_samples(detection_of_genes, alpha)
        converged = True
    pass


def main(file_path, alpha=0.5, beta=1, gamma=3):
    data, samples = get_data_from_txt_file(file_path)
    print('Loaded data from these samples: %s' % samples)
    mean_coverage_in_samples = get_mean_coverage_in_samples(data, samples)
    print(mean_coverage_in_samples)
    std_in_samples = get_std_in_samples(data,samples)
    print(std_in_samples)
    detection_of_genes = get_detection_of_genes(data, samples, mean_coverage_in_samples, std_in_samples, gamma)
    print(detection_of_genes)
    detection_of_genome_in_samples = get_detection_of_genome_in_samples(detection_of_genes, samples, alpha)
    print(detection_of_genome_in_samples)


if __name__ == '__main__':
    input_name = 'test_20'
    file_path = '/Users/alonshaiber/PycharmProjects/MACg/tests/sandbox/' + input_name + '.txt'
    main(file_path)