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
        detection_of_genome_in_samples[sample_id] = number_of_detected_genes_in_sample > alpha * len(
            detection_of_genes)
    return detection_of_genome_in_samples


def get_adjusted_std_for_gene_id(data, gene_id, samples, mean_coverage_in_samples, detection_of_genes):
    """Returns the adjusted standard deviation for a gene_id """
    # Note: originally I thought I would only consider samples in which the genome was detected, but in fact,
    # if a gene is detected in a sample in which the genome is not detected then that is a good sign that this is
    #  a NTS gene. But I still kept here the original definition of adjusted_std
    # adjusted_std = np.std([d[gene_id, sample_id] / mean_coverage_in_samples[sample_id] for sample_id in samples if (
        #         detection_of_genes[gene_id][sample_id] and detection_of_genome_in_samples[sample_id])])
    adjusted_std = np.std([data[gene_id][sample_id]/mean_coverage_in_samples[sample_id] for sample_id in samples if (
                detection_of_genes[gene_id][sample_id])])
    return adjusted_std


def get_adjusted_stds(data, samples, mean_coverage_in_samples, detection_of_genes):
    adjusted_std = {}
    for gene_id in data:
        adjusted_std[gene_id] = get_adjusted_std_for_gene_id(data, gene_id, samples, mean_coverage_in_samples,
                                                             detection_of_genes)
    return adjusted_std


def get_taxon_specificity(adjusted_stds, beta):
    """For each gene if the adjusted standard deviation (to understand what this is refer to Alon Shaiber) is smaller
    than beta the the gene is a taxon-specific gene (TS), otherwise, it is a non-taxon-specific gene (NTS)"""
    taxon_specificity = {}
    for gene_id in adjusted_stds:
        if adjusted_stds[gene_id] < beta:
            taxon_specificity[gene_id] = 'TS'
        else:
            taxon_specificity[gene_id] = 'NTS'
    return taxon_specificity


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
        adjusted_stds = get_adjusted_stds(data,samples,mean_coverage_in_samples,detection_of_genes)
        taxon_specificity = get_taxon_specificity(adjusted_stds,beta)
        loss = get_loss_function()
        converged = True
    pass

def get_specificity_from_class_id(class_id):
    if class_id in [1, 2, 3,'1', '2', '3']:
        return 'TS'
    elif class_id in [4,5,'4','5']:
        return 'NTS'
    elif class_id in [0,'0']:
        return 'NaN'
    else:
        print('class_id %s is not a valid class_id' % class_id)
        exit(1)

def check_results_for_mock_data(input_name, alpha=0.5, beta=1, gamma=3):
    pwd = '/Users/alonshaiber/PycharmProjects/MACg/tests/sandbox/'
    file_path = pwd + input_name + '.txt'
    data, samples = get_data_from_txt_file(file_path)
    print('Loaded data from these samples: %s' % samples)
    mean_coverage_in_samples = get_mean_coverage_in_samples(data, samples)
    print(mean_coverage_in_samples)
    std_in_samples = get_std_in_samples(data, samples)
    print(std_in_samples)
    detection_of_genes = get_detection_of_genes(data, samples, mean_coverage_in_samples, std_in_samples, gamma)
    print(detection_of_genes)
    detection_of_genome_in_samples = get_detection_of_genome_in_samples(detection_of_genes, samples, alpha)
    print(detection_of_genome_in_samples)
    for gene_id in data:
        print('gene_id %s: adjasted_std: %s' % (gene_id, get_adjusted_std_for_gene_id(data, gene_id, samples,
                                                                                      mean_coverage_in_samples,
                                                                                      detection_of_genes)))
    adjusted_stds = get_adjusted_stds(data, samples, mean_coverage_in_samples, detection_of_genes)
    taxon_specificity = get_taxon_specificity(adjusted_stds, beta)
    print(taxon_specificity)

    ## comparing the specificity classification to the actual classification
    # loading actual classification from the table that was generated by gen_mock_data.py
    mock_gene_group_table_file = pwd + input_name + '_gene_group_table.txt'
    mock_gene_group_table = utils.get_TAB_delimited_file_as_dictionary(mock_gene_group_table_file, column_mapping=[int]*3)
    list_of_true_classifications=[]
    list_of_wrong_classifications = []
    for gene_id in mock_gene_group_table:
        if taxon_specificity[gene_id] == get_specificity_from_class_id(mock_gene_group_table[gene_id]['class_id']):
            list_of_true_classifications.append(gene_id)
        else:
            list_of_wrong_classifications.append(gene_id)
    if list_of_wrong_classifications == []:
        print('Success! All genes were classified properly as TS or NTS (-:')
    else:
        print('There are %s genes that were wrongly classified, and here they are: %s' % (len(
            list_of_wrong_classifications), list_of_wrong_classifications))


def main(file_path, alpha=0.5, beta=1, gamma=3):
    pass


if __name__ == '__main__':
    input_name = 'test_20'
    check_results_for_mock_data(input_name)