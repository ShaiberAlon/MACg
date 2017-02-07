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


def get_detection_of_genome_in_samples(detection_of_genes, samples, alpha, genes_to_consider=None):
    if genes_to_consider is None:
        # if no list of genes is supplied then considering all genes
        genes_to_consider = detection_of_genes.keys()
    detection_of_genome_in_samples = {}
    for sample_id in samples:
        detection_of_genome_in_samples[sample_id] = {}
        number_of_detected_genes_in_sample = len([gene_id for gene_id in genes_to_consider if detection_of_genes[
            gene_id][sample_id]])
        detection_of_genome_in_samples[sample_id]['detection'] = number_of_detected_genes_in_sample > alpha * len(
            genes_to_consider)
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


def get_loss_function_value(taxon_specificity, adjusted_stds, beta):
    loss = 0
    for gene_id in taxon_specificity:
        if taxon_specificity[gene_id] == 'TS':
            # Notice: here adjusted std includes the samples that don't have the genome detected in them (it kind of
            # makes sense, because if the gene is detected even though the genome is not, then maybe it's not
            # taxon-specific
            loss += adjusted_stds[gene_id]
        else:
            loss += beta
    return loss


def get_number_of_detections_for_gene(detection_of_genes, gene_id, samples):
    detections = 0
    for sample_id in samples:
        detections += detection_of_genes[gene_id][sample_id]
    return detections


def get_core_accessory_info(detection_of_genes, gene_id, samples_with_genome, eta):
    """ Returns 'core'/'accessory' classification for each gene. This is done using only the samples in which the
    genome is detected """
    if get_number_of_detections_for_gene(detection_of_genes, gene_id, samples_with_genome) < eta * len(samples_with_genome):
        return 'accessory'
    else:
        return 'core'


def get_gene_class(taxon_specificity, core_or_accessory):
    if taxon_specificity == 'TS':
        if core_or_accessory == 'core':
            return 'TSC'
        elif core_or_accessory == 'accessory':
            return 'TSA'
        else:
            print('%s is not valid. Value should be \'core\' or \'accessory\'' % core_or_accessory)
            exit(1)
    elif taxon_specificity == 'NTS':
        if core_or_accessory == 'core':
            return 'TNC'
        elif core_or_accessory == 'accessory':
            return 'TNA'
        else:
            print('%s is not valid. Value should be \'core\' or \'accessory\'' % core_or_accessory)
            exit(1)
    else:
        print('%s is not valid. Value should be \'TS\' or \'NTS\'' % taxon_specificity)
        exit(1)


def get_gene_classes(data, samples, alpha, beta, gamma, eta):
    """ returning the classification per gene along with detection in samples (i.e. for each sample, whether the
    genome has been detected in the sample or not """
    taxon_specific_genes = None
    converged = False
    loss = None
    TSC_genes = None
    while not converged:
        # mean of coverage of all TS genes in each sample
        mean_coverage_of_TS_in_samples = get_mean_coverage_in_samples(data,samples,taxon_specific_genes)
        # Get the standard deviation of the taxon-specific genes in a sample
        # TODO: right now, single copy, and multi-copy genes would be treated identically. Hence, multi-copy genes
        # would skew both the mean and the std of the taxon-specific genes.
        std_of_TS_in_samples = get_std_in_samples(data, samples, taxon_specific_genes)
        detection_of_genes = get_detection_of_genes(data, samples, mean_coverage_of_TS_in_samples, std_of_TS_in_samples, gamma)
        detection_of_genome_in_samples = get_detection_of_genome_in_samples(detection_of_genes, samples, alpha, TSC_genes)
        samples_with_genome = [sample_id for sample_id in samples if detection_of_genome_in_samples[sample_id][
            'detection']]
        adjusted_stds = get_adjusted_stds(data,samples,mean_coverage_of_TS_in_samples,detection_of_genes)
        taxon_specificity = get_taxon_specificity(adjusted_stds,beta)
        new_loss = get_loss_function_value(taxon_specificity, adjusted_stds, beta)
        epsilon = 0.5 * beta
        if loss is not None:
            if abs(new_loss - loss) < epsilon:
                converged = True
        loss = new_loss
        print('current value of loss function: %s ' % loss)

        gene_class_information = {}
        for gene_id in data:
            gene_class_information[gene_id] = {}
            gene_class_information[gene_id]['gene_specificity'] = taxon_specificity[gene_id]
            gene_class_information[gene_id]['number_of_detections'] = get_number_of_detections_for_gene(
                detection_of_genes, gene_id, samples)
            gene_class_information[gene_id]['core_or_accessory'] = get_core_accessory_info(detection_of_genes, gene_id,
                                                                                           samples_with_genome, eta)
            gene_class_information[gene_id]['gene_class'] = get_gene_class(gene_class_information[gene_id][
                                               'gene_specificity'], gene_class_information[gene_id]['core_or_accessory'])

        TSC_genes = [gene_id for gene_id in gene_class_information if gene_class_information[gene_id][
            'gene_class']=='TSC']
    final_detection_of_genome_in_samples = get_detection_of_genome_in_samples(detection_of_genes, samples, alpha,
                                                                        genes_to_consider=TSC_genes)
    return gene_class_information, final_detection_of_genome_in_samples

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
    taxon_specificity = get_gene_classes(data, samples, alpha, beta, gamma)
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


def main(file_path, additional_layers_file, sample_information_txt, alpha, beta, gamma, eta, additional_layers_to_append):
    data, samples = get_data_from_txt_file(file_path)
    gene_class_information, detection_of_genome_in_samples = get_gene_classes(data, samples, alpha, beta, gamma, eta)
    a = lambda dictionary, field, value : len([dict_id for dict_id in dictionary if dictionary[
        dict_id][field]==value])
    number_of_TS = a(gene_class_information, 'gene_specificity','TS')
    number_of_TSC = a(gene_class_information, 'gene_class','TSC')
    number_of_TSA = a(gene_class_information, 'gene_class','TSA')
    number_of_TNC = a(gene_class_information, 'gene_class','TNC')
    number_of_TNA = a(gene_class_information, 'gene_class','TNA')
    number_of_positive_samples = a(detection_of_genome_in_samples, 'detection', True)

    print('The number of TS is %s' %number_of_TS )
    print('The number of TSC is %s' % number_of_TSC)
    print('The number of TSA is %s' % number_of_TSA)
    print('The number of TNC is %s' % number_of_TNC)
    print('The number of TNA is %s' % number_of_TNA)
    print('The number of samples with the genome is %s' % number_of_positive_samples)
    if additional_layers_to_append is None:
        additional_column_titles = []
        additional_layers_dict = gene_class_information
    else:
        additional_column_titles = utils.get_columns_of_TAB_delim_file(additional_layers_to_append)
        additional_layers_dict = utils.get_TAB_delimited_file_as_dictionary(additional_layers_to_append,
                                                                            dict_to_append=gene_class_information,
                                                                            assign_none_for_missing=True)
    utils.store_dict_as_TAB_delimited_file(additional_layers_dict, additional_layers_file,headers=['gene_callers_id',
                                                                                                   'gene_class',
                                                                                                   'number_of_detections'] + additional_column_titles)
    utils.store_dict_as_TAB_delimited_file(detection_of_genome_in_samples, sample_information_txt,
                                               headers=['samples','detection'])


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--data', metavar='FILE', dest='input_file', help='input table of coverage values of '
                                                                               'genes in samples')
    parser.add_argument('-a', '--alpha', metavar='NUM', dest='alpha', type=float, default=0.5,
                        help='portion of TSC genes required to decide a genome is detected in a sample')
    parser.add_argument('-b', '--beta', metavar='NUM', dest='beta', type=float, default=1,
                        help='Weight of the number of non-taxon-specific genes in the loss function (default = 1). '
                             'This means that if the adjusted standard deviation of a gene is greater than beta, '
                             'then the gene will be classified as non-taxon-specific')
    parser.add_argument('-g', '--gamma', metavar='NUM', dest='gamma', type=float, default=3,
                        help='This is used to define detection of genes. It is the number of standard deviation below '
                             'the mean value of the taxon specific genes for which a gene is still considered detected')
    parser.add_argument('-e', '--eta', metavar='NUM', dest='eta', type=float, default=0.95,
                        help='Treshold for deciding whether a gene is core or accessory (default is 0.95, hence if a '
                             'gene is detected in at least 95% of the samples (in which the genome is detected) then '
                             'it is considered core')
    parser.add_argument('-o', '--out', metavar='FILE', dest='output', help='Output file for classes information')
    parser.add_argument('-s', '--sample-detection', metavar='FILE', dest='sample_detection_output', help='Output file '
                                                                                   'for sample detection information')
    parser.add_argument('-A', '--additional-layers', metavar='FILE', dest='additional_layers_to_append', default=None,
                        help='An additional layer file to append to the one created by the algorithm')
    # parser.add_argument('--test', action='store_true', dest='test', help='test that everything is ok and exit')
    args = parser.parse_args()

    main(args.input_file, args.output, args.sample_detection_output, args.alpha, args.beta, args.gamma, args.eta,
         args.additional_layers_to_append)
    # # check_results_for_mock_data(input_name)
    # additional_layers_file = '/Users/alonshaiber/PycharmProjects/MACg/tests/sandbox/my_best_delete_me_so_far'
    # input_file = '/Users/alonshaiber/PycharmProjects/MACg/tests/sandbox/p214_Bfrag_positive_with_M_GG_gene_coverage.txt'
    # old_sample_information_txt = '/Users/alonshaiber/PycharmProjects/MACg/tests/sandbox' \
    #                           '/Bfrag_positive_samples_information.txt'
    # new_sample_information_txt = '/Users/alonshaiber/PycharmProjects/MACg/tests/sandbox' \
    #                           '/p214_Bfrag_positive_with_M_GG_gene_coverage_samples_information.txt'
    #
    # main(input_file,additional_layers_file,new_sample_information_txt)
