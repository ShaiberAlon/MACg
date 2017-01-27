# -*- coding: utf-8
# pylint: disable=line-too-long
import csv
import random
import numpy as np


def get_data_from_txt_file(file_name):
    data_list = list(csv.reader(open(file_name), delimiter='\t'))
    Nsamples = len(data_list[0]) - 1
    Ngenes = len(data_list) - 1

    sample_name_dictionary = {}
    for i in range(1,Nsamples):
        sample_name_dictionary[i-1] = data_list[0][i]

    gene_callers_id_dictionary = {}
    for i in range(1,Ngenes):
        gene_callers_id_dictionary[i-1] = data_list[i][0]

    data = np.loadtxt(file_name, delimiter='\t', skiprows=1, usecols=range(1, Nsamples + 1))
    return(data, sample_name_dictionary, gene_callers_id_dictionary)



def get_positive_samples(data,alpha=0.05,beta=0.5):
    # Identify samples that contain your organism of interest (from here on: positive samples):
    # input:
    #     data - data matrix
    #     alpha - cutoff for gene detection (portion of median value)
    #     beta - cutoff for positive sample (portion of detected genes)
    # output:
    #     positive_samples_list - list of positive samples

    gene_detection_matrix = np.zeros_like(data)
    Ngenes = len(data)
    positive_samples_list = []
    for sample_number in range(len(data[0])):
        # getting the median of the non-zero coverage values for
        median_value = np.median(data[np.nonzero(data[:,sample_number]),sample_number])
        gene_detection_matrix[:,sample_number] = data[:,sample_number] > alpha * median_value
        if sum(gene_detection_matrix[:,sample_number]) > beta * Ngenes:
            positive_samples_list.append(sample_number)

    return positive_samples_list, gene_detection_matrix


def get_taxon_specific_candidates(data, positive_samples_list, gamma=10):
    # Find the taxon specific candidate genes in each sample
    # input:
    #     data - data matrix
    #     positive_samples_list - list of positive samples
    # output:
    #     taxon_specific_candidates_matrix - a matrix in which for each gene there are 1's in samples in which it is
    #       identified as a taxon specific candidate and 0 otherwise

    Ngenes = len(data)
    taxon_specific_candidates_matrix = np.zeros_like(data, dtype=bool)

    for sample_number in positive_samples_list:
        converged = False
        gene_number = 1
        # approximation of the median value:
        median_coverage_index = np.argsort(data[:,sample_number])[len(data[:,sample_number])//2]
        median_coverage = data[median_coverage_index,sample_number]
        sorted_indexes = np.argsort(np.absolute(data[:, sample_number] - median_coverage))
        print('sample number %s, the median value is %s in index number %s' % (sample_number, median_coverage, median_coverage_index))
        print(sorted_indexes)
        var = 0
        mean = median_coverage
        cluster = {'gene_ids': [median_coverage_index], 'gene_coverages': [median_coverage]}
        while not converged and gene_number < Ngenes:
            new_gene_number = sorted_indexes[gene_number]
            new_gene_coverage = data[new_gene_number, sample_number]
            if var == 0:
                cutoff = 0.1 * mean
            else:
                cutoff = gamma * np.sqrt(var)
            if abs(mean - new_gene_coverage) > cutoff:
                converged = True
            else:
                cluster['gene_ids'].append(new_gene_number)
                cluster['gene_coverages'].append(new_gene_coverage)
                # updating the new mean
                new_mean = (mean * gene_number + new_gene_coverage) / (gene_number + 1)
                # updating the new variance
                var = (gene_number * var + (new_gene_coverage - new_mean) * (new_gene_coverage - mean)) / (gene_number + 1)
                mean = new_mean
            gene_number += 1
        # setting the value of the genes in the cluster as 1
        taxon_specific_candidates_matrix[cluster['gene_ids'], sample_number] = True

    return taxon_specific_candidates_matrix


def get_taxon_specific_labels_from_taxon_specific_candidates_matrix(taxon_specific_candidates_matrix,
                                                                    gene_detection_matrix, eta=0.1):
    # Decide which genes are taxon specific according to a majority vote
    # input:
    #     taxon_specific_candidates_matrix
    # output:
    #     taxon_specific_labels - dictionary with gene_callers_id as keys and values of 'TS' for Taxon-Specific and
    # 'NTS' for Non Taxon-Specific
    Ngenes = len(taxon_specific_candidates_matrix)
    taxon_specific_genes = []
    for gene_number in range(Ngenes):
        if sum(np.multiply(taxon_specific_candidates_matrix[gene_number,:],gene_detection_matrix[gene_number,
                                                                           :])) / sum(gene_detection_matrix[
                                                                                      gene_number,:]) > eta:
            taxon_specific_genes.append(gene_number)
    return taxon_specific_genes


def save_taxon_specific_labels_to_txt(taxon_specific_genes, txt_output):
    with open(txt_output, 'w') as txt_file:
        writer = csv.writer(txt_file, delimiter='\t')
        # writing the title row
        first_column_title = ['gene_callers_id', 'taxon_specific_label']
        writer.writerow(first_column_title)
        for key, value in taxon_specific_genes.items():
            writer.writerow([key] + list(value.values()))


def tests():
    input_data = '/Users/alonshaiber/PycharmProjects/MACg/tests/sandbox/p214_Bfrag_positive_with_M_GG_gene_coverage.txt'
    data, sample_name_dictionary, gene_callers_id_dictionary = get_data_from_txt_file(input_data)
    print(data[0:4,0])
    positive_samples_list , gene_detection_matrix = get_positive_samples(data)
    print('The number of positive samples is %s'%len(positive_samples_list))
    # testing get_taxon_specific_candidates
    taxon_specific_candidates_matrix = get_taxon_specific_candidates(data, positive_samples_list)
    print(len(taxon_specific_candidates_matrix))
    print(len(np.nonzero(taxon_specific_candidates_matrix)[0]))

    # testing get_taxon_specific_labels_from_taxon_specific_candidates_matrix
    taxon_specific_genes = get_taxon_specific_labels_from_taxon_specific_candidates_matrix(
        taxon_specific_candidates_matrix, gene_detection_matrix, eta=0.8)
    print(taxon_specific_genes)
    print(len(taxon_specific_genes))
if __name__ == '__main__':
    tests()