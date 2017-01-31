# -*- coding: utf-8
# pylint: disable=line-too-long
import csv
import numpy as np


def get_data_from_txt_file(file_name):
    data_list = list(csv.reader(open(file_name), delimiter='\t'))
    Nsamples = len(data_list[0]) - 1
    Ngenes = len(data_list) - 1

    sample_name_dictionary = {}
    for i in range(0,Nsamples):
        sample_name_dictionary[i] = data_list[0][i+1]

    gene_callers_id_dictionary = {}
    for i in range(0,Ngenes):
        gene_callers_id_dictionary[i] = data_list[i+1][0]

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


def alternative_algorithm(data, alpha=0.5, beta=1):
    Ns = len(data[0])
    Ngenes = len(data)
    # Initialize list of TS (Taxon Specific genes)
    taxon_specific_genes = range(Ngenes)
    # Initialize list of positive/negative samples
    sample_detection = np.ones(Ns)
    converged = False
    loss = None
    while not converged:
        # mean of coverage of all TS genes in each sample
        mean = np.mean(data[taxon_specific_genes, :], axis=0)  # calculating the mean along the columns

        # determining the detection of the Genome in each sample
        detection_portion = sum(np.abs(np.abs(data[taxon_specific_genes, :]-mean)-3*np.sqrt(np.var(data[
                                                                                                      taxon_specific_genes, :],axis=0))))
        for sample in range(Ns):
            if detection_portion[sample] >= alpha * Ngenes:
                sample_detection[sample] = 1
            else:
                sample_detection[sample] = 0

        # calculate adjusted variance of each gene (adjusted variance is just a name I made-up for this term
        positive_samples = np.nonzero(sample_detection)[0]
        v = np.var(data[:, positive_samples] / mean[positive_samples], axis=1)

        # classifying genes (TS or NTS)
        taxon_specific_genes = []
        for gene_id in range(Ngenes):
            if v[gene_id] <= beta:
                taxon_specific_genes.append(gene_id)

        # calculating the loss function
        number_of_NTS = Ngenes - len(taxon_specific_genes)
        new_loss = beta * number_of_NTS + sum(v[taxon_specific_genes])

        # Check convergence
        if loss is not None:
            if new_loss >= loss:
                converged = True
        loss = new_loss

    return taxon_specific_genes, positive_samples

def get_taxon_specific_candidates(data, positive_samples_list, gene_detection_matrix, gamma=10):
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
        # approximation of the median value (only of detected genes):
        detected_genes = np.nonzero(gene_detection_matrix[:,sample_number])
        median_coverage_index = np.argsort(data[detected_genes,sample_number][0])[len(
            data[detected_genes,sample_number])//2]
        median_coverage = data[median_coverage_index,sample_number]
        sorted_indexes = np.argsort(np.absolute(data[:, sample_number] - median_coverage))
        print('sample number %s, the median value is %s in index number %s' % (sample_number, median_coverage, median_coverage_index))
        print('the sorted indexes: %s' % sorted_indexes)
        var = 0
        mean = median_coverage
        cluster = {'gene_ids': [median_coverage_index], 'gene_coverages': [median_coverage]}
        while not converged and gene_number < Ngenes:
            new_gene_number = sorted_indexes[gene_number]
            new_gene_coverage = data[new_gene_number, sample_number]
            if var == 0:
                cutoff = 0.1 * mean
            else:
                cutoff = max(gamma * np.sqrt(var), 0.1 * mean)
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
                                                                    gene_detection_matrix, eta=0.8):
    # Decide which genes are taxon specific according to a majority vote
    # input:
    #     taxon_specific_candidates_matrix
    # output:
    #     taxon_specific_labels - dictionary with gene_callers_id as keys and values of 'TS' for Taxon-Specific and
    # 'NTS' for Non Taxon-Specific
    Ngenes = len(taxon_specific_candidates_matrix)
    taxon_specific_genes = []
    for gene_number in range(Ngenes):
        print('gene %s is detected in %s samples and TS in %s samples, and overlap in %s samples' % (gene_number,
                sum(gene_detection_matrix[gene_number,:]),sum(taxon_specific_candidates_matrix[gene_number,:]), sum(np.multiply(
                taxon_specific_candidates_matrix[gene_number,:],gene_detection_matrix[gene_number,:]))))
        if sum(np.multiply(taxon_specific_candidates_matrix[gene_number,:],gene_detection_matrix[gene_number,
                                                                           :])) / sum(gene_detection_matrix[
                                                                                      gene_number,:]) > eta:
            taxon_specific_genes.append(gene_number)
    return taxon_specific_genes


def get_accessory_genes(data, taxon_specific_genes, positive_samples, a=3):
    Ngenes = len(data)
    Ns = len(data[0])
    # indices of the sub-matrix containing all Taxon-specific genes only in positive samples:
    # calculating the mean of the taxon specific genes in each positive sample
    gene_detection = np.zeros((Ngenes, Ns), dtype=bool)
    for sample in positive_samples:
        mean_of_TS_coverage = np.mean(data[taxon_specific_genes, sample])
        std_of_TS_coverage = np.std(data[taxon_specific_genes, sample])
        print('sample number %s, mean is %s, std is %s' % (sample, mean_of_TS_coverage, std_of_TS_coverage))
        # if the gene coverage is more than 'a' times the std smaller than average then it is considered to be not
        # detected (or disconnected) from the sample
        # gene_detection[:, sample] = data[:, sample] - mean_of_TS_coverage > -a * std_of_TS_coverage
        gene_detection[:, sample] = np.logical_and(data[:, sample] - mean_of_TS_coverage > -a * std_of_TS_coverage,
                                                   data[:,sample] > 0)
        print(np.sum(gene_detection[:,sample]))

    # array showing in how many samples the gene is detected
    gene_detection_layer = np.sum(gene_detection,axis=1)

    accessory_genes = np.zeros(Ngenes)
    for gene_id in range(Ngenes):
        if gene_detection_layer[gene_id] > 0.9 * len(positive_samples):
            accessory_genes[gene_id] = 0
        else:
            accessory_genes[gene_id] = 1

    return gene_detection_layer, accessory_genes


def get_gene_classes_dictionary(taxon_specific_dictionary, accessory_genes, gene_callers_id_dictionary):
    from gen_mock_data import gene_class_id_dictionary_reverese
    gene_classes_dictionary = {}
    for gene_id in gene_callers_id_dictionary.keys():
        if taxon_specific_dictionary[gene_callers_id_dictionary[gene_id]] == 'TS' and accessory_genes[gene_id] == 0:
            # Taxon specific core
            gene_classes_dictionary[gene_callers_id_dictionary[gene_id]] = gene_class_id_dictionary_reverese[1]
        elif taxon_specific_dictionary[gene_callers_id_dictionary[gene_id]] == 'TS' and accessory_genes[gene_id] == 1:
            # Taxon specific accessory
            gene_classes_dictionary[gene_callers_id_dictionary[gene_id]] = gene_class_id_dictionary_reverese[3]
        elif taxon_specific_dictionary[gene_callers_id_dictionary[gene_id]] == 'NTS' and accessory_genes[gene_id] == 0:
            # Non taxon specific core
            gene_classes_dictionary[gene_callers_id_dictionary[gene_id]] = gene_class_id_dictionary_reverese[4]
        elif taxon_specific_dictionary[gene_callers_id_dictionary[gene_id]] == 'NTS' and accessory_genes[gene_id] == 1:
            # Non taxon specific accessory
            gene_classes_dictionary[gene_callers_id_dictionary[gene_id]] = gene_class_id_dictionary_reverese[5]

    return gene_classes_dictionary


def gen_taxon_specific_dictionary_from_list(taxon_specific_genes,gene_callers_id_dictionary):
    taxon_specific_dictionary = dict(zip(gene_callers_id_dictionary.values(),['NTS'] * len(
        gene_callers_id_dictionary)))
    print(gene_callers_id_dictionary)
    for gene_id in taxon_specific_genes:
        taxon_specific_dictionary[gene_callers_id_dictionary[gene_id]] = 'TS'
    return taxon_specific_dictionary


def save_tabular_to_txt(dictionary, new_txt_output, first_column_title, additional_columns_title,
    old_txt=None):
    if old_txt is None:
        with open(new_txt_output, 'w') as txt_file:
            writer = csv.writer(txt_file, delimiter='\t')
            # writing the title row
            first_row = [first_column_title] + additional_columns_title
            writer.writerow(first_row)
            for key, value in dictionary.items():
                writer.writerow([key, value])
    else:
        with open(old_txt, 'r') as old_file:
            reader = csv.reader(old_file, delimiter='\t')
            with open(new_txt_output, 'w') as txt_file:
                writer = csv.writer(txt_file, delimiter='\t')
                first_row = list(next(reader)) + additional_columns_title
                writer.writerow(first_row)
                for row in reader:
                    print(row + [dictionary[row[0]]])
                    writer.writerow(row + [dictionary[row[0]]])


def save_taxon_specific_labels_to_txt(taxon_specific_dictionary, txt_output, additional_layers_txt=None):
    save_tabular_to_txt(dictionary=taxon_specific_dictionary,new_txt_output=txt_output,
                        first_column_title='gene_callers_id',additional_columns_title=['taxon_specific_label'],
                        old_txt=additional_layers_txt)


def get_samples_detection_dictionay(positive_samples_list, sample_name_dictionary):
    samples_detection_dictionay = {}
    for sample in sample_name_dictionary:
        if sample in positive_samples_list:
            samples_detection_dictionay[sample_name_dictionary[sample]] = 'P'
        else:
            samples_detection_dictionay[sample_name_dictionary[sample]] = 'N'
    return samples_detection_dictionay

def get_negative_samples(positive_samples_list,sample_name_dictionary):
    negative_samples = list(set(sample_name_dictionary.keys()) - set(positive_samples_list))
    return negative_samples


def save_sample_detection_information_to_sample_information_file(positive_samples_list, sample_name_dictionary,
                                                                 txt_output, sample_information_txt=None):
    samples_detection_dictionay = get_samples_detection_dictionay(positive_samples_list, sample_name_dictionary)
    save_tabular_to_txt(dictionary=samples_detection_dictionay, new_txt_output=txt_output,
                        first_column_title='samples', additional_columns_title=['detection_of_genome'],
                                                                                old_txt=sample_information_txt)


def tests():
    # input_name = 'test_200'
    input_name = 'p214_Bfrag_positive_with_M_GG_gene_coverage'
    input_data = '/Users/alonshaiber/PycharmProjects/MACg/tests/sandbox/' + input_name + '.txt'
    data, sample_name_dictionary, gene_callers_id_dictionary = get_data_from_txt_file(input_data)
    print(data.shape)
    # get the positive samples
    positive_samples_list , gene_detection_matrix = get_positive_samples(data)
    print('The number of positive samples is %s'%len(positive_samples_list))
    negative_samples = get_negative_samples(positive_samples_list,sample_name_dictionary)
    print('The following samples are negative: %s' % ([sample_name_dictionary[key] for key in negative_samples]))
    # testing get_taxon_specific_candidates
    taxon_specific_candidates_matrix = get_taxon_specific_candidates(data, positive_samples_list, gene_detection_matrix)
    print(len(taxon_specific_candidates_matrix))
    print(len(np.nonzero(taxon_specific_candidates_matrix)[0]))

    # testing get_taxon_specific_labels_from_taxon_specific_candidates_matrix
    taxon_specific_genes = get_taxon_specific_labels_from_taxon_specific_candidates_matrix(
        taxon_specific_candidates_matrix, gene_detection_matrix, eta=0.8)
    print('The number of taxon specific genes is: %s ' % len(taxon_specific_genes))

    # save the results
    taxon_specific_dictionary = gen_taxon_specific_dictionary_from_list(taxon_specific_genes,
                                                                        gene_callers_id_dictionary)
    txt_output = '/Users/alonshaiber/PycharmProjects/MACg/tests/sandbox/' + input_name + '_taxon_specific_genes.txt'
    additional_layers_txt = None # '/Users/alonshaiber/PycharmProjects/MACg/tests/sandbox/test_additional_layers.txt'
    save_taxon_specific_labels_to_txt(taxon_specific_dictionary, txt_output, additional_layers_txt)

    sample_information_txt = '/Users/alonshaiber/PycharmProjects/MACg/tests/sandbox/' + input_name + \
                             '_sample_information.txt'
    old_sample_information =  None
    save_sample_detection_information_to_sample_information_file(positive_samples_list, sample_name_dictionary,
                                                                 sample_information_txt,old_sample_information)

    # test get_accessory_genes(data, taxon_specific_genes, positive_samples)
    gene_detection_layer, accessory_genes = get_accessory_genes(data, taxon_specific_genes, positive_samples_list)

    # get_gene_classes
    gene_classes_dictionary = get_gene_classes_dictionary(taxon_specific_dictionary, accessory_genes, gene_callers_id_dictionary)
    txt_output = '/Users/alonshaiber/PycharmProjects/MACg/tests/sandbox/' + input_name + '_additional_layers.txt'
    additional_layers_txt = '/Users/alonshaiber/PycharmProjects/MACg/tests/sandbox/' + input_name + '_taxon_specific_genes.txt'
    save_tabular_to_txt(gene_classes_dictionary, txt_output, 'gene_callers_id', ['gene_class'], additional_layers_txt)

    # save_

    # # running the alternative algorithm
    # taxon_specific_genes_alt, positive_samples_list_alt = alternative_algorithm(data, alpha=0.5, beta=1)
    #
    # # save the results
    # taxon_specific_dictionary_alt = gen_taxon_specific_dictionary_from_list(taxon_specific_genes_alt,
    #                                                                     gene_callers_id_dictionary)
    # print(taxon_specific_dictionary_alt)
    # txt_output = '/Users/alonshaiber/PycharmProjects/MACg/tests/sandbox/' + input_name + '_taxon_specific_genes_alt.txt'
    # additional_layers_txt = None
    # save_taxon_specific_labels_to_txt(taxon_specific_dictionary_alt, txt_output, additional_layers_txt)
    #
    # sample_information_txt = '/Users/alonshaiber/PycharmProjects/MACg/tests/sandbox/' + input_name + \
    #                          '_sample_information_alt.txt'
    # old_sample_information =  None
    # save_sample_detection_information_to_sample_information_file(positive_samples_list_alt, sample_name_dictionary,
    #                                                              sample_information_txt,old_sample_information)


if __name__ == '__main__':
    tests()