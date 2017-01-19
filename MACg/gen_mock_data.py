# -*- coding: utf-8
# pylint: disable=line-too-long
import csv
import random
import numpy
# dictionary to translate from class name to class_caller_id
gene_class_id_dictionary = dict(NaN=0, STC=1, MTC=2, TA=3, NTC=4, NTA=5)


def read_genome_classes_table(genome_classes_table_filename):
    """ read  the genome class table
     input: genome_classes_table_filename"""
    with open(genome_classes_table_filename) as f:
        reader = csv.DictReader(f, delimiter='\t')
        genome_classes_table = {}
        for row in reader:
            genome_classes_table[row['class_name']] = {}
            for key in row.keys():
                if key != 'class_name':
                    genome_classes_table[row['class_name']][key] = row[key]
    return genome_classes_table


def gen_group_class_table(genome_classes_table):
    """ generates a table in which for each group the class number is specified """
    group_class_table = {}
    group_caller_id = 0
    for gene_class in genome_classes_table.keys():
        if genome_classes_table[gene_class]['number_of_groups_per_class'] is '':
            group_class_table[group_caller_id] = {}
            group_class_table[group_caller_id]['class_id'] = gene_class_id_dictionary[gene_class]
            group_class_table[group_caller_id]['number_of_genes'] = int(genome_classes_table[gene_class][
                'number_of_genes'])
            group_caller_id += 1
        else:
            number_of_groups = int(genome_classes_table[gene_class]['number_of_groups_per_class'])
            if int(genome_classes_table[gene_class]['number_of_genes']) < number_of_groups:
                raise BaseException('number of groups in a class must')
            # the maximum number of genes in a group is set so every group in the class would have at least one member
            max_number_of_genes = int(genome_classes_table[gene_class]['number_of_genes']) - number_of_groups
            for i in range(number_of_groups):
                group_class_table[group_caller_id] = {}
                if i == number_of_groups - 1:
                    number_of_genes = max(max_number_of_genes + number_of_groups, 1)
                else:
                    number_of_genes = random.randint(1, max(max_number_of_genes, 1))
                max_number_of_genes -= number_of_genes
                group_class_table[group_caller_id]['number_of_genes'] = number_of_genes
                group_class_table[group_caller_id]['class_id'] = gene_class_id_dictionary[gene_class]
                group_caller_id += 1
    return group_class_table


def gen_gene_group_table_from_genome_classes_table(group_class_table):
    """ Generate a gene group table from a class table. In this table, for each gene, the group_caller_id is
    input: genome_classes_table
    output: gene_group_table
    specified and the class_id is specified

    """
    gene_group_table = {}
    number_of_genes = 0
    for group in group_class_table.keys():
        genes_callers_ids = range(number_of_genes, number_of_genes + group_class_table[group]['number_of_genes'])
        genes_group_and_class_ids = [{'group_callers_id': group, 'class_id': group_class_table[group]['class_id']}] \
                                    * group_class_table[group]['number_of_genes']
        gene_group_table.update(dict(zip(genes_callers_ids, genes_group_and_class_ids)))
        number_of_genes += group_class_table[group]['number_of_genes']
    print('gene_group_table was created with %s genes' % number_of_genes)
    return gene_group_table


def save_dict_as_tab_delimited_txt(d, txt_output, first_column_title=''):
    with open(txt_output, 'w') as txt_file:
        writer = csv.writer(txt_file, delimiter='\t')
        # writing the title row
        writer.writerow([first_column_title] + list(d[next(iter(d.keys()))].keys()))
        for key, value in d.items():
            writer.writerow([key] + list(value.values()))


def save_gene_group_table_as_txt(gene_group_table, txt_output):
    save_dict_as_tab_delimited_txt(gene_group_table, txt_output, 'gene_callers_id')


def gen_abundance_of_groups_in_sample(group_class_table):
    """ takes a group_class_table and generates a random abundance number for each class
     input: group_class_table
     output: group_abundance_table (a dictionary)
     The sum of all abundances equals 1
    """
    number_of_groups = len(group_class_table)
    max_ration_between_max_and_min_abundant_group = 1000000000000
    group_abundance_table = dict(zip(list(group_class_table.keys()), random.sample(range(
        max_ration_between_max_and_min_abundant_group), number_of_groups)))

    return group_abundance_table


def gen_gene_abundance_from_group_abundance_table(group_abundance_table, gene_group_table):
    """ Takes a group abundance table and generate a stochastic abundance for each gene of every group
    input: group_abundance_table, gene_group_table
    output: gene_abundance_table
    """
    gene_abundance_list = []
    for gene_callers_id in gene_group_table.keys():
        group_callers_id = gene_group_table[gene_callers_id]['group_callers_id']
        mu = group_abundance_table[group_callers_id]
        sigma = 0.5 * mu
        gene_abundance = max(0, numpy.around(numpy.random.normal(loc=mu, scale=sigma), decimals=0))
        gene_abundance_list.append(gene_abundance)
    gene_abundance_table = dict(zip(gene_group_table.keys(), gene_abundance_list))
    return gene_abundance_table


def gen_mock_merged_coverage_table(group_class_table, gene_group_table, number_of_samples=16):
    """ Takes a file name (of a genome_classes_table), and number of samples and generates the coverage table
    input: genome_classes_table_filename, number of required samples
    output: mock_merged_coverage_table - table with rows as abundance of a gene accross metagenomes
    """
    mock_merged_coverage_table = {}
    for sample_id in range(number_of_samples):
        group_abundance_table=gen_abundance_of_groups_in_sample(group_class_table)
        mock_merged_coverage_table['sample_' + str(sample_id)] = gen_gene_abundance_from_group_abundance_table(
            group_abundance_table, gene_group_table)
    return mock_merged_coverage_table


def transpose_txt(input_file, output_file):
    # copied from: https://gist.github.com/mikeboers/4997162
    # Read all Tab-delimited rows from stdin.
    import itertools
    with open(input_file,'r') as i:
        tsv_reader = csv.reader(i, delimiter='\t')
        all_data = list(tsv_reader)

        # Transpose it.
        all_data = list(itertools.zip_longest(*all_data, fillvalue=''))

        # Write it back out.
        with open(output_file,'w') as o:
            tsv_writer = csv.writer(o, delimiter='\t')
            for row in all_data:
                tsv_writer.writerow(row)

def save_mock_merged_coverage_table_as_txt(mock_merged_coverage_table, output):
    save_dict_as_tab_delimited_txt(mock_merged_coverage_table, output, first_column_title='gene_callers_id')


def tests():
    genome_classes_table_filename = '../tests/sandbox/example_genome_classes_table.txt'
    number_of_samples = 100

    # test read_genome_classes_table
    genome_classes_table = read_genome_classes_table(genome_classes_table_filename)
    print('genome_classes_table is: %s' % genome_classes_table)

    # test gen_group_class_table
    group_class_table = gen_group_class_table(genome_classes_table)
    print('group_class_table is: %s' % group_class_table)

    # test gen_gene_group_table_from_genome_classes_table
    gene_group_table = gen_gene_group_table_from_genome_classes_table(group_class_table)
    print('gene_group_table is: %s' % gene_group_table)

    # test save_gene_group_table_as_txt
    gene_group_table_txt = '../tests/sandbox/test_gene_group_table.txt'
    save_gene_group_table_as_txt(gene_group_table, gene_group_table_txt)

    # test gen_abundance_of_groups_in_sample
    group_abundance_table = gen_abundance_of_groups_in_sample(group_class_table)
    print('group_abundance_table is: %s' % group_abundance_table)

    # test gen_gene_abundance_from_group_abundance_table
    gene_abundance_table=gen_gene_abundance_from_group_abundance_table(group_abundance_table, gene_group_table)
    print('gene_abundance_table is: %s' % gene_abundance_table)

    # test gen_mock_merged_coverage_table
    mock_merged_coverage_table = gen_mock_merged_coverage_table(group_class_table, gene_group_table,
                                                                number_of_samples)
    print('mock_merged_coverage_table is: %s' % mock_merged_coverage_table)

    # test save_mock_merged_coverage_table_as_txt
    test_output='../tests/sandbox/test_output.txt'
    save_mock_merged_coverage_table_as_txt(mock_merged_coverage_table, test_output)
    test_output_transposed = '../tests/sandbox/test_output_transposed.txt'
    transpose_txt(test_output, test_output_transposed)


# Testing all methods
tests()

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-g', '--genome_classes_table', metavar='FILE', dest='genome_classes_table_filename',
                        help='input genome class table')
    parser.add_argument('-N', '--number_of_samples', metavar='INT', dest='number_of_samples', type=int, default=16,
                        help='Number of samples in the mock data')
    parser.add_argument('-o', '--out', metavar='FILE', dest='output', default='mock_coverage_data.txt',
                        help='Output file')
    parser.add_argument('--test', action='store_true', dest='test', help='test that everything is ok and exit')
    args = parser.parse_args()

    genome_classes_table_filename = args.genome_classes_table_filename
    number_of_samples = args.number_of_samples
    output = args.output
    # merged_coverage_table = gen_mock_merged_coverage_table(genome_classes_table_filename, number_of_samples)