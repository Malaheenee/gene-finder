#! /usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division, print_function
import os
import sys
import re
import argparse

def arg_receive():
    "This script ..."
    parser = argparse.ArgumentParser(
        description = 'For use of this script it is necessary have \
        "abstract.txt" file with abstacts (text format of PubMed \
        (http://www.ncbi.nlm.nih.gov/pubmed/)). There are must be 2 \
        empty strings between abstracts. The result of searching for \
        contents of the following data: PMID (publication medicine \
        identification number) of the article,\ miRNA name, disease \
        localization (organ or tissue), keywords (methods, change fold,\
        cellular processes, functions, animal species, types of cells,\
        biological liquids, etc.) and genes. The results will be \
        selected according to each found gene in the abstract.',
        epilog = 'Send the bugs to devolia18@mail.ru',
        argument_default=argparse.SUPPRESS)
    
    parser.add_argument(
        '-g',
        type=str,
        required=False,
        nargs='?',
        metavar='The "gene_dict.txt" file',
        help='The filename with genes (default "gene_dict.txt")')

    parser.add_argument(
        '-a',
        type=str,
        required=False,
        nargs='?',
        metavar='The "abstract.txt" file',
        help='The file name with abstracts (default "abstract.txt")')

    parser.add_argument(
        '-p',
        type=str,
        required=False,
        nargs='?',
        metavar='The "patterns.txt" file',
        help='The file name with search patterns (default "patterns.txt")')
        
    args = parser.parse_args()
    
    return args.g, args.a, args.p

# Prepare search patterns
def pattern_prepare(pattern_file):
    pattern_dict = {}
    PATT_OPEN = open(pattern_file, 'r')
    while True:
        string = PATT_OPEN.readline()
        if not string:
            PATT_OPEN.close()
            break
        string = string.rstrip()
        if string.startswith('#'):
            tmp_pat_name = string[1:]
            tmp_pat = ''
        elif string.endswith('|'):
            tmp_pat += string
        else:
            tmp_pat += string
            pattern_dict[tmp_pat_name] = re.compile(tmp_pat, re.IGNORECASE)
    return pattern_dict

# Prepare gene dictionary
def gene_prepare(gene_dict_file):
    gene_dict = {}
    for string in open(gene_dict_file):
        if string.startswith('#'):
            continue
        string = string.rstrip().split(';')
        gene_dict[string[0]] = [string[0], string[0].upper(), \
                                string[0].lower(), string[0].lower().capitalize()]
        if not string[0][-1].isdigit():
            gene_dict[string[0]] += [string[0][0:-1].lower() + string[0][-1].upper()]
        for tmp in string[1:]:
            if tmp == '':
                continue
            gene_dict[string[0]] += re.split(r',\s+', tmp)
    return gene_dict

# Run program if it call here
if __name__ == '__main__':
    gene_dict_file, abstract_file, pattern_file = arg_receive()

    pattern_dict = pattern_prepare(os.path.abspath(pattern_file))

    gene_dict = gene_prepare(os.path.abspath(gene_dict_file))

    RES_OPEN = open('result.txt', 'w')
    RES_OPEN.write(', '.join(['PMID','Journal','Gene','Exact match',\
                              'miR','Features 1','Features 2','Features 3']))

    ABS_OPEN = open(abstract_file, 'r')
    abstract = ''
    while True:
        string = ABS_OPEN.readline()
        if not string:
            ABS_OPEN.close()
            break
        if not string.startswith('PMID'):
            abstract += string
        else:
            abstract_pmid = re.split(r'\s+', string.rstrip())[1]
            abstract = re.split(r'\r\n\r\n', abstract)
            if len(abstract) >= 7:
                abstract_journ = re.sub(r'^\d+. ', '', abstract[1])
                abstract_text = re.sub(r'\r\n', ' ', abstract[-2])

                for key in gene_dict.keys():
                    for gene in gene_dict[key]:
                        if abstract_text.find(gene) != -1:
                            RES_OPEN.write('\n')
                            RES_OPEN.write(abstract_pmid)
                            RES_OPEN.write('\t')
                            RES_OPEN.write(abstract_journ)
                            RES_OPEN.write('\t')
                            RES_OPEN.write(key)
                            RES_OPEN.write('\t')
                            RES_OPEN.write(gene)
                            RES_OPEN.write('\t')

                            result = dict.fromkeys(pattern_dict.keys())
                            for pattern in sorted(pattern_dict.keys()):
                                result[pattern] = []
                                for match in pattern_dict[pattern].finditer(abstract_text, re.MULTILINE):
                                    match = str(match.group(0))
                                    if match not in result[pattern]:
                                        result[pattern].append(match)
                                RES_OPEN.write(", ".join(result[pattern]))
                                RES_OPEN.write('\t')
            abstract = ''
    RES_OPEN.close()


