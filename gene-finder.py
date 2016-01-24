#! /usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division, print_function
import os
import sys
import re
import argparse

# Arguments parser
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
        default='gene_dict.txt',
        metavar='FILE',
        help='The file with genes (default "gene_dict.txt")')

    parser.add_argument(
        '-a',
        type=str,
        required=False,
        nargs='?',
        default='abstract.txt',
        metavar='FILE',
        help='The file with abstracts (default "abstract.txt")')

    parser.add_argument(
        '-p',
        type=str,
        required=False,
        nargs='?',
        default='patterns2.txt',
        metavar='FILE',
        help='The file with search patterns (default "patterns.txt")')

    parser.add_argument(
        '-k',
        type=str,
        required=False,
        nargs='?',
        default='kegg_dict.txt',
        metavar='FILE',
        help='The file with genes in KEGG pathway (default "kegg_dict.txt")')

    parser.add_argument(
        '-m',
        type=str,
        required=False,
        nargs='?',
        default='host_dict.txt',
        metavar='FILE',
        help='The file with genes-hosts for miRNAs (default "host_dict.txt")')

    return vars(parser.parse_args())

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
def gene_prepare(gene_dict_file, search_bit):
    gene_dict = {}
    for string in open(gene_dict_file):
        if string.startswith('#'):
            continue
        string = string.rstrip().split(';')
        gname = string[0]
        if gname not in gene_dict:
            gene_dict[gname] = []
        for tmp in string[1:]:
            if tmp == '':
                continue
            gene_dict[gname] += re.split(r',\s+', tmp)
        # Make dictionary for search if set search bit
        if search_bit != 0:
            tmp_list = []
            gene_dict[gname] += [gname]
            for tmp in gene_dict[gname]:
                tmp_list.append(re.compile(''.join(['\se?', tmp, '(\s|\(|\-)?']), re.IGNORECASE))
            gene_dict[gname] = tmp_list
    return gene_dict

# Run program if it call here
if __name__ == '__main__':

    arg_files = arg_receive()

    for key in arg_files.keys():
        if os.path.exists(arg_files[key]):
            arg_files[key] = os.path.abspath(arg_files[key])
        else:
            arg_files[key] = None

    pattern_dict = pattern_prepare(arg_files['p'])

    gene_dict = gene_prepare(arg_files['g'], 1)

    if arg_files['k']:
        kegg_dict = gene_prepare(arg_files['k'], 0)
    else:
        kegg_dict = None
    if arg_files['m']:
        host_dict = gene_prepare(arg_files['m'], 0)
    else:
        host_dict = None

    RES_OPEN = open('result.txt', 'w')
    RES_OPEN.write(';'.join(['PMID', 'Gene', 'Exact match',\
                              'miR', 'Features 1', 'Features 2', 'Features 3',\
                              'KEGG Pathway', 'Host miR', 'Journal']))

    ABS_OPEN = open(arg_files['a'], 'r')
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
                        match = gene.search(abstract_text, re.MULTILINE)
                        if match:
                            match = abstract_text[match.start(0)-20:match.end(0)+20]
                            RES_OPEN.write('\n')
                            RES_OPEN.write(';'.join([abstract_pmid, key, match, '']))

                            result = dict.fromkeys(pattern_dict.keys())
                            for pattern in sorted(pattern_dict.keys()):
                                result[pattern] = []
                                for match in pattern_dict[pattern].finditer(abstract_text, re.MULTILINE):
                                    match = str(match.group(0))
                                    if match not in result[pattern]:
                                        result[pattern].append(match)
                                RES_OPEN.write(';'.join([', '.join(result[pattern]), '']))
                            if kegg_dict and key in kegg_dict:
                                RES_OPEN.write(', '.join(kegg_dict[key]))
                            RES_OPEN.write(';')
                            if host_dict and key in host_dict:
                                RES_OPEN.write(', '.join(host_dict[key]))
                            RES_OPEN.write(';')
                            RES_OPEN.write(abstract_journ)
            abstract = ''
    RES_OPEN.close()
