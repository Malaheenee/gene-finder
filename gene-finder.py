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
    tmp_pat = ''
    tmp_pat_name = ''
    PATT_OPEN = open(pattern_file, 'r')
    while True:
        string = PATT_OPEN.readline()
        if not string:
            PATT_OPEN.close()
            break
        string = string.rstrip()
        if string.startswith('#'):
            tmp_pat_name = string[1:]
        elif string.endswith('|'):
            tmp_pat += string
        else:
            tmp_pat += string
            pattern_dict[tmp_pat_name] = re.compile(tmp_pat, re.IGNORECASE)
            tmp_pat = ''
            tmp_pat_name = ''
    return pattern_dict
