#! /usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division, print_function
import os
import sys
import re
import argparse
from Bio import Medline

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
        default='patterns.txt',
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

    parser.add_argument(
        '-o',
        type=str,
        required=False,
        nargs='?',
        default='result.txt',
        metavar='FILE',
        help='The file with results (default "result.txt")')

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
        gene_dict[gname] += re.split(r',\s+', string[1])
        # Make dictionary for search if set search bit
        if search_bit != 0:
            tmp_list = []
            gene_dict[gname] += [gname]
            for tmp in gene_dict[gname]:
                if str(type(tmp)).find('_sre.SRE_Pattern') != -1 or \
                tmp == '' or re.match(r"""\b\d?\-?[a-z]\d?\-?\b|
                                        ai?(d|r|m|s|my)?|
                                        ar(c|g)h?|
                                        acts?|
                                        ag(o|e)|
                                        ap(ex|ril)|
                                        adipose|
                                        b(e|ig|ut)?|
                                        b(ank|ad|ase)|
                                        c(an|m|o|s?t)|
                                        coil|
                                        cl(i|am)p|
                                        d(amage|iabetes)|
                                        (d|t)o|
                                        e(asy|x|cho|ndo?|stimate|yes)|
                                        ex(is|ac)t|
                                        fa(r|t|c(e|t))|
                                        f(ind|lap|or)|
                                        g(ap|as|enesis|et|reat|ut|)|
                                        h(e|(i|a)s)|
                                        h(int|eld|ole)|
                                        i(d|f|n|t)|
                                        impact|
                                        jaundice|
                                        k(d|ca|it)|
                                        l(ag|obe|ow|ed)|
                                        li(ttle|ght)|
                                        la(rge|st)|
                                        I(I|V)I?|
                                        m(a|e)n|
                                        mi((n(ute)?)|nor)|
                                        m(ass|g|l|m|y|et|cm|u)|
                                        mal(-|e)|
                                        n(e|g|m|o)t?|
                                        o(bese|f|r|n|d|s|ut)|
                                        per|
                                        pi(lot|gs)|
                                        (pa|re)(d|st)|
                                        ra(m|n|w)|
                                        r(im|ro)|
                                        s(k|p)in|
                                        s(ac|alt|ex|he|i|kip|lim|o|oft|pastic|tar)|
                                        se(al|t)|
                                        st(o|e)p|
                                        simple|
                                        spasm|
                                        ta(ctile|sk)|
                                        t(i|a)p?|
                                        t(en|ube|his)|
                                        v(ia|s)|
                                        w(e|(a?s)|t)|
                                        white|
                                        wave|
                                        wi(re|sh)|u(p|s)|
                                        urea|
                                        zeta""", tmp, re.IGNORECASE | re.VERBOSE):
                    continue
                tmp_list.append(re.compile(''.join(['\se?', re.escape(tmp), '(\s|\(|\-)']), re.IGNORECASE))
            gene_dict[gname] = tmp_list
    return gene_dict

# Run program if it call here
if __name__ == '__main__':

    arg_files = arg_receive()

    for key in arg_files.keys():
        if os.path.exists(arg_files[key]) or key == 'o':
            arg_files[key] = os.path.abspath(arg_files[key])
        else:
            arg_files[key] = None

    print('Prepare patterns... ', end='')
    pattern_dict = pattern_prepare(arg_files['p'])
    print(len(pattern_dict), 'patterns.', sep=' ')

    print('Prepare gene dictionary... ', end='')
    gene_dict = gene_prepare(arg_files['g'], 1)
    print(len(gene_dict), 'genes.', sep=' ')

    if arg_files['k']:
        print('Prepare KEGG dictionary... ', end='')
        kegg_dict = gene_prepare(arg_files['k'], 0)
        print(len(kegg_dict), 'genes.', sep=' ')
    else:
        kegg_dict = None
    if arg_files['m']:
        print('Prepare miRNA-hosts dictionary... ', end='')
        host_dict = gene_prepare(arg_files['m'], 0)
        print(len(host_dict), 'genes.', sep=' ')
    else:
        host_dict = None

    RES_OPEN = open(arg_files['o'], 'w')
    RES_OPEN.write(';'.join(['PMID', 'Gene', 'Exact match', 'Wide string', \
                              'miR', 'Features 1', 'Features 2', 'Features 3',\
                              'KEGG Pathway', 'Host miR', 'Journal']))

    print('Searching...', end='')
    ABS_OPEN = open(arg_files['a'], 'r')
    abstract_size = round(os.path.getsize(arg_files['a'])/1024/1024)
    all_abstracts = Medline.parse(ABS_OPEN)
    for abstract in all_abstracts:
        if 'AB' in abstract:
            abstract_text = abstract['AB']
            abstract_pmid = 'Unknown'
            abstract_journ = 'Unknown'
            if 'PMID' in abstract:
                abstract_pmid = abstract['PMID']
            if 'SO' in abstract:
                abstract_journ = abstract['SO']
            print('\rSearching in PMID', abstract_pmid, '(', \
                   round(ABS_OPEN.tell()/1024/1024, 1), 'MB read from', abstract_size, \
                   'MB )', sep=' ', end='')
            abstract_text = re.sub(r'\r\n|\n|\s+|;', ' ', abstract_text)
            for key in gene_dict.keys():
                for gene in gene_dict[key]:
                    match = gene.search(abstract_text, re.MULTILINE)
                    if match:
                        RES_OPEN.write('\n')
                        RES_OPEN.write(';'.join([abstract_pmid, key, \
                                        match.group(0), \
                                        abstract_text[match.start(0)-20:match.end(0)+20], '']))

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
    print('\n')
    ABS_OPEN.close()
    RES_OPEN.close()
