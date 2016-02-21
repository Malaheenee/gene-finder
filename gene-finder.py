#! /usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division, print_function
import os
import sys
import re
import argparse
from Bio import Medline
from multiprocessing import cpu_count, Queue, Process, TimeoutError

# Arguments parser
def arg_receive():
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

    parser.add_argument(
        '-n',
        required=False,
        default=False,
        action='store_true',
        help='Search without gene synonyms (default "False")')

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
def gene_prepare(gene_dict_file, search_bit=0, nosyn_bit=False):
    gene_dict = {}

    if search_bit != 0:
        bad_words = re.compile(r"""
        \b\d?\-?[a-z]\d?\-?\b|
        \b(a((cts?|dipose)|g(o|e)|i?(d|r|m|s)|p(ex|ril)|r(c|g)h?)|
        b(ank|ad|ase|e|ig|ut)|
        c((a|o)(n|m|o|s?t)?\d?|l(i|am)p|oil)|
        d(amage|iabetes|o)|
        e(asy|x\-?|cho|ndo?|stimate|yes|x(is|ac)t)|
        f(a(r|t|c(e|t))|ind|lap|or)|
        g(ap|as|enesis|et|reat|ut|)|
        h(e|(i|a)s|int|(e|o)ld|oles?)|
        i(d|f|mpact|n|t|(I|V)I?)|
        jaundice|
        k(d|ca|it)|
        l(ag|obe|ow|ed|arge|(a|i)(st|ttle|ght))|
        m(al(\-|e)|(a|e)n|i(n(ute)?|nor)|ass|g|l|m|y|et|cm|u)|
        n(e|g|m|o)t?|
        o(bese|f|r|n|d|s|ut)|
        p(a(d|st)|i(lot|gs?)|er)|
        r(a(m|n|w)|e(d|st)|im|ro)|
        s((k|p)in|e(al|t)|t(o|e)p|ac|alt|ex|he|i|imple|kip|lim|o|oft|pasm|pastic|tar)|
        t(a(ctile|sk)|(i|a)p|en|o|ube|his)|
        u(p|s|rea)|
        v(ia|s)|
        w(e|(a?s)|t|h(at|ite|o)|ave|i(re|sh))|
        zeta)\b
        """, re.IGNORECASE | re.VERBOSE)

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
            if nosyn_bit:
                gene_dict[gname] = [gname]
            else:
                gene_dict[gname] += [gname]
            for tmp in gene_dict[gname]:
                if str(type(tmp)).find('_sre.SRE_Pattern') != -1 or \
                tmp == '' or bad_words.match(tmp):
                    continue
                tmp_list.append(re.compile(''.join(['\se?', re.escape(tmp), '\W']), re.IGNORECASE))
            gene_dict[gname] = tmp_list

    # Prepare gene dictionary for multirpocessing
    cpus = cpu_count()
    if cpus >= 1 and search_bit != 0 :
        keys = sorted(gene_dict.keys())
        chunk_len = len(keys) // cpus
        rest_count = len(keys) % cpus

        for x in range(cpus):
            chunk = keys[:chunk_len]
            keys = keys[chunk_len:]
            if rest_count and keys:
                chunk.append(keys.pop(0))
                rest_count -= 1
            gene_dict[x] = {}
            for y in chunk:
                gene_dict[x][y] = gene_dict[y]
                del(gene_dict[y])

    return gene_dict

# Search in abstract
def abs_search(gene_dict, pattern_dict, abstract_file, out_queue):
    result_dict = {}
    ABS_OPEN = open(abstract_file, 'r')
    all_abstracts = Medline.parse(ABS_OPEN)
    for abstract in all_abstracts:
        if 'AB' in abstract:
            abstract_text = re.sub(r'\r\n|\n|\s+|;', ' ', abstract['AB'])
            abstract_pmid = 'Unknown'
            abstract_journ = 'Unknown'
            if 'PMID' in abstract:
                abstract_pmid = abstract['PMID']
            if 'SO' in abstract:
                abstract_journ = abstract['SO']

            for key in gene_dict.keys():
                for gene in gene_dict[key]:
                    match = gene.search(abstract_text, re.MULTILINE)
                    if match:
                        if key not in result_dict:
                            result_dict[key] = []
                        result_dict[key].append([abstract_pmid, match.group(0), \
                                                 abstract_text[match.start(0)-30:match.end(0)+30]])
                        result = dict.fromkeys(pattern_dict.keys())
                        for pattern in sorted(pattern_dict.keys()):
                            result[pattern] = []
                            for match in pattern_dict[pattern].finditer(abstract_text, re.MULTILINE):
                                match = str(match.group(0))
                                if match not in result[pattern]:
                                    result[pattern].append(match)
                            result_dict[key][-1].append(', '.join(result[pattern]))
                        result_dict[key][-1].append(abstract_journ)
    ABS_OPEN.close()
    out_queue.put(result_dict)

# Run program if it call here
if __name__ == '__main__':

    arg_files = arg_receive()

    for key in arg_files.keys():
        if key == 'n':
            continue
        if os.path.exists(arg_files[key]) or key == 'o':
            arg_files[key] = os.path.abspath(arg_files[key])
        else:
            arg_files[key] = None

    print('Prepare patterns... ', end='')
    pattern_dict = pattern_prepare(arg_files['p'])
    print(len(pattern_dict), 'patterns.', sep=' ')

    print('Prepare gene dictionary... ', end='')
    gene_dict = gene_prepare(arg_files['g'], 1, arg_files['n'])
    genes = 0
    for chunk in gene_dict.keys():
        genes += len(gene_dict[chunk])
    print(genes, 'genes.', sep=' ')

    if arg_files['k']:
        print('Prepare KEGG dictionary... ', end='')
        kegg_dict = gene_prepare(arg_files['k'])
        print(len(kegg_dict), 'genes.', sep=' ')
    else:
        kegg_dict = None
    if arg_files['m']:
        print('Prepare miRNA-hosts dictionary... ', end='')
        host_dict = gene_prepare(arg_files['m'])
        print(len(host_dict), 'genes.', sep=' ')
    else:
        host_dict = None

    print('Searching...')

    result_dict = {}
    out_queue = Queue()
    searchers = [Process(target=abs_search, args=(gene_dict[cpu], pattern_dict, arg_files['a'], out_queue)) for cpu in range(cpu_count())]
    for search in searchers:
        search.start()
    for search in searchers:
        search.join()
        result_dict[search] = out_queue.get(timeout=1)

    RES_OPEN = open(arg_files['o'], 'w')
    RES_OPEN.write(';'.join(['PMID', 'Gene', 'Exact match', 'Wide string', \
                              'miR', 'Features 1', 'Features 2', 'Features 3',\
                              'KEGG Pathway', 'Host miR', 'Journal']))
    for proc in result_dict.keys():
        for gene in result_dict[proc]:
            for match in result_dict[proc][gene]:
                RES_OPEN.write('\n')
                RES_OPEN.write(';'.join([match[0], gene, ';'.join([x for x in match[1:-1]]), '']))
                if kegg_dict and gene in kegg_dict:
                    RES_OPEN.write(', '.join(kegg_dict[gene]))
                RES_OPEN.write(';')
                if host_dict and gene in host_dict:
                    RES_OPEN.write(', '.join(host_dict[gene]))
                RES_OPEN.write(';')
                RES_OPEN.write(match[-1])
    RES_OPEN.close()
    out_queue.close()
