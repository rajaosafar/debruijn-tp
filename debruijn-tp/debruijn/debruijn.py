#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import matplotlib.pyplot as plt

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()

def read_fastq(fastq_file):
    with open(fastq_file, "rt") as myfile:
        for line in myfile:
            yield next(myfile).replace('\n', '')
            next(myfile)
            next(myfile)

def cut_kmer(read, kmer_size):
    for i, letter in enumerate(read):
        if i<=len(read)-kmer_size:
            yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    lines = read_fastq(fastq_file)
    kmers = []
    for line in lines:
        for kmer in cut_kmer(line, kmer_size):
            kmers.append(kmer)

    dict_kmers = dict()
    for kmer in kmers:
        if kmer not in dict_kmers.keys():
            dict_kmers[kmer]=1
        else:
            dict_kmers[kmer]+=1
    return dict_kmers

def build_graph(kmer_dict):
    graph = nx.DiGraph()
    for kmer, weight in kmer_dict.items():
        graph.add_edge(kmer[:(len(kmer)-1)], kmer[1:], weight=weight)
    return graph

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    pass

def std(data):
    pass


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    pass

def path_average_weight(graph, path):
    pass

def solve_bubble(graph, ancestor_node, descendant_node):
    pass

def simplify_bubbles(graph):
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    nodes = []
    for node in graph.nodes:
        count=0
        for p in graph.predecessors(node):
            count +=1
        if count==0:
            nodes.append(node)
    return nodes

def get_sink_nodes(graph):
    nodes = []
    for node in graph.nodes:
        count=0
        for p in graph.successors(node):
            count +=1
        if count==0:
            nodes.append(node)
    return nodes

def get_contigs(graph, starting_nodes, ending_nodes):
    contigs = []
    for s in starting_nodes:
        for e in ending_nodes:
            #all paths between s and e
            all_paths = nx.all_simple_paths(graph,s,e)
            all_paths = list(all_paths)

            for path in all_paths:
                contig = path[0]
                for st in path[1:]:
                    contig += st[len(st)-1:]
                contigs.append((contig,len(contig)))
    return contigs

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(contigs_list, output_file):
    with open(output_file, "wt") as myfile:
        for i, (contig, len_contig) in enumerate(contigs_list):
            myfile.write(">contig_"+str(i)+ " len="+str(len_contig)+"\n")
            myfile.write(fill(contig, width=80)+"\n")


def draw_graph(graph, graphimg_file):
    """Draw the graph
    """
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5,
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    fastq_file = args.fastq_file
    kmer_dict = build_kmer_dict(fastq_file, 21)

    graph = build_graph(kmer_dict)

    starting_nodes=get_starting_nodes(graph)
    ending_nodes=get_sink_nodes(graph)

    contigs = get_contigs(graph, starting_nodes, ending_nodes)
    save_contigs(contigs,"test.fna")

if __name__ == '__main__':
    main()
