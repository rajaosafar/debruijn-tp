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
import itertools

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
            # Remove the '\n' at the end of each line
            yield next(myfile).replace('\n', '')
            next(myfile)
            next(myfile)

def cut_kmer(read, kmer_size):
    for i, letter in enumerate(read):
        # Only take kmer of size kmer_size
        if i<=len(read)-kmer_size:
            yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    kmers = []
    # Get the list of kmers
    lines = read_fastq(fastq_file)
    for line in lines:
        for kmer in cut_kmer(line, kmer_size):
            kmers.append(kmer)

    dict_kmers = dict()
    # For each kmer, if it is not in the dict, add occurence of 1, else increment the occurence
    for kmer in kmers:
        if kmer not in dict_kmers.keys():
            dict_kmers[kmer]=1
        else:
            dict_kmers[kmer]+=1
    return dict_kmers

def build_graph(kmer_dict):
    graph = nx.DiGraph()
    # Divide each kmer in two parts and add a weighted edge
    for kmer, weight in kmer_dict.items():
        graph.add_edge(kmer[:(len(kmer)-1)], kmer[1:], weight=weight)
    return graph

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    # For each path, delete the nodes between the entry and sink nodes if they exist in the graph
    for path in path_list:
        if not delete_entry_node and not delete_sink_node:
            for i in range(1,len(path)-1):
                if path[i] in graph.nodes:
                    graph.remove_node(path[i])
        elif delete_entry_node and not delete_sink_node:
            for i in range(0,len(path)-1):
                if path[i] in graph.nodes:
                    graph.remove_node(path[i])
        elif not delete_entry_node and delete_sink_node:
            for i in range(1,len(path)):
                if path[i] in graph.nodes:
                    graph.remove_node(path[i])
        else:
            for i in range(0,len(path)):
                if path[i] in graph.nodes:
                    graph.remove_node(path[i])
    return graph

def std(data):
    # Return the standard deviation
    return statistics.stdev(data)

def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):

    # Get the list of the paths with an average weight equal to the max average weight
    # Get the list of their corresponding lengths
    heavier_val = max(weight_avg_list)
    heavier_list = []
    new_lengths = []
    for i in range(len(path_list)):
        if weight_avg_list[i] == heavier_val:
            heavier_list.append(path_list[i])
            new_lengths.append(path_length[i])

    # Get the list of the paths with a length equal to the max length among the paths in the new list
    longest_val = max(new_lengths)
    longest_list = []
    for i in range(len(heavier_list)):
        if new_lengths[i] == longest_val:
            longest_list.append(heavier_list[i])

    # Select a random path among the paths in the new list
    random_index = random.randint(0,len(longest_list)-1)

    # Remove it from the initial path_list
    path_list.remove(longest_list[random_index])

    # Remove all the paths in path_list
    return remove_paths(graph, path_list, delete_entry_node, delete_sink_node)


def path_average_weight(graph, path):
    s=0
    # Add all the weights
    for i in range(len(path)-1):
        s += graph[path[i]][path[i+1]]['weight']
    # Divide by the number of weights
    return s/(len(path)-1)

def solve_bubble(graph, ancestor_node, descendant_node):
    # Get all the paths between the ancestor node and the descendant node
    all_paths = nx.all_simple_paths(graph,ancestor_node,descendant_node)
    all_paths = list(all_paths)

    # Get a list with all the lengths
    # Get a list with all the average weights
    lengths = []
    weights = []
    for path in all_paths:
        lengths.append(len(path))
        weights.append(path_average_weight(graph,path))

    # Select the best path among these paths
    return select_best_path(graph, all_paths, lengths, weights)

def simplify_bubbles(graph):
    couples = []
    # For each node
    for e in graph.nodes:
        # Get a list of all the predecessors
        preds = list(graph.predecessors(e))
        ancestors = []
        # If there is more than one predecessor
        if len(preds)>1:
            # Look at the lowest common ancestor between each pair of predecessors
            for node1, node2 in itertools.combinations(preds,2):
                # If there is one, add it to the list of ancestors
                # (there is a bubble between s and e)
                s = nx.lowest_common_ancestor(graph, node1, node2, None)
                if s is not None:
                    ancestors.append(s)
            # Get a list of couples representing a bubble
            for s in list(set(ancestors)):
                couples.append((s,e))

    # Remove each bubble
    for s,e in couples:
        # If the node has not been removed yet
        if s in graph.nodes and e in graph.nodes:
            graph = solve_bubble(graph,s,e)

    return graph

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    nodes = []
    # For each node, count the number of predecessors
    for node in graph.nodes:
        count=0
        for p in graph.predecessors(node):
            count +=1
        # If it has no predecessor, it is a starting node
        if count==0:
            nodes.append(node)
    return nodes

def get_sink_nodes(graph):
    nodes = []
    # For each node, count the number of successors
    for node in graph.nodes:
        count=0
        for p in graph.successors(node):
            count +=1
        # If it has no successor, it is a sink node
        if count==0:
            nodes.append(node)
    return nodes

def get_contigs(graph, starting_nodes, ending_nodes):
    contigs = []
    for s in starting_nodes:
        for e in ending_nodes:
            # Get all the paths between s and e
            all_paths = nx.all_simple_paths(graph,s,e)
            all_paths = list(all_paths)

            for path in all_paths:
                # Get the content of the first node
                contig = path[0]
                # Add the last letter of all the other nodes to the contig
                for st in path[1:]:
                    contig += st[len(st)-1:]
                contigs.append((contig,len(contig)))
    return contigs

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(contigs_list, output_file):
    # Write each contig to the file
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

    #path = ['TCAGAGCTCTAGAGTTGGTT', 'CAGAGCTCTAGAGTTGGTTC']
    #print(path_average_weight(graph, path))

    graph_1 = nx.DiGraph()
    graph_1.add_weighted_edges_from([(3,2,1),(2,8,2),(8,6,3),(6,7,4),(2,4,5),(4,7,6),
                                    (10,7,7),(3,12,8),(12,11,9),(11,7,10)])
    graph_1 = simplify_bubbles(graph_1)

    print(graph_1.edges)

    graph_2 = nx.DiGraph()
    graph_2.add_weighted_edges_from([(3,2,1),(2,8,2),(8,9,3),(9,5,4),(5,6,5),(2,4,6),
                                    (4,9,7),(4,5,8),(2,10,9),(10,5,10),(5,7,11)])
    graph_2 = simplify_bubbles(graph_2)

    print(graph_2.edges)


if __name__ == '__main__':
    main()
