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
import statistics
import random
from random import randint
random.seed(9001)

__author__ = "Colin Davidson"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Colin Davidson"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Colin Davidson"
__email__ = "davidsonco@eisti.eu"
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
    """ Reads a fastq file to return a generator
      :Parameters:
        fastq_file : Path to the file
      Returns: yields the sequences
    """
    with open(fastq_file, 'rt') as f:
        for line in f:
            # Remove \n
            out = line.rstrip()
            # Could improve by checking all characters from sentence, but longer
            if out[0] in ['A','T','C','G']:
                yield out


def cut_kmer(read, kmer_size):
    """ Reads a sequence and yields the k-mer
      :Parameters:
        read : Sequence (str)
        kmer_size : Size of kmer (int)
      Returns: yields the kmers
    """
    for i in range(len(read) - kmer_size + 1):
        yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    """ Builds a dictionnary of kmers and occurences
      :Parameters:
        fastq_file : Path to the file
        kmer_size : Size of kmer (int)
      Returns:
        out : dictionnary of kmers and occurences (dict)
    """
    out = dict()

    for seq in read_fastq(fastq_file):
        for kmer in cut_kmer(seq, kmer_size):
            if kmer in out:
                out[kmer] += 1
            else:
                out[kmer] = 1
    return out


def build_graph(kmer_dict):
    """ Builds the graph using Networkx
      :Parameters:
        kmer_dict : dictionnary of kmers and occurences
      Returns:
        graph : graph of kmers (nx.Digraph)
    """
    graph = nx.DiGraph()
    for kmer in kmer_dict:
        pre = kmer[:-1]
        suf = kmer[1:]
        graph.add_weighted_edges_from([(pre, suf, kmer_dict[kmer])])
    return graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """ Removes paths from a graph
      :Parameters:
        graph : (nx.Digraph)
        path_list : paths to remove (list)
        delete_entry_node : delete entry node of the path (bool)
        delete_sink_node : delete sink node of the path (bool)
      Returns:
        graph : updated graph (nx.Digraph)
    """
    for path in path_list:
        for i in range(len(path)):
            if path[i] in graph:
                if i == 0 and delete_entry_node:
                    graph.remove_node(path[i])
                elif i == (len(path)-1) and delete_sink_node:
                    graph.remove_node(path[i])
                elif i not in [0, len(path) - 1]:
                    graph.remove_node(path[i])
    return graph


def std(data):
    """ Returns standard deviation
    """
    return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list,
                    delete_entry_node=False, delete_sink_node=False):
    """ Cleans the graph from its unecessary paths
      :Parameters:
        graph : (nx.Digraph)
        path_list : paths that we want to compare (list)
        path_length : length of each path (list (int))
        weight_avg_list : weight average of each path (list (int))
        delete_entry_node : delete entry node of the path (bool)
        delete_sink_node : delete sink node of the path (bool)
      Returns:
        graph : updated graph (nx.Digraph)
    """
    if path_list == []:
        return graph
    # Temporary best path index
    best_path_index = 0
    for i in range(1, len(path_list)):
        w1 = weight_avg_list[best_path_index]
        w2 = weight_avg_list[i]
        l1 = path_length[best_path_index]
        l2 = path_length[i]
        # Highest weight
        if w1 < w2:
            best_path_index = i
        elif w1 == w2:
            # Longest path
            if l1 < l2:
                best_path_index = i
            # Random
            elif l1 == l2:
                best_path_index = random.choice([best_path_index, i])

    # Update graph
    path_list.remove(path_list[best_path_index])
    graph = remove_paths(graph, path_list, delete_entry_node, delete_sink_node)

    return graph

def path_average_weight(graph, path):
    """ Computes the average weight of a path
      :Parameters:
        graph : (nx.Digraph)
        path : list of nodes
      Returns:
        out : average weight (int)
    """
    out = 0
    if len(path) != 0:
        for i in range(len(path)-1):
            for n, nbrs in graph.adj.items():
                if n == path[i]:
                    for nbr, eattr in nbrs.items():
                        out += eattr['weight']
                        break
                    break
        out /= (len(path)-1)
    return out

def solve_bubble(graph, ancestor_node, descendant_node):
    """ Removes bubbles from parts of a graph
      :Parameters:
        graph : (nx.Digraph)
        ancestor_node : node
        descendant_node : node
      Returns:
        graph : (nx.Digraph)
    """
    # Get all possibles paths from ancestor to descendant
    paths = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))

    # Find best path in graph
    path_length = []
    weight_avg_list = []
    for path in paths:
        path_length.append(len(path))
        weight_avg_list.append(int(path_average_weight(graph, path)))

    graph = select_best_path(graph, paths, path_length, weight_avg_list)

    return graph

def simplify_bubbles(graph):
    """ Removes every bubbles from a graph
      :Parameters:
        graph : (nx.Digraph)
      Returns:
        graph : (nx.Digraph)
    """
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)

    for start in starting_nodes:
        for end in ending_nodes:
            graph = solve_bubble(graph, start, end)

    return graph

def solve_entry_tips(graph, starting_nodes):
    """ TODO
    """
    pass

def solve_out_tips(graph, ending_nodes):
    """ TODO
    """
    pass

def get_starting_nodes(graph):
    """ Finds every starting nodes in a graph
      :Parameters:
        graph : (nx.Digraph)
      Returns:
        nodes : nodes that do not have predecessors (list)
    """
    nodes = []
    for node in graph.nodes:
        if list(graph.predecessors(node)) == []:
            nodes.append(node)
    return nodes

def get_sink_nodes(graph):
    """ Finds every sink nodes in a graph
      :Parameters:
        graph : (nx.Digraph)
      Returns:
        nodes : nodes that do not have successors (list)
    """
    nodes = []
    for node in graph.nodes:
        if list(graph.successors(node)) == []:
            nodes.append(node)
    return nodes

def get_contigs(graph, starting_nodes, ending_nodes):
    """ Finds contigs in a graph
      :Parameters:
        graph : (nx.Digraph)
        starting_nodes : list of nodes (list)
        ending_nodes : list of nodes (list)
      Returns:
        out : list of tuples with a path and its length (list)
    """
    out = []
    for start in starting_nodes:
        for end in ending_nodes:
            paths = nx.all_simple_paths(graph, start, end)
            paths = list(paths)
            if len(paths) != 0:
                path = paths[0]
                str_path = path[0]
                for i in range(1,len(path)):
                    str_path += path[i][-1]
                out.append((str_path, len(path) + 1))
    return out


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def save_contigs(contigs_list, output_file):
    """ Save a contig list following a specific format
      :Parameters:
        contigs_list : list of contigs (list)
        output_file : Path of the file
      Returns:
        void
    """
    with open(output_file, 'w') as f:
        i = 0
        for contig in contigs_list:
            text = ">contig_{} len={}\n{}\n".format(i, contig[1], contig[0])
            text = fill(text)
            f.write(text)
            i += 1

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    """
    # Get arguments
    args = get_arguments()

    # Build Graph
    d = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(d)
    starting_nodes = get_starting_nodes(graph)
    sink_nodes = get_sink_nodes(graph)

    # Get contigs
    contigs_list = get_contigs(graph, starting_nodes, sink_nodes)

    # Save contigs
    save_contigs(contigs_list, args.output_file)
    """
    pass

if __name__ == '__main__':
    main()
