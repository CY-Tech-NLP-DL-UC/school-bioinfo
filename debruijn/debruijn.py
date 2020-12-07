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
    """ TODO
    """
    with open(fastq_file, 'rt') as f:
        for line in f:
            # Remove \n
            out = line.rstrip()
            # Could improve by checking all characters from sentence, but longer
            if out[0] in ['A','T','C','G']:
                yield out


def cut_kmer(read, kmer_size):
    """ TODO
    """
    for i in range(len(read) - kmer_size + 1):
        yield read[i:i+kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    """ TODO
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
    """ TODO
    """
    graph = nx.DiGraph()
    for kmer in kmer_dict:
        pre = kmer[:-1]
        suf = kmer[1:]
        graph.add_weighted_edges_from([(pre, suf, kmer_dict[kmer])])
    return graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """ TODO
    """
    pass

def std(data):
    """ TODO
    """
    pass


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                    delete_entry_node=False, delete_sink_node=False):
    pass

def path_average_weight(graph, path):
    """ TODO
    """
    pass

def solve_bubble(graph, ancestor_node, descendant_node):
    """ TODO
    """
    pass

def simplify_bubbles(graph):
    """ TODO
    """
    pass

def solve_entry_tips(graph, starting_nodes):
    """ TODO
    """
    pass

def solve_out_tips(graph, ending_nodes):
    """ TODO
    """
    pass

def get_starting_nodes(graph):
    """ TODO
    """
    nodes = []
    for node in graph.nodes:
        if list(graph.predecessors(node)) == []:
            nodes.append(node)
    return nodes

def get_sink_nodes(graph):
    """ TODO
    """
    nodes = []
    for node in graph.nodes:
        if list(graph.successors(node)) == []:
            nodes.append(node)
    return nodes

def get_contigs(graph, starting_nodes, ending_nodes):
    """ TODO
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
    """ TODO
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


if __name__ == '__main__':
    main()
