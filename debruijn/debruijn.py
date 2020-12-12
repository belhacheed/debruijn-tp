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
import statistics
import random
import operator
import networkx as nx
import matplotlib


#from operator import itemgetter
random.seed(9001)
#from random import randint


__author__ = "Edouard Belhache"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Edouard Belhache"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Edouard Belhache"
__email__ = "belhacheed@eisti.eu"
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


def read_fastq(fastq_file):#retourne un générateur de séquences
    with open(fastq_file, "r") as f:
        for line in f:
            yield next(f).strip()
            next(f)
            next(f)


def cut_kmer(read, kmer_size):#retourne un générateur de k-mer
    for i in range(0,len(read)+1-kmer_size):
        kmer_parts=read[i:i+kmer_size]
        yield(kmer_parts) #obliger de retourner les kmers (warning en conséquence)


def build_kmer_dict(fastq_file, kmer_size):
#retourne un dictionnaire ayant pour clé le k-mer et pour valeur le nombre d’occurrence de ce k-mer
    kmer_compteur = {}
    for read in read_fastq(fastq_file):#appel fonction
        kmer_parts = cut_kmer(read, kmer_size)#appel fonction
        for kmer in kmer_parts:
            if kmer in kmer_compteur:
                kmer_compteur[kmer] = kmer_compteur[kmer]+1
            else:
                kmer_compteur[kmer] = 1
    return kmer_compteur


def build_graph(kmer_dict):# créera l’arbre de k-mers préfixes et suffixes
    graphe = nx.DiGraph()
    for kmer in kmer_dict.keys():
        prefixe = kmer[:-1]
        suffixe = kmer[1:]
        graphe.add_edge(prefixe, suffixe, weight = kmer_dict[kmer])
    return graphe


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        if delete_entry_node == True and delete_sink_node == False:
            graph.remove_nodes_from(path[:-1])
        elif delete_entry_node == False and delete_sink_node == True:
            graph.remove_nodes_from(path[1:])
        elif delete_entry_node == False and delete_sink_node == False:
            graph.remove_nodes_from(path[1:-1])
        else: #case: delete_entry_node == True and delete_sink_node == True
            graph.remove_nodes_from(path)
    return graph


def std(data):
    return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    pass


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    #  First get the best weight
    liste = []
    b_liste = []
    max_weight = max(weight_avg_list)
    max_length = max(path_length)
    #  List of best weight index
    for i in range(len(weight_avg_list)):
        if weight_avg_list[i] == max_weight:
            liste.append(i)
    if len(liste) !=1:
        #  Best weight paths
        for i in liste:
            if path_length[i] == max_length:
                b_liste.append(i)
            elif path_length[i] > max_length:
                max_length = path_length[i]
        if len(b_liste) > 1:
            best_path = best_path[0]
        else:
            best_path = b_liste[0]
    else:
        best_path = liste[0]
    remove = path_list[:best_path] + path_list[best_path +1:]
    graph = remove_paths(graph, remove, delete_entry_node, delete_sink_node)#appel fonction
    return graph


def path_average_weight(graph, path):
    liste = []
    subgraphe=graph.subgraph(path)
    for a, b, edge in subgraphe.edges(data=True):
        liste.append(edge["weight"])
    moy = statistics.mean(liste)
    return moy


def solve_bubble(graph, ancestor_node, descendant_node):
    # List of all simple paths
    path_lengths = []
    weight_avg_list = []
    list_paths = list(nx.shortest_simple_paths(graph, ancestor_node, descendant_node))
    for path in list_paths:
        path_lengths.append(len(path))
        weight_avg_list.append(path_average_weight(graph, path))#appel fonction
    graph = select_best_path(graph, list_paths, path_lengths, weight_avg_list)#appel fonction
    return graph


def simplify_bubbles(graph):
    print("test121")


def solve_entry_tips(graph, starting_nodes):
    pass


def solve_out_tips(graph, ending_nodes):
    pass


def get_starting_nodes(graph):#retourne une liste de noeuds d'entrée
    liste_noeuds = []
    for noeud in graph.nodes():
        if not graph.pred[noeud]:
            liste_noeuds.append(noeud)
    return liste_noeuds

def get_sink_nodes(graph):#retourne une liste de noeuds de sortie
    liste_noeuds = []
    for noeud in graph.nodes():
        if not graph.succ[noeud]:
            liste_noeuds.append(noeud)
    return liste_noeuds


def get_contigs(graph, starting_nodes, ending_nodes):#saut de ligne pour pyint
    #retourne une liste de tuple(contig, taille du contig)
    liste_contigs = []
    for deb in starting_nodes:
        for fin in ending_nodes:
            for path in nx.shortest_simple_paths(graph, deb, fin):
                contig = [c[:-1] for c in path]
                contig.append(fin[-1])
                liste_contigs.append(("".join(contig), len(contig)))
    return liste_contigs


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    #le commentaire au dessus ajoute du score (???)
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(contigs_list, output_file):#écrit un fichier de sortie contenant les contigs
    with open(output_file, "w") as f:
        for i in range(len(contigs_list)):
            f.write('>contig_{} len={}\n'.format(i, contigs_list[i][1]))
            f.write(fill(contigs_list[i][0]))
            f.write('\n')
    f.close()


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    taille = args.kmer_size
    fichier = (args.fastq_file)
    dictionnaire = build_kmer_dict(fichier,taille)
    graphe = build_graph(dictionnaire)
    #graphe = simplify_bubbles(graphe)
    #...

if __name__ == '__main__':
    main()
