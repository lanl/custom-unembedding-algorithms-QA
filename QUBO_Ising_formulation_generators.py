"""
These functions generate Ising or QUBO formulations given in input graph instance for the four NP-Hard problems
Minimum Vertex Cover
Maximum Cut
Maximum Clique
Graph Partitioning

The QUBO formulations are dictionaries, structured to be consistent with the D-Wave SDK; each QUBO is a single dictionary.
The Ising formulations are two dictionaries, one composed of Linear terms and the other of Quadratic terms

The input graph instances are assumed to be Networkx
"""

import networkx as nx
import math

def minimum_vertex_cover_QUBO(G):
    Q = {}
    for a in list(G.edges()):
        Q[a] = 2
    for i in list(G.nodes()):
        Q[(i, i)] = (-2*G.degree(i))+1
    return Q
def maximum_cut_Ising(G):
    h = {}
    j = {}
    for a in list(G.edges()):
        #Antiferromagnetic weight
        j[a] = 1
    return h, j
def maximum_clique_QUBO(G):
    Q = {}
    GC = nx.algorithms.operators.unary.complement(G)  # complement of the graph
    for i in list(GC.nodes()):
        Q[(i, i)] = -1
    for a in list(GC.edges()):
        Q[a] = 2
    return Q
def graph_partitioning_Ising(G):
        h = {}
        j = {}
        B = 1
        degree_sequence = sorted([d for n, d in G.degree()], reverse=True)
        delta = max(degree_sequence)
        bound = min([2*delta, len(G.nodes())])/float(8)
        A = math.ceil(bound)
        for edge in list(G.edges()):
                j[edge] = -0.5*B + A*2
        GC = nx.algorithms.operators.unary.complement(G)
        for edge in list(GC.edges()):
                j[edge] = A*2
        return h, j
