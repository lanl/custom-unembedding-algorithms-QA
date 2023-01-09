"""
The function custom_unembed_Graph_Partitioning is the primary method; 
it accepts raw samples, the problem graph, and the minor embedding, and executes the custom unembedding algorithm on all of the samples
"""

import random
import networkx as nx
import math

def search(subg1, subg0, broken, G, chain_freq):
        N = G.number_of_nodes()
        size0 = random.choice([math.ceil(N/2.0), math.floor(N/2.0)])
        size1 = N-size0
        for node in broken:
                if len(subg1) == size1:
                        subg0.append(node)
                        continue
                if len(subg0) == size0:
                        subg1.append(node)
                        continue
                tmp_subg0 = G.subgraph(subg1+[node])
                tmp_subg1 = G.subgraph(subg0+[node])
                deg0 = tmp_subg0.degree(node)
                deg1 = tmp_subg1.degree(node)
                if deg0 == deg1:
                        if chain_freq[node] < 0.5:
                                subg0.append(node)
                        if chain_freq[node] > 0.5:
                                subg1.append(node)
                        if chain_freq[node] == 0.5:
                                if len(subg0) < len(subg1):
                                        subg0.append(node)
                                if len(subg1) < len(subg0):
                                        subg1.append(node)
                                if len(subg1) == len(subg0):
                                        choice = random.choice([0, 1])
                                        if choice == 0:
                                                subg0.append(node)
                                        if choice == 1:
                                                subg1.append(node)
                if deg0 > deg1:
                        subg0.append(deg0)
                if deg1 > deg0:
                        subg1.append(deg1)
        return nx.cut_size(G, subg0, subg1)
def find_chain_agreements(vector, minor_embedding):
        broken = []
        subg_0 = []
        subg_1 = []
        for a in minor_embedding:
                chain = []
                for var in minor_embedding[a]:
                        chain.append(vector[var])
                if list(set(chain)) == [-1]:
                        subg_0.append(a)
                if list(set(chain)) == [1]:
                        subg_1.append(a)
                if len(list(set(chain))) == 2:
                        broken.append(a)
        return broken, subg_0, subg_1
def find_chain_freq(broken_chains, vector, G_embedding):
        freq = {}
        for node in broken_chains:
                chain = G_embedding[node]
                length = 0
                count = 0
                for v in chain:
                        length += 1
                        if vector[v] == 1:
                                count += 1
                freq[node] = count/float(length)
        return freq
def custom_unembed_Graph_Partitioning(samples, G, G_minor_embedding):
        """
        Main unembedding method
        
        Parameters
        ----------
        samples : List
            List of raw samples from the Quantum Annealer. These samples are based on the Graph Partitioning Ising problem defined on the unweighted, undirected graph G. 
            These samples are assumed to be -1/1 variable assignments, i.e. spins from sampling an Ising problem
        G : Networkx graph
            The problem graph instance.
        G_minor_embedding : Dictionary
            Minor Embedding dictionary, structured as integer keys (for the logical variables) and values of lists of physical qubits. Must match the problem graph G
        
        Returns
        -------
        results : List
            Unembedded cut sizes for all samples - the cut size is computed on the problem graph G.
    
        """
        N = G.number_of_nodes()
        cuts = []
        for vector in samples:
                broken, subg_0, subg_1 = find_chain_agreements(vector, G_minor_embedding)
                chain_freq = find_chain_freq(broken, vector, G_minor_embedding)
                if len(subg_0) > math.ceil(N/2.0) or len(subg_1) > math.ceil(N/2.0):
                        continue
                if len(broken) == 0:#If all chains agree, we do not need to process any further
                        cuts.append(nx.cut_size(G, subg_0, subg_1))
                        continue
                if len(subg_0) == math.floor(N/2.0) or len(subg_0) == math.ceil(N/2.0):#case where the solution is already half partitioned
                        subg_1 += broken
                        cuts.append(nx.cut_size(G, subg_0, subg_1))
                        continue
                if len(subg_1) == math.floor(N/2.0) or len(subg_1) == math.ceil(N/2.0):#case where the solution is already half partitioned
                        subg_0 += broken
                        cuts.append(nx.cut_size(G, subg_0, subg_1))
                        continue
                cuts.append(search(subg_1, subg_0, broken, G, chain_freq))
        return cuts

