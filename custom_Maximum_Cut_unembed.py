"""
The function custom_unembed_Maximum_Cut is the primary method; 
it accepts raw samples, the problem graph, and the minor embedding, and executes the custom unembedding algorithm on all of the samples
"""

import random
import networkx as nx

def find_broken_chains(vector, G_embedding):
        broken = []
        for a in G_embedding:
                chain = []
                for qubit in G_embedding[a]:
                        chain.append(vector[qubit])
                if len(list(set(chain))) == 2:
                        broken.append(a)
        return broken
def subgraph_agreement(vector, G_embedding):
        variables_0 = []
        variables_1 = []
        for a in G_embedding:
                chain = []
                for qubit in G_embedding[a]:
                        chain.append(vector[qubit])
                if list(set(chain)) == [-1]:
                        variables_0.append(a)
                if list(set(chain)) == [1]:
                        variables_1.append(a)
        return variables_0, variables_1
def maxcut_degree_search(G, broken, variables_0, variables_1, chain_freq):
        for node in broken:
                tmp_subg_0 = G.subgraph(variables_0+[node])
                tmp_subg_1 = G.subgraph(variables_1+[node])
                if tmp_subg_0.degree(node) > tmp_subg_1.degree(node):
                        variables_1.append(node)
                if tmp_subg_0.degree(node) < tmp_subg_1.degree(node):
                        variables_0.append(node)
                if tmp_subg_0.degree(node) == tmp_subg_1.degree(node):
                        if chain_freq[node] < 0.5:
                                variables_0.append(node)
                        if chain_freq[node] > 0.5:
                                variables_1.append(node)
                        if chain_freq[node] == 0.5:
                                choice = random.choice([0, 1])
                                if choice == 0:
                                        variables_0.append(node)
                                if choice == 1:
                                        variables_1.append(node)
        return variables_0, variables_1
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
def custom_unembed_Maximum_Cut(samples, G, G_minor_embedding):
        """
        Main unembedding method
        
        Parameters
        ----------
        samples : List
            The raw D-Wave samples list. Here the variables are assumed to be spins; meaning they are +1/-1.
        G : Networkx graph
            The problem graph instance which we are computing the Maximum Cut of.
        G_minor_embedding : Dictionary
            The minor embedding for the problem instance.
        
        Returns
        -------
        cuts : List
            List of cut values. Each sample corresponds to one cut value.
        """
        cuts = []
        for vector in samples:
                broken = find_broken_chains(vector, G_minor_embedding)
                variables_0, variables_1 = subgraph_agreement(vector, G_minor_embedding)
                chain_freq = find_chain_freq(broken, vector, G_minor_embedding)
                added_variables_0, added_variables_1 = maxcut_degree_search(G, broken, variables_0, variables_1, chain_freq)
                cuts.append(nx.cut_size(G, added_variables_0, added_variables_1))
        return cuts

