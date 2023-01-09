"""
The function custom_unembed_Minimum_Vertex_Cover is the primary method; 
it accepts raw samples, the problem graph, and the minor embedding, and executes the custom unembedding algorithm on all of the samples
"""

def minvc_core_unembed(vector, G, broken_chains, chain_freq):
        """
        vector is a dictionary where each key is a node and each value is 0 or 1 corresponding to chains that agreed completely on being either 0 or 1
        G is the problem graph
        broken_chains is a list of nodes where each chain for that node was broken
        remark: length(broken_chains)+length(vector) = size of the problem instance
        """
        Z = []#variables assigned 0
        core = []#Variables assigned 1
        for a in vector:
                if vector[a] == 0:
                        Z.append(a)
                if vector[a] == 1:
                        core.append(a)
        for node in Z:
                if len(list(set(Z) & set(list(G.neighbors(node))))) != 0:
                        return len(G)
        X = set(broken_chains) & (set([G.neighbors(z) for z in Z]) - set(Z)) #broken ch. neighbors of 0
        core += list(X) # assign 1 to x in X
        broken_chains  = set([b for b in broken_chains if b not in X]) # remove X from broken_chains
        # sort u by minimum degree first in order to satisfy condition 2
        while len(broken_chains) > 0:
            u, G_br = best_chain(broken_chains, chain_freq, G)
            if len(list(set(Z) & set(list(G_br.neighbors(u))))) != 0:
                core.append(u)
            else:
                Z.append(u)
        return len(core)
def unbroken_chain(vector, G_minor_embedding):
    unembedded = {}
    for a in G_minor_embedding:
        chain = []
        for qubit in G_minor_embedding[a]:
            chain.append(vector[qubit])
        if list(set(chain)) == [1]:
            unembedded[a] = 1
        if list(set(chain)) == [0]:
            unembedded[a] = 0
    return unembedded
def find_chain_freq(broken_chains, a, G_minor_embedding):
        freq = {}
        for node in broken_chains:
                chain = G_minor_embedding[node]
                length = 0
                count = 0
                for v in chain:
                        length += 1
                        count += a[v]
                freq[node] = count/float(length)
        for a in freq:
                if freq[a] > 0.9999:
                        freq[a] = 0.99
        return freq
def best_chain(broken_chains, chain_freq, G):#sorts by degree first, then chain agreement on taking the value of 1 in the QUBO
        G_br = G.subgraph(list(broken_chains))
        pr = {x:G_br.degree[x]+chain_freq[x]  for x in G_br}
        best_chain = max(pr.keys(), key=(lambda k: pr[k]))
        broken_chains.remove(best_chain)
        return best_chain, G_br
def find_broken_chains(vector, G_minor_embedding):
        broken = []
        for a in G_minor_embedding:
                chain = []
                for qubit in G_minor_embedding[a]:
                        chain.append(vector[qubit])
                if len(list(set(chain))) == 2:
                        broken.append(a)
        return broken
def custom_unembed_Minimum_Vertex_Cover(samples, G, G_minor_embedding):
        """
        Main unembedding method
        
        Parameters
        ----------
        samples : List
            Samples from the quantum annealer. 
            These samples are assumed to be from sampling a QUBO; so the active variables will be either 0 or 1
        G : Networkx graph
            The problem unweighted, undirected graph which are are attempring to find the minimum vertex cover of.
        G_minor_embedding : Dictionary
            The minor embedding dictionary for the problem graph G (or equivalently the QUBO).
            The keys are integers corresponding to the logical variables, and the values are the lists of physical qubits which form that chain.
    
        Returns
        -------
        covers : List
            List of vertex cover objective function values.
            The length of this list is exactly equal to the length of the samples

        """
        covers = []
        for a in samples:
                unembedded = unbroken_chain(a, G_minor_embedding)
                broken_chains = find_broken_chains(a, G_minor_embedding)
                chain_freq = find_chain_freq(broken_chains, a, G_minor_embedding)
                cover = minvc_core_unembed(unembedded, G, broken_chains, chain_freq)
                covers.append(cover)
        return covers
    
