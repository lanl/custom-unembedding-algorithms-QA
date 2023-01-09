"""
The function custom_unembed_Maximum_Clique is the primary method; 
it accepts raw samples, the problem graph, and the minor embedding, and executes the custom unembedding algorithm on all of the samples
"""

import random

def is_clique(G):
    """
    parameters: G (networkx.Graph() object)
    description: checks if the input graph is a clique
    return Boolean True/False
    """
    n = len(list(G.nodes()))
    m = len(list(G.edges()))
    if m == (n*(n-1))/2:  # A complete graph has number of edges equal to n choose 2 for n nodes
        return True
    else:
        return False
def chain_based_tie_break(nodes, vector, embedding):
    results = {}
    results_list = []
    for node in nodes:
        tmp = 0
        c = 0
        for var in embedding:
            c += 1
            if vector[var] == 1:
                tmp += 1
        results[node] = tmp/float(c) # fraction of the 1s
        results_list.append(tmp/float(c))
    if len(list(set(results_list))) == 1:
        return random.choice(nodes)
    else:
        for a in results:
            if results[a] == max(results_list):
                return a
def broken_chain_search(subg, G, nodes_to_search, vector, embedding):
    while len(nodes_to_search) != 0:
        degrees = {}
        degrees_list = []
        nodes_to_search2 = set()
        subg_nodes = list(subg.nodes())
        edges_to_search = G.subgraph(nodes_to_search).edges()
        for e in edges_to_search: #try first adding an edge
            tmp_subg = G.subgraph(subg_nodes+[e[0],e[1]])
            if is_clique(tmp_subg):
                nodes_to_search2.add(e[0])
                nodes_to_search2.add(e[1])
        if len(nodes_to_search2) == 0:
            for v in nodes_to_search: #try now with single vertices
                tmp_subg = G.subgraph(subg_nodes+[v])
                if is_clique(tmp_subg):
                    nodes_to_search2.add(v)
            if len(nodes_to_search2) == 0:
                return subg
        nodes_to_search = list(nodes_to_search2) #  nodes that may extend the clique
        tmp_subg = G.subgraph(nodes_to_search)  # graph of candidate nodes
        for node in nodes_to_search: #compute all degrees
            deg = tmp_subg.degree(node)
            degrees_list.append(deg)
            degrees[node] = deg
        nodes_to_add = [] #fmax-degree nodes
        for node in degrees: #compute nodes_to_add
            if degrees[node] == max(degrees_list):
                nodes_to_add.append(node)
        if len(nodes_to_add) == 1:
            node_to_add = nodes_to_add[0]
            subg = G.subgraph(subg_nodes+[node_to_add])
            nodes_to_search.remove(node_to_add)
        else:
            node_to_add = chain_based_tie_break(nodes_to_add, vector, embedding)
            subg = G.subgraph(subg_nodes+[node_to_add])
            nodes_to_search.remove(node_to_add)
    return subg
def clique_core_unembed(vector, G, broken_chains, embedded_vector, embedding):
    """
    parameters: vector (dictionary), G (networkx Graph()), broken_chains (list)
    description: vector is an unembedded dictionary where the keys are the graph node variables, and the values are either 0 or 1
    G is the problem graph, and broken chains is a list of all variables that had broken chains
    If the group of nodes that took the value of 1 in the unembedded sample forms a clique, then the broken_chain_search is called
    in order to attempt to find a larger clique size, otherwise no clique was found and 0 is returned
    return: clique size found (integer)
    """
    h = []
    for a in vector:
        if vector[a] == 1:  # If the unembedding has a 1 for this variable
            h.append(a)  # add the node to the list
    # create subgraph from the variables that took the value of 1
    subg = G.subgraph(h)
    if is_clique(subg) == True:  # Check if the subgraph is a clique
        # Try to improve the solution using clique-core
        return len(broken_chain_search(subg, G, broken_chains, embedded_vector, embedding))
    else:
        return 0
def find_broken_chains(vector, minor_embedding):
    """
    parameters: vector (list), minor_embedding (dictionary)
    description: The input is a single D-Wave sample of length 2048, and the problem embedding. 
    The function returns a list of variables which represent some subset of graph nodes, where each of those variables
    had a broken chain for this particular D-Wave sample
    return broken (list)
    """
    broken = []
    for a in minor_embedding:  # for each node in the embedding
        chain = []
        # Iterate through the entire chain for variable a
        for qubit in minor_embedding[a]:
            chain.append(vector[qubit])  # create list of the chain variables
        # in other words, if the chain list is made up of both 0 and 1, which means it is broken
        if len(list(set(chain))) == 2:
            broken.append(a)
    return broken
def unbroken_chain_1s(vector, minor_embedding):
    """
    parameters: vector (list), minor_embedding (dictionary)
    description: Given an input of the list vector (which is a single D-Wave sample, and has a length of 2048 for 2000Q devices), and the embedding used,
    this function finds all chains for this sample that agree (i.e. every qubit took the value of 1). The output is a dictionary where the keys are all nodes from the graph,
    and all values are either 0 or 1, where the value is 0 is the chain had any 0's in it.
    return: unembedded (dictionary)
    """
    unembedded = {}
    for a in minor_embedding:
        chain = []
        # Iterate through the entire chain for variable a
        for qubit in minor_embedding[a]:
            chain.append(vector[qubit])  # create list of the chain variables
        if list(set(chain)) == [1]:
            unembedded[a] = 1
        else:
            # If the chain has any 0 at all, then we say the embedded value is 0
            unembedded[a] = 0
    return unembedded
def custom_unembed_Maximum_Clique(samples, G, minor_embedding):
    """
    Main unembedding method
    
    Parameters
    ----------
    samples : List
        List of D-Wave quantum annealer samples; QUBO format so the variables are 0/1.
    G : Networkx graph
        The problem graph which we are attempting to find the maximum clique of.
    minor_embedding : Dictionary
        Minor embedding for the problem instance.

    Returns
    -------
    List
        List of unembedded clique numbers.

    """
    cliques = []
    for a in samples:
        # Unembed the sample just by whatever chains agree in 1's
        unembedded = unbroken_chain_1s(a, minor_embedding)
        # Find all broken chains. The variable itself is a list of graph nodes
        broken_chains = find_broken_chains(a, minor_embedding)
        # The particular clique number found for this given sample
        clique = clique_core_unembed(unembedded, G, broken_chains, a, minor_embedding)
        cliques.append(clique)
    return cliques  # return all cliques found as a list

