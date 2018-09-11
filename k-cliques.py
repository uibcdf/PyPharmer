''' Cliques algorithm implementations

This library was written by https://github.com/gallexis at https://github.com/gallexis/nsd-K-cliques '''

from utils import *
import matrix
import itertools

def get_neighbours(node,graph):
    return graph[node]

def get_clique(node,k,B,graph):
    A = set()
    A.add(node)
    while len(B) > 0 and len(A) < k:
        n = B.pop(0)
        A.add(n)
        B = list(set(B).intersection(get_neighbours(n, graph)))
        if k == len(A):
            return {" ".join(sorted(A))}

    return []

def get_k_clique(k,graph):
    clique = []

    for node in graph:
        B = get_neighbours(node,graph)
        for subset in itertools.permutations(B, k-1):
            clique += get_clique(node,k,list(subset),graph)

    return clique

def merge_cliques(cliques_set):
    list_of_sets = []
    final_list_of_sets = []
    clique_list = sorted(cliques_set, key=str.__len__)

    [list_of_sets.append(set(c.split(' '))) for c in clique_list]

    index = 0
    max_len = len(list_of_sets)
    while index < max_len:
        set1 = list_of_sets[index]

        if is_biggest_subset(set1, list_of_sets, index+1, max_len):
            final_list_of_sets.append(sorted(set1))

        index+=1

    return final_list_of_sets

def is_biggest_subset(set1, sets, index, max_len):
    while index < max_len:
        if set1.issubset(sets[index]):
            return False

        index +=1

    return True

def get_community(clique_list, matrix):
    community=set()
    columnPos=0
    for line in matrix:
        if line[columnPos] == 1:
            for element in clique_list[columnPos]:
                community.add(element)
        columnPos=columnPos+1
    return community



if "__main__" == __name__:

    with open("louvain", "r+") as file:
        graph= file.read().splitlines()

        lg = load_graph(graph)
        delete_loop(lg)

        #print("Loaded graph:", lg,"\n")

        #Test with predefined cliques
        #test = merge_cliques({"1 2","1 2 3","4 3 1 2","4","2 9"})
        #print("------> Test merged clique: ", test, "\n")


        ### Display a k-clique

        k = 3
        cliques = get_k_clique(k,lg)

        #print(k,"-clique: ",cliques)
        #print("merged clique: ",merge_cliques(cliques),"\n")


        ### Generate then merge k1 to kn cliques
        """
        cliques = []
        k1 = 2
        k2 = get_max_degree_nodes(lg)+1
        while True:
            print("Generation of a",k1,"clique")
            c= get_k_clique(k1,lg)
            if c == []:
                break
            k1 += 1
            cliques += c
        """
        #print("All cliques: ", cliques)
        merged_cliques = merge_cliques(cliques)

        print("Merged cliques: ", merged_cliques)

        overlapMatrix = matrix.createOverlapMatrix(merged_cliques)
        print(overlapMatrix)
        finlaMatrix=matrix.createCommunitiesMatrix(overlapMatrix,3 )
        print(finlaMatrix)
        print("community: "+str(get_community(merged_cliques,finlaMatrix)))


'''From matrix.py file of k-cliques library'''


# optimize comparison. Can be done by traingular part
def createOverlapMatrix(cliqueList):
    overlapMatrix = [[0 for _ in range(len(cliqueList))] for _ in range(len(cliqueList))]
    cliquePos = 0
    for clique in cliqueList:
        compCliquePos = 0
        for compClique in cliqueList:
            overlapMatrix[cliquePos][compCliquePos] = len(set(clique).intersection(compClique))
            compCliquePos = compCliquePos + 1
        cliquePos = cliquePos + 1
    # print (overlapMatrix)
    return overlapMatrix

def createCommunitiesMatrix(overlapMatrix, k):
    communititesMatrix = [[0 for _ in range(len(overlapMatrix))] for _ in range(len(overlapMatrix))]
    linePos = 0
    for line in overlapMatrix:
        columnPos = 0
        for node in line:
            if linePos == columnPos:
                if node >= k:
                    communititesMatrix[linePos][columnPos] = 1
            else:
                if node >= k - 1:
                    communititesMatrix[linePos][columnPos] = 1
            columnPos = columnPos + 1
        linePos = linePos + 1
    # print (communititesMatrix)
    return communititesMatrix


if __name__ == "__main__":
    cliqueList = [['d', 'e', 'h', 'i', 'j'], ['d', 'e', 'g', 'h'], ['b', 'd', 'e'], ['b', 'e', 'f'],
                  ['c', 'd', 'i', 'j'], ['a', 'b', 'c', 'd']]
    overlapMatrix = createOverlapMatrix(cliqueList)
    print(createCommunitiesMatrix(overlapMatrix, 4))


''' From utils.py of k-cliques library'''


def get_nodes_from_edge(edge):
    return edge.split(' ')

def get_nodes_from_graph(graph):
    nodesSet = set()

    for edge in graph:
        nodes = get_nodes_from_edge(edge)
        [nodesSet.add(n) for n in nodes]

    return list(nodesSet)

def size_of_graph(graph):
    nodes = get_nodes_from_graph(graph)
    return len(nodes)


def load_graph(graph):
    nodes = {}

    for edge in graph:
        ns = get_nodes_from_edge(edge)
        n0 = ns[0]

        # if node's degree == 0, empty list
        if len(ns) == 1:
            nodes[n0] = set()

        else:
            n1 = ns[1]
            # if nodes already in dict, append his new neighbour
            if n0 in nodes:
                nodes[n0].add(n1)

            # if node not yet in dict, initialise a new list & add his neighbour
            else:
                nodes[n0] = set(n1)

            # same as before, but with opposite nodes, in case there is
            # A B, but not B A
            if n1 in nodes:
                nodes[n1].add(n0)
            else:
                nodes[n1] = set(n0)

    return nodes

def delete_loop(loadedGraph):
    for node in loadedGraph:
        if node in loadedGraph[node]:
            loadedGraph[node].remove(node)

    return loadedGraph

def get_max_degree_nodes(loadedGraph):
    return len( max(loadedGraph.values()) )