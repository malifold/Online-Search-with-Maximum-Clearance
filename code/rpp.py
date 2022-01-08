# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 16:48:09 2020

@author: Malachi
"""

EPSILON=1e-5

import networkx as nx
from cpp import cpp as cpp
import numpy as np



def rpp(g,r,start=None):
    #build subgraph we need to visit
    rgraph = nx.Graph()
    for [x,y] in r:
        rgraph.add_edge(x,y,weight=g[x][y]['weight'])
    #unused nodes
    unused=list(g)
    for x in list(rgraph):
        unused.remove(x)
    unused = [set([x]) for x in unused]
    #print("unused",unused)
    #get connected components
    ccs = list(nx.connected_components(rgraph))
    #print("{} ccs".format(len(ccs)))
    #contract g with the ccs of rgraph
    
    def edge_data_min_weight(A,B):
        if len(A)>len(B):
            A,B=B,A
        #print('\n')
        #print('A : ',A)
        #print('B : ',B)                        
        dists, paths = nx.multi_source_dijkstra(g,A)
        #print(dists,paths)
        m = np.Infinity
        e = None
        for b in B:
            if b in dists and dists[b]<m+EPSILON:
                m = dists[b]
                e = paths[b]
        #print('\n')
        #print('A : ',A)
        #print('B : ',B)
        #print("data : ", {'weight':m,'path':e})
        #print('\n')
        if e==None:
            #print(dists[frozenset({3})])
            print(dists)
            print('A : ',A)
            print('B : ',B)
            print("data : ", {'weight':m,'path':e})
            print('\n')
        return {'weight':m,'path':e}
    cong = nx.quotient_graph(g,ccs+unused,edge_data=edge_data_min_weight)
    #print("cong",list(cong))
    #calculate MST
    t = nx.minimum_spanning_tree(cong)
    #print("MST",len(list(t)))
    #print(t.edges)
    #cut down to a tree with only the ccs as leaves
    if {start} in unused:
        unused.remove({start})
    #print("unused : ",unused)
    leaves = [x for x in t if t.degree(x)==1]
    """
    idx=0
    while idx < len(leaves):
        to_remove=[]
        for nleaf in t.neighbors(leaves[idx]):
            if nleaf in unused and t.degree(nleaf)==2:
                leaves.append(nleaf)
                to_remove.append(nleaf)
        for x in to_remove:
            t.remove_edge(leaves[idx],x)
        idx+=1
    """
    while leaves:
        leaf = leaves.pop()
        #print("leaf",leaf)
        if leaf in unused:
            nleaf = list(t.neighbors(leaf))[0]
            t.remove_node(leaf)
            if t.degree(nleaf)==1:
                leaves.append(nleaf)
    #print("pruned tree : {} nodes : ".format(len(list(t))),list(t))
    #print(t.edges)
    #add edges to connect rgraph
    for (x,y) in t.edges:
        #print(x,y,t[x][y])
        pxy = t[x][y]['path']
        for idx in range(0,len(pxy)-1):
            rgraph.add_edge(pxy[idx],pxy[idx+1],weight=g[pxy[idx]][pxy[idx+1]]['weight'])
    
    #print("{} rgraphnodes".format(len(list(rgraph))),list(rgraph))
    #print("{} rgraphedges".format(len(rgraph.edges)),rgraph.edges)
    
    
    #Now do cpp-like thingy
    #tour = cpp(rgraph,start)
    odds=[]
    for n in rgraph:
        if rgraph.degree(n)%2==1:
            odds.append(n)
    spd = nx.floyd_warshall(g)
    
    tomatch = nx.Graph()
    for idx,i in enumerate(odds):
        for jdx,j in enumerate(odds,idx+1):
            tomatch.add_edge(i,j,weight=-spd[i][j])
    m = nx.max_weight_matching(tomatch,True)
    multig = nx.MultiGraph(rgraph)
    for (i,j) in m:
        #recalculate shortest path from i to j,
        # add all edges in path to the graph
        pij = nx.shortest_path(g,i,j,weight='weight')
        for idx in range(0,len(pij)-1):
            multig.add_edge(pij[idx],pij[idx+1],weight=g[pij[idx]][pij[idx+1]])
    euler = nx.eulerian_circuit(multig,start)
    return euler


def cppaeourfysegf(g,start=None):
    #find odd-parity nodes
    odds=[]
    for n in g:
        if g.degree(n)%2==1:
            odds.append(n)
    #get shortest paths between them
    spd = nx.floyd_warshall(g)
    #print(spd)
    #construct new graph for the matching
    tomatch = nx.Graph()
    for idx,i in enumerate(odds):
        for jdx,j in enumerate(odds,idx+1):
            tomatch.add_edge(i,j,weight=-spd[i][j])
    #get matching
    m = nx.max_weight_matching(tomatch,True)
    #add edges back to g
    multig = nx.MultiGraph(g)    
    for (i,j) in m:
        #recalculate shortest path from i to j,
        # add all edges in path to the graph
        pij = nx.shortest_path(g,i,j,weight='weight')
        for idx in range(0,len(pij)-1):
            multig.add_edge(pij[idx],pij[idx+1],weight=g[pij[idx]][pij[idx+1]])
    #get euler tour
    euler = nx.eulerian_circuit(multig,start)
    return euler