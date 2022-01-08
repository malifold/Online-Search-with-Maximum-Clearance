# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 12:31:40 2020

@author: Malachi
"""

import networkx as nx
import numpy as np

def cpp(g,start=None):
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
    
    for u,v in g.edges:
        if g[u][v]['weight']<0:
            print(u,v,g[u][v])
    
    for (i,j) in m:
        #recalculate shortest path from i to j,
        # add all edges in path to the graph
        pij = nx.shortest_path(g,i,j,weight='weight')
        for idx in range(0,len(pij)-1):
            multig.add_edge(pij[idx],pij[idx+1],weight=g[pij[idx]][pij[idx+1]])
    #get euler tour
    euler = nx.eulerian_circuit(multig,start)
    return euler