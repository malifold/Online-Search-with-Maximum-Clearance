# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 11:42:39 2020

@author: Malachi
"""

import networkx as nx
import os

def getgraph(city,positioning=True):
    g = nx.Graph()
    #takes city as the name of the file
    file = open("..\\data\\"+city+"_net.tntp.txt")
    file.readline() #nzones
    nnodes = int(file.readline().split()[3]) #nnodes
    ftn = int(file.readline().split()[3]) #ftn
    nedges = int(file.readline().split()[3])
    file.readline() #oh
    file.readline() #eom
    file.readline()
    file.readline()
    file.readline() #header
    
    d = {}
    g.graph["ftn"]=ftn
    g.graph["nnodes"]=nnodes
    
    
    minw = 1e10
    for i in range(nedges):
        line = file.readline().split()
        #print(line)
        x=int(line[0])
        y=int(line[1])
        w = float(line[3])
        if x<ftn:
            d[y]=x
        if y>=ftn:
            if x in d:
                x = d[x]
            if y in d:
                y = d[y]
            if w>0:
                g.add_edge(x,y,weight=w)
                minw = min(minw,w)
    for u,v in g.edges:
        g[u][v]["weight"]/=minw
    
    pos=None
    if positioning:
        pos={}
        file = open("..\\data\\"+city+"_node.tntp.txt")
        file.readline()
        for i in range(nnodes):
            line = file.readline().split()
            #print(line)
            pos[int(line[0])] = (float(line[1]),float(line[2]))
    return g,pos
    