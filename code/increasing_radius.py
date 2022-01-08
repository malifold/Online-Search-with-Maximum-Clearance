# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 09:12:33 2020

@author: Malachi
"""

EPSILON=1e-9
SCALE = 0.25

import networkx as nx
import numpy as np
from cpp import cpp as cpp
from rpp import rpp as rpp

def getrad(g,c,r):
    #get ball of radius r with center c in graph g
    #color g as unvisited, distance to c large
    for x in g:
        g.nodes[x]['visit']=False
        g.nodes[x]['distc']=np.Infinity
        #print(x)
    g.nodes[c]['distc']=0
    #do the BFS
    toconsider = [c]
    weirdedges=[]
    normaledges=[]
    while toconsider:
        x = toconsider.pop(0)
        #print("considering",x)
        if not g.nodes[x]['visit']:  
            for y in g.neighbors(x):
                #print('y is ',y)
                if not g.nodes[y]['visit']:
                    if g[x][y]['weight']+g.nodes[x]['distc']<=r+EPSILON:
                        toconsider.append(y)
                        normaledges.append((x,y))
                        if g.nodes[y]['distc']>g.nodes[x]['distc']+g[x][y]['weight']+EPSILON:
                            g.nodes[y]['distc'] = g.nodes[x]['distc']+g[x][y]['weight']
                    else:
                        weirdedges.append((x,y))
            g.nodes[x]['visit']=True
    #get the ball
    ball = nx.Graph()
    while weirdedges:
        u,v = weirdedges.pop()
        if g.has_edge(u,v):
            if g.nodes[u]['visit'] and g.nodes[u]['distc']+EPSILON<r:
                if g.nodes[v]['visit'] and g.nodes[v]['distc']+EPSILON<r:
                    if g.nodes[u]['distc']+g.nodes[v]['distc']+g[u][v]['weight']>2*r+EPSILON:
                        g.add_edge(u,(frozenset({u,v}),'left'),weight=r-g.nodes[u]['distc'])
                        g.add_edge(v,(frozenset({u,v}),'right'),weight=r-g.nodes[v]['distc'])
                        g.add_edge((frozenset({u,v}),'left'),(frozenset({u,v}),'right'),
                                   weight=g[u][v]['weight']+g.nodes[u]['distc']+g.nodes[v]['distc']-2*r)
                        g.remove_edge(u,v)
                        normaledges.append((u,(frozenset({u,v}),'left')))
                        normaledges.append((v,(frozenset({u,v}),'right')))
                        g.nodes[(frozenset({u,v}),'right')]['visit']=True
                        g.nodes[(frozenset({u,v}),'right')]['distc']=r
                        g.nodes[(frozenset({u,v}),'left')]['visit']=True
                        g.nodes[(frozenset({u,v}),'left')]['distc']=r
                    else:
                        normaledges.append((u,v))
                else:
                    g.add_edge(u,(frozenset({u,v}),'center'),weight=r-g.nodes[u]['distc'])
                    g.add_edge(v,(frozenset({u,v}),'center'),
                               weight=g[u][v]['weight']+g.nodes[u]['distc']-r)
                    g.remove_edge(u,v)
                    normaledges.append((u,(frozenset({u,v}),'center')))
                    g.nodes[(frozenset({u,v}),'center')]['visit']=True
                    g.nodes[(frozenset({u,v}),'center')]['distc']=r
            elif g.nodes[v]['distc']+EPSILON<r:
                g.add_edge(v,(frozenset({v,u}),'center'),weight=r-g.nodes[v]['distc'])
                g.add_edge(u,(frozenset({v,u}),'center'),
                           weight=g[v][u]['weight']+g.nodes[v]['distc']-r)
                g.remove_edge(v,u)
                g.nodes[(frozenset({v,u}),'center')]['visit']=True
                g.nodes[(frozenset({v,u}),'center')]['distc']=r
                normaledges.append((v,(frozenset({v,u}),'center')))
    for (u,v) in normaledges:
        #print(u,v)
        ball.add_edge(u,v,weight=g[u][v]['weight'])
    #cleanup
    for x in g:
        del g.nodes[x]['visit']
        del g.nodes[x]['distc']
    return ball

def cppsearching(g,c,r=2):
    strategy=[]
    currradius=r*SCALE
    ball = getrad(g,c,currradius)
    strategy.append(cpp(ball,c))
    while len(ball.edges)<len(g.edges):
        currradius*=r
        ball = getrad(g,c,currradius)
        #print(list(ball),r)
        strategy.append(cpp(ball,c))
    return strategy

def rppsearching(g,c,r=2):
    strategy=[]
    currradius=r*SCALE
    currroot = c
    lastball = nx.Graph()
    ball = getrad(g,c,currradius)
    strategy.append(clean(list(rpp(ball,[e for e in ball.edges if not e in lastball.edges],currroot)),nx.Graph()))
    _,currroot = strategy[-1][-1]
    while len(ball.edges)<len(g.edges):
        currradius*=r
        lastball=ball
        ball = getrad(g,c,currradius)
        for u,v in ball.edges:
            if ball[u][v]["weight"]<=0:
                print(u,v,ball[u][v]["weight"])
        #print("ball is : ",list(ball),"radius is : ",currradius)
        strategy.append(clean(list(rpp(ball,[e for e in ball.edges if not e in lastball.edges],currroot)),lastball))
        _,currroot = strategy[-1][-1]
    return strategy

def strategylength(g,s):
    return sum(g[u][v]['weight'] for (u,v) in s)

def discoverylength(g,s,l,discovered=None,d=0,t=0):
    if discovered == None:
        discovered = set()
    #calculates time to discover length l using strategy s in graph g
    for (u,v) in s:
        if not frozenset({u,v}) in discovered:
            if g[u][v]['weight'] + d >l:
                return True,None,None,t+l-d
            else:
                discovered.add(frozenset({u,v}))
                d+= g[u][v]['weight']
        t+= g[u][v]['weight']
    return False,discovered,d,t

def discoverytime(g,s,T,discovered=None,d=0,t=0):
    if discovered == None:
        discovered = set()
    #calculates time to discover length l using strategy s in graph g
    for (u,v) in s:
        if g[u][v]['weight'] + t >T:
            if frozenset({u,v}) in discovered:
                return True,None,d,None
            else:
                return True,None,d+T-t,None
        else:
            if not frozenset({u,v}) in discovered:
                discovered.add(frozenset({u,v}))
                d+= g[u][v]['weight']
        t+= g[u][v]['weight']
    return False,discovered,d,t

def cppl(g,c,r,l):
    #search graph g from center c using radius r,
    #until a length l is discovered
    #returns the time taken
    strategy=[]
    currradius=r
    ball = getrad(g,c,currradius)
    strategy.append(list(cpp(ball,c)))
    b,discovered,d,t = discoverylength(g,strategy[-1],l)
    if b:
        return t
    while len(ball.edges)<len(g.edges):
        currradius*=r
        ball = getrad(g,c,currradius)
        #print(list(ball),r)
        strategy.append(list(cpp(ball,c)))
        b,discovered,d,t = discoverylength(g,strategy[-1],l,discovered,d,t)
        if b:
            return t
    return False

def rppl(g,c,r,l):
    #search graph g from center c using radius r,
    #until a length l is discovered
    #returns the time taken
    strategy=[]
    currradius=r
    lastball = nx.Graph()
    ball = getrad(g,c,currradius)
    strategy.append(rpp(ball,[e for e in ball.edges if not e in lastball.edges],c))
    b,discovered,d,t = discoverylength(g,strategy[-1],l)
    if b:
        return t
    while len(ball.edges)<len(g.edges):
        #print(b,d,t)
        currradius*=r
        lastball=ball
        ball = getrad(g,c,currradius)
        #print(list(ball),r)
        strategy.append(rpp(ball,[e for e in ball.edges if not e in lastball.edges],c))
        b,discovered,d,t = discoverylength(g,strategy[-1],l,discovered,d,t)
        if b:
            return t
    return False

def rppt(g,c,r,T):
    #search graph g from center c using radius r,
    #until a length l is discovered
    #returns the time taken
    strategy=[]
    currradius=r
    lastball = nx.Graph()
    ball = getrad(g,c,currradius)
    strategy.append(rpp(ball,[e for e in ball.edges if not e in lastball.edges],c))
    b,discovered,d,t = discoverytime(g,strategy[-1],T)
    if b:
        return d
    while len(ball.edges)<len(g.edges):
        #print(b,d,t)
        currradius*=r
        #print("radius is " currradius)
        lastball=ball
        ball = getrad(g,c,currradius)
        #print(list(ball),r)
        strategy.append(rpp(ball,[e for e in ball.edges if not e in lastball.edges],c))
        b,discovered,d,t = discoverytime(g,strategy[-1],T,discovered,d,t)
        if b:
            return d
    return d

def cppt(g,c,r,T):
    #search graph g from center c using radius r,
    #until a length l is discovered
    #returns the time taken
    strategy=[]
    currradius=r
    ball = getrad(g,c,currradius)
    strategy.append(list(cpp(ball,c)))
    b,discovered,d,t = discoverytime(g,strategy[-1],T)
    if b:
        return d
    while len(ball.edges)<len(g.edges):
        currradius*=r
        ball = getrad(g,c,currradius)
        #print(list(ball),r)
        strategy.append(list(cpp(ball,c)))
        b,discovered,d,t = discoverytime(g,strategy[-1],T,discovered,d,t)
        if b:
            return d
    return d

def clean(strategy,lastball):
    lastidx = 0
    visited={}
    for u,v in lastball.edges:
        visited[frozenset({u,v})]=True
    for i,(u,v) in enumerate(strategy):
        if not frozenset({u,v}) in visited:
            lastidx=i
            visited[frozenset({u,v})]=True
    return strategy[:lastidx+1]
            