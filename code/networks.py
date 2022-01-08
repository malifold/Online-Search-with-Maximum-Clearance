# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 15:01:19 2020

@author: Malachi
"""


import copy
import matplotlib.pyplot as plt
import networkx as nx
import random
import increasing_radius
import numpy as np
import cityparser
import time
import pickle

cities = ["SiouxFalls","EMA","friedrichshain-center","berlin","ChicagoSketch"]
poss = [True,False,True,True,True]

def perf(g,strategy):
    l=0
    t=0
    ts=[]
    ls=[]
    for u,v in g.edges:
        g[u][v]['visited']=False
    for u,v in strategy:
        if not g[u][v]['visited']:
            g[u][v]['visited']=True
            l+=g[u][v]['weight']
        t+=g[u][v]['weight']
        ts.append(t)
        ls.append(l)
    return ls,ts

def calcrhos(g,strategy):
    rhos={}
    root,_ = strategy[0]
    dists = nx.shortest_path_length(g,root,weight='weight')
    bonus = list(g.neighbors(root))
    bonusrho={}
    for x in bonus:
        bonusrho[x] = 0
    
    for u,v in g.edges:
        g[u][v]['rhovisited']=False
    
    t=0
    ts=[]
    rho=0
    rholist=[]
    for u,v in strategy:
        #deal with bonus
        if not g[u][v]['rhovisited']:
            g[u][v]['rhovisited']=True
            if dists[u] != 0:
                rhos[u] = t/dists[u]
                rho = max(rho,t/dists[u])
            t+=g[u][v]['weight']
            if dists[v] != 0:
                rhos[v] = t/dists[v]
                rho = max(rho,t/dists[v])
        else:
            t+=g[u][v]['weight']
        rholist.append(rho)
        ts.append(t)
    return rholist,ts

def clearance(g,s,T):
    l=0
    t=0
    for u,v in g.edges:
        g[u][v]['visited']=False
    for u,v in s:
        t+=g[u][v]['weight']
        if t>T:
            return l
        if not g[u][v]['visited']:
            g[u][v]['visited']=True
            l+=g[u][v]['weight']
    return l

def timetoclear(g,s,L):
    l=0
    t=0
    for u,v in g.edges:
        g[u][v]['visited']=False
    for u,v in s:
        if l>=L:
            return t
        if not g[u][v]['visited']:
            g[u][v]['visited']=True
            l+=g[u][v]['weight']
        t+=g[u][v]['weight']
    return np.Infinity

def fixs(sl):
    s=[]
    for sx in sl:
        s+=sx
    return s


def dummycity(idx):
    g,pos = cityparser.getgraph(cities[idx],poss[idx])
    startnode = random.randint(1,g.graph["nnodes"])
    while not g.has_node(startnode):
        startnode = random.randint(1,g.graph["nnodes"])
    return startnode

def docity(idx, draw=True, r=[2], startnode='random'):
    
    g,pos = cityparser.getgraph(cities[idx],poss[idx])
    if draw:
        if pos==None:
            nx.draw(g)
        else:
            nx.draw(g,pos,node_size=30)
        print("{} nodes and {} edges".format(len(g.nodes),len(g.edges)))
    
        #plt.savefig("citygraphs/{}pic.png".format(cities[idx]))
        plt.show()
    
    if startnode=='random':
        startnode = random.randint(1,g.graph["nnodes"])
        while not g.has_node(startnode):
            startnode = random.randint(1,g.graph["nnodes"])   
    
    rhos1=[]
    rhos2 = []
    finalts1=[]
    finalts2 = []
    
    def plotstuff(scpp,srpp):
        ls,ts = perf(g2,srpp)
        plt.plot(ts,ls,label="RPT(r)")
        ls,ts = perf(g1,scpp)
        plt.plot(ts,ls,label="CPT(r)")
        
        plt.ylabel("Clearance")
        plt.xlabel("Time")
        plt.legend()
        plt.savefig("..\\figs\\cities\\{}_{}_{}.png".format(cities[idx],startnode,int(r[-1]*1000)),dpi=300)
        plt.show()
    
    def factorplot(scpp,srpp):
        ls1,ts1 = perf(g1,scpp)
        ls2,ts2 = perf(g2,srpp)
        ts = np.linspace(0,ts2[-1])
        def getti(t,tlist):
            lowi=0
            highi=len(tlist)-1
            while lowi<highi-1:
                ci = (lowi+highi)//2
                if tlist[ci]>t:
                    highi=ci
                else:
                    lowi=ci
            return lowi
        ratios = [ls2[getti(t,ts2)]/ls1[getti(t,ts1)] for t in ts]
        plt.plot(ts,ratios)
        plt.ylabel("Advantage")
        plt.xlabel("Time")
        plt.savefig("..\\figs\\cities\\{}_{}_{}factor.png".format(cities[idx],startnode,int(r[-1]*1000)),dpi=300)
        plt.show()
    
    print("startnode = {}".format(startnode))
    for rx in r:
        starttime = time.time()
        g1 = copy.deepcopy(g)
        scpp = fixs(increasing_radius.cppsearching(g1,startnode,rx))
        print("cppdone in {}".format(int(time.time()-starttime)))
        with open("..\\figs\\cities\\{}_{}_{}_scpp.txt".format(cities[idx],startnode,int(r[-1]*1000)),"wb") as fp:
            pickle.dump(scpp,fp)
        starttime2 = time.time()
        g2 = copy.deepcopy(g)
        srpp = fixs(increasing_radius.rppsearching(g2,startnode,rx))
        print("rppdone in {}".format(int(time.time()-starttime2)))
        with open("..\\figs\\cities\\{}_{}_{}_srpp.txt".format(cities[idx],startnode,int(r[-1]*1000)),"wb") as fp:
            pickle.dump(srpp,fp)
        
        plotstuff(scpp,srpp)
        factorplot(scpp,srpp)
       
    
        r1,t1 = calcrhos(g1,scpp)
        r2,t2 = calcrhos(g2,srpp)
        rhos1.append(r1[-1])
        rhos2.append(r2[-1])
        finalts1.append(t1[-1])
        finalts2.append(t2[-1])
    
    if len(r)>1:
        plt.plot(r,rhos2,label="RPT(r)")
        plt.plot(r,rhos1,label="CPT(r)")
        plt.xlabel("r")
        plt.ylabel("Competitive ratio")
        plt.legend()
        plt.savefig("..\\figs\\cities\\compr.png",dpi=300)
        plt.show()
        plt.plot(r,finalts2,label="RPT(r)")
        plt.plot(r,finalts1,label="CPT(r)")
        plt.xlabel("r")
        plt.ylabel("T")
        plt.legend()
        plt.savefig("..\\figs\\cities\\comprt.png",dpi=300)
        plt.show()
    
    return None
    

def testruns():
    random.seed(42)
    for i in range(5):
        docity(i)


def berlinruns():
    random.seed(42)
    for i in range(10):
        docity(3,draw=False)

def chicagoruns():
    random.seed(42)
    for i in range(49):
        docity(4,draw=False)

def ballify(g,c,r=2):
    currradius = r*0.25
    ball = increasing_radius.getrad(g,c,currradius)
    while len(ball.edges)<len(g.edges):
        currradius*=r
        ball = increasing_radius.getrad(g,c,currradius)
    return None

def clear(T,city):
    if city=='berlin':
        gmaster,pos = cityparser.getgraph("berlin",False)
        starts = [655,115,26,760,282,251,229,755,914,90]

    elif city=='ChicagoSketch':
        gmaster,pos = cityparser.getgraph("ChicagoSketch",False)
        starts = [7,26,28,31,33,90,96,105,115,143,160,164,
                  204,221,224,226,229,239,251,282,285,349,
                  430,433,460,518,559,575,604,605,617,655,
                  666,693,715,719,734,755,759,760,778,826,
                  829,891,914]
    else:
        return "I don't know this city"
    nruns=len(starts)
    clearcpps=[]
    clearrpps=[]
    for x in starts:
        g = copy.deepcopy(gmaster)
        ballify(g,x)
        with open("..\\figs\\cities\\{}_{}_2000_scpp.txt".format(city,x),"rb") as fp:
            scpp = pickle.load(fp)
        clearcpps.append(clearance(g,scpp,T))
        with open("..\\figs\\cities\\{}_{}_2000_srpp.txt".format(city,x),"rb") as fp:
            srpp = pickle.load(fp)
        clearrpps.append(clearance(g,srpp,T))
    plt.plot(range(nruns),clearrpps,label="RPT(r)")
    plt.plot(range(nruns),clearcpps,label="CPT(r)")
    plt.legend()
    plt.show()
    return np.mean([clearrpps[i]/clearcpps[i] for i in range(nruns)])

def timer(L,city):
    if city=='berlin':
        gmaster,pos = cityparser.getgraph("berlin",False)
        starts = [655,115,26,760,282,251,229,755,914,90]

    elif city=='ChicagoSketch':
        gmaster,pos = cityparser.getgraph("ChicagoSketch",False)
        starts = [7,26,28,31,33,90,96,105,115,143,160,164,
                  204,221,224,226,229,239,251,282,285,349,
                  430,433,460,518,559,575,604,605,617,655,
                  666,693,715,719,734,755,759,760,778,826,
                  829,891,914]
    else:
        return "I don't know this city"
    nruns=len(starts)
    clearcpps=[]
    clearrpps=[]
    for x in starts:
        g = copy.deepcopy(gmaster)
        ballify(g,x)
        with open("..\\figs\\cities\\{}_{}_2000_scpp.txt".format(city,x),"rb") as fp:
            scpp = pickle.load(fp)
        clearcpps.append(timetoclear(g,scpp,L))
        with open("..\\figs\\cities\\{}_{}_2000_srpp.txt".format(city,x),"rb") as fp:
            srpp = pickle.load(fp)
        clearrpps.append(timetoclear(g,srpp,L))
    plt.plot(range(nruns),clearrpps,label="RPT(r)")
    plt.plot(range(nruns),clearcpps,label="CPT(r)")
    plt.legend()
    plt.show()
    return np.mean([clearcpps[i]/clearrpps[i] for i in range(nruns)])

def compt(city):
    if city=='berlin':
        maxb=250000
    elif city=='ChicagoSketch':
        maxb=100000
    else:
        return "I don't know this city"
    budgets = np.linspace(10,maxb)
    cleargains = [clear(b,city) for b in budgets]
    plt.plot(budgets[1:],cleargains[1:])
    plt.xlabel("T")
    plt.ylabel("Clearance ratio")
    plt.savefig("..\\figs\\cities\\compall_{}.png".format(city),dpi=300)
    plt.show()

def compl(city):
    g,pos = cityparser.getgraph(city,False)
    sizecity=0
    for u,v in g.edges:
        sizecity+=g[u][v]["weight"]
    budgets = np.linspace(10,sizecity)
    cleargains = [timer(b,city) for b in budgets]
    plt.plot(budgets[1:-1],cleargains[1:-1])
    plt.xlabel("L")
    plt.ylabel("Avg-time ratios for clearing L")
    plt.savefig("..\\figs\\cities\\earliestclearance_{}.png".format(city),dpi=300)
    plt.show()



def rhos(city):
    if city=='berlin':
        gmaster,pos = cityparser.getgraph("berlin",False)
        starts = [655,115,26,760,282,251,229,755,914,90]

    elif city=='ChicagoSketch':
        gmaster,pos = cityparser.getgraph("ChicagoSketch",False)
        starts = [7,26,28,31,33,90,96,105,115,143,160,164,
                  204,221,224,226,229,239,251,282,285,349,
                  430,433,460,518,559,575,604,605,617,655,
                  666,693,715,719,734,755,759,760,778,826,
                  829,891,914]
    else:
        return "I don't know this city"
    nruns=len(starts)
    rhocpps=[]
    rhorpps=[]
    for x in starts:
        g = copy.deepcopy(gmaster)
        ballify(g,x)
        print(x)
        with open("..\\figs\\cities\\{}_{}_2000_scpp.txt".format(city,x),"rb") as fp:
            scpp = pickle.load(fp)
        rhocpps.append(calcrhos(g,scpp)[0][-1])
        with open("..\\figs\\cities\\{}_{}_2000_srpp.txt".format(city,x),"rb") as fp:
            srpp = pickle.load(fp)
        rhorpps.append(calcrhos(g,srpp)[0][-1])
    plt.plot(range(nruns),rhorpps,label="Rural")
    plt.plot(range(nruns),rhocpps,label="CPT")
    plt.legend()
    plt.show()
    return rhocpps,rhorpps
def compr(city):
    rhocpps,rhorpps = rhos(city)
    nruns = len(rhocpps)
    allrhos = [(rhorpps[i],rhocpps[i]) for i in range(nruns)]
    allrhos.sort()
    rhocpps=[]
    rhorpps=[]
    for rr,rc in allrhos:
        rhocpps.append(rc)
        rhorpps.append(rr)
    plt.scatter(range(nruns),rhorpps,label="RPT(r)")
    plt.scatter(range(nruns),rhocpps,label="CPT(r)")
    plt.ylabel("Competitive factor")
    plt.legend()
    plt.savefig("..\\figs\\cities\\r{}.png".format(city),dpi=300)
    plt.show()

def competitiveness():
    r=np.linspace(1.2,4.5,150)
    startnode=10 
    docity(0,False,r,startnode)
