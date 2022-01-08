# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 16:19:47 2020

@author: Malachi
"""


import matplotlib.pyplot as plt
import numpy as np

def getab(m,rho,alpha,beta):
    #aA+bB=rho
    #aC+BD=0
    #solve for a,b
    am = alpha**m
    bm = beta**m
    A=(am-1)/(alpha-1)-1
    B=(bm-1)/(beta-1)-1
    C=(alpha*am-1)/(alpha-1)-rho*alpha-1
    D=(beta*bm-1)/(beta-1)-rho*beta-1
    
    det=A*D-B*C
    a = rho*D/det
    b = -rho*C/det
    return a,b

def rootfinder(rho,m):
    prec=1e-10
    def p(t):
        return t**m-rho*t+rho
    def rooter(f,a,b):
        while b-a>prec:
            c = (a+b)/2
            fc = f(c)
            if fc*f(a)>=0:
                a=c
            else:
                b=c
        return c
    return rooter(p,1,m/(m-1)),rooter(p,m/(m-1),rho**(1/(m-1)))

def staraggressivestop(T,m,rho=None):
    if rho==None:
        rho = m*(m/(m-1))**(m-1)
    if abs(rho-m*(m/(m-1))**(m-1))>1e-7:
        p=np.zeros(m+1)
        p[0]=1
        p[-1]=rho
        p[-2]=-rho
        r = list(np.roots(p))
        
        #############3
        #
        # NEEDS MORE STABILITY
        #
        ################3
                
        r.sort()
        alpha=abs(r.pop())
        beta=abs(r.pop())
        
        beta,alpha = rootfinder(rho,m)
        
        a = beta*(alpha-1)/(alpha-beta)
        #print(a,b)
        x=[0]
        ai=1
        bi=1
        i=0
        while 2*sum(x[:-1])+x[-1]<T:
            i+=1
            ai*=alpha
            bi*=beta
            x.append((1+a)*ai-a*bi)
        
    else:
        r=m/(m-1)
        ri=1
        x=[0]
        i=0
        while 2*sum(x[:-1])+x[-1]<T:
            i+=1
            ri*=r
            x.append((m+i-1)/(m-1)*ri)
    #return(x)
    #print(x)
    x.pop()
    return (sum(x[-m:]))

def staraggressiveshrink(T,m,rho=None):
    if rho==None:
        rho = m*(m/(m-1))**(m-1)
    if abs(rho-m*(m/(m-1))**(m-1))>1e-7:
        p=np.zeros(m+1)
        p[0]=1
        p[-1]=rho
        p[-2]=-rho
        r = list(np.roots(p))
        r.sort()
        alpha=abs(r.pop())
        beta=abs(r.pop())
        #print(alpha,beta)
        
        beta,alpha = rootfinder(rho,m)
        
        a = beta*(alpha-1)/(alpha-beta)
        x=[0]
        ai=1
        bi=1
        i=0
        while 2*sum(x[:-1])+x[-1]<T:
            i+=1
            ai*=alpha
            bi*=beta
            x.append((1+a)*ai-a*bi)
        
    else:
        r=m/(m-1)
        ri=1
        x=[0]
        i=0
        while 2*sum(x[:-1])+x[-1]<T:
            i+=1
            ri*=r
            x.append((m+i-1)/(m-1)*ri)
    #return(x)
    #print(x)
    return (sum(x[-m:]))*T/(2*sum(x[:-1])+x[-1])

def staraggressiveoldstop(T,m,rho=None):
    if rho==None:
        rho = m*(m/(m-1))**(m-1)
    elif rho!=m*(m/(m-1))**(m-1):
        return "ERROR"
    
    a,b = getab(m)
    r=m/(m-1)
    ri=1
    x=[0]
    i=0
    while 2*sum(x[:-1])+x[-1]<T:
        i+=1
        ri*=r
        x.append(a*ri+b*i*ri)
    #return(x)
    #print(x)
    x.pop()
    return (sum(x[-m:]))

def staraggressiveoldshrink(T,m,rho=None):
    if rho==None:
        rho = m*(m/(m-1))**(m-1)
    elif rho!=m*(m/(m-1))**(m-1):
        return "ERROR"
    
    a,b = getab(m)
    r=m/(m-1)
    ri=1
    x=[0]
    i=0
    while 2*sum(x[:-1])+x[-1]<T:
        i+=1
        ri*=r
        x.append(a*ri+b*i*ri)
    #return(x)
    #print(x)
    return (sum(x[-m:]))*T/(2*sum(x[:-1])+x[-1])

def staraggrtime(L,m,rho=None):
    return max(staraggressivestop(L,m,rho),staraggressiveshrink(L,m,rho))


def stargeomshrink(T,m,rho=None):
    if rho==None:
        rho = m*(m/(m-1))**(m-1)
        r = m/(m-1)
    if abs(rho-m*(m/(m-1))**(m-1))>1e-7:
        p=np.zeros(m+1)
        p[0]=1
        p[-1]=rho
        p[-2]=-rho
        r = np.real(max(np.roots(p)))
        _,r = rootfinder(rho,m)
    else:
        r = m/(m-1)
    x=[0]*m+[r]
    t=r
    for i in range(m-1):
        t+=x[-1]
        x.append(x[-1]*r)
        t+=x[-1]
    while t<T:
        t+=x[-1]
        x.append(x[-1]*r)
        t+=x[-1]
    
    #print(x,t)
    #print([xi*L/(sum(x[-m:])) for xi in x])
    #return [xi*L/(sum(x[-m:])) for xi in x]
    return sum(x[-m:])*T/t+max(T-t,0)

def stargeomstop(T,m,rho=None):
    if rho==None:
        rho = m*(m/(m-1))**(m-1)
        r = m/(m-1)
    if abs(rho-m*(m/(m-1))**(m-1))>1e-7:
        _,r = rootfinder(rho,m)
    else:
        r = m/(m-1)
    if T<=r:
        return min(T,r)
    x=[0]*m+[r]
    t=r
    while t<T:
        #go back to the origin
        t+=x[-1]
        #print(x,t)
        if t+x[-m]>=T:
            #if we run out of time doing nothing interesting
            return sum(x[-m:])
        if t+x[-1]*r>=T:
            #if we run out of time before the next turn-around
            return sum(x[-m:])+T-t-x[-m]
        x.append(x[-1]*r)
        t+=x[-1]
    return t
    
def stargeomtime(T,m,rho=None):
    return stargeomshrink(T,m,rho)
    #return max(stargeomshrink(T,m,rho),stargeomstop(T,m,rho))

def starpolyk(T,m,k,rho=None):
    #FOR NOW I don't want to look at this because it seems like it's
    #going to need a lot more work to find a good value for k
    #and all that
    if rho==None:
        rho = m*(m/(m-1))**(m-1)
    x = [T/(2*m-1) for i in range(m)]
    while sum(x[:m-1])>rho:
        x.insert(0,x[0]-x[m-1]/rho)
    #print(x)
    #return x
    return (2*sum(x[:-1])+x[-1])

def starrhotopdown(T,m,rho=None):
    if rho==None:
        rho = m*(m/(m-1))**(m-1)
    low = 2
    high = int(m*np.log(T)/np.log(rho))+2
    hx = gausssolveexacttopdown(high,T,m,rho)
    while not hx:
        low=high
        high*=2
        hx = gausssolveexacttopdown(high,T,m,rho)
    #print(low,high)
    lx = gausssolveexacttopdown(low,T,m,rho)
    #print(lx)
    while (not lx>0) and high>low+1:
        c = (high+low)//2
        #print(c)
        cx = gausssolveexacttopdown(c,T,m,rho)
        if cx:
            high=c
            hx = cx
        else:
            low=c
            lx = cx
    #print(high)
    return hx

def starrho(T,m,rho=None):
    if rho==None:
        rho = m*(m/(m-1))**(m-1)
    low = 2
    high = int(m*np.log(T)/np.log(rho))+2
    hx = gausssolveexacttopdown(high,T,m,rho)
    while not hx:
        low=high
        high*=2
        hx = gausssolveexacttopdown(high,T,m,rho)
    #print(low,high)
    lx = gausssolveexacttopdown(low,T,m,rho)
    #print(lx)
    while (not lx>0) and high>low+1:
        c = (high+low)//2
        #print(c)
        cx = gausssolveexacttopdown(c,T,m,rho)
        if cx:
            high=c
            hx = cx
        else:
            low=c
            lx = cx
    #print(high)
    bu = gausssolveexactbottomup(high-1,T,m,rho)
    if bu > hx:
        return bu
    return hx


def solvetopdown(k,T,m,rho):
    mat = np.zeros((k,k))
    for i in range(k-2):
        mat[i][i] = rho
        mat[i][i+1] = -rho
    for i in range(k-m):
        mat[i][i+m] = 1
    mat[k-2] = np.ones(k)
    mat[k-2][k-2] = 1-rho
    mat[k-1][k-2] = 2*rho
    mat[k-1][k-1] = -1
    #print(mat)
    if np.linalg.det(mat)==0:
        print(mat,k,T,m,rho)
    inv = np.linalg.inv(mat)
    #print(inv)
    b = np.zeros(k)
    b[k-1] = T
    #print(mat)
    #print(b)
    ans=inv.dot(b)
    print(ans)
    #print(sum(ans[:m-1]))
    if sum(ans[:m-1])>rho:
        return False
    return sum(ans[-m:])
    
    
    return(ans[-2])
    print(ans)
    print([ans[i+1]/ans[i] for i in range(k-1)])
    print(ans[:m-1],sum(ans[:m-1]))
    #return ans
    return 2*sum(ans[:k-1])+ans[k-1]

def solvebottomup(k,T,m,rho,verbose=False):
    if k==1:
        if T>rho:
            return False
        return T
    mat = np.zeros((k,k))
    for i in range(k-2):
        mat[i][i] = rho
        mat[i][i+1] = -rho
    for i in range(k-m):
        mat[i][i+m] = 1
    mat[k-2] = np.ones(k)
    mat[k-2][k-2] = 1-rho
    for i in range(m-1):
        mat[k-1][i] = 1
    #print(mat)
    if np.linalg.det(mat)==0:
        print(mat,k,T,m,rho)
    inv = np.linalg.inv(mat)
    #print(inv)
    b = np.zeros(k)
    b[k-1] = rho
    #print(mat)
    #print(b)
    ans=inv.dot(b)
    #print(ans)
    
    if verbose:
        return ans
    
    if 2*sum(ans[:-1])+ans[-1] > T:
        return False
    return sum(ans[-m:])
    
    return(ans[-2])
    print(ans)
    print([ans[i+1]/ans[i] for i in range(k-1)])
    print(ans[:m-1],sum(ans[:m-1]))
    #return ans
    return 2*sum(ans[:k-1])+ans[k-1]

def gausssolvetopdown(k,T,m,rho):
    mat = np.zeros((k,k))
    for i in range(k-2):
        mat[i][i] = rho
        mat[i][i+1] = -rho
    for i in range(k-m):
        mat[i][i+m] = 1
    mat[k-1] = np.ones(k)
    mat[k-1][k-2] = 1-rho
    mat[k-2][k-2] = 2*rho
    mat[k-2][k-1] = -1
    #print(mat)
    if np.linalg.det(mat)==0:
        print(mat,k,T,m,rho)
    for i in range(k-m):
        factor = mat[k-1][i]/mat[i][i]
        mat[k-1][i]-=factor*mat[i][i]
        mat[k-1][i+1]-=factor*mat[i][i+1]
        mat[k-1][i+m]-=factor*mat[i][i+m]
        #print(mat)
    for i in range(k-m,k-1):
        factor = mat[k-1][i]/mat[i][i]
        mat[k-1][i]-=factor*mat[i][i]
        mat[k-1][i+1]-=factor*mat[i][i+1]
        #print(mat)
    ans = np.zeros(k)
    ans[k-1] = -factor*T/mat[k-1][k-1]
    ans[k-2] = (T-ans[k-1]*mat[k-2][k-1])/mat[k-2][k-2]
    
    for i in range(k-3,k-m-1,-1):
        ans[i] = -ans[i+1]*mat[i][i+1]/mat[i][i]
    
    for i in range(k-m-1,-1,-1):
        ans[i] = -(ans[i+m]*mat[i][i+m]+ans[i+1]*mat[i][i+1])/mat[i][i]
    
    print(ans)
    #print(sum(ans[:m-1]))
    if sum(ans[:m-1])>rho:
        return False
    return sum(ans[-m:])
    
    
    return(ans[-2])
    print(ans)
    print([ans[i+1]/ans[i] for i in range(k-1)])
    print(ans[:m-1],sum(ans[:m-1]))
    #return ans
    return 2*sum(ans[:k-1])+ans[k-1]

def gausssolveexacttopdown(k,T,m,rho):
    
    lastline = np.ones(k)
    lastline[k-2] = 1-rho
    
    for i in range(k-m):
        lastline[i+1]+=lastline[i]
        lastline[i+m]-=lastline[i]/rho
    for i in range(max(k-m,0),k-2):
        lastline[i+1]+=lastline[i]
    lastline[k-1] += lastline[k-2]/(2*rho)
    
    #print(lastline)
    
    ans = np.zeros(k)
    ans[k-1] = -lastline[k-2]*T/(lastline[k-1]*2*rho)
    ans[k-2] = (T+ans[k-1])/(2*rho)
    for i in range(k-3,max(k-m-1,-1),-1):
        ans[i] = ans[i+1]
    for i in range(k-m-1,-1,-1):
        ans[i] = ans[i+1]-ans[i+m]/rho
    
    #print(ans)
    #print(sum(ans[:m-1]))
    if sum(ans[:m-1])>rho:
        return False
    return sum(ans[-m:])
    
    
    return(ans[-2])
    print(ans)
    print([ans[i+1]/ans[i] for i in range(k-1)])
    print(ans[:m-1],sum(ans[:m-1]))
    #return ans
    return 2*sum(ans[:k-1])+ans[k-1]

def gausssolveexactbottomup(k,T,m,rho):
    
    lastline = np.ones(k)
    lastline[k-2] = 1-rho
    
    for i in range(k-m):
        lastline[i+1]+=lastline[i]
        lastline[i+m]-=lastline[i]/rho
    for i in range(max(k-m,0),k-2):
        lastline[i+1]+=lastline[i]
    lastline[k-1] += lastline[k-2]/(2*rho)
    
    #print(lastline)
    
    ans = np.zeros(k)
    ans[k-1] = -lastline[k-2]*T/(lastline[k-1]*2*rho)
    ans[k-2] = (T+ans[k-1])/(2*rho)
    for i in range(k-3,max(k-m-1,-1),-1):
        ans[i] = ans[i+1]
    for i in range(k-m-1,-1,-1):
        ans[i] = ans[i+1]-ans[i+m]/rho
    
    #print(ans)
    #print(sum(ans[:m-1]))
    
    factor = rho/sum(ans[:m-1])
    ans = [ansi*factor for ansi in ans]
    
    if 2*sum(ans[:k-1])+ans[k-1]>T:
        return False
    return sum(ans[-m:])
    
    
    return(ans[-2])
    print(ans)
    print([ans[i+1]/ans[i] for i in range(k-1)])
    print(ans[:m-1],sum(ans[:m-1]))
    #return ans
    return 2*sum(ans[:k-1])+ans[k-1]


def testrho(x,m,rhofactor,npoints):
    
    rhom = m**m/(m-1)**(m-1)
    rhos = np.linspace(rhom,rhofactor*rhom,npoints)
    yfact = [stargeomtime(x,m,rho) for rho in rhos]
    yfact2 = [starrho(x,m,rho) for rho in rhos]
    yfact3 = [staraggrtime(x,m,rho) for rho in rhos]
    
    plt.plot(rhos*2+1,yfact2,label="Optimal strategy")
    plt.plot(rhos*2+1,yfact3,label="Mixed aggressive")
    plt.plot(rhos*2+1,yfact,label="Scaled geometric")
    plt.xlabel("R")
    plt.ylabel("Clearance")

    plt.legend()



def testfactor(x0,xn,m,rhofactor,npoints):
    rho = m**m/(m-1)**(m-1)*rhofactor
    
    xs = np.logspace(x0,xn,npoints)
    
    yfact = [starrho(x,m,rho)/stargeomtime(x,m,rho) for x in xs]
    yfact2 = [starrho(x,m,rho)/staraggrtime(x,m,rho) for x in xs]
    yfact3 = [staraggrtime(x,m,rho)/stargeomtime(x,m,rho) for x in xs]
    
    plt.plot(xs,yfact,label="Optimal/scaled geometric")
    plt.plot(xs,yfact2,label="Optimal/mixed aggressive")
    plt.plot(xs,yfact3,label="Mixed agressive/scaled geometric")
    
    plt.xlabel("T")
    plt.ylabel("Clearance ratio")
    plt.xscale('log')
    plt.legend()

def testx(x0,xn,m,rhofactor,npoints):
    xs = np.logspace(x0,xn,npoints)
    rho = m**m/(m-1)**(m-1)*rhofactor
    
    ygeom=[stargeomtime(x,m,rho) for x in xs]
    yrho=[starrho(x,m,rho) for x in xs]
    yaggr=[staraggrtime(x,m,rho) for x in xs]
    
    plt.plot(xs,ygeom,label="Geometric strategy")
    plt.plot(xs,yaggr,label="Aggressive strategy")
    plt.plot(xs,yrho,label="Optimal strategy")
    
    plt.xlabel("T")
    plt.legend()

def testm(x,m0,mn,rhofactor,npoints):
    ms = list(range(m0,mn))
    ygeom=[stargeomtime(x,m,rhofactor*m**m/(m-1)**(m-1)) for m in ms]
    yrho=[starrho(x,m,rhofactor*m**m/(m-1)**(m-1)) for m in ms]
    yaggr=[staraggrtime(x,m,rhofactor*m**m/(m-1)**(m-1)) for m in ms]
    
    plt.plot(ms,yrho,label="Optimal strategy")
    plt.plot(ms,yaggr,label="Aggressive strategy")
    plt.plot(ms,ygeom,label="Geometric strategy")

    plt.xlabel("m")
    plt.ylabel("Clearance")
    plt.legend()

def gain(m,x,rhofact=1):
    rho = m**m/(m-1)**(m-1)*rhofact
    return starrho(x,m,rho)/staraggrtime(x,m,rho)


def test(testid,filename=None,x=1e5,x0=1,xn=15,m=4,m0=3,mn=20,rhofactor=1,npoints=1000):
    if testid=='factor':
        testfactor(x0,xn,m,rhofactor,npoints)
    elif testid=='x':
        testx(x0,xn,m,rhofactor,npoints)
    elif testid=='rho':
        testrho(x,m,rhofactor,npoints)
    elif testid=='m':
        testm(x,m0,mn,rhofactor,npoints)
    elif testid=='table':
        for m in [3,4,5,10,20,50,100]:
            for rhofactor in [1,2,5,10]:
                print(gain(m,1e16,rhofactor))
        return None
    else:
        return "Not recognized test type"
    if filename:
        plt.savefig("../figs/star/"+filename,dpi=300)
    plt.show()
    return None

test('factor',filename='factor.png')
test('m',filename='m.png')
test('rho',filename='rho.png',rhofactor=3)
test('table')