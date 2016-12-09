#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 00:21:36 2016

@author: thaut
"""

"""
import the packages we need next"""

from multiprocessing import Pool
from functools import partial
import gc  
import numpy as np
import copy
from matplotlib import pyplot as plt
rnd=np.random.random
try:
    pool
except:
    pool=Pool(4)
def dot(x,I,f1,f2,w=None):  
    """inner product between primary function1 with primary function2
    """
    if type(w)==type(None):
        w=np.ones((x.shape[0],))
    return np.sum(w*f1(x,I)*f2(x,I))

def dot_y(x,y,fi,I,w=None):    
    """inner product between primary function_i with the function to fit
        return a scalar"""
    if type(w)==type(None):
        w=np.ones((x.shape[0],))
    return np.sum(w*fi(x,I)*y)
def dot_matrix(x,f,I,w=None):
    """ inner product between primary function_i with primary function_j 
        i=1,2,3,...,n,j=1,2,3,...,n
        return a matrix"""
    n=len(f)
    frow=lambda j:[dot(x,I,f[j],f[i],w) for i in range(n)]
    result=np.array(list(map(frow,range(n))))
    return result
def get_A(x,y,f,I,w=None):
    """ get the parameters of all the primary functions 
        return a vector"""
    n=len(f)
    B=np.array([dot_y(x,y,f[i],I,w) for i in range(n)])
    M=dot_matrix(x,f,I,w)
    A=np.dot(B,np.linalg.inv(M))
    return A
def plot_A(x,y,f,I,A,mode=True):
    """get the score of the functions which perfrom at fitting the function y=F(x)""" 
    """f: a (numpy)array,each element is a primary function"""
    """if use 'mode = True', the Curves of testing function 
    ( y=F(x) )  and fitted function 
    ( y=f(x) )
    will be plotted """  
    n=len(f)
    fall=lambda x:np.sum([A[i]*f[i](x,I) for i in range(n)],axis=0)
    y_pre=fall(x)
    if mode:
        plt.title('compare')
        plt.plot(x,y,'red')
        plt.plot(x,y_pre,'blue')
    return np.sum(np.abs(y_pre-y))
def score(x,y,f,I,w=None):
    """package the function get_A and plot_A(mode=False)  """
    A=get_A(x,y,f,I,w)
    scores=plot_A(x,y,f,I,A,mode=False)
    return scores


    
"""Genetic Algorithm"""

""" support functions """
def makeC():
    def f(t,I):
        return t**0
    return f

def makeL(a,b):
    def f(t,I):
        return np.log(np.abs(a-I*t)/np.abs(b+I*t))
        
    return f
def GD(X,y,I,Group,delta=0.001,epoch=10):
            """Calculate the gredient of function with the current parameters"""
            deltas=np.zeros((Group.shape[0],))
            fs=np.array(makeup(Group))
            for i in range(Group.shape[0]):
                Group_r=Group
                Group_r[i]=Group[i]+delta
                fs_r=np.array(makeup(Group_r))
                deltas[i]=(score(X,y,fs_r,I)-score(X,y,fs,I))/delta
            return deltas

def makeup(Group):
    return [makeC(),makeL(Group[0],Group[1])]
def relu(li,reas):
    for i in range(li.shape[0]):
        li[i]=max(reas[i][0],min(li[i],reas[i][1]))
    return li
def lo2col(lo):
    ind=list(range(len(lo)))
    ind=np.array(ind)[lo]
    return ind
def fscore(x,y,I,w,Group): 
    f=makeup(Group)
    A= get_A(x,y,f,I,w)
    score=plot_A(x,y,f,I,A,mode=False)
    return score

class GA:
    """Genetic Algorithm Class"""
    

    def __init__(self,GroupN,SurviveN,W=None):
        """initial the properties of a population"""
        self.GroupN=GroupN;"""define the number of survivors in each generation"""
        self.SurviveN=SurviveN; """define the number of survivors in each generation"""
        self.W=W
    def born(self,bornlist,limitlist):
        """born the original ones"""
        self.Groups=np.array(bornlist) 
        self.L=limitlist
        return self
    def datasets(self,X,y,I):
        """get the datasets to fit, in the view of GA, 
        it means the Natural Enviroment to screening"""
        self.X=X
        self.y=y
        self.I=I
    def fit(self):
        """search for the parameter with the Linear Combination of the current parameters"""
        X=self.X
        y=self.y
        I=self.I
        n_g=self.Groups.shape[0]
        news=[]
        for j in range(n_g*self.GroupN):
            news.append(2*np.sum(rnd(self.Groups.shape)*np.array(self.Groups)/n_g,axis=0))
        for i in news:
            i=relu(i,self.L)
        news=np.vstack((news,self.Groups))
        fs=np.array(list(map(makeup,news)))
        scores=np.array(list(map(lambda f:score(X,y,f,I,self.W),fs)))
        index=np.isnan(scores)
        index=~index
        index=np.array(lo2col(index))
        news=news[index]
        scores=scores[index]
        index=np.argsort(scores)
        news=news[index]
        scores=scores[index]
        self.Groups=news
        self.Groups=self.Groups[:self.SurviveN]
        self.Scores=scores[:self.SurviveN]
        gc.collect()
        return scores
    def random_fit(self): 
        """search for the parameters in the range setted advanced randomly,
            using Genetic method"""
        X=self.X
        y=self.y
        I=self.I
        n_g=self.Groups.shape[0]
        news=[]
        vec0=np.array([i[0] for i in self.L])
        ranges=np.array([i[1]-i[0] for i in self.L])
        for j in range(n_g*self.GroupN):
            news.append(vec0+rnd()*ranges)
        news=np.array(news)
        news=np.vstack((news,self.Groups))
        fs=np.array(list(map(makeup,news)))
        scores=np.array(list(map(lambda f:score(X,y,f,I,self.W),fs)))
        index=np.isnan(scores)
        index=~index
        index=np.array(lo2col(index))
        news=news[index]
        scores=scores[index]
        index=np.argsort(scores)
        news=news[index]
        scores=scores[index]
        self.Groups=news
        self.Groups=self.Groups[:self.SurviveN]  
        self.Scores=scores[:self.SurviveN]
        del scores,index,news
        gc.collect()
        return self.Scores
    
    def SGD(self,epoch=10,step=10):
        """ using the gredient method to search for the local opyimum solutions of
        each individual of the population"""
        

        while epoch>0:   
            gd=partial(GD,self.X,self.y,self.I)
            Gredient=pool.map(gd,self.Groups)
            self.Groups-=rnd(self.Groups.shape)*step*np.array(Gredient)
            for i in self.Groups:
                i=relu(i,self.L)
            epoch-=1

        fs0=np.array(list(map(makeup,self.Groups)))
        #fs1=np.array(list(map(makeup,Groups_afterGD)))
        self.Scores=np.array(list(map(lambda f:score(self.X,self.y,f,self.I,self.W),fs0)))
        index=np.argsort(self.Scores)
        self.Scores=self.Scores[index]
        self.Groups=self.Groups[index]
        del gd,fs0,Gredient,index
        gc.collect()
        return self.Scores
    def searchrandom(self,epoch=50,step=0.1,N=30,n_jobs=False):
        """search for the local optimum solutions of each individual of the population."""     
        if self.GroupN==1:
            self.Groups=np.array([copy.deepcopy(self.Groups) for i in range(N)])
        
        while epoch>0:   
            Groups_r=np.array(self.Groups+rnd(self.Groups.shape)*step)
            for i in self.Groups:
                i=relu(i,self.L)
            funcs=partial(fscore,self.X,self.y,self.I,self.W)
            if n_jobs:    
                scores_r=np.array(pool.map(funcs,Groups_r))#pool.map(GD,Groups)
            else:
                fs0=np.array(list(map(makeup,self.Groups)))
                scores_r=np.array(list(map(lambda f:score(self.X,self.y,f,self.I,self.W),fs0)))
            
            index=scores_r<self.Scores
            index=np.array(lo2col(index))

            self.Scores[index]=scores_r[index]
            self.Groups[index]=Groups_r[index]
            epoch-=1
        gc.collect()
        return self.Scores

        
        
"""Lagrange Interpolation"""
def Extreme(y,x):
    def GetExtreme(y,step=10):
        for i in range(len(y)-2*step):
            if y[i:i+step].mean()<y[i+step//2:i+3*step//2].mean()\
                 and y[i+step:i+2*step].mean()<y[i+step//2:i+3*step//2].mean():

                return i+10
        return -1
    ranges=[1,30]
    st=ranges[0]
    while st<ranges[1]:
        ret=GetExtreme(y,st)
        if ret!=-1:
            return y[ret:],x[ret:]
        st+=1
    return y[0:],x[0:]

def LagI(Ilist,index,I):
    N=len(Ilist)
    """get the denominator of  i_th of the Lagrange paramaters ,i = index"""
    Ili=copy.deepcopy(Ilist)
    del Ili[index]
    Ili=[(I-i)/N for i in Ili]
    return np.prod(Ili)
def Lags(Ilist,fs):
    """get the Lagrange Polynimial Function"""
    Ili0=copy.deepcopy(Ilist)
    Ili=[LagI(Ili0,ind,i) for ind,i in enumerate(Ilist)]
    def f(X,I):
        return np.sum([fs[ind](X,I)*LagI(Ilist,ind,I)/i \
                       for ind,i in enumerate(Ili)])
    return f
def makeitok(a,b,c,d):
     """get the inverse function with the structure of our mathmatical model"""
     def f(y,I):
        
        return (c-d*np.exp((y-a)/b))/I/(1+np.exp((y-a)/b))
     return f
def makeitok2(a,b,c,d):
    """get the function with the structure of our mathmatical model"""
    def f(t,I):
        return a+b*np.log((c-I*t)/np.abs(d+I*t))
    return f
def LagP_inv(A,B,C,D,I):
    """get a Lagrange Polynomial by using the parameter arrays A,B,C and 
    D, while each array contains a kind of parameter in the mathmatical model T(I,U) we've reached"""
    fs=[]
    for i in range(len(A)):
        fs.append(makeitok(A[i],B[i],C[i],D[i]))
    return Lags(I,fs)
def LagP(A,B,C,D,I):
    """get a Lagrange Polynomial by using the parameter arrays A,B,C and 
    D, while each array contains a kind of parameter in the mathmatical model U(I,T) we've reached"""
    fs=[]
    for i in range(len(A)):
        fs.append(makeitok2(A[i],B[i],C[i],D[i]))
    return Lags(I,fs)
def getMRE(f,y,x,I):
    """calculate MRE"""
    yy,xx=Extreme(y,x)
    x_pre=[f(i,I)for i in yy]
    rate=np.abs((xx-x_pre)/xx)
    return np.mean(rate[::-1][:231])
def predict(f,y,I,isy=False):
    """use the established model to predict
    return a scalar"""
    x=range(len(y))
    yy,xx=Extreme(y,x)
    x_pre=[f(i,I)for i in yy]
    if isy==True:
        ret=[]
        for i in x_pre:
            if i!=i or i<9.0:
                return ret
            ret.append(i)
    return x_pre
def Vectorize(f,y,I):
    """use the established model to predict
    return a vector"""
    try:
        ret=[f(i,I) for i in y]
    except:
        ret=[f(y,i) for i in I]
    return ret
def Matrize(f,yvec,Ivec):
    """use the established model to predict
    return a matrix"""
    ymat,Imat=np.meshgrid(yvec,Ivec)
    subfuncs=np.frompyfunc(f,2,1)
    X_pre=subfuncs(ymat,Imat)
    return X_pre
def Matrize1(f,ymat,Imat):
    """use the established model to predict
    return a matrix and plot the 3-D figure of f(y,I)"""
    subfuncs=np.frompyfunc(f,2,1)
    X_pre=subfuncs(ymat,Imat)
    plt.subplot(projection='3d').plot_surface(ymat,Imat,X_pre,alpha=0.2)
    return X_pre



