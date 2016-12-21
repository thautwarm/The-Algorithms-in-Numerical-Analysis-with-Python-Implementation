#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 19:01:34 2016

@author: thaut
"""
from __future__ import division
def std_vec(x):
    if x==[]:
        return 0
    length=len(x)
    Mean=sum(x)/length
    stds=[(i-Mean)**2 for i in x]
    ret=( sum(stds)/length )**(1/2)
    return ret
        
class Polyn:
    def __init__(self,pow_index,param):
        N1=len(param)
        N2=len(pow_index)
        if N1==N2:
            N=N1
            del N1,N2
        else:
            print "length not match for pows and params!"
            return ValueError
        def f(x):
            y=[param[i] * x**pow_index[i]for i in range(N)]
            return sum(y)
        self.body=f
    def __call__(self,param):
        return self.body(param)
    def type(self):
        return 'Function'
    #def __add__(self,f):

def Linear(**kw):
    a=kw['a']
    b=kw['b']
    def f(x):
        ret= (a+b*x)
        return ret
    return f
    
    
def FixedIter(Func ,Param,x0=1,epoch=20,depth=5,limit=10e-5,verbose=False):
    stack=[];r=0
    while r<epoch:
        
        x_=Func(**Param)(x0)
        stack.append(x_-x0)
        x0,x_=x_,x0
        if r>depth:
            stack=stack[-depth:]
            if verbose:
                print std_vec(stack)
        r+=1
    if  std_vec(stack)>limit :
        print 'Not Convergenced!'
    
    elif verbose:
            print 'Convergenced!'
    _x=Func(**Param)(x0)
    if _x-2*x0+x_!=0:
        ret=(x_*_x-x0**2)/(_x-2*x0+x_)
    else:
        ret=x0
    return ret
    
def NewtonIter(Func,dFunc,Param,x0=1,epoch=10,depth=5,limit=10e-5,verbose=False):
    try:
        
        p2=Param[1]
        p1=Param[0]
    except:
        p2=p1=Param
    stack=[];r=0

    while r<epoch:
        div=dFunc(**p2)(x0)
        if div==0:
            print 'with division 0!'
            return ValueError

        x_=Func(**p1)(x0)/div
        stack.append(x_-x0)
        x0,x_=x_,x0
        if r>depth:
            stack=stack[-depth:]
            if verbose:
                print std_vec(stack)
        r+=1
    if  std_vec(stack)>limit :
        print 'Not Convergenced!'
    
    elif verbose:
            print 'Convergenced!'
    div=dFunc(**p2)(x0)
    if div:
        _x=Func(**p1)(x0)/div
        if _x-2*x0+x_!=0:
            ret=(x_*_x-x0**2)/(_x-2*x0+x_)
    else:
        ret=x0
    return ret
        
        
    
        
        