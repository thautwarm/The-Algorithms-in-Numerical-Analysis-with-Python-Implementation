#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 19:01:34 2016

@author: thaut
"""
from __future__ import division
from FunctionalPro import Function

def std_vec(x):
    if x==[]:
        return 0
    length=len(x)
    Mean=sum(x)/length
    stds=[(i-Mean)**2 for i in x]
    ret=( sum(stds)/length )**(1/2)
    return ret
        

    

from matplotlib import pyplot as plt
def FixedIter(Func ,Param,x0=1,epoch=20,depth=10,limit=10e-5,verbose=False):
    stack=[];r=0
    fc=lambda x: Func(**Param)(1/x)
    fc=Function(fc)+Function(lambda x:x)
    while r<epoch:
        x_=fc(x0)
        
        stack.append(x_)
        x0,x_=x_,x0
        if verbose:
            print stack[-1]
        if r>depth:
            stack=stack[-depth:]
            
        r+=1
    if  std_vec(stack)>limit :
        print 'Not Convergenced!'
        

    elif verbose:
            print 'Convergenced!'
    plt.plot(stack)
    
    _x=fc(x0)
    if _x-2*x0+x_!=0:
        ret=(x_*_x-x0**2)/(_x-2*x0+x_)
    else:
        ret=x0
    return ret
    
def StrangeIter(Func,Param,dFunc=None,x0=1,epoch=10,depth=5,limit=10e-7,verbose=False):
    print 'Specific'
    try:
        p2=Param[1]
        p1=Param[0]
    except:
        p2=p1=Param
    
 
    fc=Func(**p1)
    if not dFunc:
        def _dFunc(x):
            return (fc(x+limit)-fc(x))/limit
        df=_dFunc
    stack=[];r=0
    flag=True
    Xstack=[x0,x0]
    try:
        df
    except:
        df=dFunc(**p2)
    
    while True:
        div=df(Xstack[-1]) if r%2==0 else df(Xstack[-2])
        if div==0:
            print 'with division 0!'
            return ValueError

        x_=Xstack[-1]-fc(Xstack[-1])/div
        stack.append(x_)
        if r>1 and r%20==0 :
            flag=False
            print r,'    ',x_
        if verbose and flag:
                print stack[-1]
        if r>depth:
            if  std_vec(stack[-depth:])<limit :
                print str(r)+' times Convergenced!'
                break
          
        r+=1

        
        if r>epoch:
            break
        del Xstack[0]
        Xstack.append(x_)
    Xstack.append(x_)    
    plt.xlabel('Times')
    plt.ylabel('f(x_now)')
    plt.plot(stack,label='Specific Method')
    plt.legend()
    if  std_vec(stack[-depth:])>limit :
        print 'Not Convergenced!'
        plt.plot(stack[-depth:])
    elif verbose:
            print 'Convergenced!'
        
    return x_ 
    
    
def NewtonIter(Func,Param,dFunc=None,x0=1,epoch=10,depth=5,limit=10e-7,verbose=False):
    print 'Newton'
    try:
        p2=Param[1]
        p1=Param[0]
    except:
        p2=p1=Param

    fc=Func(**p1)

    if not dFunc:
        def _dFunc(x):
            return (fc(x+limit)-fc(x))/limit
        df=_dFunc
    stack=[];r=0
    flag=True
    try:
        dFunc
    except:
        df=dFunc(**p2)
    stack=[];r=0
    while r<epoch:
        div=df(x0)
        if div==0:
            print 'with division 0!'
            return ValueError
        x_=x0-fc(x0)/div
        if r>1 and  r%20==0:
            flag=False   
            print r,'  ',x_
        stack.append(x_)
        if verbose and flag :
            print stack[-1]
        x0,x_=x_,x0
        if r>depth:
            if  std_vec(stack[-depth:])<limit :
                print str(r)+' times Convergenced!'
                break
          
            
        r+=1
 
        
    plt.plot(stack,color='green',label='Newton')
    plt.legend()
    if  std_vec(stack[-depth:])>limit :
        print 'Not Convergenced!'
    
    elif verbose:
            print 'Convergenced!'
            

    
    return x0

    
def NewtonIter2(Func,Param,dFunc=None,x0=1,epoch=10,depth=5,limit=10e-7,verbose=False):
    print 'Newton2'
    try:
        p2=Param[1]
        p1=Param[0]
    except:
        p2=p1=Param

    fc=Func(**p1)

    if not dFunc:
        def _dFunc(x):
            return (fc(x+limit)-fc(x))/limit
        df=_dFunc
    stack=[];r=0
    flag=True
    try:
        dFunc
    except:
        df=dFunc(**p2)
    stack=[];r=0
    while r<epoch:
        div=df(x0)
        if div==0:
            print 'with division 0!'
            return ValueError
        x_=x0-fc(x0)/div
        if abs(fc(x_))>abs(fc(x0)):
            x_=(x_+x0)/2
        if r>1 and  r%20==0:
            flag=False   
            print r,'  ',x_
        stack.append(x_)
        if verbose and flag :
            print stack[-1]
        x0,x_=x_,x0
        if r>depth:
            if  std_vec(stack[-depth:])<limit :
                print str(r)+' times Convergenced!'
                break
          
            
        r+=1
 
        
    plt.plot(stack,color='red',label='Newton_adjusted')
    plt.legend()
    if  std_vec(stack[-depth:])>limit :
        print 'Not Convergenced!'
    
    elif verbose:
            print 'Convergenced!'
            

    
    return x0
    
def UnknownIter(Func,Param,dFunc=None,ddFunc=None,x0=1,epoch=10,depth=5,limit=10e-9,verbose=False):
    try:
        p3=Param[2]
        p2=Param[1]
        p1=Param[0]

    except:
        p3=p2=p1=Param
    
    reduces=[]
    fc=Func(**p1)
    if not dFunc:
        def _dFunc(x):
            return (fc(x+limit)-fc(x))/limit
        df=_dFunc
    else:
        df=dFunc(**p2)
    if not ddFunc:
        def _ddFunc(x):
            return (df(x+limit)-df(x))/limit
        ddf=_ddFunc
    else:
        ddf=dFunc(**p3)
    stack=[];r=0
    flag=True
    while r<epoch:
        div=2*df(x0)**2-ddf(x0)*fc(x0)
        if div==0:
            print 'with division 0!'
            return ValueError
        x_=x0-2*fc(x0)*df(x0)/div
        if r>1 and  r%20==0:
            flag=False   
            print r,'x_now',x_,' |x_now - x_last|', abs(x_-x0)
            
        stack.append(x_)
        if verbose and flag :
            print r,'x_now',x_,' |x_now - x_last|', abs(x_-x0)
            reduces.append(abs(x_-x0))
        x0,x_=x_,x0
        if r>depth:
            if  std_vec(stack[-depth:])<limit :
                print str(r)+' times Convergenced!'
                break
          
            
        r+=1
 
    
    plt.plot(stack,color='green',label='This method')
    plt.xlabel('Times')
    plt.ylabel('f(x_now)')
    plt.legend()
    if  std_vec(stack[-depth:])>limit :
        print 'Not Convergenced!'
    
    elif verbose:
            print 'Convergenced!'
    print 'convergence ratio:( |delta x_n | / |delta x_n+1|)'
    function=lambda x:x[0]/x[1] if x[1]!=0 else 'inf'
    ratio=map(function,zip(reduces[:-1],reduces[1:]))
    for i in ratio:
        print '\t',i

    
    return x0
    
    

    
 

