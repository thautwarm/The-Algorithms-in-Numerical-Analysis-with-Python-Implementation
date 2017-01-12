#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 10:47:32 2017

@author: root
"""



from itertools import repeat
def sign(x):
    if x>0:
        return  1
    elif x<0:
        return -1
    else:
        return 0
class StrangeNum:
    def __init__(self,t):
      self.t=t
    def __lt__(self,x):
         if isinstance(x,StrangeNum):
             return self.t<x.t
         return (2*self.t)<sign(x)
    def __gt__(self,x):
         if isinstance(x,StrangeNum):
             return self.t>x.t
         return (2*self.t)>sign(x)
    def __neg__(self):
        return StrangeNum(-1*self.t)
    
inf=StrangeNum(1)
nan=StrangeNum(0)

class Function:
    def __init__(self,basefunc):
        self.f=basefunc
    def __call__(self,x):
        if type(x)==int:
            x=float(x)
        if type(x)==list:
            return [self.f(x_i) for x_i in x]
        return self.f(x)
    def __add__(self,b):
        if  type(b) != type(self):
            def f(x):
                return self.f(x)+b
            return Function(f)
        def f(x):
        
            return self.f(x)+b(x)
        return Function(f)
    def __rsub__(self,b):
        def f(x):
            return self.f(x)-b
            return Function(f)
        return Function(f)
    
    def __sub__(self,b):
        if  type(b) != type(self):
            def f(x):
                return self.f(x)-b
            return Function(f)
        def f(x):
            return self.f(x)-b(x)
        return Function(f)
    def __div__(self,b):
        if  type(b) != type(self):
            def f(x):
                if b==0.0:
                    return inf
                return self.f(x)/b
            return Function(f)
        def f(x):
            if b(x)==0.0:
                return inf
            return self.f(x)/b(x)
        return Function(f)
    def __rtruediv__(self,b):
        print 'tp'
        if  type(b) != type(self):
            def f(x):
                if b==0.0:
                    return inf
                return self.f(x)/b
            return Function(f)
        def f(x):
            if b(x)==0.0:
                return inf
            return self.f(x)/b(x)
        return Function(f)
        
    def __mul__(self,b):
        if  type(b) != type(self):
            def f(x):
                return self.f(x)*b
            return Function(f)
        def f(x):
            return self.f(x)*b(x)
        return Function(f)


            
          
        
def LagBase(Xvec,x0):
    vec=[x_i for x_i in Xvec if x_i!=x0]
    f=reduce(lambda x,y:x*y ,map(lambda x_i:Function(lambda x: float(x-x_i) ) , vec))
    f=f/f(x0)
    return f    

def Lagrange_Interpolation(Xvec,Yvec):
    f=reduce(lambda x,y:x+y,map(lambda x,y,Xvec:LagBase(Xvec,x)*y, Xvec,Yvec,repeat(Xvec,len(Yvec))))
    def cf(x):
        x=float(x)
        return f(x)
    return f
def Lins(Xrange):
        b=Xrange[1]
        a=Xrange[0]
        N=Xrange[2]
        interval=1.0*(b-a)/Xrange[2]
        try:
            Xvec=[a+i*interval  for i in range(N+1)]
        except:
            print Xrange
        return Xvec
def Function_Interpolation(f,Xrange=None,Xvec=None):
    if not Xvec:
        Xvec=Lins(Xrange)
    Yvec=f(Xvec)
    return Lagrange_Interpolation(Xvec,Yvec)

        