#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 14:09:08 2016

@author: root
"""
from __future__ import division
import copy
import math
from itertools import cycle as repeat
def add_vec(x,y,w=1):
    x_=repeat(x)
    do=zip(x_,y)
    ret=[i[0]+w*i[1] for i in do]
    return ret
def dot_vec(x,y):
    
    x_=copy.deepcopy(x)
    if type(x)!=list:
        x_=[x_]
    elif len(x)==0 or len(y)==0: return [0]
    x_=repeat(x_)
    do=zip(x_,y)
    ret=[i[0]*i[1] for i in do]
    return ret
def argabsmax(x):
    max=x[0]
    arg=0
    for ind,i in enumerate(x[0:]):
        if abs(max)<abs(i):
            max=i
            arg=ind
    return arg
def absmax(x):
    max=x[0]
    for ind,i in enumerate(x[0:]):
        if abs(max)<abs(i):
            max=i
    return max
def jacobi_sign(x):
    if x>=0:
        ret=1
    else :
        ret=-1
    return ret
def UnitMat_rev(body,**kw):
    a=kw['a']
    b=kw['b']
    body[a],body[b]=body[b],body[a]
def UnitMat_add(body,**kw):
    a=kw['a']
    b=kw['b']
    w=kw['w']
    body[a]=add_vec(body[a],body[b],w)
def UnitMat_dot(body,**kw):
    a=kw['a']
    w=kw['w']
    body[a]=[w*i for i in body[a]]
def calculate_givensSy_cos_sin(body,i,j):
    a=body[i][i]
    b=body[j][j]
    ab=body[j][i]
    A=a-b
    if ab==0:
        d='inf'
    else:
        d=A/2/ab
    if d=='inf':
        c=0;s=1
    else:
        tg=jacobi_sign(d)/( abs(d)+math.sqrt(d**2 + 1) )
        c=1/math.sqrt(1+ tg**2)
        s=tg/math.sqrt(1+ tg**2)
    return c,s
def calculate_givensQR_cos_sin(body,i,j):
    a=body[i][i]
    b=body[j][i]
    if abs(b)>=abs(a):
        tg=a/b
        s=1/math.sqrt(1+ tg**2)
        c=tg/math.sqrt(1+ tg**2)
    else:
        tg=b/a
        c=1/math.sqrt(1+ tg**2)
        s=tg/math.sqrt(1+ tg**2)
    return c,s
def cg_qr(*tup):
    return calculate_givensQR_cos_sin(*tup)
def Givens(body,**kw):
    a=kw['a']
    b=kw['b']
    c=kw['c']
    s=kw['s']
    if c==0:
        kw0={'a':a,'b':b}
        UnitMat_rev(body,**kw0)
    else:
        
        kw1={'a':a,'w':1/c}
        kw2={'a':b,'w':c}
        kw3={'a':a,'b':b,'w':c*s}
        kw4={'a':b,'b':a,'w':-s/c}
        UnitMat_add(body,**kw4)
        UnitMat_add(body,**kw3)
        UnitMat_dot(body,**kw2)
        UnitMat_dot(body,**kw1)
    if c==0:
        def f(A):
            UnitMat_rev(A,**kw0)
    else:    
        def f(A):
            UnitMat_add(A,**kw4)
            UnitMat_add(A,**kw3)
            UnitMat_dot(A,**kw2)
            UnitMat_dot(A,**kw1)
    return f

UnitMethod={'rev':UnitMat_rev,
            'add':UnitMat_add,
            'dot':UnitMat_dot,
            'givens':Givens
            }
class functionStack:
    def __init__(self):
        self.stack=[]
    def getin(self,func_item):
        self.stack.append(func_item)
    def getout(self):
        return self.stack.pop()
    def pullout(self,mat):
        n=len(self.stack)
        while n:
            self.getout()(mat)
            n-=1
        return mat
class UnitMatrix:
    def __init__(self):
        self.stack=[]
    def getin(self,method,**kw):
        self.stack.append((kw,method))
    def func(self):
        kw,method=self.stack.pop()
        return lambda A:UnitMethod[method](body=A,**kw)
class matrix:
    def __init__(self,lists,eye=True):
        
        if type(lists)==list and type(lists[0])!=str:
            if type(lists[0])!=list:
                lists=[lists]
            self.data=copy.deepcopy(lists)
            self.size=(len(lists),len(lists[0]))
        elif type(lists)==tuple:
            m=lists[0];n=lists[1]
            self.size=(m,n)
            if eye==True:
                self.data=[[1*(i==j) for j in xrange(n)]for i in xrange(m)]
            else:
                self.data=[[0 for j in xrange(n)]for i in xrange(m)]
        elif type(lists)==int:
            n=lists
            self.size=(n,n)
            if eye==True:
                self.data=[[1*(i==j) for j in xrange(n)]for i in xrange(n)]
            else:
                self.data=[[0 for j in xrange(n)]for i in xrange(n)]
        else:
            self=self.__init__(lists.data)
    
    def __setitem__(self,i,x):
        self.data[i]=x
    def c(self,index1,index2,mat=False):
        data=self.data
        if type(data)==type(None):
            print('error: not initialized.\n')
        if type(index2)!=list:
            index2=[index2]
        if type(index1)!=list:
            index1=[index1]
        mat=[[data[i][j] for j in index2] for i in index1]
        ret=[]
        for i in mat:
            ret+=i
        return ret
    def disp(self,dig=2):
        data=self.data
        for i in data:
            for j in i:
                print " %f " % (j),
            print '\n'
    def __getitem__(self,I):
        if type(I) in [slice,int]:
            ret=self.data[I]
        return ret
    def T(self,retdata=False):
        data=self.data
        size=self.size
        datap=[[data[j][i] for j in xrange(size[0])] for i in xrange(size[1])]
        if retdata:
            return datap
        else:
            return matrix(datap)
    
    def __rmul__(self,x):
            try:
                data=copy.deepcopy(self.data)
                data=[[j*x for j in i]for i in data]
                return matrix(data)
            except:
                return 0
    def __mul__(self,x):         
            data2=x.T().data
            data1=self.data
            mat=[[sum(dot_vec(i,j)) for j in data2 ]for i in data1]
            return matrix(mat)
    def __add__(self,x):
        if not isinstance(x,matrix):
            try:
                data=copy.deepcopy(self.data)
                data=[[j+x for j in i]for i in data]
                return matrix(data)
            except:
                None
        
        size=x.size
        try:
            data1=x.data
            data2=self.data
            mat=[[data1[i][j]+data2[i][j] for j in xrange(size[1])] for i in xrange(size[0])]
        except:
            print('error: unexpected size of matrix!\n')
        return matrix(mat)
    
    def __sub__(self,x):
        size=x.size
        try:
            data1=x.data
            data2=self.data
            mat=[[data2[i][j]-data1[i][j] for j in xrange(size[1])] for i in xrange(size[0])]
        except:
            print('error: unexpected size of matrix!\n')
        return matrix(mat)
    def __neg__(self):
        data=self.data
        ret=[[-j for j in i]for i in data]
        return matrix(ret)
    
    def inv(self,overall=False,mat=False):
        size=self.size
        data=copy.deepcopy(self.data)
        row_arg=range(size[0])
        if size[0]!=size[1]:
            print('error\n')
            return ValueError
        inverse=matrix((size)).data
        size=size[0]
        for i in row_arg[:-1]:
            vec=matrix(data).c(range(i,size),i)
            max_index=i+argabsmax(vec)
            data[max_index],data[i]=data[i],data[max_index]
            inverse[max_index],inverse[i]=inverse[i],inverse[max_index]
            for index in range(i+1,size):
                param=data[index][i]/data[i][i]
                data[index]=add_vec(data[index],data[i],-param)
                inverse[index]=add_vec(inverse[index],inverse[i],-param)
        row_arg.reverse()
        for i in row_arg[:-1]:
                for index in range(0,i):
                    param=data[index][i]/data[i][i]
                    data[index]=add_vec(data[index],data[i],-param)
                    inverse[index]=add_vec(inverse[index],inverse[i],-param)
        for i in row_arg:
             inverse[i]=dot_vec(1/data[i][i],inverse[i])
             data[i]=dot_vec(1/data[i][i],data[i])
        if mat==False:
            return inverse
        elif mat==True:
            return matrix(inverse)
        

    
    