#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 19:56:57 2016

@author: root
"""

from __future__ import division
from Structure import matrix,argabsmax,absmax,dot_vec,add_vec,cg_qr,functionStack,Givens
import copy
import math
def cumsum():pass
def GaussEq(A,B):
        if isinstance(B,matrix):
            B=B.c(range(B.size[0]),0)
        B=copy.deepcopy(B)
        size=A.size
        data=copy.deepcopy(A.data)
        row_arg=range(size[0])
        if size[0]!=size[1]:
            print('error\n')
            return ValueError
        size=size[0]
        for i in row_arg[:-1]:
            vec=matrix(data).c(range(i,size),i)
            max_index=i+argabsmax(vec)
            data[max_index],data[i]=data[i],data[max_index]
            B[max_index],B[i]=B[i],B[max_index]
            for index in range(i+1,size):
                try:
                    param=data[index][i]/data[i][i]
                    data[index]=add_vec(data[index],data[i],-param)
                    B[index]+=B[i]*(-param)
                except:
                    print('error: a singular matrix')
                    return ValueError
        row_arg.reverse()
        x=[0 for i in range(size)]
        for i in row_arg:
            x[i]=(B[i]-sum(dot_vec(x[i:],data[i][i:])))/data[i][i]
        return x
def LU(A):
    size=A.size
    if size[0]!=size[1]:
        return ValueError
    L=matrix(size,eye=True)
    U=matrix(size,eye=False)
    size=size[0]
    aij=A.data
    for i in range(0,size):
        rang=range(0,i)
        for r in range(i,size):
            U.data[i][r]=aij[i][r]-sum(dot_vec(L.c([i],rang),U.c(rang,[r]) ))
        for r in range(i+1,size):
            L.data[r][i] =  ( aij[r][i]-sum(dot_vec( L.c([r],rang),U.c(rang,[i]) )) )/U.data[i][i]

    return L.data,U.data
def check_symmetrical(A):
    n=A.size[0]
    data=A.data
    for i in range(1,n):
        for j in range(i):
            if data[i][j]!=data[j][i]:
                return False
    return True
    
"""about Eigenvalue """    

class PowIter:
    def __init__(self,A):
        self.Mat=A
        size=A.size
        if size[0]!=size[1]:
            return ValueError
        self.X=[1 for i in range(size[0])]        
    def fit(self,epoch=10):
        data=self.Mat.data
        X=self.X
        while epoch:
            
            m=absmax(X)
            X=[i/m for i in X]
            X=[sum([ A_i[j]*X_j for j,X_j in enumerate(X) ]) for A_i in data]
            epoch-=1
        #get rayleigh div
        AX=[sum([ A_i[j]*X_j for j,X_j in enumerate(X) ]) for A_i in data] 
        m1=sum(dot_vec(AX,X))/sum(dot_vec(X,X))
        #m2=absmax(X)
        self.X=X
        return m1#,m2
def QR_split(B,NumLimit=10e-9,P_is_a_matrix=False):
    A=copy.deepcopy(B)
    if isinstance(A,matrix):
        size=A.size[0]
        A=A.data
    else:
        size=len(A)
    P=functionStack()
    for i in xrange(size):
        for j in xrange(i+1,size):
            c,s=cg_qr(A,i,j)
            kw={'a':i,'b':j,'c':c,'s':s}
            P.getin(Givens(A,**kw))
            if abs(A[j][i])<NumLimit:
                    A[j][i]=0
    if P_is_a_matrix:
        e=matrix(size).data
        e=P.pullout(e)
        del P
        P=matrix(e)
    return P,matrix(A)
def QR_Iter(A,epoch=10):
    while epoch:
        P,A=QR_split(A,P_is_a_matrix=True)
        A=A*P.T()
        epoch-=1
    return A
            
        

                
        

def Norm(A,p='inf',epoch_if_norm2=200):
    if isinstance(A,matrix):
        data=A.data
        tdata=A.T().data
    else :
        data=A
        A=matrix(data)
        tdata=A.T().data
    switch={1:lambda A:max([sum(abs(j) for j in i) for i in tdata]),
            2:lambda A:(PowIter(A.T()*A).fit(epoch_if_norm2)**(1/2)),
            'inf':lambda A:max([sum(abs(j) for j in i) for i in data]),
            'F':lambda A: math.sqrt(sum([data[i][i]**2 for i in range(len(data))]) )
            }
    return switch[p](A)
def Cond(A,p=2):
    if isinstance(A,matrix):
        inv=A.inv()
    else:
        inv=matrix(A).inv()
    return Norm(A,p)*Norm(inv,p)
    
    
def Chol(A):
    size=A.size
    if size[0]!=size[1]:
        return ValueError
    if not check_symmetrical(A) :
        print 'Not symmetrical!\n'
        return ValueError
    aij=A.data
    L=matrix(size,eye=False)
    size=size[0]
    for j in range(size):
        try:
            L[j][j] = math.sqrt(aij[j][j]-sum( dot_vec(L[j][:j],L[j][:j]) ) )
        except:
            print 'error:Not positive!\n'
            
        for i in range(j+1,size):
            try:
                L[i][j]=(aij[i][j] - sum( dot_vec(L[i][:j],L[j][:j] ) ) ) / L[j][j]
            except :
                print 'error:Not positive!\n'
    return L   

def Thomas(a,b,c):
    n=len(b)
    if n<=1:
        return ValueError 
        
    l=[0 for i in xrange(n-1)]
    u=[0 for i in range(n)]
    u[0]=b[0]
    for i in range(n-1):
       l[i]=a[i]/u[i]
       u[i+1]=b[i+1]-l[i]*c[i]
    return l,u,c
class ThomasEq:
    def __init__(self,a,b,c):
        self.tool=Thomas(a,b,c)
        
    def fit(self,D):
        l,u,c=self.tool
        n=len(u)
        y=[0 for i in range(n)]
        x=[0 for i in range(n)]
        y[0]=D[0]
        for i,l_i in enumerate(l):
            y[i+1]=D[i+1]-l_i*y[i]
        y.reverse()
        u.reverse()
        c.reverse()
        x[0]=y[0]/u[0]
        for i in range(1,n):
            x[i]=(y[i]-c[i-1]*x[i-1])/u[i]
        x.reverse()
        return x
class LUEq:
    def __init__(self,A):
        L,U=LU(A)
        self.L=matrix(L)
        self.U=matrix(U)
        print U
    def fit(self,B):
        if isinstance(B,matrix):
            B=B.c(range(B.size[0]),0)
        n=len(B)
        L=self.L
        U=self.U
        y=[0 for i in xrange(n)]
        x=[0 for i in xrange(n)]
        for i in  xrange(n):
            rang=range(i)
            y[i]=B[i]-sum(dot_vec( L.c([i],rang), y[:i]))
        for i in xrange(n-1,-1,-1):
            rang=range(i+1,n)
            x[i]= (y[i] - sum(dot_vec(U.c([i],rang),x[(i+1):]) ) )/U.data[i][i]
        return x
def Jacobi(self,b,mat=True):
        D=self.D
        A=self.A.data
        B=[[-A_ij/D_i if i!=j \
           else 1-A_ij/D_i for j,A_ij in enumerate(A[i])]for i,D_i in enumerate(D)]
        if mat==True:
            B=matrix(B)
        f=[b_i/D[i] for i,b_i in enumerate(b)]
        f=matrix(f).T()
        return B,f
def Gauss_Seidel(self,b,mat=True):
        D=self.D
        L=self.L
        A=self.A
        size=A.size
        Dinv=[[-L[i-1][j] if j<i  else 0 for j in xrange(size[1])] for i in xrange(size[0])]
        for i in range(size[0]):
            Dinv[i][i]=D[i]
        Dinv=matrix(Dinv).inv(mat=True)
        B=matrix(size,eye=True)-Dinv*A
        if mat==False:
            B=B.data
        f=Dinv*matrix(b).T()
        return B,f
def _SOR(self,b,w,mat=True):
        D=self.D
        L=self.L
        A=self.A
        size=A.size
        data=A.data
        Dinv=[[-w*L[i-1][j] if j<i  else 0 for j in xrange(size[1])] for i in xrange(size[0])]
        LA=[[-w*data[i][j] if j>=i  else 0 for j in xrange(size[1])] for i in xrange(size[0])]
        
        for i in range(size[0]):
            Dinv[i][i]=D[i]   
            LA[i][i]+=D[i]
        Dinv=matrix(Dinv).inv(mat=True)
        LA=matrix(LA)
        SORmat=Dinv*LA
        b=matrix(b).T()
        f=w*Dinv*b
        return SORmat,f
class IterEq:
    def __init__(self,A):
        data=A.data
        self.size=A.size
        size=A.size[0]
        self.A=copy.deepcopy(A)       
        self.L=[[-data[i][j] for j in range(i)]for i in range(1,size)]
        self.D=[data[i][i] for i in range(size)]
        self.X=matrix([1 for i in self.D]).T()
    def IterParamSet(self,**kw):
        b=kw['b']
        method=kw['method']
        try:
            w=kw['w']
            B,f=method(self,b,w,mat=True)
        except:
            B,f=method(self,b,mat=True)
        self.B=B
        self.f=f
    def fit(self,epoch=10):
        try:
            B=self.B
            f=self.f
        except:
            print 'Not initialized B and f!\n'
        x=self.X
        while epoch:
            x=B*x+f
            epoch-=1
        return x
    def GS(self,b,epoch=10):
        kw={'method':Gauss_Seidel,
            'b':b}
        self.IterParamSet(kw)
        x=self.fit(epoch)
        return x
    def Jaco(self,b,epoch=10):
        kw={'method':Jacobi,
            'b':b}
        self.IterParamSet(kw)
        x=self.fit(epoch)
        return x
    def SOR(self,b,w=1.18,epoch=10):
        kw={'method':_SOR,
            'b':b,
            'w':w}
        self.IterParamSet(**kw)
        x=self.fit(epoch)
        return x

        
        
        
        
#def Jacobi_Eq(A,B):
    