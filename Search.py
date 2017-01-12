#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 12:49:39 2017

@author: root
"""
import random
import copy
from Structure import argabsmax
from Structure import argabsmin
from Method import vecnorm,vecstd,vecmean
from FunctionalPro import Lins,Function

def Searchabsmax(f,Xrange,retv=False,ab=[]):
    

    Xvec=Lins(Xrange)

    if not isinstance(f,Function):
        f=Function(f)
    for ab_i in ab:
        try:
            Xvec.remove(ab_i)
        except:
            None
    Y=f(Xvec)
    if retv:
        return Y[argabsmax(Y)]
        
    return Xvec[argabsmax(Y)]
def Searchabsmin(f,Xrange,retv=False,ab=[]):
    Xvec=Lins(Xrange)
    if not isinstance(f,Function):
        f=Function(f)
    for ab_i in ab:
        try:
            Xvec.remove(ab_i)
        except:
            None
    Y=f(Xvec)
    if retv:
        return Y[argabsmin(Y)]
        
    return Xvec[argabsmin(Y)]
        
