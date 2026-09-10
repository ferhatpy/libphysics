#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
libdiscretemath.py
Created on Sep 1, 2017

@author: Ferhat Nutku
'''
import numpy as np
import scipy
from sympy import *
from libsympy import *


class cdisretemath(object):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''
    
    def calc(self):
        print("optics")
        
        
### Aperiodic Sequences        
"""
Aperiodic Sequences
"""
[a,b] = symbols('a b')

def CheckerBoard(n, conj=False, text=False, verbose=False):
    m = int(n/2)
    if text == False:
        if conj==False:
            # CheckerBoard sequence
            a = np.array(([A(d_A), B(d_B)]*m + [B(d_B), A(d_A)]*m)*m).reshape((n,n))
        else:
            # Conjugated CheckerBoard sequence
            a = np.array(([B(d_B), A(d_A)]*m + [A(d_A), B(d_B)]*m)*m).reshape((n,n))
    
    elif text == True:
        if conj==False:
            # CheckerBoard sequence
            a = np.array((["A", "B"]*m + ["B", "A"]*m)*m).reshape((n,n))
        else:
            # Conjugated CheckerBoard sequence
            a = np.array((["B", "A"]*m + ["A", "B"]*m)*m).reshape((n,n))
        if verbose: print("CB({0}) = {1}".format(n,a))
    
    return flatten(a)
        

def Fibonacci(n, conj=False, text=False, verbose=False):
    """
    Generates a Fibonacci replacement sequence.
    A -> AB, B -> A.
    
    Fibonacci(0, text=True) -> B
    Fibonacci(1, text=True) -> A
    Fibonacci(2, text=True) -> AB
    Fibonacci(3, text=True) -> ABA
    """
    if text == False:
        if conj==False:
            # Fib sequence
            a, b = B(d_B), A(d_A)
        else:
            # Conjugated Fib sequence
            a, b = A(d_A), B(d_B)
        for i in range(0, n):
            a, b = b, b+a
    
    elif text == True:
        if conj==False:
            # Fib sequence
            a, b = "B","A"
        else:
            # Conjugated Fib sequence
            a, b = "A","B"
        for i in range(0, n):
            if verbose: print("F({0}) = {1}".format(i,a))
            a, b = b, b+a
        if verbose: print("F({0}) = {1}".format(n,a))
    
    return a


def Periodic(n, conj=False, text=False, verbose=False):
    """
    Generates a Periodic sequence.
    A -> AB, B -> AB
    """
    if text == False:
        if conj==False:
            # Periodic sequence
            a, b = A(d_A), B(d_B)
        else:
            # Conjugated Periodic sequence
            a, b = B(d_B), A(d_A)
        for i in range(0, n):
            a, b = a+b, a+b
    
    elif text == True:
        if conj==False:
            # Periodic sequence
            a, b = "A","B"
        else:
            # Conjugated Periodic sequence
            a, b = "B","A"
        for i in range(0, n):
            if verbose: print("PR({0}) = {1}".format(i,a))
            a, b = a+b, a+b
        if verbose: print("PR({0}) = {1}".format(n,a))
    
    return a

def Periodic2D(n, conj=False, text=False, verbose=False):
    """
    Generates a Periodic sequence. # todo same as periodic modify
    A -> AB, B -> AB
    """
    if text == False:
        if conj==False:
            # 2D Periodic sequence
            a, b = A(d_A), B(d_B)
        else:
            # 2D Conjugated Periodic sequence
            a, b = B(d_B), A(d_A)
        for i in range(0, n):
            a, b = a+b, a+b
    
    elif text == True:
        if conj==False:
            # 2D Periodic sequence
            a, b = "A","B"
        else:
            # 2D Conjugated Periodic sequence
            a, b = "B","A"
        for i in range(0, n):
            if verbose: print("PR2D({0}) = {1}".format(i,a))
            a, b = a+b, a+b
        if verbose: print("PR2D({0}) = {1}".format(n,a))
    
    return a


def PeriodDoubling(n, conj=False, text=False, verbose=False):
    """
    Generates a period-doubling replacement sequence.
    A -> AB, B -> AA.
    """
    if text == False:
        if conj==False:
            # Double-period sequence.
            a, b = A(d_A), B(d_B)
        else:
            # Conjugate double-period sequence.
            a, b = B(d_B), A(d_A)
        for i in range(0, n):
            a, b = a+b, a+a
    
    elif text == True:
        if conj==False:
            # Double-period sequence.
            a, b = "A","B"
        else:
            # Conjugate double-period sequence.
            a, b = "B","A"
        for i in range(0, n):
            if verbose:  print("DP({0}) = {1}".format(i,a))
            a, b = a+b, a+a
        if verbose: print("DP({0}) = {1}".format(n,a))
        
    return a

            
def RudinShapiro(n, conj=False, text=False,verbose=False):
    """    
    Generates a Rudin-Shapiro replacement sequence.
    AA -> AAAB, AB -> AABA
    BA -> BBAB, BB -> BBBA
    
    RudinShapiro(0, text=True) -> AA
    RudinShapiro(1, text=True) -> AA
    RudinShapiro(2, text=True) -> AAAB
    RudinShapiro(3, text=True) -> AAABAABA
    """
    if text == False:
        if conj==False:
            # RudinShapiro sequence
            aa, ab = A(d_A)+A(d_A), A(d_A)+B(d_B)
            ba, bb = B(d_B)+A(d_A), B(d_B)+B(d_B)
        else:
            # Conjugated RudinShapiro sequence
            aa, ab = B(d_B)+B(d_B), B(d_B)+A(d_A)
            ba, bb = A(d_A)+B(d_B), A(d_A)+A(d_A)
        for i in range(1, n):
            aa, ab, ba, bb = aa+ab, aa+ba, bb+ab, bb+ba
    
    elif text == True:
        if conj==False:
            # RudinShapiro sequence
            aa, ab = "AA", "AB"
            ba, bb = "BA", "BB"
        else:
            # Conjugated RudinShapiro sequence
            aa, ab = "BB", "BA"
            ba, bb = "AB", "AA"
        for i in range(1, n):
            if verbose: print("RS({0}) = {1}".format(i, aa))
            aa, ab, ba, bb = aa+ab, aa+ba, bb+ab, bb+ba
        if verbose: print("RS({0}) = {1}".format(n, aa))
        
    return aa

            
def ThueMorse(n, conj=False, text=False, verbose=False):
    """    
    Generates a Thue-Morse replacement sequence.
    A -> AB, B -> BA.
    """
    """
    # length of each unit
    nh = A_m.n()
    nl = B_m.n()
    n0 = nh-(nh-nl)/2.0       # no=(nh+nl)/2, Deltan=nh-n0=n0-nl
    d_A = (get_lambda()/4./n0).real
    d_B = (get_lambda()/4./n0).real
    
    ThueMorse(0, text=True) -> A
    ThueMorse(1, text=True) -> AB
    ThueMorse(2, text=True) -> ABBA
    ThueMorse(3, text=True) -> ABBABAAB
    """
    if text == False:
        if conj==False:
            # ThueMorse sequence.
            a, b = A(d_A), B(d_B)
        else:
            # Conjugated ThueMorse sequence.
            a, b = B(d_B), A(d_A)
        for i in range(0, n):
            a, b = a+b, b+a            
    
    elif text == True:
        if conj==False:
            # ThueMorse sequence.
            a, b = "A","B"
        else:
            # Conjugated ThueMorse sequence.
            a, b = "B","A"
        for i in range(0, n):
            if verbose: print("TM({0}) = {1}".format(i,a))
            a, b = a+b, b+a
        if verbose: print("TM({0}) = {1}".format(n,a))
        
    return a

    
def binarize(inpt, subs, verbose=False):
    """
    inpt = "ABBA"
    subs= {"A":1, "B":0}
    out = "1001"
    """
    
    if type(inpt) is type(""):
        out = inpt
        for ikey in subs:
            out = out.replace(ikey, subs[ikey].__str__())
        
        if verbose:
            print("Substitutions:")
            for ikey in subs:
                print(ikey, " -> " , subs[ikey])
            print("Input=\t", inpt)
            print("Output=\t", out)
    
    return(out)
    
def get_matrix(seq, nm=(5,5), subs = {"A":1, "B":0}):
    """
    Converts a sequence to a binary nm=(nrow x ncolum) matrix.
    get_matrix(libdiscretemath.ThueMorse(8, text=True), (16,16))
    """
    seq = ''.join(map(str, seq))
    abin = np.array([int(i) for i in list(binarize(seq, subs))])
    res = abin.reshape(nm)
    return(res)
