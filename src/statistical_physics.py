#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
statistical_physics.py
Created on Fri Mar 11 12:53:36 2022

"""
import mpmath as mp
import numpy as np
import scipy.constants as pc
from sympy import*
from libreflection import *
from sympy.physics.quantum.constants import * # hbar etc.

class statistical_mechanics(branch):
    """

    """
    _name = "statistical_mechanics"
    global numeric
    
    def define_symbols(self):
        """
        Common global symbols, functions.
        i:   state index.
        Zsp: Single particle partition function.
        """
        # Integer symbols
        global i,j
#        i,j  = symbols('i j', cls=Idx)
        i,j = symbols('i j', integer=True)
        global B,V
        B,V = symbols('B V', real=True)
        global eng,n,p_i,q_i,r
        eng,n,p_i,q_i,r = symbols('varepsilon n p_i q_i r', real=True)
        global alpha,beta,gamma,phi,theta
        alpha,beta,gamma,phi,theta = symbols('alpha beta gamma phi theta', real=True)
        global h,kB,m,M,N,T
        h,kB,m,M,N,T = symbols('h k_B m M N T', real=True, positive=True)
        
        # Global Functions
        global g, engF
        if self.class_type in ["micro_canonical_discrete_distinguihable",
                               "micro_canonical_discrete_indistinguihable",
                               "micro_canonical_continuous_indistinguihable"]:
            g   = Function('g')(i)           # Degeneracy function.
            engF= Function('varepsilon')(i)  # Energy function.
        
        elif self.class_type in ["canonical"]:
            g   = Function('g')(i)           # Degeneracy function.
            engF= Function('varepsilon')(i)  # Energy function.

    def __init__(self, class_type='micro_canonical_discrete_indistinguihable', numeric=False):
        """
        class_type = \
        {1:"micro_canonical_discrete_distinguihable",
         2:"micro_canonical_discrete_indistinguihable",
         3:"micro_canonical_continuous_indistinguihable",
         4:"canonical",
         5:"grand_canonical"}
        """
        super().__init__()
        self.class_type = class_type
        self.numeric = numeric
        self.define_symbols()
        
        # File settings
        self.input_dir  = "input/statistical_physics"
        self.output_dir = "output/statistical_physics"
        
        class subformulary:
            """
            Sub formulary class.
            
            Define global symbols in the outer class.
            """
            def __init__(self):
                self.Z_Ideal_Gas = V**N/factorial(N)*(2*pi*m*kB*T/h**2)**(3*N/2)
                self.Z_Ideal_Gas_beta = V**N/(factorial(N)*h**(3*N))*(2*m*pi/beta)**(3*N/2)
        self.subformulary = subformulary()


#### 1) micro_canonical_discrete_distinguihable
        if self.class_type == "micro_canonical_discrete_distinguihable":
            if self.numeric:
                """
                Numeric Formulary
                =================                
                def nZsp(self, engF=lambda _:_, g=lambda _:1, T=1, j=1, n=mp.inf):
                    Usage:
                        Zsp(lambda i:i**2)             # g(i)=1
                        Zsp(lambda i:i**2, lambda i:1) # g(i)=1
                        
                    Zsp() === Zsp(engF=i, g=1, T=1, j=1, n=mp.inf)
                    Zsp(lambda _:_**2) === Zsp(lambda i:i**2)
                    
                    Equivalent single-line function:
                    Zsp = lambda engF=lambda _:_, g=lambda _:1, T=1, j=1, n=mp.inf: mp.nsum(lambda i: g(i)*mp.exp(-engF(i)/(T)), [1, n])
                    return mp.nsum(lambda i: g(i)*mp.exp(-engF(i)/(pc.k*T)), [j, n])
                """
                # mp.dps = 15 # Set numerical precission
                print("mpmath library is being used !")
                
                self.Zsp = lambda engF=lambda _:_, g=lambda _:1, T=1, j=1, n=mp.inf, kB=1: (
                    mp.nsum(lambda i: g(i)*mp.exp(-engF(i)/(kB*T)), [j, n]) )
                
                self.ZN  = lambda engF=lambda _:_, g=lambda _:1, T=1, j=1, n=mp.inf, kB=1, N=2: (
                    self.Zsp(engF, g, T, j, n, kB)**N )
                
                self.U   = lambda engF=lambda _:_, g=lambda _:1, T=1, j=1, n=mp.inf, kB=1, N=1, point=1.0: (
                    N*kB*T**2*mp.diff(lambda T: mp.log(self.Zsp(engF, g, T, j, n, kB)), point) )
                
                self.S   = lambda engF=lambda _:_, g=lambda _:1, T=1, j=1, n=mp.inf, kB=1, N=1, point=1.0: (
                    N*kB*mp.log(self.Zsp(engF, g, T, j, n, kB)) + N*kB*T*mp.diff(lambda T: mp.log(self.Zsp(engF, g, T, j, n, kB)), point) )
                
                self.F   = lambda engF=lambda _:_, g=lambda _:1, T=1, j=1, n=mp.inf, kB=1, N=1: (
                    -N*kB*T*mp.log(self.Zsp(engF, g, T, j, n, kB)))
                
                self.Cv  = lambda engF=lambda _:_, g=lambda _:1, T=1, j=1, n=mp.inf, kB=1, N=1, point=1.0: (
                    mp.diff(lambda T: self.U(engF, g, T, j, n, kB, N, point), point))
                
                # kaldik todo check whether B must be an argument of the lambda function
                self.M   = lambda engF=lambda _:_, g=lambda _:1, T=1, j=1, n=mp.inf, kB=1, N=1, point=1.0: (
                    -mp.diff(lambda B: self.F(engF, g, T, j, n, kB, N), point))

            else:
                """
                Analytic Formulary
                ==================  
                """
                self.Zsp = Eq(symbols('Z_sp'), Sum(g*exp(-engF/(kB*T)), (i, j, n)))
                self.ZN  = Eq(symbols('Z_N'),  self.Zsp.rhs**N)
                self.U   = Eq(symbols('U'),   N*kB*T**2*diff(log(self.Zsp.rhs), T, evaluate=False))
                self.S   = Eq(symbols('S'),   N*kB*log(self.Zsp.rhs) + N*kB*T*diff(log(self.Zsp.rhs), T, evaluate=False))
                self.F   = Eq(symbols('F'),  -N*kB*T*log(self.Zsp.rhs))
                self.Cv  = Eq(symbols('C_v'), diff(self.U.rhs, T, evaluate=False))
                self.M   = Eq(symbols('M'),  -diff(self.F.rhs, B, evaluate=False))
            

#### 2) micro_canonical_discrete_indistinguihable
        if self.class_type == "micro_canonical_discrete_indistinguihable":
            self.Zsp = Eq(symbols('Z_sp'), Sum(g*exp(-engF/(kB*T)), (i,j,n)))
            self.ZN  = Eq(symbols('Z_N'),  self.Zsp.rhs**N/factorial(N))
            self.U   = Eq(symbols('U'),    N*kB*T**2*diff(log(self.Zsp.rhs), T, evaluate=False))
            self.S   = Eq(symbols('S'),    N*kB*log(self.Zsp.rhs) + N*kB*T*diff(log(self.Zsp.rhs), T, evaluate=False) - kB*log(factorial(N)))
            self.F   = Eq(symbols('F'),   -N*kB*T*log(self.Zsp.rhs) + kB*T*log(factorial(N)))
            self.Cv  = Eq(symbols('C_v'),  diff(self.U.rhs, T, evaluate=False))


#### 3) micro_canonical_continuous_indistinguihable
        if self.class_type == "micro_canonical_continuous_indistinguihable":
            self.Zsp = Eq(symbols('Z_sp'), Integral(g*exp(-engF/(kB*T)), (eng,0,oo)))
            self.ZN  = Eq(symbols('Z_N'), self.Zsp.rhs**N/factorial(N))
            self.U   = Eq(symbols('U'),  N*kB*T**2*diff(log(self.Zsp.rhs), T, evaluate=False))
            self.S   = Eq(symbols('S'),  N*kB*log(self.Zsp.rhs) + N*kB*T*diff(log(self.Zsp.rhs), T, evaluate=False) - kB*log(factorial(N)) )
            self.F   = Eq(symbols('F'), -N*kB*T*log(self.Zsp.rhs) + kB*T*log(factorial(N)))
            self.Cv  = Eq(symbols('C_v'), diff(self.U.rhs, T, evaluate=False))


#### 4) canonical            
        if self.class_type == "canonical":
            self.ZN_ideal_gas  = Eq(symbols('Z_N(Ideal gas)'), (1/(factorial(N)*h**(3*N)))*Integral(exp(-engF/(kB*T)), (p_i,-oo,oo), (q_i,-oo,oo)) )
            self.ZN  = Eq(symbols('Z_N'), self.subformulary.Z_Ideal_Gas) # todo sil
            self.F   = Eq(symbols('F'), -kB*T*log(self.ZN.rhs))


#### 5) grand_canonical        
        if self.class_type == "grand_canonical":
            pass

        # Common text definitions.
        # self.internal_energy = self.U
        # self.entropy = self.S
        # self.Helmholtz_free_energy = self.F
        if hasattr(self, "Zsp"): self.partition_function_sp = self.Zsp
        if hasattr(self, "ZN"):  self.partition_function_Np = self.ZN
        if hasattr(self, "M"):   self.magnetization = self.M
        
    @staticmethod
    def __doc__():
        return("Document of statistical_mechanics class.")
        
ostat = statistical_mechanics()