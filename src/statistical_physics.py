#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
statistical_physics.py
Created on Fri Mar 11 12:53:36 2022

"""
from sympy import*
from libreflection import *
import libphyscon as pc

exec(open("../src/libreflection.py").read())

class statistical_mechanics(branch):
    """

    """
    _name = "statistical_mechanics"
    class_type = {1:"micro_canonical_discrete_distinguihable",
                  2:"micro_canonical_discrete_indistinguihable",
                  3:"micro_canonical_continuous_indistinguihable",
                  4:"canonical",
                  5:"grand_canonical"}[1]
    
    def define_symbols(self):
        """
        Common global symbols, functions.
        i:   state index.
        Zsp: Single particle partition function.
        """
        global alpha,beta,gamma,phi,theta
        global i,j
        global eng,n,p_i,q_i,r
        global h,kB,m,M,N,T
        global g,engF,_Zsp,_ZN,_Cv,_F,_M,_S,_U
        global B,V
        
        # Global Symbols
#        i,j  = symbols('i j', cls=Idx)
        i,j = symbols('i j', integer=True)
        B,V = symbols('B V', real=True)
        eng,n,p_i,q_i,r = symbols('varepsilon n p_i q_i r', real=True)
        alpha,beta,gamma,phi,theta = symbols('alpha beta gamma phi theta', real=True)
        h,kB,m,M,N,T = symbols('h k_B m M N T', real=True, positive=True)
        
        # Global Functions
        
        if self.class_type in ["micro_canonical_discrete_distinguihable",
                               "micro_canonical_discrete_indistinguihable",
                               "micro_canonical_continuous_indistinguihable"]:
            g   = Function('g')(i)           # Degeneracy function.
            engF= Function('varepsilon')(i)  # Energy function.
            _Zsp= Function('Z_sp')(engF, T)  # Partition function of a single particle.
            _ZN = Function('Z_N')(engF, T)   # Partition function of N particles.
            _Cv = Function('C_v')(T)         # Heat capacity at constant volume.
            _F  = Function('F')(T,B)         # Helmholtz free energy.
            _M  = Function('M')(T,B)         # Magnetization.
            _S  = Function('S')(T)           # Entropy.
            _U  = Function('U')(T)           # Internal energy.
        
        elif self.class_type in ["canonical"]:
            g   = Function('g')(i)           # Degeneracy function.
            engF= Function('varepsilon')(i)  # Energy function.
            _Zsp= Function('Z_sp')(N,engF,T) # Partition function of a single particle.
            _ZN = Function('Z_N')(N,engF,T)  # Partition function of N particles.
            _Cv = Function('C_v')(N,T)       # Heat capacity at constant volume.
            _F  = Function('F')(N,T,B)       # Helmholtz free energy.
            _M  = Function('M')(N,T,B)       # Magnetization.
            _S  = Function('S')(N,T)         # Entropy.
            _U  = Function('U')(N,T)         # Internal energy.

    def __init__(self):
        super().__init__()
        self.define_symbols()
        
        class subformulary:
            """
            Sub formulary class.
            
            Define global symbols in the outer class.
            """
            def __init__(self):
                # List of Moment of Inertia
                self.Z_Ideal_Gas = V**N/factorial(N)*(2*pi*m*kB*T/h**2)**(3*N/2)
                self.Z_Ideal_Gas_beta = V**N/(factorial(N)*h**(3*N))*(2*m*pi/beta)**(3*N/2)
        self.subformulary = subformulary()

        if self.class_type == "micro_canonical_discrete_distinguihable":
            self.Zsp = Eq( _Zsp, Sum(g*exp(-engF/(kB*T)), (i,j,n)) )
            self.U   = Eq( _U,  N*kB*T**2*diff(log(self.Zsp.rhs), T, evaluate=False) )
            self.S   = Eq( _S,  N*kB*log(self.Zsp.rhs) + N*kB*T*diff(log(self.Zsp.rhs), T, evaluate=False) )
            self.F   = Eq( _F, -N*kB*T*log(self.Zsp.rhs) )
            self.Cv  = Eq( _Cv,  diff(self.U.rhs, T, evaluate=False) )
            self.M   = Eq( _M,  -diff(self.F.rhs, B, evaluate=False) )
           
        if self.class_type == "micro_canonical_discrete_indistinguihable":
            self.Zsp = Eq( _Zsp, Sum(g*exp(-engF/(kB*T)), (i,j,n)) )
            self.ZN  = Eq( _ZN, self.Zsp.rhs**N/factorial(N) )
            self.U   = Eq( _U,  N*kB*T**2*diff(log(self.Zsp.rhs), T, evaluate=False) )
            self.S   = Eq( _S,  N*kB*log(self.Zsp.rhs) + N*kB*T*diff(log(self.Zsp.rhs), T, evaluate=False) - kB*log(factorial(N)) )
            self.F   = Eq( _F, -N*kB*T*log(self.Zsp.rhs) + kB*T*log(factorial(N)) )
            self.Cv  = Eq( _Cv, diff(self.U.rhs, T, evaluate=False) )

        if self.class_type == "micro_canonical_continuous_indistinguihable":
            self.Zsp = Eq( _Zsp, Integral(g*exp(-engF/(kB*T)), (eng,0,oo)) )
            self.ZN  = Eq( _ZN, self.Zsp.rhs**N/factorial(N) )
            self.U   = Eq( _U,  N*kB*T**2*diff(log(self.Zsp.rhs), T, evaluate=False) )
            self.S   = Eq( _S,  N*kB*log(self.Zsp.rhs) + N*kB*T*diff(log(self.Zsp.rhs), T, evaluate=False) - kB*log(factorial(N)) )
            self.F   = Eq( _F, -N*kB*T*log(self.Zsp.rhs) + kB*T*log(factorial(N)) )
            self.Cv  = Eq( _Cv, diff(self.U.rhs, T, evaluate=False) )
            
        if self.class_type == "canonical":
#            self.ZN  = Eq( _ZN, (1/(factorial(N)*h**(3*N)))*Integral(exp(-engF/(kB*T)), (p_i,-oo,oo), (q_i,-oo,oo)) )
            self.ZN  = Eq( _ZN, self.subformulary.Z_Ideal_Gas) # todo sil
            self.F   = Eq( _F, -kB*T*log(self.ZN.rhs))
        
        if self.class_type == "grand_canonical":
            pass

        # Common text definitions.
        self.internal_energy = self.U
        self.entropy = self.S
        self.Helmholtz_free_energy = self.F
        if hasattr(self, "Zsp"): self.partition_function_sp = self.Zsp
        if hasattr(self, "ZN"):  self.partition_function_Np = self.ZN
        if hasattr(self, "M"):   self.magnetization = self.M
        
    @staticmethod
    def __doc__():
        return("Document of statistical_mechanics class.")
        
ostat = statistical_mechanics()