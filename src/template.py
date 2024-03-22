#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
template.py todo copy from mechanics!!!
Created on Fri Mar 11 12:53:36 2022

Find and replace template with desired class name.

Mathematica Equivalance
=======================
Newton2 := EqualTo[F][m a]
TraditionalForm[Sort[Newton2]]
Solve[Newton2, a]

todo correct

"""
from sympy import*
from sympy.abc import*
from sympy.diffgeom import *
from sympy.diffgeom.rn import *
from sympy.diffgeom.rn import R3_r, R3_s
from sympy.physics.quantum.constants import * # hbar etc.
from sympy.physics.vector import *
from sympy.plotting import plot_parametric
from sympy.vector import CoordSys3D

from libreflection import *
import libphyscon as pc

#exec(open("../src/libreflection.py").read())

class template(branch):
    """

    """
    _name = "template"

    global g,engF

    def define_symbols(self):
        """
        Common global symbols, functions.
        a: 
        F: 
        """
        global C
        global C1,C2,C3 # Integration constants.
        global alpha,beta,gamma,phi,theta
        global a,b,c,d,r
        global k,m,t,tau,w
        global F,M,T
        global vA,vB,vC,vD
        global Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz
        global fV, fE, psi
        global x
        global U
        global _U
        
        # Global Indeces
        i,j  = symbols('i j', cls=Idx)
        u = IndexedBase('u')     # Generates an error.
        
        # Global Integer symbols
        [i, j] = symbols('i j', positive=True, integer=True, nonzero=True)
        
        # Global Real symbols
        alpha,beta,gamma,phi,theta = symbols('alpha beta gamma phi theta', real=True)
        a,b,c,d,r = symbols('a b c d r', real=True)
        k,m,t,tau,w = symbols('k m t tau w', real=True, positive=True)
        F,M,T = symbols('F M T', real=True)
        vA,vB,vC,vD = symbols('A B C D', vector=True)
        Ax,Bx,Cx,Dx = symbols('A_x B_x C_x D_x', real=True)
        Ay,By,Cy,Dy = symbols('A_y B_y C_y D_y', real=True)
        Az,Bz,Cz,Dz = symbols('A_z B_z C_z D_z', real=True)
        tau         = Symbol('tau',  real=True, positive=True)
        #w = Symbol("w", positive=True, integer=True, nonzero=True)
        
        # Global Functions
        u  = dynamicsymbols('u') # gives time dependent function
        x  = Function('x')(t)
        U  = Function('U')(T)    # Function is accesible out of the module.
        _U = Function('U')(T)    # Function is not accesible out of the module.
        Gnon  = Function('G')(t,tau)                    # A noncallable function.
        Gcal1 = Lambda((t,tau), Function('G')(t,tau))   # A callable function.
        Gcal2 = lambda t,tau: Function('G')(t,tau)      # A callable function.
        lst_functions = ['fV','fE','psi']
        [fV, fE, psi] = [Function(ifun)(x) for ifun in lst_functions]
    
    
        # Global Functions
        if self.class_type in ["micro_canonical_discrete_distinguihable",
                               "micro_canonical_discrete_indistinguihable",
                               "micro_canonical_continuous_indistinguihable"]:
            g   = Function('g')(i)           # Degeneracy function.
            engF= Function('varepsilon')(i)  # Energy function.
        
        elif self.class_type in ["canonical"]:
            g   = Function('g')(i)           # Degeneracy function.
            engF= Function('varepsilon')(i)  # Energy function.
    
    	# Common definitions.
        if self.class_type in ["scalar", "vectorial"]:
            _H = Function('H')(t)           # Total energy.
            
    def __init__(self, class_type='scalar'):
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
        self.define_symbols()
        
        # File settings
        self.input_dir  = "input/template"
        self.output_dir = "output/template"
        
        class subformulary:
            """
            Sub formulary class.
            
            Define global symbols in the outer class.
            """
            def __init__(self):
                # List of Moment of Inertia
                self.Icm_sphere = S(2)/5*M*r**2
                
        self.subformulary = subformulary()
        
        #### Quantum Harmonic Oscillator
        class sho(branch):
            """
            Sub class for quantum simple harmonic oscillator.
            """
            def __init__(self):
                super().__init__()
                self.name = "Quantum Harmonic Oscillator"
                self.xi = Eq(S('xi'), sqrt(m*w/hbar)*x)
                self.psi = lambda n=n: Eq(S(f'psi_{n}'), (m*w/(pi*hbar))**(S(1)/4)*(1/sqrt((2**n)*factorial(n)))*hermite(n, xi)*exp(-xi**2/2))
        self.sho = sho()
        
        if self.class_type == "scalar":
            # Construct a cascaded formulary structure.
            self.NewtonsLaw2_1 = Eq(var('F'), m*a) # Default one is var.
            self.NewtonsLaw2_2 = Eq(symbols('F'), m*a)
            self.NewtonsLaw2_3 = Eq(S('F'), m*a)
            self.HookesLaw   = Eq(F, -k*x)
            self.S           = Eq(F, -k*x)
         
        # Common text definitions.
        self.entropy = self.S
        if hasattr(self, "Zsp"): self.partition_function_sp = self.Zsp
    
    @staticmethod
    # Abstract method
    def __doc__():
        return("Document of template class.")
        
otemp = template() # Create an otemp object from template class.