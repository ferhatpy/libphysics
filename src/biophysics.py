#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
biophysics.py
Created on Fri Mar 11 12:53:36 2022

"""
from sympy import*
from libreflection import *
import libphyscon as pc

class biophysics(branch):
    """

    """
    _name = "biophysics"

    def define_symbols(self):
        """
        Common global symbols, functions.
        a: 
        F: 
        """
        global alpha,beta,gamma,phi,theta
        global a,b,c,d,r
        global k,m,t,w
        global F,M,T
        global vA,vB,vC,vD
        global Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz
        global x
        global U
        global _U
        
        alpha,beta,gamma,phi,theta = symbols('alpha beta gamma phi theta', real=True)
        a,b,c,d,r   = symbols('a b c d r', real=True)
        k,m,t,w = symbols('k m t w', real=True, positive=True)
        F,M,T = symbols('F M T', real=True)
        vA,vB,vC,vD = symbols('A B C D', vector=True)
        Ax,Bx,Cx,Dx = symbols('A_x B_x C_x D_x', real=True)
        Ay,By,Cy,Dy = symbols('A_y B_y C_y D_y', real=True)
        Az,Bz,Cz,Dz = symbols('A_z B_z C_z D_z', real=True)
        
        x = Function('x')(t)
        U  = Function('U')(T)  # Function is accesible out of the module.
        _U = Function('U')(T)  # Function is not accesible out of the module.
    
    def __init__(self, class_type='default'):
        super().__init__()
        self.class_type = class_type
        self.define_symbols()
        
        # File settings
        self.input_dir  = "input/biophysics"
        self.output_dir = "output/biophysics"
        
        class subformulary:
            """
            Sub formulary class.
            
            Define global symbols in the outer class.
            """
            def __init__(self):
                # List of Moment of Inertia
                self.Icm_sphere = S(2)/5*M*r**2
        self.subformulary = subformulary()
        
        if self.class_type == "default":
            self.NewtonsLaw2 = Eq(F, m*a)
            self.HookesLaw   = Eq(F, -k*x)

        
    @staticmethod
    def __doc__():
        return("Document of <template> class.")
        
obiop = biophysics()
