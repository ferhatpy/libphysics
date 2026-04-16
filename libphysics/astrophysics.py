# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
astrophysics.py
Created on Fri Mar 11 12:53:36 2022

"""
from sympy import symbols, Function, S, Eq, diff, integrate, solve, sqrt, sin, cos, pi, exp, Rational
from .libsympy import *
from libphysics.libreflection import branch
import libphysics.libphyscon as pc

# Global symbols
alpha, beta, gamma, phi, theta = symbols('alpha beta gamma phi theta', real=True)
a, b, c, d, r   = symbols('a b c d r', real=True)
k, m, t, w = symbols('k m t w', real=True, positive=True)
F, M, T = symbols('F M T', real=True)
vA, vB, vC, vD = symbols('A B C D', vector=True)
Ax, Bx, Cx, Dx = symbols('A_x B_x C_x D_x', real=True)
Ay, By, Cy, Dy = symbols('A_y B_y C_y D_y', real=True)
Az, Bz, Cz, Dz = symbols('A_z B_z C_z D_z', real=True)

x = Function('x')(t)
U  = Function('U')(T)
_U = Function('U')(T)

__all__ = [
    'astrophysics', 'oastr',
    'alpha', 'beta', 'gamma', 'phi', 'theta', 'a', 'b', 'c', 'd', 'r',
    'k', 'm', 't', 'w', 'F', 'M', 'T', 'vA', 'vB', 'vC', 'vD',
    'Ax', 'Bx', 'Cx', 'Dx', 'Ay', 'By', 'Cy', 'Dy', 'Az', 'Bz', 'Cz', 'Dz',
    'x', 'U', '_U'
]

class astrophysics(branch):
    """

    """
    _name = "astrophysics"
    
    def define_symbols(self):
        """
        Common global symbols, functions.
        """
        # Symbols are now defined at the module level for better IDE support.
        pass
    
    def __init__(self, class_type='default'):
        super().__init__()
        self.class_type = class_type
        self.define_symbols()
        
        # File settings
        self.input_dir  = "input/astrophysics"
        self.output_dir = "output/astrophysics"
        
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
        
oastr = astrophysics()
