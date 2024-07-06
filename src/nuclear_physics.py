#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nuclear_physics.py
Created on Fri Mar 11 12:53:36 2022

"""
from sympy import*
from sympy.abc import*
from sympy.assumptions.assume import global_assumptions
from sympy.assumptions import assuming, Q, ask
from sympy.integrals.manualintegrate import manualintegrate
from sympy.integrals.manualintegrate import integral_steps
#from sympy.integrals.rubi.utility_function import Simplify
from sympy.integrals.transforms import inverse_fourier_transform
from sympy.diffgeom.rn import *
from sympy.vector import CoordSys3D

from sympy.physics.quantum import *
from sympy.physics.quantum.cartesian import *
from sympy.physics.quantum.constants import * # hbar etc.
from sympy.physics.quantum.cg import *
from sympy.physics.quantum.dagger import Dagger
from sympy.physics.quantum.operator import *
from sympy.physics.quantum.qapply import *
from sympy.physics.quantum.represent import *
from sympy.physics.quantum.sho1d import *
from sympy.physics.quantum.spin import *
from sympy.physics.quantum.state import *
from sympy.physics.paulialgebra import *
from sympy.physics.paulialgebra import Pauli, evaluate_pauli_product

from libreflection import *

class nuclear_physics(branch):
    """

    """
    _name = "nuclear_physics"
    
    def define_symbols(self):
        """
        Common global symbols, functions.
        a: 
        F: 
        """
        global alpha,beta,gamma,phi,theta
        global vA,vB,vC,vD
        global Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz
        global x
        global U
        global _U
        
        # Real symbols
        # Elementary particles
        global e, pe, nu_e, anu_e
        global m_e, m_n, m_p, m_nu_e, m_nubar_e
        e, pe, nu_e, anu_e = symbols('e^- e^+ nu_e nubar_e', real=True)
        m_e, m_n, m_p, m_nu_e, m_nubar_e = symbols('m_e m_n m_p m_nu_e m_nubar_e', real=True, positive=True)
        
        global m,w,t
        m,w,t = symbols('m w t', real=True)
        global x,y,z, xmin, xmax
        x,y,z = symbols('x y z', real=True)
        xmin,xmax = symbols('x_{min} x_{max}', real=True)
        global V0
        V0 = symbols('V_0', real=True)
        global vA,vB,vC,vD
        vA,vB,vC,vD = symbols('A B C D', vector=True)
        global Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz
        Ax,Bx,Cx,Dx = symbols('A_x B_x C_x D_x', real=True)
        Ay,By,Cy,Dy = symbols('A_y B_y C_y D_y', real=True)
        Az,Bz,Cz,Dz = symbols('A_z B_z C_z D_z', real=True)

        # Numerical values
        global total
        
        # Function symbols
        global f
        f = Function('f')
        U = Function('U')(T)  # Function is accesible out of the module.
        _U = Function('U')(T) # Function is not accesible out of the module.
        global En, En0,  psi0b, psi0k
        psi0b, psi0k = [lambda n=n:Bra(f"psi_{n}^0"), lambda n=n:Ket(f"psi_{n}^0")]
        En  = lambda n=n: Function(var(f'E_{n}'))(n) # A callable function.
        En0 = lambda n=n: Eq(S(f'E_{n}'), Function(var(f'E_{n}^0'))(n)) # A callable function.
        # En0 = lambda n=n: Eq(S(f'E_{n}'), Function(var('E_n^0'))(n)) # A callable function.
    
    def __init__(self, class_type='default'):
        super().__init__()
        self.class_type = class_type
        self.define_symbols()
        
        # File settings
        self.input_dir  = "input/nuclear_physics"
        self.output_dir = "output/nuclear_physics"
        
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
            self.beta_negative_decay = Eq(n, p+e+anu_e)
            self.beta_positive_decay = Eq(p, n+pe+nu_e)

#### # Global Methods
    #### Ket state to wave function convertion.
    def Q(self, react, substitutions):
        """
        Calculate Q of a reaction.
        
        USAGE
        =====
        Returns:
        """
        m_reac = react.xreplace(substitutions)
        Q = Eq(var('Q'), (m_reac.lhs - m_reac.rhs)*c**2)
        
        if self.verbose:
            libsympy.pprints(react,
                             m_reac,
                             Q,
                             self.output_style)
        return Q
        
    @staticmethod
    def __doc__():
        return("Document of <template> class.")
        
onucl = nuclear_physics()