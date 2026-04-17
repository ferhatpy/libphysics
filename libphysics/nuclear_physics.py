# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
nuclear_physics.py
Created on Fri Mar 11 12:53:36 2022

"""
from sympy import symbols, Function, S, Eq, diff, integrate, solve, sqrt, sin, cos, pi, exp, Rational, Matrix
from .libsympy import *
from libphysics.libreflection import branch

# Global symbols
alpha, beta, gamma, phi, theta = symbols('alpha beta gamma phi theta', real=True)
e, pe, nu_e, anu_e = symbols('e^- e^+ nu_e nubar_e', real=True)
p = symbols('p', real=True)  # proton symbol
m_e, m_n, m_p, m_nu_e, m_nubar_e = symbols('m_e m_n m_p m_nu_e m_nubar_e', real=True, positive=True)
m, w, t, T, c, M = symbols('m w t T c M', real=True)
x, y, z = symbols('x y z', real=True)
xmin, xmax = symbols('x_{min} x_{max}', real=True)
V0 = symbols('V_0', real=True)
vA, vB, vC, vD = symbols('A B C D', vector=True)
Ax, Bx, Cx, Dx = symbols('A_x B_x C_x D_x', real=True)
Ay, By, Cy, Dy = symbols('A_y B_y C_y D_y', real=True)
Az, Bz, Cz, Dz = symbols('A_z B_z C_z D_z', real=True)

f = Function('f')
U = Function('U')(T)
_U = Function('U')(T)
psi0b, psi0k = [lambda n=n:Bra(f"psi_{n}^0"), lambda n=n:Ket(f"psi_{n}^0")]
En  = lambda n=n: Function(var(f'E_{n}'))(n)
En0 = lambda n=n: Eq(S(f'E_{n}'), Function(var(f'E_{n}^0'))(n))

__all__ = [
    'nuclear_physics', 'onucl',
    'alpha', 'beta', 'gamma', 'phi', 'theta', 'e', 'pe', 'nu_e', 'anu_e',
    'm_e', 'm_n', 'm_p', 'm_nu_e', 'm_nubar_e', 'm', 'w', 't', 'T', 'c', 'M',
    'x', 'y', 'z', 'xmin', 'xmax', 'V0', 'vA', 'vB', 'vC', 'vD',
    'Ax', 'Bx', 'Cx', 'Dx', 'Ay', 'By', 'Cy', 'Dy', 'Az', 'Bz', 'Cz', 'Dz',
    'f', 'U', '_U', 'psi0b', 'psi0k', 'En', 'En0'
]

class nuclear_physics(branch):
    """

    """
    _name = "nuclear_physics"
    
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
