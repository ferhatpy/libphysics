#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quantum_mechanics.py
Created on Fri Mar 11 12:53:36 2022

"""
from sympy import*
from sympy.integrals.manualintegrate import manualintegrate
from sympy.integrals.manualintegrate import integral_steps
#from sympy.integrals.rubi.utility_function import Simplify
from sympy.integrals.transforms import inverse_fourier_transform
from sympy.physics.quantum import *
from sympy.physics.quantum.state import *
from sympy.physics.quantum.operator import *
from sympy.physics.quantum.qapply import *

from libreflection import *
import libphyscon as pc

exec(open("../src/libreflection.py").read())

class quantum_mechanics(branch):
    """
    solve(oqmec.expX, _expX)
    """
    _name = "quantum_mechanics"
    class_type = {1:"position_space", 2:"momentum_space"}[1]
    
    def define_symbols(self):
        """
        Common global symbols, functions.
        a: 
        F: 
        """
        global x, xmin, xmax
        global _exp_x, _exp_x2
        global _exp_px, _exp_px2
        global _delta_x2, _delta_px2
        global _px, _py, _pz, _px2, _py2, _pz2
        
        global k,m,t,w
        global Psi, Psix
        global phi
        
        global f
        global alpha,beta,gamma,phi,theta
        global a,b,c,d,r
        global F,M,T
        global vA,vB,vC,vD
        global Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz
        global U
        global _U
        
        # Global Symbols
        x,y,z = symbols('x y z', real=True)
        xmin,xmax = symbols('x_{min} x_{max}', real=True)
        k,m,t = symbols('k m t', real=True)
        
        alpha,beta,gamma,phi,theta = symbols('alpha beta gamma phi theta', real=True)
        a,b,c,d,r   = symbols('a b c d r', real=True)

        F,M,T = symbols('F M T', real=True)
        vA,vB,vC,vD = symbols('A B C D', vector=True)
        Ax,Bx,Cx,Dx = symbols('A_x B_x C_x D_x', real=True)
        Ay,By,Cy,Dy = symbols('A_y B_y C_y D_y', real=True)
        Az,Bz,Cz,Dz = symbols('A_z B_z C_z D_z', real=True)
        
        # Global Functions
        f = Function('f')
        U  = Function('U')(T)  # Function is accesible out of the module.
        _U = Function('U')(T)  # Function is not accesible out of the module.

        if self.class_type in ["position_space"]:
            Psix = Function('Psi')(x,t)
            Psi = Psix
            phi  = Function('phi')(k)
            _exp_x,   _exp_x2  = symbols(r'\langle{x}\rangle, \langle{x^2}\rangle')
            _exp_px, _exp_px2  = symbols(r'\langle{p_x}\rangle, \langle{p_x^2}\rangle')
            _delta_x2, _delta_px2 = symbols(r'\Delta{x^2}, \Delta{p_x^2}')
            _px, _py, _pz = symbols('p_x p_y p_z')
            _px2, _py2, _pz2 = symbols('p^2_x p^2_y p^2_z')
            
        if self.class_type in ["momentum_space"]:
            pass

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
                self.Icm_sphere = S(2)/5*M*r**2
                
        self.subformulary = subformulary()
        
        if self.class_type in ["position_space"]:
            self.Psix = Psix
            self.Psi = Psix
            self.normalization = Eq(Integral(conjugate(self.Psi)*self.Psi, (x,xmin,xmax)), 1)
            # Position related
            self.exp_x    = Eq( _exp_x, Integral(conjugate(self.Psi)*x*self.Psi, (x,xmin,xmax)))
            self.exp_x2   = Eq(_exp_x2, Integral(conjugate(self.Psi)*x**2*self.Psi, (x,xmin,xmax)))
            self.delta_x2 = Eq(_delta_x2, self.exp_x2.rhs - self.exp_x.rhs**2)
            # Momentum related
#            self.px      = Eq( _px,  -I*hbar*Derivative(self.Psi, x))
            self.px      = Eq( _px,  -I*hbar*Derivative(self.Psi, x))
            self.px2     = Eq( _px2, (-I*hbar)**2*Derivative(self.Psi, x, 2))
            self.exp_px  = Eq( _exp_px, Integral(conjugate(self.Psi)*self.px.rhs, (x,xmin,xmax)))
            self.exp_px2 = Eq( _exp_px2, Integral(conjugate(self.Psi)*self.px2.rhs, (x,xmin,xmax)))
            self.delta_px2 = Eq(_delta_px2, self.exp_px2.rhs - self.exp_px.rhs**2)
            # Operator Definitions via SymPy sympy.physics.quantum class 
            self.op_px   = DifferentialOperator(-I*hbar*Derivative(self.Psi, x), self.Psi)
            self.Px      = Eq( _px, self.op_px)
            self.exp_Px  = Eq( _exp_px, Integral(conjugate(self.Psi)*qapply(self.Px.rhs*self.Psi), (x,xmin,xmax)))
            self.exp_Px2 = Eq( _exp_px, Integral(conjugate(self.Psi)*qapply(self.Px.rhs*self.Px.rhs*self.Psi), (x,xmin,xmax)))
            
            self.V = self.potential_energy = Function('V')(x)
            self.SchrodingerEq = Eq(-(hbar**2)/(2*m)*diff(self.Psi, x, 2) + self.V*self.Psi, I*hbar*diff(self.Psi, t, 1))
        
            self.FourierTransform_phi_t0 = Eq( Psi.subs(t,0), 1/(sqrt(2*pi))*Integral(phi*exp(I*k*x) , (k,-oo,oo)) )
            self.FourierTransform_Psi_t0 = Eq( phi, 1/(sqrt(2*pi))*Integral(Psi.subs(t,0)*exp(-I*k*x) , (k,-oo,oo)) )
        
        if self.class_type in ["momentum_space"]:
            pass

        
    @staticmethod
    def __doc__():
        return("Document of <template> class.")
        
oqmec = quantum_mechanics()
