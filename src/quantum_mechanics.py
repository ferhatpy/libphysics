#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quantum_mechanics.py
Created on Fri Mar 11 12:53:36 2022

"""
from sympy import*
from sympy.abc import*
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
        global x,y,z, xmin, xmax
        global phi, Psi, psix
        global vA,vB,vC,vD
        global Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz
        global _U
        
        # Global Integer symbols
        global n,m
        n,m = symbols('n m', positive=True, integer=True, nonzero=True)
        
        # Global Real symbols
        global a,b
        A,a,b = symbols('A a b', real=True, positive=True)
        x,y,z = symbols('x y z', real=True)
        xmin,xmax = symbols('x_{min} x_{max}', real=True)
        k,m,t = symbols('k m t', real=True)
        vA,vB,vC,vD = symbols('A B C D', vector=True)
        Ax,Bx,Cx,Dx = symbols('A_x B_x C_x D_x', real=True)
        Ay,By,Cy,Dy = symbols('A_y B_y C_y D_y', real=True)
        Az,Bz,Cz,Dz = symbols('A_z B_z C_z D_z', real=True)
        
        # Global Functions
        global f
        f = Function('f')
        U  = Function('U')(T)  # Function is accesible out of the module.
        _U = Function('U')(T)  # Function is not accesible out of the module.

        if self.class_type in ["position_space"]:
            Psi  = Function('Psi')(x,y,z,t)
            psix = Function('psi')(x)
            phi  = Function('phi')(k)
            
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
            def __init__(self, class_type):
                # List of Moment of Inertia
                self.Icm_sphere = S(2)/5*M*r**2
                
                if class_type in ["position_space"]:
                    self.psi_infinite_qw = Eq(var(r'\psi_{n}(x)'), sqrt(2/a)*sin(n*pi*x/a))
                    self.E_inifinite_wq = Eq(var('E_n'), (n**2*pi**2*hbar**2)/(2*m*a**2))
                
        self.subformulary = subformulary(self.class_type)
        
        if self.class_type in ["position_space"]:
            self.Psi = Psi
            self.psix = psix
            self.normalization = Eq(Integral(conjugate(self.Psi)*self.Psi, (x,xmin,xmax)), 1)
            
            # Position related formulary.
            self.exp_x    = Eq(var(r'\langle{x}\rangle'),   Integral(conjugate(self.Psi)*x*self.Psi, (x,xmin,xmax)))
            self.exp_x2   = Eq(var(r'\langle{x^2}\rangle'), Integral(conjugate(self.Psi)*x**2*self.Psi, (x,xmin,xmax)))
            self.delta_x  = Eq(var(r'\Delta{x}'), sqrt(self.exp_x2.rhs - self.exp_x.rhs**2))
            self.delta_x2 = Eq(var(r'(\Delta{x})^2'), self.exp_x2.rhs - self.exp_x.rhs**2)
            
            # Momentum related formulary.
            self.px      = Eq(var(r'\hat{p}_x'),  -I*hbar*Derivative(self.Psi, x))
            self.py      = Eq(var(r'\hat{p}_y'),  -I*hbar*Derivative(self.Psi, y))
            self.pz      = Eq(var(r'\hat{p}_z'),  -I*hbar*Derivative(self.Psi, z))
            self.px2     = Eq(var(r'\hat{p}_x^2'), (-I*hbar)**2*Derivative(self.Psi, x, 2))
            self.exp_px  = Eq(var(r'\langle{p_x}\rangle'),   Integral(conjugate(self.Psi)*self.px.rhs,  (x,xmin,xmax)))
            self.exp_px2 = Eq(var(r'\langle{p_x^2}\rangle'), Integral(conjugate(self.Psi)*self.px2.rhs, (x,xmin,xmax)))
            self.delta_px  = Eq(var(r'\Delta{p_x}'), sqrt(self.exp_px2.rhs - self.exp_px.rhs**2))
            self.delta_px2 = Eq(var(r'(\Delta{p_x})^2'), self.exp_px2.rhs - self.exp_px.rhs**2)
            self.delta_XP  = self.uncertainityXP  = Eq(var(r'\Delta{x}\Delta{p_x}'), self.delta_x.rhs*self.delta_px.rhs)
            
            # Orbital angular momentum formulary.
#            self.Lx      = Eq(var('L_x'), var('L_x'))
#            self.Ly      = Eq(var('L_y'), var('L_y'))
            self.Lplus   = Eq(var('L_+'), var('L_x') + I*var('L_y'))
            self.Lminus  = Eq(var('L_-'), var('L_x') - I*var('L_y'))
            
            # Operator Definitions via SymPy sympy.physics.quantum class todo ERROR in DifferantialOperator application
            # Position related formulary.
            self.exp_opX   = Eq(var(r'\langle{x}\rangle'),   Integral(conjugate(self.Psi)*qapply(x*self.Psi), (x,xmin,xmax)))
            self.exp_opX2  = Eq(var(r'\langle{x^2}\rangle'), Integral(conjugate(self.Psi)*qapply(x**2*self.Psi), (x,xmin,xmax)))
            self.delta_opX = Eq(var(r'\Delta{x}'), sqrt(self.exp_opX2.rhs - self.exp_opX.rhs**2))
            self.delta_opX2= Eq(var(r'(\Delta{x})^2'), self.exp_opX2.rhs - self.exp_opX.rhs**2)
            # Momentum related formulary.
            self.opPx     =  Eq(var(r'\hat{p}_x'), DifferentialOperator(-I*hbar*Derivative(self.Psi, x), self.Psi))
            self.exp_opPx  = Eq(var(r'\langle{\hat{p}_x}\rangle'),   Integral(conjugate(self.Psi)*qapply(self.opPx.rhs*self.Psi), (x,xmin,xmax)))
            self.exp_opPx2 = Eq(var(r'\langle{\hat{p}_x^2}\rangle'), Integral(conjugate(self.Psi)*qapply(self.opPx.rhs*self.opPx.rhs*self.Psi), (x,xmin,xmax)))
            self.delta_opPx  = Eq(var(r'\Delta{p_x}'), sqrt(self.exp_opPx2.rhs - self.exp_opPx.rhs**2))
            self.delta_opPx2 = Eq(var(r'(\Delta{p_x})^2'), self.exp_opPx2.rhs - self.exp_opPx.rhs**2)
            self.delta_opXopPx = self.uncertainityXP  = Eq(var(r'\Delta{x}\Delta{p_x}'), self.delta_opX.rhs*self.delta_opPx2.rhs)
            
            self.V = self.potential_energy = Function('V')(x)
            self.H = Eq(S('H'), -(hbar**2)/(2*m)*diff(self.Psi, x, 2) + self.V*self.Psi)
            self.exp_H = Eq(var(r'\langle{H}\rangle'), Integral(conjugate(self.Psi)*self.H.rhs, (x,xmin,xmax)))
            self.opH = Eq(S('H'), hbar**2/(2*m)*(DifferentialOperator(Derivative(self.Psi, x), self.Psi))**2 + Operator(self.V*self.Psi))
            self.exp_opH = Eq(var(r'\langle{H}\rangle'), integrate(conjugate(self.Psi)*qapply(self.opH.rhs*self.Psi), (x,xmin,xmax)))
            
            self.SchrodingerEq = Eq(-(hbar**2)/(2*m)*diff(self.Psi, x, 2) + self.V*self.Psi, I*hbar*diff(self.Psi, t, 1))
            # Fourier Transforms formulary.
            self.FourierTransform_phi_t0 = Eq( Psi.subs(t,0), 1/(sqrt(2*pi))*Integral(phi*exp(I*k*x) , (k,-oo,oo)) )
            self.FourierTransform_Psi_t0 = Eq( phi, 1/(sqrt(2*pi))*Integral(Psi.subs(t,0)*exp(-I*k*x) , (k,-oo,oo)) )
        
        if self.class_type in ["momentum_space"]:
            pass
        
        # Common text definitions.
        self.Hamiltonian = self.H

        
    @staticmethod
    def __doc__():
        return("Document of <template> class.")
        
oqmec = quantum_mechanics()