# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
quantum_mechanics.py
Created on Fri Mar 11 12:53:36 2022

"""
from __future__ import annotations
from typing import Any, Callable, List, Tuple, Union, Optional, Dict
import mpmath as mp
from IPython.display import display
# from abc import ABC, abstractmethod
from sympy.abc import x, y, z, t, k, n
from sympy import (
    symbols, Function, Eq, Integral, Sum, Product, Derivative, diff, integrate, simplify, sqrt, exp, sin, cos, pi, oo, I,
    factorial, Range, Matrix, UnevaluatedExpr, Piecewise, S, var, cot, hermite, assoc_laguerre, Ynm, zeros, log, sinh, conjugate, DiracDelta
)
from sympy import Derivative as D
from sympy.assumptions.assume import global_assumptions
from sympy.assumptions import assuming, Q, ask
from sympy.integrals.manualintegrate import manualintegrate
from sympy.integrals.manualintegrate import integral_steps
#from sympy.integrals.rubi.utility_function import Simplify
from sympy.integrals.transforms import inverse_fourier_transform, fourier_transform
# from sympy.diffgeom import Manifold, Patch, CoordSystem, BaseVectorField
# from sympy.diffgeom.rn import R2, R3
from sympy.vector import CoordSys3D

# Hydrogen Atom
from sympy.physics.hydrogen import Psi_nlm, R_nl, E_nl, E_nl_dirac

from sympy.physics.quantum import (
    Operator, Ket, Bra, Commutator, Dagger, qapply, DifferentialOperator, TensorProduct
)
# No unique direct imports needed from quantum.cartesian
from sympy.physics.quantum.constants import hbar
# No unique direct imports needed from quantum.cg
from sympy.physics.quantum.dagger import Dagger
# from sympy.physics.quantum.operator import Operator, DifferentialOperator
# from sympy.physics.quantum.qapply import qapply
# No unique direct imports needed from quantum.represent
from sympy.physics.quantum.sho1d import SHOKet, SHOBra, RaisingOp, LoweringOp
from sympy.physics.quantum.spin import Jx, Jy, Jz, J2, Jplus, Jminus, JzKet
from sympy.physics.quantum.state import Wavefunction
from sympy.physics.paulialgebra import Pauli, evaluate_pauli_product


from .libsympy import *
from libphysics.libreflection import branch
# import libphyscon as pc

# Global symbols
xi = symbols('xi')
xA, xB, yA, yB, pA, pB = symbols('x_A x_B y_A y_B p_A p_B', real=True)
xmin, xmax, pmin, pmax = symbols('x_min x_max p_min p_max', real=True)
Hp, psib, psik, nb, nk = symbols('Hp psib psik nb nk')
Psi, psi, psix, psixt, PsiSph, psiSph, phik = [Function(s) for s in ['Psi', 'psi', 'psi_x', 'psi_xt', 'Psi_sph', 'psi_sph', 'phi_k']]
En, En0 = symbols('En En0', cls=Function)

__all__ = [
    'quantum_mechanics', 'oqmec', 'Wavefunction',
    'xi', 'xA', 'xB', 'yA', 'yB', 'pA', 'pB', 'xmin', 'xmax', 'pmin', 'pmax',
    'Hp', 'psib', 'psik', 'nb', 'nk', 'Psi', 'psi', 'psix', 'psixt', 'PsiSph', 'psiSph', 'phik',
    'En', 'En0', 'eps0', 'r', 'x', 'y', 'z', 'p', 't', 'hbar', 'l', 'k', 'm', 'w', 'w0', 'M', 'm_e', 'm_p', 'm_n', 'm_N', 'V', 'psi0b', 'psi0k', 'C'
]

# Real symbols
alpha, beta, gamma, theta, phi = symbols('alpha beta gamma theta phi', real=True)
eps0 = symbols('varepsilon_0', real=True, positive=True)
r, x, y, z, p = symbols('r x y z p', real=True)
t = symbols('t', real=True, positive=True)
hbar = symbols('hbar', real=True, positive=True)
l = symbols('lambda', real=True)
k, m, w, w0, M = symbols('k m w w0 M', real=True, positive=True)
m_e, m_p, m_n, m_N = symbols('m_e m_p m_n m_N', real=True, positive=True)

# Function symbols
V = Function('V')(x,y,z)
psi0b, psi0k = [lambda n=n:Bra(f"psi_{n}^0"), lambda n=n:Ket(f"psi_{n}^0")]
En  = lambda n=n: Function(var(f'E_{n}'))(n) # A callable function.
En0 = lambda n=n: Eq(S(f'E_{n}'), Function(var(f'E_{n}^0'))(n)) # A callable function.

# Object symbols
C = CoordSys3D('C') # Cartesian coordinate system.

Hp = Operator('Hprime') # Perturbation onto Hamiltonian
psib, psik, nb, nk = [Bra('psi'), Ket('psi'), Bra('n'), Ket('n')]
Psi    = Function('Psi')(x,y,z,t)
psi    = Function('psi')
psix   = Function('psi')(x)
psixt  = Function('psi')(x,t)
phik   = Function('phi')(k)
PsiSph = Function('Psi')(r,theta,phi,t)
psiSph = Function('psi')(r,theta,phi)

class quantum_mechanics(branch):
    """
    solve(oqmec.expX, _expX)
    """

    _name: str = "quantum_mechanics"

    # Integer symbols
    global k,k1,k2,k3,n,i,j,m
    k, k1,k2,k3 = symbols('k k1 k2 k3', integer=True)
    n, i, j, m = symbols('n i j m', positive=True, integer=True)

    # Real symbols
    global xA,xB,yA,yB,pA,pB
    xA, xB = symbols('x_A x_B', real=True)
    yA, yB = symbols('y_A y_B', real=True)
    pA, pB = symbols('p_A p_B', real=True)


    def define_symbols(self) -> None:
        """
        # Real symbols
        global eps0
        eps0 = symbols('varepsilon_0', real=True, positive=True)
        global A,B,a,b,p,M,q,l,omega,alpha,Z,e
        A,B,a,b,p,M,q,l,omega,alpha,Z,e = symbols('A B a b p M q l omega alpha Z e', real=True, positive=True)
        global m,r,w,t,T,H,V
        m,r,w,t,T,H,V = symbols('m r w t T H V', real=True, positive=True, nonzero=True)
        # Mass of electron, proton, neutron, nucleus.
        global m_e,m_p,m_n,m_N
        m_e,m_p,m_n,m_N = symbols('m_e m_p m_n m_N', real=True, positive=True)
        global x,y,z, xmin, xmax
        x,y,z = symbols('x y z', real=True)
        xmin,xmax = symbols('x_{min} x_{max}', real=True)
        global pmin, pmax
        pmin,pmax = symbols('p_{min} p_{max}', real=True)
        global B0,V0, theta, phi
        B0,V0, theta, phi = symbols('B_0 V_0 theta phi', real=True)
        global vA,vB,vC,vD
        vA,vB,vC,vD = symbols('A B C D', vector=True)
        global Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz
        Ax,Bx,Cx,Dx = symbols('A_x B_x C_x D_x', real=True)
        Ay,By,Cy,Dy = symbols('A_y B_y C_y D_y', real=True)
        Az,Bz,Cz,Dz = symbols('A_z B_z C_z D_z', real=True)

        # Numerical values
        Common global symbols, functions, objects.
        """
        # Symbols are now defined at code level for better IDE support.
        pass

    def __init__(self, class_type: str = 'position_space') -> None:
        """
        class_type = \
        {1:"position_space",
         2:"momentum_space"}
        """
        super().__init__()
        self.class_type = class_type
        self.define_symbols()

        # File settings
        self.input_dir: str = "input/quantum_mechanics"
        self.output_dir: str = "output/quantum_mechanics"

        class subformulary:
            """
            Sub formulary class.

            Define global symbols in the outer class.
            """
            def __init__(self, class_type: str) -> None:
                if class_type in ["position_space"]:
                    self.Icm_sphere = S(2)/5*M*r**2
        self.subformulary = subformulary(self.class_type)
        
#### POSITION SPACE
        if self.class_type in ["position_space"]:
            """
            Append to class
            2nd edition
            1.40
            2.5 OK
            2.17
            
            ch2.2
            2.27 OK
            2.28 OK
            2.34 ???
            """
            self.psib, self.psik = [psib, psik] 
            self.Psi = Psi
            self.psi = psi
            self.psix = psix
            self.V = V
            self.normalization   = Eq(Integral(conjugate(self.Psi)*self.Psi, (x,xmin,xmax)), 1)
            self.normalizationXP = Eq(Integral(conjugate(self.Psi)*self.Psi, (x,xmin,xmax), (p,pmin,pmax)), 1)
            self.cn = Eq(var(r'c_n'), Integral(conjugate(self.psix)*f(x), (x,xmin,xmax)))
            self.psicn = Eq(S(r'psi(x,t)'), Sum(self.cn.rhs*f(x)*exp(-I*En()*t), (n,1,oo)))
            
#### ----> Position Related Formulary
            # Do not use Lambda of sympy it requires definifion of fx. (NameError: name 'fx' is not defined)
            """
            Examples:
            =========
            oqmec.Psi = oqmec.hydrogen.psi_sy(1,0,0).rhs
            oqmec.exp_fxSph(1/r).doit()        
            """
            # Generic expectation value calculation routine for a given operator fx.
            self.exp_fx    = lambda fx:Eq(var(r'\langle{'+str(fx)+r'}\rangle'), Integral(conjugate(self.Psi)*fx*self.Psi, (x,xmin,xmax)))
            self.exp_fx2   = lambda fx:Eq(var(r'\langle{'+str(fx)+r'}\rangle'), Integral(conjugate(self.Psi)*fx**2*self.Psi, (x,xmin,xmax)))
            self.exp_fxSph = lambda fx:Eq(var(r'\langle{'+str(fx)+r'}\rangle'), Integral(conjugate(self.Psi)*fx*self.Psi*r**2*sin(theta), (r,0,oo), (phi,0,2*pi), (theta,0,pi)))
            self.exp_fxSphR= lambda fx:Eq(var(r'\langle{'+str(fx)+r'}\rangle'), Integral(conjugate(self.Psi)*fx*self.Psi*r**2*4*pi, (r,0,oo))) # 4*pi comes from solid angle.
            self.delta_fx  = lambda fx:Eq(var(r'\Delta{'+str(fx)+'}='+str(fx)+r'-\langle{'+str(fx)+'}\rangle'), sqrt(self.exp_fx2.rhs - self.exp_fx.rhs**2))
            
            self.exp_x    = Eq(var(r'\langle{x}\rangle'),   integrate(conjugate(self.Psi)*x*self.Psi, (x,xmin,xmax)))
            self.exp_x2   = Eq(var(r'\langle{x^2}\rangle'), Integral(conjugate(self.Psi)*x**2*self.Psi, (x,xmin,xmax)))
            self.delta_x  = Eq(var(r'\Delta{x}=x-\langle{x}\rangle'), sqrt(self.exp_x2.rhs - self.exp_x.rhs**2))
            self.delta_x2 = Eq(var(r'(\Delta{x})^2=\langle{x^2}\rangle-\langle{x}\rangle^2'), self.exp_x2.rhs - self.exp_x.rhs**2)


#### Entanglement
#### ----> Bell States
            """
            The Bell's states or EPR pairs are specific quantum states of 
            two qubits that represent the simplest examples of quantum entanglement. 
            The Bell's states are a form of entangled and normalized basis vectors.
            """
            self.Bell_psi_p = Eq( var(r'Psi^+'), 1/sqrt(2)*(TensorProduct(Ket(0,A), Ket(1,B)) + TensorProduct(Ket(1,A), Ket(0,B))) )
            self.Bell_psi_m = Eq( var(r'Psi^-'), 1/sqrt(2)*(TensorProduct(Ket(0,A), Ket(1,B)) - TensorProduct(Ket(1,A), Ket(0,B))) )
            self.Bell_phi_p = Eq( var(r'Phi^+'), 1/sqrt(2)*(TensorProduct(Ket(0,A), Ket(0,B)) + TensorProduct(Ket(1,A), Ket(1,B))) )
            self.Bell_phi_m = Eq( var(r'Phi^+'), 1/sqrt(2)*(TensorProduct(Ket(0,A), Ket(0,B)) - TensorProduct(Ket(1,A), Ket(1,B))) )
            
            # oqmec.Bell_psi_p_func(x,y,x**2,y**2)
            # oqmec.Bell_psi_p_func(Ket(0,A), Ket(0,B), Ket(1,A), Ket(1,B))
            self.Bell_psi_p_func = lambda func_0A=lambda _:_, func_0B=lambda _:_, func_1A=lambda _:_, func_1B=lambda _:_ : (1 / sqrt(2)) * (func_0A * func_1B + func_1A * func_0B)
                        

#### Operator Definitions via SymPy sympy.physics.quantum class 
            # todo ERROR in DifferantialOperator application
            self.exp_xop   = Eq(var(r'\langle{x}\rangle'),   Integral(conjugate(self.Psi)*qapply(x*self.Psi), (x,xmin,xmax)))
            self.exp_x2op  = Eq(var(r'\langle{x^2}\rangle'), Integral(conjugate(self.Psi)*qapply(x**2*self.Psi), (x,xmin,xmax)))
            self.delta_xop = Eq(var(r'\Delta{x}'), sqrt(self.exp_x2op.rhs - self.exp_xop.rhs**2))
            self.delta_x2op= Eq(var(r'(\Delta{x})^2'), self.exp_x2op.rhs - self.exp_xop.rhs**2)
            
#### ----> Momentum Related Formulary.
            """
            Examples:
            =========    
            oqmec.px.subs(Psi, x**2).doit()
            oqmec.px.subs({oqmec.Psi:x**2}).doit()
            oqmec.px.replace(Psi, x**2).doit()
            oqmec.px.xreplace({oqmec.Psi:x**2}).doit()  psi=x**2 --> px=-2*hb*I*x
            
            oqmec.pxop.xreplace({oqmec.Psi:x**2}).doit()
            oqmec.pxop.subs({oqmec.Psi:x**2}).doit()
            oqmec.pxop.rhs.subs({Psi:x**2}).doit()
            qapply(oqmec.pxop).subs({Psi:x**2}).doit()
            
            oqmec.JxSph.replace(psiSph, Ynm(1,1,theta,phi)).doit().expand(func=True)
            """
            self.px      = Eq(var(r'\hat{p}_x'),  -I*hbar*D(self.Psi, x)) 
            self.py      = Eq(var(r'\hat{p}_y'),  -I*hbar*D(self.Psi, y))
            self.pz      = Eq(var(r'\hat{p}_z'),  -I*hbar*Derivative(self.Psi, z))
            self.px2     = Eq(var(r'\hat{p}_x^2'), (-I*hbar)**2*D(self.Psi, x, 2))
            self.exp_px  = Eq(var(r'\langle{p_x}\rangle'),   Integral(conjugate(self.Psi)*self.px.rhs,  (x,xmin,xmax)))
            self.exp_px2 = Eq(var(r'\langle{p_x^2}\rangle'), Integral(conjugate(self.Psi)*self.px2.rhs, (x,xmin,xmax)))
            self.delta_px  = Eq(var(r'\Delta{p_x}'), sqrt(self.exp_px2.rhs - self.exp_px.rhs**2))
            self.delta_px2 = Eq(var(r'(\Delta{p_x})^2'), self.exp_px2.rhs - self.exp_px.rhs**2)
            self.delta_xp  = self.uncertainity_xp = Eq(var(r'\Delta{x}\Delta{p_x}'), self.delta_x.rhs*self.delta_px.rhs)

#### Operator Definitions via SymPy sympy.physics.quantum class 
            # DifferantialOperator applications.
            self.pxop     =  Eq(var(r'\hat{p}_x'), DifferentialOperator(-I*hbar*D(self.Psi, x), self.Psi).expr)
            self.px2op    =  Eq(var(r'\hat{p}^2_x'), DifferentialOperator(-I*hbar*D(self.pxop.rhs, x), self.Psi).expr)
            self.exp_pxop  = Eq(var(r'\langle{\hat{p}_x}\rangle'),   Integral(conjugate(self.Psi)*self.pxop.rhs,  (x,xmin,xmax)))
            self.exp_px2op = Eq(var(r'\langle{\hat{p}_x^2}\rangle'), Integral(conjugate(self.Psi)*self.px2op.rhs, (x,xmin,xmax)))
            self.delta_pxop  = Eq(var(r'\Delta{p_x}'), sqrt(self.exp_px2op.rhs - self.exp_pxop.rhs**2))
            self.delta_px2op = Eq(var(r'(\Delta{p_x})^2'), self.exp_px2op.rhs - self.exp_pxop.rhs**2)
            self.delta_xop_pxop = self.uncertainity_xop_pxop  = Eq(var(r'\Delta{x}\Delta{p_x}'), self.delta_xop.rhs*self.delta_px2op.rhs)
            
#### ----> Orbital Angular Momentum in Cartesian Coordinates
            """
            
            """
#            self.Lx      = Eq(var('L_x'), var('L_x'))
#            self.Ly      = Eq(var('L_y'), var('L_y'))
            self.Lx,self.Ly,self.Lz = [Operator(r'\hat{L}_x'), Operator(r'\hat{L}_y'), Operator(r'\hat{L}_z')]
            self.L2      = Eq(Operator(r'\hat{L}^2'), self.Lx**2 + self.Ly**2 + self.Lz**2)
            self.comLxLy = Eq(Commutator(self.Lx, self.Ly), I*hbar*self.Lz)
            self.comLzLx = Eq(Commutator(self.Lz, self.Lx), I*hbar*self.Ly)
            self.comLyLz = Eq(Commutator(self.Ly, self.Lz), I*hbar*self.Lx)
            self.Lplus   = Eq(Operator(r'\hat{L}_+'), self.Lx + I*self.Ly)
            self.Lminus  = Eq(Operator(r'\hat{L}_-'), self.Ly - I*self.Ly)
            
            self.comJxJy = Eq(Commutator(Jx, Jy), I*hbar*Jz)
            self.comJzJx = Eq(Commutator(Jz, Jx), I*hbar*Jy)
            self.comJyJz = Eq(Commutator(Jy, Jz), I*hbar*Jx)
            self.comJ2Ji = lambda i=x:Eq(Commutator(J2, {x:Jx, y:Jy, z:Jz}[i]), 0)
            self.comJzJp = Eq(Commutator(Jz, Jplus),   hbar*Jplus)
            self.comJzJm = Eq(Commutator(Jz, Jminus), -hbar*Jminus)
            self.comJpJm = Eq(Commutator(Jplus, Jminus), 2*hbar*Jz)
            self.J2ket   = Eq(J2*JzKet(j,m), simplify(qapply(J2*JzKet(j,m))))
            self.Jzket   = Eq(Jz*JzKet(j,m), simplify(qapply(Jz*JzKet(j,m))))
            self.Jpket   = Eq(Jplus*JzKet(j,m), simplify(qapply(Jplus*JzKet(j,m))))
            self.Jmket   = Eq(Jminus*JzKet(j,m), simplify(qapply(Jminus*JzKet(j,m))))

#### ----> Total Angular Momentum in Spherical Coordinates.
            self.JxSph    = Eq(Operator(r'\hat{J}_x'),  hbar/I*(-sin(phi)*D(psiSph, theta) - cos(phi)*cot(theta)*D(psiSph, phi)) )
            self.JySph    = Eq(Operator(r'\hat{J}_y'),  hbar/I*( cos(phi)*D(psiSph, theta) - sin(phi)*cot(theta)*D(psiSph, phi)) )
            self.JzSph    = Eq(Operator(r'\hat{J}_z'),  hbar/I*(Derivative(psiSph, phi)) )
            self.JpSph    = Eq(Operator(r'\hat{J}_+'),  hbar*exp( I*phi)*( D(psiSph, theta) + I*cot(theta)*D(psiSph, phi)) )
            self.JmSph    = Eq(Operator(r'\hat{J}_-'),  hbar*exp(-I*phi)*(-D(psiSph, theta) + I*cot(theta)*D(psiSph, phi)) )
            self.J2Sph    = Eq(Operator(r'\hat{J}^2'), -hbar**2*(1/sin(theta)*D(sin(theta)*D(psiSph, theta), theta) + 1/(sin(theta))**2*D(psiSph, phi, 2)) ) # Schwabl2007 (5.18e)
            self.L2Sph    = self.J2Sph

#### ----> Orbital Angular Momentum in Spherical Coordinates
            # Define the Laplacian in spherical coordinates todo copy to omethod
            self.laplacian_spherical = Eq(Operator(r'\nabla^2'), \
                                       (1/r**2)*D(r**2*D(psiSph, r), r) + \
                                       (1/(r**2*sin(theta)))*D(sin(theta)*D(psiSph, theta), theta) + \
                                       (1/(r**2*sin(theta)**2))*D(psiSph, phi, phi))
            self.nabla2Sph = self.laplacian_spherical

            self.prSph    = Eq(Operator(r'\hat{p}_r'), hbar/I*(1/r)*D(r*psiSph,r)) # Schwabl2007 (6.6)
            self.pr2Sph   = Eq(Operator(r'\hat{p}_r^2'), qapply(self.prSph.rhs * qapply(self.prSph.rhs*psiSph)) )
            self.p2Sph    = Eq(Operator(r'\hat{p^2}'), self.pr2Sph.rhs + (1/r**2)*self.L2Sph.rhs ) # Schwabl2007 (6.4')
            self.p2Sph2   = Eq(Operator(r'\hat{p^2}'), -hbar**2*self.nabla2Sph.rhs)       
            
#### ----> Spin Angular Momentum
            """
            oqmec.Sr.rhs.doit().eigenvals()
            oqmec.Sr.rhs.doit().eigenvects()
            """
            # self.sigmax   = Eq(var(r'\sigma_x'), UnevaluatedExpr(Matrix([[0,1], [1,0]])))
            self.sigmax   = Eq(var(r'\sigma_x'), UnevaluatedExpr(Matrix([[0,1], [1,0]])))
            self.sigmay   = Eq(var(r'\sigma_y'), UnevaluatedExpr(Matrix([[0,-I], [I,0]])))
            self.sigmaz   = Eq(var(r'\sigma_z'), UnevaluatedExpr(Matrix([[1,0], [0,-1]])))
#            self.Sx       = Eq(Operator(r'\hat{S}_x'), Operator(UnevaluatedExpr(Matrix([[0,1], [1,0]]))))
            self.Sx       = Eq(Operator(r'\hat{S}_x'), hbar/2*self.sigmax.rhs)
            self.Sy       = Eq(Operator(r'\hat{S}_y'), hbar/2*self.sigmay.rhs)
            self.Sz       = Eq(Operator(r'\hat{S}_z'), hbar/2*self.sigmaz.rhs)
            self.Sr       = Eq(Operator(r'\hat{S}_r'), 
                               UnevaluatedExpr(hbar/2*Matrix(
                                   [[cos(theta),            exp(-I*phi)*sin(theta)],
                                    [exp(I*phi)*sin(theta), -cos(theta)]])))
            
            
            # Spin 1/2 particle eigenvectors.
            self.sx_up    = Eq(var(r'|+x\rangle'), 1/sqrt(2)*UnevaluatedExpr(Matrix([[1],[1]])))
            self.sx_down  = Eq(var(r'|-x\rangle'), 1/sqrt(2)*UnevaluatedExpr(Matrix([[1],[-1]])))
            self.sy_up    = Eq(var(r'|+y\rangle'), 1/sqrt(2)*UnevaluatedExpr(Matrix([[1],[I]])))
            self.sy_down  = Eq(var(r'|-y\rangle'), 1/sqrt(2)*UnevaluatedExpr(Matrix([[1],[-I]])))
            self.sz_up    = Eq(var(r'|+z\rangle'), UnevaluatedExpr(Matrix([[1],[0]])))
            self.sz_down  = Eq(var(r'|-z\rangle'), UnevaluatedExpr(Matrix([[0],[1]])))

            
#### ----> Schrödinger Equation
            # self.V = self.potential_energy = Function('V')(x) todo erase after tests
            self.H = Eq(S('H'), -(hbar**2)/(2*m)*diff(self.Psi, x, 2) + self.V*self.Psi)
            self.exp_V = Eq(var(r'\langle{V}\rangle'), Integral(conjugate(self.Psi)*qapply(self.V*self.Psi), (x,xmin,xmax)))
            self.exp_T = Eq(var(r'\langle{T}\rangle'), self.exp_px2.rhs/(2*m))
            self.exp_H = Eq(var(r'\langle{H}\rangle'), Integral(conjugate(self.Psi)*self.H.rhs, (x,xmin,xmax)))
            self.HOp = Eq(S('H'), hbar**2/(2*m)*(DifferentialOperator(diff(self.Psi, x), self.Psi))**2 + Operator(self.V*self.Psi))
            self.SchrodingerEq = Eq(-(hbar**2)/(2*m)*diff(self.Psi, x, 2) + self.V*self.Psi, I*hbar*diff(self.Psi, t, 1))
            self.SchrodingerEq_TI = self.SchrodingerEq_Time_Independent = Eq(-(hbar**2)/(2*m)*diff(self.Psi,x,2) + self.V*self.Psi , En()*self.Psi)

            
####----> Fourier Transforms
            # todo look Schwabl not Griffiths
            """
            #===========================================================================
            # # todo use in Fourier Transform
            # x = XKet()
            # pprints("|x>=",x,
            #         "|x> in x-space=",rep_innerproduct(XKet(), basis = XOp()),
            #         "|x> in p-space=",rep_innerproduct(XKet(), basis = PxOp()),
            #         )
            #===========================================================================
            """
            self.FourierTransform_phi_t0 = Eq( Psi.subs(t,0), 1/(sqrt(2*pi))*Integral(phi*exp(I*k*x) , (k,-oo,oo)) )
            self.FourierTransform_Psi_t0 = Eq( phi, 1/(sqrt(2*pi))*Integral(Psi.subs(t,0)*exp(-I*k*x) , (k,-oo,oo)) )
            # fourier_transform(1/sqrt(2*pi)*npsi.expr, x, k).subs({k:k/(2*pi)})
            self.fourier_transform_psix         = Eq( phik, fourier_transform(1/sqrt(2*pi)*psix, x, k).subs({k:k/(2*pi)}))
            self.inverse_fourier_transform_phik = Eq( psix, inverse_fourier_transform(1/sqrt(2*pi)*phik*exp(-I*hbar*(k**2)*t/(2*m)), k, x).subs({x:x/(2*pi)}))


####----> Time Evoluation
            self.time_evolution_psixt_from_phik = Eq(psixt, 1/(sqrt(2*pi))*Integral(phik*exp(I*(k*x-hbar*k**2/(2*m))) , (k,-oo,oo)))
            self.time_evolution_psixt_from_psix = Eq(psixt, 1/(sqrt(2*pi))*Integral( fourier_transform(1/sqrt(2*pi)*psix, x, k).subs({k:k/(2*pi)})*exp(I*(k*x-hbar*k**2/(2*m))) , (k,-oo,oo) ))

            
#### Perturbation Theory
#### ----> Time-independent Perturbation Theory
            """
            ND: Nondegenerate perturbation
            PT: Perturbation Theory
            
            qapply(Bra(i)*A*Ket(j))
            """
            self.Hp = Hp
            # 1st order nondegenerate perturbation for energy. E_n^1.
            self.En1_ND_PT_   = Eq( var(r"E_n^1=\langle{\psi_n^0}|H^\prime|\psi_n^0\rangle"), Integral(conjugate(self.Psi)*qapply(self.Hp*self.Psi), (x,xmin,xmax)) ) # todo add summation 2nd order pert.
            self.En1_ND_PT_bk = Eq( var(r"E_n^1=\langle{n}|H^\prime|n\rangle"), qapply(nb*self.Hp*nk) ) # <n|H'|n>
            self.En_ND_PT1 = En0().rhs + self.En1_ND_PT_.rhs
            self.nondegenerate_perturbation_theory_En_1ord = self.En1_ND_PT_
            self.nondegenerate_perturbation_theory_En_1ord_bk = self.En1_ND_PT_bk
            self.psi1n = lambda imin=-2, imax=2, total=0: [total := total + qapply(psi0b(k)*Hp*psi0k(n)).subs({k:i})*psi0k(i)/(En0(n).rhs-En0(i).rhs) for i in list(Range(n+imin, n+imax+1)) if i != n][-1]


#### PHASE SPACE
#### ----> Wigner Quasi-Probability Distribution Function
            """
            oqmec.Wigner1D(x**2)
            oqmec.Wigner2D(xA**2+xB**3)
            """
            self.Wigner1D = lambda psi = psi(x): \
                Eq( S(r'W(x,p)'),
                   1/(2*pi*hbar)*Integral( conjugate(psi.subs({x:x+y/2}))*psi.subs({x:x-y/2})*exp(I*p*y/hbar), (y,-oo,oo) ) ) # Agarwal2004, Eq.3
            self.Wigner2D = lambda psi = psi(xA,xB): \
                Eq( S(r'W(x_A,p_A,x_B,p_B)'),
                   1/(2*pi*hbar)**2*Integral( conjugate(psi.subs({xA:xA+yA/2, xB:xB+yB/2}))*psi.subs({xA:xA-yA/2, xB:xB-yB/2})*exp(I*(pA*yA+pB*yB)/hbar), (yA,-oo,oo), (yB,-oo,oo) ) ) # Bhatt2008, Eq.7

            
#### --- CLASSES ---

#### Hamiltonians
            class Hamiltonians(branch):
                """
                Sub class for Hamiltonians that are used in quantum mechanics.
                
                Examples:
                ---------    
                oqmec.Hamiltonian.B = B0*C.k
                oqmec.Hamiltonian.Sp = Sx*C.i + Sy*C.j + Sz*C.k
                
                substitutions = {oqmec.Hamiltonian.B:B0*C.k, oqmec.Hamiltonian.Sp:Sx*C.i + Sy*C.j + Sz*C.k}
                oqmec.Hamiltonian.e_in_B.xreplace(substitutions).doit()
                """
                global Sx,Sy,Sz,Sp
                Sx,Sy,Sz,Sp = symbols('S_x S_y S_z Sp')
                
                def __init__(self):
                    super().__init__()
                    self.name = "Hamiltonians"
                    self.gamma = var('gamma', real=True)
                    self.B  = Bx*C.i + By*C.j + Bz*C.k
                    self.Sp = Sx*C.i + Sy*C.j + Sz*C.k
                    self.e_in_B = self.e_in_magnetic_field = lambda B=self.B, Sp=self.Sp: Eq(S('H'), -self.gamma*B.dot(Sp))
            self.Hamiltonians = Hamiltonians()

            
#### Quantum Harmonic Oscillator
            class qho(branch):
                """
                Sub class for 1D quantum simple harmonic oscillator.
                """
                def __init__(self):
                    super().__init__()
                    self.name = "Quantum Harmonic Oscillator"
                    self.xi = Eq(S('xi'), sqrt(m*w/hbar)*x)
                    self.psi = lambda n=n: Eq(S(f'psi_{n}'), (m*w/(pi*hbar))**(S(1)/4)*(1/sqrt((2**n)*factorial(n)))*hermite(n, xi)*exp(-xi**2/2))
                    self.psix= lambda n=n: Eq(S(f'psi_{n}'), self.psi(n).rhs.subs({xi:self.xi.rhs}))
                    self.ad = RaisingOp('a') # Raising ladder operator.  Creation operator.
                    self.a  = LoweringOp('a')# Lowering ladder operator. Annihilation operator.
                    self.nk = lambda n=n:SHOKet(n) # qapply(oqmec.qho.nb(n)*oqmec.qho.nk(n)).doit() -> 1
                    self.nb = lambda n=n:SHOBra(n) # simplify(qapply(oqmec.qho.nb(j)*oqmec.qho.x2op.rhs*oqmec.qho.nk(k)))
                    self.px = Eq(var(r'\hat{p}_x'), DifferentialOperator(-I*hbar*D(f(x), x, 1), f(x)) )
                    self.xop = Eq(S('xhat'), sqrt(hbar/(2*m*w))*(self.ad+self.a))
                    self.pop = Eq(S('phat'), I*sqrt(hbar*m*w/2)*(self.ad-self.a))
                    self.x2op = Eq(S('xhat^2'), simplify(self.xop.rhs*self.xop.rhs))
                    self.p2op = Eq(S('phat^2'), simplify(self.pop.rhs*self.pop.rhs))
                    self.V  = Eq(V, S(1)/2*m*w**2*self.x2op.rhs)
                    self.Vx = Eq(V, S(1)/2*m*w**2*x**2)
                    self.H = Eq(H, simplify(self.p2op.rhs/(2*m) + self.V.rhs))
                    self.Hx = Eq(H, DifferentialOperator(-hbar**2/(2*m)*Derivative(f(x),x,2) + S(1)/2*m*w**2*x**2, f(x)))
                    self.En = lambda n=n:Eq(S(f'E_{n}'), (n+S(1)/2)*hbar*w)
                    # todo 2.51
                
                """
                def psi(self, n=n):
                    # Same as self.psix
                    xi = sqrt(m*w/hbar)*x
                    return Eq(S(f'psi_{n}'), Wavefunction((m*w/(pi*hbar))**(1/4)*(1/sqrt((2**n)*factorial(n)))*hermite(n, xi)*exp(-xi**2/2), (x,-oo, oo)).expr)
                """
            self.qho = qho()
            
            
#### q-Deformed Quantum Harmonic Oscillator
            class qdefho(branch):
                """
                Sub class for 1D q-deformed quantum harmonic oscillator.

                oqmec.qdefho.__init__(numeric=False)
                plot_sympfunc([abs(1/sqrt(2)*oqmec.qdefho.psix(x,2,0.001,1/2).doit().n().rhs + 1/sqrt(2)*oqmec.qdefho.psix(x,6,0.001,1/2).doit().n().rhs),], (-4,4,2000) )
                plot_sympfunc([re(oqmec.qdefho.psix(x,6,0.001,1/2).doit().n().rhs),], (-4,4,400) )
                
                oqmec.qdefho.__init__(numeric=True)
                mp.plot(lambda x: mp.re(oqmec.qdefho.psix(x,6,0.001,1/2)), [-4, 4], points=400)
                
                References:
                    Jafarov2010 todoefe
                """
                def __init__(self, numeric=False):
                    super().__init__()
                    self.numeric = numeric
                    self.name = "q-Deformed Quantum Harmonic Oscillator"
                    self.ad = RaisingOp('a') # Raising ladder operator.  Creation operator.
                    self.a  = LoweringOp('a')# Lowering ladder operator. Annihilation operator.
                    
                    if not self.numeric:
                        self.qPochhammer = lambda a=a, q=q, n=n: Eq( var(f'({a};{q})_{n}'), Product(1-a*q**k, (k, 0, n-1)) )
                        self.cn = lambda n=n, q=q, l=l: Eq( var(f'c_{n}'), (2*l/pi)**(S(1)/4)*(I)**n*q**(n/2)*self.qPochhammer(q,q,n).rhs**(-S(1)/2) )
                        self.psix = lambda x=x, n=n, q=q, l=l: (
                            Eq( var(f'psi_{n}'), self.cn(n,q,l).rhs*exp(-l*x**2)*
                               Sum( ( self.qPochhammer(q**(-n),q,k).rhs/self.qPochhammer(q,q,k).rhs )*q**(n*k-(k**2)/2)*exp(-2*I*l*(sqrt(-log(q)/l))*x*k), (k,0,n)) ) )
                        self.En = lambda n=n, q=q: Eq( S(f'En(n, q)'), hbar*omega/2*sinh(log(q)*(n + S(1)/2)) / sinh(log(q)/2) )
                                                
                    else:
                        self.cn = lambda n=1, q=0.001, l=1/2: ( (2*l/mp.pi)**(1/4)*(mp.j)**n*q**(n/2)*mp.qp(q,q,n)**(-1/2) )
                        self.psix = lambda x, n=1, q=0.001, l=1/2: (self.cn(n,q,l) * mp.exp(-l*x**2) * 
                            mp.nsum( lambda k: ( mp.qp(q**(-n),q,k) / mp.qp(q,q,k) )*mp.power(q, n*k-(k**2)/2)*mp.expj(-2*l*(mp.sqrt(-mp.ln(q)/l))*x*k), [0,n]) )
                        
                    """
                    self.nk = lambda n=n:SHOKet(n) # qapply(oqmec.qho.nb(n)*oqmec.qho.nk(n)).doit() -> 1
                    self.nb = lambda n=n:SHOBra(n) # simplify(qapply(oqmec.qho.nb(j)*oqmec.qho.x2op.rhs*oqmec.qho.nk(k)))
                    self.px = Eq(var(r'\hat{p}_x'), DifferentialOperator(-I*hbar*D(f(x), x, 1), f(x)) )
                    self.xop = Eq(S('xhat'), sqrt(hbar/(2*m*w))*(self.ad+self.a))
                    self.pop = Eq(S('phat'), I*sqrt(hbar*m*w/2)*(self.ad-self.a))
                    self.x2op = Eq(S('xhat^2'), simplify(self.xop.rhs*self.xop.rhs))
                    self.p2op = Eq(S('phat^2'), simplify(self.pop.rhs*self.pop.rhs))
                    self.V  = Eq(V, S(1)/2*m*w**2*self.x2op.rhs)
                    self.Vx = Eq(V, S(1)/2*m*w**2*x**2)
                    self.H = Eq(H, simplify(self.p2op.rhs/(2*m) + self.V.rhs))
                    self.Hx = Eq(H, DifferentialOperator(-hbar**2/(2*m)*Derivative(f(x),x,2) + S(1)/2*m*w**2*x**2, f(x)))
                    self.En = lambda n=n:Eq(S(f'E_{n}'), (n+S(1)/2)*hbar*w)
                    """
            self.qdefho = qdefho()            

            
#### Delta Function Quantum Well
            class dqw(branch):
                """
                Sub class for Delta Function Quantum Well. todo check
                
                2.114 use sympy delta function
                2.129
                2.141
                """
                def __init__(self):
                    super().__init__()
                    self.name = "Delta Function Quantum Well"
#                    2.114 use sympy delta function
                    self.V = Eq(var(r'V(x)'), -alpha*DiracDelta(x))
#                    2.129
                    self.psi= Eq(var(r'psi(x)'), sqrt(m*alpha)/hbar * exp(-m*alpha*abs(x)/hbar**2) )
                    self.En = Eq(var('E'), -(m*alpha**2)/(2*hbar**2))
#                    2.141
                    self.R = self.reflection_coefﬁcient   = Eq(var('R'), 1/(1 + ((2*hbar**2*self.En.rhs) / (m*alpha**2))))
                    self.T = self.transmission_coefﬁcient = Eq(var('T'), 1/(1 + ((m*alpha**2)/(2*hbar**2*self.En.rhs))))
            self.dqw = dqw()


#### Infinite Square Quantum Well
            class iqw(branch):
                """
                Sub class for Infinite Square Quantum Well
                """
                def __init__(self):
                    super().__init__()
                    self.name = "Infinite Square Quantum Well"
                    self.psix =lambda n=n: Eq(var(rf'psi_{n}(x)'), sqrt(2/a)*sin(n*pi*x/a))
#                    self.V = Eq(V, S(1)/2*m*w**2*self.x2op.rhs)
#                    self.H = Eq(H, simplify(self.p2op.rhs/(2*m) + self.V.rhs))
                    self.En = lambda n=n:Eq(var(f'E_{n}'), (n**2*pi**2*hbar**2)/(2*m*a**2))
            self.iqw = iqw()
            

#### Finite Square Quantum Well            
            class fqw(branch):
                """
                Sub class for Finite Square Quantum Well
                
                class finite_square_well(branch):
                        2.145  use piecewise
                        2.157
                        2.171
                        2.169
                """
                def __init__(self):
                    super().__init__()
                    self.name = "Finite Square Quantum Well"
                    self.psix =lambda n=n: Eq(S(f'psi_{n}'), sqrt(2/a)*sin(n*pi*x/a))
#                    self.V = Piecewise((r'-V_0' , -a <= x <= a), (0, abs(x) > a)) todo
#                    self.H = Eq(H, simplify(self.p2op.rhs/(2*m) + self.V.rhs)) todo
                    self.En = lambda n=n:Eq(S(f'E_{n}'), (n**2*pi**2*hbar**2)/(2*m*a**2))
#                    2.157
#                    E_bwd = E_n + V_0 for wide deep well todo
#                    self.E_bwd = Eq(var(r'E_n + V_0'), (n**2 * pi**2 * hbar**2)/(2*m*(2*a)**2)) todo
#                    2.171
#                    E_bsn = E_n + V_0 for shallow narrow well
#                    self.E_bsn = Eq(var(r'E_n + V_0'), (n**2 * pi**2 * hbar**2)/(2*m*(2*a)**2))
#                    2.169
                    self.T = self.transmission_coefﬁcient = Eq(var(r'T^{-1}'), 1 + ((V0**2)/(4*En()*(En()+V0)) * (sin((2*a)/(hbar) * sqrt(2*m*(En()+V0))))**2)) # todo look desai
            self.fqw = fqw()

            
#### Hydrogen Atom
            class hydrogen(branch):
                """
                Sub class for Hydrogen-like atom.
                References:
                - Wikipedia Hydrogen-like atom.
                """
                global e
                e = symbols('e', positive=True, real=True)
                
                def __init__(self):
                    super().__init__()
                    self.name = "Hydrogen-like atom"
                    self.a0 = self.bohr_radius = S('a_0')
                    self.mu = m_N*m_e/(m_N+m_e)
                    self.a_mu = Eq(S('a_mu'), Eq(4*pi*eps0**2*hbar**2/(self.mu*e**2), m_e/self.mu*self.a0))
                    self.Rnl = lambda n=n, l=l, r=r, Z=1:Eq( S(f'R_{n}{l}(r)'), sqrt((2*Z/(n*self.a_mu.lhs))**3*factorial(n-l-1)/(2*n*factorial(n+l)))*exp(-Z*r/(n*self.a_mu.lhs))*(2*Z*r/(n*self.a_mu.lhs))**l*assoc_laguerre(n-l-1, 2*l+1, 2*Z*r/(n*self.a_mu.lhs)) )
                    self.psi = lambda n=n, l=l, m=m, r=r, Z=1, phi=phi, theta=theta:Eq(S(f'psi_{n}{l}{m}({r},{phi},{theta})'), self.Rnl(n,l,r,Z).rhs*Ynm(l,m,theta,phi))
                    self.V  = Eq(S('V(r)'), -1/(4*pi*eps0)*Z*e**2/r)
                    self.En = lambda n=n, Z=Z:Eq(S(f'E_{n}'), (Z**2*self.mu*e**4/(32*pi**2*eps0**2*hbar**2))*(1/n**2))
                    """
                    todo write from Griffiths
                    self.exp_invr 
                    self.exp_invr2
                    self.KramersRelation
                    """
                    
                    # todo make a table Rnl functions.
                    """
                    def table_func(func, kwargs*=[]):
                        # table_func(oqmec.hydrogen.Rnl, kwargs*=[]):
                        res = []
                        return res
                    """
                    
                    # SymPy based properties.
                    self.psi_sy = lambda n=n, l=l, m=m, r=r, phi=phi, theta=theta:Eq(S(f'psi_{n}{l}{m}({r},{phi},{theta})'), Psi_nlm(n,l,m,r,phi,theta))
#                    self.V = Eq(V, S(1)/2*m*w**2*self.x2op.rhs) todo
#                    self.H = Eq(H, simplify(self.p2op.rhs/(2*m) + self.V.rhs)) todo
                    self.R_nl_sy = lambda n=n, l=l, r=r, Z=1:Eq(S(f'R_{n}{l}(r)'), R_nl(n, l, r, Z))
                    self.E_nZ_sy = lambda n=n, Z=Z:Eq(S(f'En({n},{Z})'), E_nl(n,Z))
                    self.E_nl_dirac_sy = lambda n=n, l=l, s=True, Z=1:Eq(S(f'En_dirac({n},{l},{s},{Z})'), E_nl_dirac(n, l, spin_up=s, Z=Z, c=137.035999037000))
            self.hydrogen = hydrogen()

#### MOMENTUM SPACE
        if self.class_type in ["momentum_space"]:
            pass
        
        # Common text definitions.
        

#### GLOBAL METHODS

#### Wavefunction Methods

#### ---->  Ket state to wave function convertion.
    def ket_to_wavefunction(
        self,
        n_: int,
        j_: int,
        psi0: Callable[[Any], Any],
        wfpsi0: Callable[..., Any],
        psipert: Any
    ) -> Any:
        """
        Replaces |kets> with wavefunctions.
        
        USAGE
        =====
        wfpsi1 = lambda n: oqmec.ket_to_wavefunction(n, j, psi0, wfpsi0, psi1n.rhs)
        wfpsi1 = lambda n,j: oqmec.ket_to_wavefunction(n, j, psi0, wfpsi0, psi1n.rhs.xreplace(substitutions))
        
        n_:          Integer state quantum number.
        j_:          Integer cutoff number. For example in x^4 potential j=4.   
        psi0(n):     Unperturbed ket state with quantum number n, |n>.
        wfpsi0(n,x): Unperturbed wavefunction with quantum number n, psi(n,x),
                     exp(-x**2/2)*hermite(n, x)/(2**(n/2)*pi**(1/4)*sqrt(factorial(n)))
        psipert:     Perturbed ket state psi1n, psi2n, etc.        
        Returns:
            Perturbed wave function for given perturbed ket state.
        """
        substitutions = [(psi0(i), wfpsi0(i)) for i in Range(n_+j_+1)]
        res = simplify(psipert.xreplace({n:n_}).subs(substitutions))
        
        if self.verbose:
            libsympy.pprints(rf"|ket> <-> Wave Function Substitutions", substitutions, 
                             self.output_style)
        return res
    

#### Phase Space Methods
    """
    def Wigner2D_num(psi, xA_val, pA, xB_val, pB, hbar=1):
        prefactor = 1 / (2 * mp.pi * hbar)**2
        # Explicitly capture xA_val and xB_val to avoid symbol/function conflicts
        integrand = lambda yA, yB: (
            mp.conj(psi(xA_val + yA/2, xB_val + yB/2)) *
            psi(xA_val - yA/2, xB_val - yB/2) *
            mp.exp(1j * (pA * yA + pB * yB) / hbar)
        )
        integral = mp.quad(integrand, [-mp.inf, mp.inf], [-mp.inf, mp.inf]) # kaldik double quadrature gerekiyor.
        return prefactor * integral
    """

    @staticmethod
    def Wigner2D_num(
        psi_num: Callable[[Any, Any], Any], 
        xA_val: float, 
        pA: float, 
        xB_val: float, 
        pB: float, 
        hbar: float = 1
    ) -> float:
        """
        Numerical 2D Wigner function:
        W(xA,pA,xB,pB)
        
        self.Wigner2D_num = lambda psi, xA, pA, xB, pB, hbar=1: (
            (1 / (2 * mp.pi * hbar)**2) * 
            mp.quad(lambda yA, yB: (
                mp.conj(psi(xA + yA / 2, xB + yB / 2)) * 
                psi(xA - yA / 2, xB - yB / 2) * 
                mp.exp(mp.j * (pA * yA + pB * yB) / hbar) ), # todo mp.expj maybe.
                [-mp.inf, mp.inf], 
                [-mp.inf, mp.inf])) # Bhatt2008, Eq.7
        """
    
        prefactor = 1 / (2 * mp.pi * hbar)**2
    
        def integrand(yA, yB):
            return (
                mp.conj(psi_num(xA_val + yA/2, xB_val + yB/2)) *
                psi_num(xA_val - yA/2, xB_val - yB/2) *
                mp.exp(1j * (pA*yA + pB*yB) / hbar)
            )
    
        integral = mp.quad(
            lambda yA: mp.quad(
                lambda yB: integrand(yA, yB),
                [-mp.inf, mp.inf]
            ),
            [-mp.inf, mp.inf]
        )
    
        return prefactor * integral


#### Nondegenerate Perturbation Theory
#### ---->  n^th order nondegenerate perturbation for energy.
    def En_ND_PT(
        self,
        order: int,
        psi0c: Callable[..., Any],
        psi0: Callable[..., Any],
        Hp: Any,
        En0: Callable[..., Any],
        k2min: Any = k2,
        k2max: Any = k2
    ) -> Any:
        """
        psi0  = |n0>, Unperturbed states, wavefunction or eigenket.
        psi0c = <m0|, Unperturbed states, dual (complex conjugate) of the 
                      wavefunction or eigenbra.
        Hp = perturbation term added to the Hamiltonian.
             Hp is assumed as Hermitian. Dagger(V(i,j)) = V(j,i)
        En0 = Unperturbed eigen-energy value function.
        imin = 
        imax =
        
        total_even = 0
        nums = [1,2,3,4,5]
        s = [total_even := total_even+x for x in nums if x % 2 == 0] # > [2,6]
        
        f  = oqmec.qho.nb(n)*oqmec.qho.x2op.rhs*oqmec.qho.nk(k)
        s1 = Dagger(f)*f
        s2 = qapply(s1)
        s3 = Sum(s2.subs({k:i}), (i,n-2,n+2))
        
        total = 0
        indices = Range(n-2,n+3)
        s = [total := total+s2.subs({k:i}) for i in indices if i != n]
        s = [total := total+norm2_matrix_s2.subs({k:i}) for i in indices if i != n]
        
        
        for i in Range(imin, imax):
            if i != n:
                # if psi0 is not a wavefunction
                S = S + norm2_matrix_s2.subs({k:i}) / (En0(n).rhs-En0(i).rhs)
        
        OLD CODE
        ========
        def En_ND_PT(self, order, psi0c, psi0, Hp, En0, **kwargs):
            kstates = Range(kwargs["kmin"], kwargs["kmax"]+1)
        
            total = 0
            matrix = psi0c(n)*Hp*psi0(k)            # Braket case. <n|H'|m> 
            norm2_matrix_s1 = Dagger(matrix)*matrix # |<n|H'|m>|^2 = <m|H'|n><n|H'|m>
            norm2_matrix_s2 = simplify(qapply(norm2_matrix_s1))
            indices = Range(imin, imax+1)
            # sum(|<n|V|m>|^2 / (En0-Em0), m!=n)
            # Iterate over i; symbolicly i=k.
            sum_list = [total := total+norm2_matrix_s2.subs({k:i})/(En0(n).rhs-En0(i).rhs) for i in indices if i != n]
            res = Eq( var(r"E_n^2"), total )
        
        Usage
        =====
            res =  oqmec.En_ND_PT(2, psi0b, psi0k, Hp, En0, k2min=k, k2max=k)
            
            kmin,kmax = [n-4,n+4]
            lmin,lmax = [kmin,kmax]
            psi0c = oqmec.qho.nb
            psi0  = oqmec.qho.nk
            En0   = oqmec.qho.En
            Hp = S(1)/2*m*alpha**2*oqmec.qho.x2op.rhs
            Hp = S(1)/4*hbar*w*(oqmec.qho.a + oqmec.qho.ad)**4
            res = oqmec.En_ND_PT(3, oqmec.qho.nb, oqmec.qho.nk, oqmec.Hp, oqmec.qho.En, k2min=n-2, k2max=n+2)
        """
        V_ij  = lambda i=i,j=j: qapply(psi0c(i)*qapply(Hp*psi0(j))) # <i0|V|j0> Application of qapply step by step gives result faster than multiplication.
        dE_ni = lambda n=n,i=i: (En0(n).rhs-En0(i).rhs)  # E_n^(0)-E_i^(0)
        
        # 1st order perturbation to energy.
        if order == 1:
            if type(psi0(n)) in (Ket, SHOKet):
                total = V_ij(n,n)
                total_qa = qapply(total)
            else:
                integral = Integral(psi0c*qapply(Hp*psi0), (x,xmin,xmax))
                res = integral.doit()
        
        # 2nd order perturbation to energy.
        if order == 2:
            sum_Vkn = 0
            k2states = Range(k2min, k2max+1)
        
            # Single summation          
            for ik2 in k2states:
                if ik2 != n:
    #                sum_Vkn += V_ij(ik2,n)*Dagger(V_ij(ik2,n)) / (dE_ni(n,ik))
                    sum_Vkn += V_ij(ik2,n)*V_ij(n,ik2) * (1/dE_ni(n,ik2))
            
            total = sum_Vkn
            total_qa = qapply(total)
          
        # 3rd order perturbation to energy.    
        if order == 3:
            sum_Vnl_Vlk_Vkn = 0
            sum_Vnl = 0
            k2states = Range(k2min, k2max+1)
            lstates = k2states
            
            # First term: Double summation
            for ik in k2states:
                if ik != n:
                    for il in lstates:
                        if il != n:
                            sum_Vnl_Vlk_Vkn += V_ij(n,il)*V_ij(il,ik)*V_ij(ik,n) * (1 / (dE_ni(n,il)*dE_ni(n,ik)))
            sum_Vnl_Vlk_Vkn_qa = qapply(sum_Vnl_Vlk_Vkn)
    
            # Second term: Single summation
            for il in lstates:
                if il != n:
                    sum_Vnl += -V_ij(n,il)*V_ij(il,n) * (1 / (dE_ni(n,il)**2))
            sum_Vnl = V_ij(n,n)*sum_Vnl
            sum_Vnl_qa = qapply(sum_Vnl)
            
            total = sum_Vnl_Vlk_Vkn + sum_Vnl
            total_qa = sum_Vnl_Vlk_Vkn_qa + sum_Vnl_qa
        
        if self.verbose:
            if order == 1:
                display(libsympy.Math(r"E_n^1=\left\langle\psi_n^0\left|H^{\prime}\right| \psi_n^0\right\rangle"))
                if type(psi0(n)) in (Ket, SHOKet):
                    libsympy.pprints(rf"E_n^({order})", total,
                                     self.output_style)
                else:
                    libsympy.pprints(rf"E_n^({order})", integral, 
                                     self.output_style)
            if order == 2:
                display(libsympy.Math(r"E_n^2=\sum_{m \neq n} \frac{\left|\left\langle\psi_m^0\left|H^{\prime}\right| \psi_n^0\right\rangle\right|^2}{E_n^0-E_m^0}"))
                libsympy.pprints(rf"E_n^({order})", total,
                                 self.output_style)
            if order == 3:
                display(libsympy.Math(r"todo write E_n^3=\sum_{m \neq n} \frac{\left|\left\langle\psi_m^0\left|H^{\prime}\right| \psi_n^0\right\rangle\right|^2}{E_n^0-E_m^0}"))
        
        res = Eq( var(rf"E_n^({order})"), total_qa)
        return(res)

#### ---->  n^th order nondegenerate perturbation for wave function.
    def psin_ND_PT(
        self,
        order: int,
        psi0c: Callable[..., Any],
        psi0: Callable[..., Any],
        Hp: Any,
        En0: Callable[..., Any],
        k1min: Any = k1,
        k1max: Any = k1,
        k2min: Any = k2,
        k2max: Any = k2,
        k3min: Any = k3,
        k3max: Any = k3
    ) -> Any:
        """
        Parameters
        ----------
        order : TYPE
            DESCRIPTION.
        psi0 : TYPE
            DESCRIPTION.
        psi0c : TYPE
            DESCRIPTION.
        Hp : TYPE
            DESCRIPTION.
        En0 : TYPE
            DESCRIPTION.
        imin : TYPE
            DESCRIPTION.
        imax : TYPE
            DESCRIPTION.

        Returns
        -------
        |perturbed ket>

        res = oqmec.psin_ND_PT(3, psi0b, psi0k, Hp, En0, k1min=k1, k1max=k1, k2min=k2, k2max=k2, k3min=k3, k3max=k3)
        
        OLD CODE:
        =========            
        psi0  = |n0>, Unperturbed states, wavefunction or eigenket.
        psi0c = <m0|, Unperturbed states, dual (complex conjugate) of the 
                      wavefunction or eigenbra.
        Hp = perturbation term added to the Hamiltonian.
        En0 = Unperturbed eigen-energy value function.
        imin = 
        imax =
        
        total_even = 0
        nums = [1,2,3,4,5]
        s = [total_even := total_even+x for x in nums if x % 2 == 0] # > [2,6]
        
        total = 0
        indices = Range(n-2,n+3)
        s = [total := total+s2.subs({k:i}) for i in indices if i != n]
        s = [total := total+norm2_matrix_s2.subs({k:i}) for i in indices if i != n]
        
        
        for i in Range(imin, imax):
            if i != n:
                # if psi0 is not a wavefunction
                S = S + norm2_matrix_s2.subs({k:i}) / (En0(n).rhs-En0(i).rhs)
        """
        V_ij  = lambda i=i,j=j: qapply(psi0c(i)*qapply(Hp*psi0(j))) # <i0|V|j0> Application of qapply step by step gives result faster than multiplication.
        dE_ni = lambda n=n,i=i: (En0(n).rhs-En0(i).rhs)  # E_n^(0)-E_i^(0)
        
        # 1st order perturbation to wave function.
        if order == 1:
            total = 0
            k1states = Range(k1min, k1max+1)
            if type(psi0(n)) in (Ket, SHOKet):
                sum_list = [total := total + V_ij(ik1,n)*psi0(ik1) * (1/dE_ni(n,ik1)) for ik1 in k1states if ik1 != n]
                total_qa = total
            
            """
            OLD CODE
            ========
            total = 0
            matrix_s1 = psi0c(k)*Hp*psi0(n) # Braket case. <m|H'|n>
            matrix_s2 = qapply(matrix_s1)
            matrix_s3 = simplify(matrix_s2)
            indices   = Range(imin, imax+1)
            # sum(<m|H'|n>*|m0> / (En0-Em0), m!=n)
            # Iterate over i; symbolicly i=k.
            sum_list = [total := total + matrix_s3.subs({k:i})*psi0(i)/(En0(n).rhs-En0(i).rhs) for i in indices if i != n]
            res = Eq( var(r"\psi_n^1"), total )
            """
            
            """
            total = 0
            k1states = Range(k1min, k1max+1)
            
            if type(psi0(n)) in (Ket, SHOKet):
                sum_Vkn_psi0k = sum([V_ij(k,n)*psi0(k) * (1/dE_ni(n,k)) for k in Range(k1min, k1max+1) if k != n])
                total_qa = total
            """
        
        # 2nd order perturbation to wave function.            
        if order == 2:
            sum_Vkl_Vln = 0
            sum_Vkn_Vnn = 0
            sum_Vkn     = 0
            k1states = Range(k1min, k1max+1)
            k2states = Range(k2min, k2max+1)
            
            # First term: Double summation
            for ik1 in k1states:
                if ik1 != n:
                    for ik2 in k2states:
                        if ik2 != n:
                            sum_Vkl_Vln += V_ij(ik1,ik2)*V_ij(ik2,n)*psi0(ik1) / (dE_ni(n,ik1)*dE_ni(n,ik2))
                            sum_Vkl_Vln_qa = qapply(sum_Vkl_Vln)
        
            # Second term: Single summation
            for ik1 in k1states:
                if ik1 != n:
                    sum_Vkn_Vnn += -V_ij(ik1,n)*V_ij(n,n)*psi0(ik1) / (dE_ni(n,ik1)**2)
                    sum_Vkn_Vnn_qa = qapply(sum_Vkn_Vnn)
        
            # Third term: Single summation          
            for ik1 in k1states:
                if ik1 != n:
                    sum_Vkn += -S(1)/2*V_ij(n,ik1)*V_ij(ik1,n)*psi0(n) / (dE_ni(ik1,n)**2)
                    sum_Vkn_qa = qapply(sum_Vkn)
                   
            total = sum_Vkl_Vln + sum_Vkn_Vnn + sum_Vkn
            total_qa = sum_Vkl_Vln_qa + sum_Vkn_Vnn_qa + sum_Vkn_qa
        
        # 3rd order perturbation to wave function.  
        if order == 3:
            sum_Vk1k2_Vk2k3_Vk3n, sum_Vnn_Vk1k2_Vk2n, sum_Vnn_Vnn_Vk1n   = zeros(1,3)
            sum_Vnk2_Vnk2_Vk1n, sum_Vnk2_Vk2k1_Vk1n, sum_Vk2n_Vk1k2_Vnk1 = zeros(1,3)
            sum_Vk1_Vnk1_Vnn = 0
            k1states = Range(k1min, k1max+1)
            k2states = Range(k2min, k2max+1)
            k3states = Range(k3min, k3max+1)
            
        # First term 
            for ik1 in k1states:
                if ik1 != n:
                    for ik2 in k2states:
                        if ik2 != n:
                            for ik3 in k3states:
                                if ik3 != n:
                                    sum_Vk1k2_Vk2k3_Vk3n += -V_ij(ik1,ik2)*V_ij(ik2,ik3)*V_ij(ik3,n)*psi0(ik1) / (dE_ni(ik1,n)*dE_ni(n,ik2)*dE_ni(n,ik3))
                                    sum_Vk1k2_Vk2k3_Vk3n_qa = qapply(sum_Vk1k2_Vk2k3_Vk3n)
                    
        # Second term
            for ik1 in k1states:
                if ik1 != n:
                    for ik2 in k2states:
                        if ik2 != n:
                            sum_Vnn_Vk1k2_Vk2n += V_ij(n,n)*V_ij(ik1,ik2)*V_ij(ik2,n)*(1/dE_ni(n,ik1) + 1/dE_ni(n,ik2))*psi0(ik1) / (dE_ni(ik1,n)*dE_ni(n,ik2))
                            sum_Vnn_Vk1k2_Vk2n_qa = qapply(sum_Vnn_Vk1k2_Vk2n)
    
        # Third term
            for ik1 in k1states:
                if ik1 != n:
                    sum_Vnn_Vnn_Vk1n += -V_ij(n,n)*V_ij(n,n)*V_ij(ik1,n)*psi0(ik1) / (dE_ni(ik1,n)**3)
                    sum_Vnn_Vnn_Vk1n_qa = qapply(sum_Vnn_Vnn_Vk1n)
                    
        # Fourth term
            for ik1 in k1states:
                if ik1 != n:
                    for ik2 in k2states:
                        if ik2 != n:
                            sum_Vnk2_Vnk2_Vk1n += V_ij(n,ik2)*V_ij(ik2,n)*V_ij(ik1,n)*(1/dE_ni(n,ik1)+S(1)/2*dE_ni(n,ik2))*psi0(ik1) / (dE_ni(ik1,n)*dE_ni(n,ik2))
                            sum_Vnk2_Vnk2_Vk1n_qa = qapply(sum_Vnk2_Vnk2_Vk1n)
                            
        # Fifth term
            for ik1 in k1states:
                if ik1 != n:
                    for ik2 in k2states:
                        if ik2 != n:
                            sum_Vnk2_Vk2k1_Vk1n += -V_ij(n,ik2)*V_ij(ik2,ik1)*V_ij(ik1,n)*psi0(n) / (2*dE_ni(n,ik2)**2*dE_ni(n,ik1))
                            sum_Vnk2_Vk2k1_Vk1n_qa = qapply(sum_Vnk2_Vk2k1_Vk1n)
                            
        # Sixth term
            for ik1 in k1states:
                if ik1 != n:
                    for ik2 in k2states:
                        if ik2 != n:
                            sum_Vk2n_Vk1k2_Vnk1 += -V_ij(ik2,n)*V_ij(ik1,ik2)*V_ij(n,ik1)*psi0(n) / (2*dE_ni(n,ik2)**2*dE_ni(n,ik1))
                            sum_Vk2n_Vk1k2_Vnk1_qa = qapply(sum_Vk2n_Vk1k2_Vnk1)
                            
        # Seventh term
            for ik1 in k1states:
                if ik1 != n:
                    sum_Vk1_Vnk1_Vnn += V_ij(n,ik1)*V_ij(ik1,n)*V_ij(n,n)*psi0(n) / (dE_ni(n,ik1)**3)
                    sum_Vk1_Vnk1_Vnn_qa = qapply(sum_Vk1_Vnk1_Vnn)
                    
            total = sum_Vk1k2_Vk2k3_Vk3n + sum_Vnn_Vk1k2_Vk2n + sum_Vnn_Vnn_Vk1n + sum_Vnk2_Vnk2_Vk1n + sum_Vnk2_Vk2k1_Vk1n + sum_Vk2n_Vk1k2_Vnk1 + sum_Vk1_Vnk1_Vnn
            total_qa = sum_Vk1k2_Vk2k3_Vk3n_qa + sum_Vnn_Vk1k2_Vk2n_qa + sum_Vnn_Vnn_Vk1n_qa + sum_Vnk2_Vnk2_Vk1n_qa + sum_Vnk2_Vk2k1_Vk1n_qa + sum_Vk2n_Vk1k2_Vnk1_qa + sum_Vk1_Vnk1_Vnn_qa
        
            
        if self.verbose:
            if order == 1:
                if type(psi0(n)) in (Ket, SHOKet):
                    display(libsympy.Math(r"\psi_n^1=\sum_{m \neq n} \frac{\left\langle\psi_m^0\left|H^{\prime}\right| \psi_n^0\right\rangle}{\left(E_n^0-E_m^0\right)} \psi_m^0"))
                else:
                    libsympy.pprints(rf"E_n^({order})", integral, 
                                     self.output_style)
            
            if order == 2:
                display(libsympy.Math(r"""
                                      \psi_n^2  =  \sum_{k \neq n} \sum_{\ell \neq n}\left|k^{(0)}\right\rangle \frac{\left\langle k^{(0)}|V| \ell^{(0)}\right\rangle\left\langle\ell^{(0)}|V| n^{(0)}\right\rangle}{\left(E_n^{(0)}-E_k^{(0)}\right)\left(E_n^{(0)}-E_{\ell}^{(0)}\right)}
                                      -\sum_{k \neq n}\left|k^{(0)}\right\rangle \frac{\left\langle k^{(0)}|V| n^{(0)}\right\rangle\left\langle n^{(0)}|V| n^{(0)}\right\rangle}{\left(E_n^{(0)}-E_k^{(0)}\right)^2}
                                      - \frac{1}{2} \left|n^{(0)}\right\rangle \sum_{k \neq n} \frac{\left|\left\langle k^{(0)}|V| n^{(0)}\right\rangle\right|^2}{\left(E_n^{(0)}-E_k^{(0)}\right)^2}"""
                                      ))
                """                    
                libsympy.pprints("sum_Vkl_Vln", sum_Vkl_Vln, 
                                 "sum_Vkl_Vln_qa", sum_Vkl_Vln_qa,
                                 self.output_style)
                """
            
            if order == 3:
                display(libsympy.Math(r"todo E_n^3=\sum_{m \neq n} \frac{\left|\left\langle\psi_m^0\left|H^{\prime}\right| \psi_n^0\right\rangle\right|^2}{E_n^0-E_m^0}"))
            
            libsympy.pprints("Total", total,
                             self.output_style)
            
        res = Eq( var(rf"\psi_n^({order})"), total_qa)
        return(res)

    @staticmethod
    def __doc__():
        return("Document of <template> class.")
        
oqmec = quantum_mechanics()
oqmec.verbose = True
