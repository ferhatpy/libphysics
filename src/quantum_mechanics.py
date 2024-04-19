#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
quantum_mechanics.py
Created on Fri Mar 11 12:53:36 2022

"""
# from abc import ABC, abstractmethod
from sympy.abc import*
from sympy import*
from sympy import Derivative as D
from sympy.assumptions.assume import global_assumptions
from sympy.assumptions import assuming, Q, ask
from sympy.integrals.manualintegrate import manualintegrate
from sympy.integrals.manualintegrate import integral_steps
#from sympy.integrals.rubi.utility_function import Simplify
from sympy.integrals.transforms import inverse_fourier_transform
from sympy.diffgeom.rn import *
from sympy.vector import CoordSys3D

# Hydrogen Atom
from sympy.physics.hydrogen import *

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
import libphyscon as pc

class quantum_mechanics(branch):
    """
    solve(oqmec.expX, _expX)
    """

    _name = "quantum_mechanics"

    # Integer symbols
    global k,k1,k2,k3,n
    k, k1,k2,k3 = symbols('k k1 k2 k3', integer=True)
    n = symbols('n', positive=True, integer=True)
    
    def define_symbols(self):
        """
        Common global symbols, functions.
        a: 
        F: 
        """
        global _U
        
        # Real symbols
        global eps0
        eps0 = symbols('varepsilon_0', real=True, positive=True)
        global A,a,b
        A,a,b = symbols('A a b', real=True, positive=True)
        global m,w,t
        m,w,t = symbols('m w t', real=True)
        # Mass of electron, proton, neutron, nucleus.
        global m_e,m_p,m_n,m_N
        m_e,m_p,m_n,m_N = symbols('m_e m_p m_n m_N', real=True, positive=True)
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

        """
        class_type specific global symbols, functions, objects.
        """ 
        if self.class_type in ["position_space"]:
            global Hp
            Hp = Operator('Hprime') # Perturbation onto Hamiltonian
            
            global psib, psik, nb, nk, phi
            global Psi, psi, psix, PsiSph, psiSph
            psib, psik, nb, nk = [Bra('psi'), Ket('psi'), Bra('n'), Ket('n')]
            Psi    = Function('Psi')(x,y,z,t)
            psi    = Function('psi')(x,y,z)
            psix   = Function('psi')(x)
            
            PsiSph = Function('Psi')(r,theta,phi,t)
            psiSph = Function('psi')(r,theta,phi)
            
        if self.class_type in ["momentum_space"]:
            phi  = Function('phi')(k)

    def __init__(self, class_type='position_space'):
        """
        class_type = \
        {1:"position_space",
         2:"momentum_space"}
        """
        super().__init__()
        self.class_type = class_type
        self.define_symbols()
        
        # File settings
        self.input_dir  = "input/quantum_mechanics"
        self.output_dir = "output/quantum_mechanics"
        
        class subformulary:
            """
            Sub formulary class.
            
            Define global symbols in the outer class.
            """
            def __init__(self, class_type):
                if class_type in ["position_space"]:
                    self.Icm_sphere = S(2)/5*M*r**2
        self.subformulary = subformulary(self.class_type)
        
#### # Position Space
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
            self.psix = psix
            self.normalization = Eq(Integral(conjugate(self.Psi)*self.Psi, (x,xmin,xmax)), 1)
            
            #----> Position related formulary.
            # Do not use Lambda of sympy it requires definifion of fx. (NameError: name 'fx' is not defined)
            # Generic expectation value calculation routine for a given operator fx.
            self.exp_fx   = lambda fx:Eq(var(r'\langle{'+str(latex(fx))+r'}\rangle'), Integral(conjugate(self.Psi)*fx*self.Psi, (x,xmin,xmax)))
            self.exp_x    = Eq(var(r'\langle{x}\rangle'),   Integral(conjugate(self.Psi)*x*self.Psi, (x,xmin,xmax)))
            self.exp_x2   = Eq(var(r'\langle{x^2}\rangle'), Integral(conjugate(self.Psi)*x**2*self.Psi, (x,xmin,xmax)))
            self.delta_x  = Eq(var(r'\Delta{x}=x-\langle{x}\rangle'), sqrt(self.exp_x2.rhs - self.exp_x.rhs**2))
            self.delta_x2 = Eq(var(r'(\Delta{x})^2=\langle{x^2}\rangle-\langle{x}\rangle^2'), self.exp_x2.rhs - self.exp_x.rhs**2)
            
#### # ----> Momentum related formulary.
            """
            oqmec.px.subs(Psi,x**2).doit()
            oqmec.px.subs({oqmec.Psi:x**2}).doit()
            oqmec.px.replace(Psi,x**2).doit()
            oqmec.px.xreplace({oqmec.Psi:x**2}).doit()  psi=x**2 --> px=-2*hb*I*x
            oqmec.pxop.xreplace({oqmec.Psi:x**2}).doit().rhs.expr
            oqmec.pxop.subs({oqmec.Psi:x**2}).doit().rhs.expr
            
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
            self.delta_XP  = self.uncertainityXP  = Eq(var(r'\Delta{x}\Delta{p_x}'), self.delta_x.rhs*self.delta_px.rhs)
            
#### # ----> Orbital Angular Momentum
            """
            
            """
#            self.Lx      = Eq(var('L_x'), var('L_x'))
#            self.Ly      = Eq(var('L_y'), var('L_y'))
            self.Lx,self.Ly,self.Lz = [Operator('\hat{L}_x'), Operator('\hat{L}_y'), Operator('\hat{L}_z')]
            self.L2      = Eq(Operator('\hat{L}^2'), self.Lx**2 + self.Ly**2 + self.Lz**2)
            self.comLxLy = Eq(Commutator(self.Lx, self.Ly), i*hbar*self.Lz)
            self.comLzLx = Eq(Commutator(self.Lz, self.Lx), i*hbar*self.Ly)
            self.comLyLz = Eq(Commutator(self.Ly, self.Lz), i*hbar*self.Lx)
            self.Lplus   = Eq(Operator('\hat{L}_+'), self.Lx + I*self.Ly)
            self.Lminus  = Eq(Operator('\hat{L}_-'), self.Ly - I*self.Ly)
            
            self.comJxJy = Eq(Commutator(Jx, Jy), i*hbar*Jz)
            self.comJzJx = Eq(Commutator(Jz, Jx), i*hbar*Jy)
            self.comJyJz = Eq(Commutator(Jy, Jz), i*hbar*Jx)
            self.comJ2Ji = lambda i=x:Eq(Commutator(J2, {x:Jx, y:Jy, z:Jz}[i]), 0)
            self.comJzJp = Eq(Commutator(Jz, Jplus),   hbar*Jplus)
            self.comJzJm = Eq(Commutator(Jz, Jminus), -hbar*Jminus)
            self.comJpJm = Eq(Commutator(Jplus, Jminus), 2*hbar*Jz)
            self.J2ket   = Eq(J2*JzKet(j,m), simplify(qapply(J2*JzKet(j,m))))
            self.Jzket   = Eq(Jz*JzKet(j,m), simplify(qapply(Jz*JzKet(j,m))))
            self.Jpket   = Eq(Jplus*JzKet(j,m), simplify(qapply(Jplus*JzKet(j,m))))
            self.Jmket   = Eq(Jminus*JzKet(j,m), simplify(qapply(Jminus*JzKet(j,m))))

            
#### # ----> Total Angular Momentum in Spherical Coordinates.
            self.JxSph    = Eq(Operator('\hat{J}_x'),  hbar/I*(-sin(phi)*D(psiSph, theta) - cos(phi)*cot(theta)*D(psiSph, phi)) )
            self.JySph    = Eq(Operator('\hat{J}_y'),  hbar/I*( cos(phi)*D(psiSph, theta) - sin(phi)*cot(theta)*D(psiSph, phi)) )
            self.JzSph    = Eq(Operator('\hat{J}_z'),  hbar/I*(Derivative(psiSph, phi)) )
            self.JpSph    = Eq(Operator('\hat{J}_+'),  hbar*exp( I*phi)*( D(psiSph, theta) + I*cot(theta)*D(psiSph, phi)) )
            self.JmSph    = Eq(Operator('\hat{J}_-'),  hbar*exp(-I*phi)*(-D(psiSph, theta) + I*cot(theta)*D(psiSph, phi)) )
            self.J2Sph    = Eq(Operator('\hat{J}^2'), -hbar**2*(1/sin(theta)*D(sin(theta)*D(psiSph, theta), theta) + 1/(sin(theta))**2*D(psiSph, phi, 2)) )

            
#### # ----> Spin Angular Momentum
            self.sigmax   = Eq(var(r'\sigma_x'), UnevaluatedExpr(Matrix([[0,1], [1,0]])))
            self.sigmay   = Eq(var(r'\sigma_y'), UnevaluatedExpr(Matrix([[0,-I], [I,0]])))
            self.sigmaz   = Eq(var(r'\sigma_z'), UnevaluatedExpr(Matrix([[1,0], [0,-1]])))
#            self.Sx       = Eq(Operator('\hat{S}_x'), Operator(UnevaluatedExpr(Matrix([[0,1], [1,0]]))))
            self.Sx       = Eq(Operator('\hat{S}_x'), hbar/2*self.sigmax.rhs)
            self.Sy       = Eq(Operator('\hat{S}_y'), hbar/2*self.sigmay.rhs)
            self.Sz       = Eq(Operator('\hat{S}_z'), hbar/2*self.sigmaz.rhs)

            
            #----> Operator Definitions via SymPy sympy.physics.quantum class todo ERROR in DifferantialOperator application
            #----> Position related formulary.
            self.exp_xop   = Eq(var(r'\langle{x}\rangle'),   Integral(conjugate(self.Psi)*qapply(x*self.Psi), (x,xmin,xmax)))
            self.exp_x2op  = Eq(var(r'\langle{x^2}\rangle'), Integral(conjugate(self.Psi)*qapply(x**2*self.Psi), (x,xmin,xmax)))
            self.delta_xop = Eq(var(r'\Delta{x}'), sqrt(self.exp_x2op.rhs - self.exp_xop.rhs**2))
            self.delta_x2op= Eq(var(r'(\Delta{x})^2'), self.exp_x2op.rhs - self.exp_xop.rhs**2)
            
            #----> Momentum
            self.pxop     =  Eq(var(r'\hat{p}_x'), DifferentialOperator(-I*hbar*D(self.Psi, x), self.Psi))
            self.exp_pxop  = Eq(var(r'\langle{\hat{p}_x}\rangle'),   Integral(conjugate(self.Psi)*qapply(self.pxop.rhs*self.Psi), (x,xmin,xmax)))
            self.exp_px2op = Eq(var(r'\langle{\hat{p}_x^2}\rangle'), Integral(conjugate(self.Psi)*qapply(self.pxop.rhs*self.pxop.rhs*self.Psi), (x,xmin,xmax)))
            self.delta_pxop  = Eq(var(r'\Delta{p_x}'), sqrt(self.exp_px2op.rhs - self.exp_pxop.rhs**2))
            self.delta_px2op = Eq(var(r'(\Delta{p_x})^2'), self.exp_px2op.rhs - self.exp_pxop.rhs**2)
            self.delta_xop_pxop = self.uncertainityXP  = Eq(var(r'\Delta{x}\Delta{p_x}'), self.delta_xop.rhs*self.delta_px2op.rhs)

            
#### # ----> Schrödinger Equation
            self.V = self.potential_energy = Function('V')(x)
            self.H = Eq(S('H'), -(hbar**2)/(2*m)*diff(self.Psi, x, 2) + self.V*self.Psi)
            self.exp_H = Eq(var(r'\langle{H}\rangle'), Integral(conjugate(self.Psi)*self.H.rhs, (x,xmin,xmax)))
            self.Hop = Eq(S('H'), hbar**2/(2*m)*(DifferentialOperator(D(self.Psi, x), self.Psi))**2 + Operator(self.V*self.Psi))
            self.exp_Hop = Eq(var(r'\langle{H}\rangle'), Integral(conjugate(self.Psi)*qapply(self.Hop.rhs*self.Psi), (x,xmin,xmax)))
            self.SchrodingerEq = Eq(-(hbar**2)/(2*m)*diff(self.Psi, x, 2) + self.V*self.Psi, I*hbar*diff(self.Psi, t, 1))
            self.SchrodingerEq_TI = self.SchrodingerEq_Time_Independent = Eq(-(hbar**2)/(2*m)*diff(self.Psi,x,2) + self.V*self.Psi , En()*self.Psi)
            
            #----> Fourier Transforms
            # todo look Schwabl not Griffiths
            self.FourierTransform_phi_t0 = Eq( Psi.subs(t,0), 1/(sqrt(2*pi))*Integral(phi*exp(I*k*x) , (k,-oo,oo)) )
            self.FourierTransform_Psi_t0 = Eq( phi, 1/(sqrt(2*pi))*Integral(Psi.subs(t,0)*exp(-I*k*x) , (k,-oo,oo)) )
            
#### Perturbation Theory
            #----> Time-independent perturbation theory
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
            
#### # Quantum Harmonic Oscillator
            class qho(branch):
                """
                Sub class for quantum simple harmonic oscillator.
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
                    self.V = Eq(V, S(1)/2*m*w**2*self.x2op.rhs)
                    self.H = Eq(H, simplify(self.p2op.rhs/(2*m) + self.V.rhs))
                    self.Hx = Eq(H, DifferentialOperator(-hbar**2/(2*m)*Derivative(f(x),x,2) + S(1)/2*m*w**2*x**2, f(x)))
                    self.En = lambda n=n:Eq(S(f'E_{n}'), (n+S(1)/2)*hbar*w)
                    # todo 2.51
            self.qho = qho()
            
#### # Delta Function Quantum Well
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

#### # Infinite Square Quantum Well
            class iqw(branch):
                """
                Sub class for Infinite Square Quantum Well
                """
                def __init__(self):
                    super().__init__()
                    self.name = "Infinite Square Quantum Well"
                    self.psix =lambda n=n: Eq(var(rf'\psi_{n}(x)'), sqrt(2/a)*sin(n*pi*x/a))
#                    self.V = Eq(V, S(1)/2*m*w**2*self.x2op.rhs)
#                    self.H = Eq(H, simplify(self.p2op.rhs/(2*m) + self.V.rhs))
                    self.En = lambda n=n:Eq(var(f'E_{n}'), (n**2*pi**2*hbar**2)/(2*m*a**2))
            self.iqw = iqw()
            

#### # Finite Square Quantum Well            
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
            
#### # Hydrogen Atom
            class hydrogen(branch):
                """
                Sub class for Hydrogen-like atom.
                References:
                - Wikipedia Hydrogen-like atom.
                """
                def __init__(self):
                    super().__init__()
                    self.name = "Hydrogen-like atom"
                    self.a0 = self.bohr_radius = S('a_0')
                    self.mu = m_N*m_e/(m_N+m_e)
                    self.a_mu = Eq(S('a_mu'), Eq(4*pi*eps0**2*hbar**2/(self.mu*e**2), m_e/self.mu*self.a0))
                    self.Rnl = lambda n=n, l=l, r=r, Z=Z:Eq( S(f'R_{n}{l}(r)'), sqrt((2*Z/(n*self.a_mu.lhs))**3*factorial(n-l-1)/(2*n*factorial(n+l)))*exp(-Z*r/(n*self.a_mu.lhs))*(2*Z*r/(n*self.a_mu.lhs))**l*assoc_laguerre(n-l-1, 2*l+1, 2*Z*r/(n*self.a_mu.lhs)) )
                    self.psi  = lambda n=n, l=l, m=m, r=r, Z=Z, phi=phi, theta=theta:Eq(S(f'psi_{n}{l}{m}({r},{phi},{theta})'), self.Rnl(n,l,r,Z).rhs*Ynm(l,m,theta,phi))
                    self.V  = Eq(S('V(r)'), -1/(4*pi*eps0)*Z*e**2/r)
                    self.En = lambda n=n, Z=Z:Eq(S(f'E_{n}'), (Z**2*self.mu*e**4/(32*pi**2*eps0**2*hbar**2))*(1/n**2))
                    
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

#### # Momentum Space
        if self.class_type in ["momentum_space"]:
            pass
        
        # Common text definitions.
        self.Hamiltonian = self.H

#### # Global Methods
    #### Ket state to wave function convertion.
    def ket_to_wavefunction(self, n_, j_, psi0, wfpsi0, psipert):
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

    #### Nondegenerate Perturbation Theory
    #---->  n^th order nondegenerate perturbation for energy.
    def En_ND_PT(self, order, psi0c, psi0, Hp, En0, k2min=k2, k2max=k2):
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
                    sum_Vkn += V_ij(ik2,n)*V_ij(n,ik2) / (dE_ni(n,ik2))
            
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
                            sum_Vnl_Vlk_Vkn += V_ij(n,il)*V_ij(il,ik)*V_ij(ik,n) / (dE_ni(n,il)*dE_ni(n,ik))
            sum_Vnl_Vlk_Vkn_qa = qapply(sum_Vnl_Vlk_Vkn)
    
            # Second term: Single summation
            for il in lstates:
                if il != n:
                    sum_Vnl += -V_ij(n,il)*V_ij(il,n) / (dE_ni(n,il)**2)
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

    #---->  n^th order nondegenerate perturbation for wave function.
    def psin_ND_PT(self, order, psi0c, psi0, Hp, En0, 
                   k1min=k1, k1max=k1,
                   k2min=k2, k2max=k2,
                   k3min=k3, k3max=k3
                   ):
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
                sum_list = [total := total + V_ij(ik1,n)*psi0(ik1)/(dE_ni(n,ik1)) for ik1 in k1states if ik1 != n]
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
                sum_list = [total := total + V_ij(ik1,n)*psi0(ik1)/(dE_ni(n,ik1)) for ik1 in k1states if i != n]
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