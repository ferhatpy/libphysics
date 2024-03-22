#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
methods.py
Created on Fri Mar 11 12:53:36 2022

Find and replace template with desired class name.

SymPy 1.10.1 documentation » Reference Documentation » Matrices » Vector

"""
from sympy import *
from sympy.abc import *
from sympy.diffgeom import *
from sympy.diffgeom.rn import *
from sympy.diffgeom.rn import R3_r, R3_s
from sympy.vector import CoordSys3D

from libreflection import *
import libphyscon as pc

class methods(branch):
    """
    Displaying Tables: display(*ometh.cYnm.Table())
    """
    _name = "methods"
    
    def define_symbols(self):
        """
        Common global symbols, functions.
        a: 
        F: 
        """
        global Gw
        global C
        
        global vA,vB,vC,vD
        vA,vB,vC,vD = symbols('Abold Bbold Cbold Dbold', vector=True)
        global Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz
        Ax,Bx,Cx,Dx = symbols('A_x B_x C_x D_x', real=True)
        Ay,By,Cy,Dy = symbols('A_y B_y C_y D_y', real=True)
        Az,Bz,Cz,Dz = symbols('A_z B_z C_z D_z', real=True)
        
        C         = CoordSys3D('C')
#        G         = Lambda((t,tau), Function('G')(t,tau)) # Callable Green's function.
#        Gw        = Lambda((w), Function('Gtilde')(w))    # Callable Green's tilde function.
        G         = Function('G')(t,tau)    # Incallable Green's function.
        Gw        = Function('Gtilde')(w)   # Incallable Green's tilde function.
        
        global dl
        dl = 1*C.i + 1*C.j + 1*C.k  # Differential position element.
    
    def __init__(self, class_type='default'):
        super().__init__()
        self.class_type = class_type
        self.define_symbols()
        
        # File settings
        self.input_dir  = "input/methods"
        self.output_dir = "output/methods"

        class cderivatives:
            """
            Sub Formulary Class for Derivatives.
            
            ometh.cderivatives.cosh.replace(u,x**2).rhs.doit()
            """
            global u 
            u = Function('u')(x)
            
            def __init__(self):
                self.name   = "Derivatives"
                
                # Trigonometric Functions.
                self.sin   = Eq(Derivative(sin(u),x), sin(u).diff(x))
                self.cos   = Eq(Derivative(cos(u),x), cos(u).diff(x))
                self.tan   = Eq(Derivative(tan(u),x), tan(u).diff(x))
                self.cot   = Eq(Derivative(cot(u),x), cot(u).diff(x))
                self.sec   = Eq(Derivative(sec(u),x), sec(u).diff(x))
                self.csc   = Eq(Derivative(csc(u),x), csc(u).diff(x))
                
                # Hyperbolic Functions.
                self.sinh  = Eq(Derivative(sinh(u),x), sinh(u).diff(x))
                self.cosh  = Eq(Derivative(cosh(u),x), cosh(u).diff(x))
                self.tanh  = Eq(Derivative(tanh(u),x), tanh(u).diff(x))
                self.coth  = Eq(Derivative(coth(u),x), coth(u).diff(x))
                self.sech  = Eq(Derivative(sech(u),x), sech(u).diff(x))
                self.csch  = Eq(Derivative(csch(u),x), csch(u).diff(x))
        
        class cintegrals:
            """
            Sub Formulary Class for Integrals.
            """
            def __init__(self):
                self.name   = "Indefinite Integrals"
                # Trigonometric Functions.
                self.sin   = Eq(Integral(sin(x),x), integrate(sin(x),x))
                self.cos   = Eq(Integral(cos(x),x), integrate(cos(x),x))
                self.tan   = Eq(Integral(tan(x),x), integrate(tan(x),x))
                self.cot   = Eq(Integral(cot(x),x), integrate(cot(x),x))
                self.sec   = Eq(Integral(sec(x),x), integrate(sec(x),x))
                self.csc   = Eq(Integral(csc(x),x), integrate(sin(x),x))
                
                # Hyperbolic Functions.
                self.sinh  = Eq(Integral(sinh(x),x), integrate(sinh(x),x))
                self.cosh  = Eq(Integral(cosh(x),x), integrate(cosh(x),x))
                self.tanh  = Eq(Integral(tanh(x),x), integrate(tanh(x),x))
                self.coth  = Eq(Integral(coth(x),x), integrate(coth(x),x))
                self.sech  = Eq(Integral(sech(x),x), integrate(sech(x),x))
                self.csch  = Eq(Integral(csch(x),x), integrate(sinh(x),x))
                
        class cYnm:
            """
            Sub Formulary Class for Spherical Harmonics.
            """
            def __init__(self):
                self.name   = "Spherical Harmonics"
                self.expand = Eq(Ynm(n,m,theta,phi), Ynm(n,m,theta,phi).expand(func=True))
            def Table(self, ns=Range(3), ms=Range(3)):
                self.table = []
                for i in ms:
                    for j in ns:
                        if j>=i:self.table.append( Eq(Ynm(j,i,theta,phi), simplify(Ynm(j,i,theta,phi).expand(func=True))) )
                return self.table                               
                
        self.cderivatives = cderivatives()        
        self.cintegrals   = cintegrals()
        self.cYnm = self.spherical_harmonics = cYnm()
        
        class subformulary:
            """
            Sub formulary class.
            
            Define global symbols in the outer class.
            """
            def __init__(self):
                # Transformations from Polar to Cartesian Coordinates.
                self.pol_to_cart_x = r*cos(theta)*C.i
                self.pol_to_cart_y = r*sin(theta)*C.j
                # Transformations from Cylindrical to Cartesian Coordinates.
                self.cyl_to_cart_x = r*sin(phi)*C.i
                self.cyl_to_cart_y = r*sin(phi)*C.j
                self.cyl_to_cart_z = z*C.k
                # Transformations from Spherical to Cartesian Coordinates.
                self.sph_to_cart_x = r*sin(theta)*cos(phi)*C.i
                self.sph_to_cart_y = r*sin(theta)*sin(phi)*C.j
                self.sph_to_cart_z = r*cos(theta)*C.k
        self.subformulary = subformulary()
        
        if self.class_type == "default":
            #----Green's Function Methods
            self.G  = G     # Green's function.
            self.Gw = Gw    # Green's tilde function.
            self.IFT_Gw = 1/sqrt(2*pi)*Integral(Gw*exp(I*w*(t-tau)), (w))
            # IFT_Gw = 1/sqrt(2*pi)*Integral(Gw*exp(I*w*(t-tau)), (w))
            self.IFT_Dirac_delta = 1/(2*pi)*Integral(exp(I*w*(t-tau)), (w))

            # Common text definitions.
            self.Greens_function  = self.G
            self.Greensw_function = self.Gw
            self.inverse_Fourier_transform_Gw = self.IFT_Gw
            self.inverse_Fourier_transform_Dirac_delta = self.IFT_Dirac_delta
            
            
            #----Integral Transforms
            #----> Fourier Transform
            # todo
            
            #----> Laplace Transform
            self.L = self.Laplace_transform = Lambda( (f,s,t), Integral(f*exp(-s*t), (t, S.Zero, S.Infinity)) )
            # ometh.L(exp(-a*t),s,t).doit().args[0][0]

            
            #----Vector Algebra
            self.vA = Ax*C.i + Ay*C.j + Az*C.k
            self.vB = Bx*C.i + By*C.j + Bz*C.k
            self.vC = Cx*C.i + Cy*C.j + Cz*C.k
            self.vD = Dx*C.i + Dy*C.j + Dz*C.k
        
            self.dot_product = Eq(S('C'), self.vA.dot(self.vB))
            self.dot_product_theta = Eq(S('C'), self.vA.magnitude()*self.vB.magnitude()*cos(theta))
            self.cross_product = Eq(self.vC, self.vA.cross(self.vB))
            self.cross_product_theta = Eq(S('C'), self.vA.magnitude()*self.vB.magnitude()*sin(theta))
            self.triple_scalar_product = Eq(S('D'), self.vA.dot(self.vB.cross(self.vC)))
            self.triple_vector_product = Eq(self.vD, self.vA.cross(self.vB.cross(self.vC)))
#            self.triple_vector_product = Eq(self.vD, self.vA^(self.vB^(self.vC)))
        
        
    @staticmethod
    def __doc__():
        return("Document of methods class.")
        
ometh = methods()