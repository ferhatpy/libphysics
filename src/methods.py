#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
template.py
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

    """
    _name = "methods"
    class_type = "default"

    
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
    
    def __init__(self):
        super().__init__()
        self.define_symbols()

        class cderivatives:
            """
            Sub formulary class for derivatives.
            """
            def __init__(self):
                # Trigonometric Functions.
                self.dsin   = Eq(Derivative(sin(x),x), sin(x).diff(x))
                self.dcos   = Eq(Derivative(cos(x),x), cos(x).diff(x))
                self.dtan   = Eq(Derivative(tan(x),x), tan(x).diff(x))
                self.dcot   = Eq(Derivative(cot(x),x), cot(x).diff(x))
                self.dsec   = Eq(Derivative(sec(x),x), sec(x).diff(x))
                self.dcsc   = Eq(Derivative(csc(x),x), csc(x).diff(x))
                
                # Hyperbolic Functions.
                self.dsinh  = Eq(Derivative(sinh(x),x), sinh(x).diff(x))
                self.dcosh  = Eq(Derivative(cosh(x),x), cosh(x).diff(x))
                self.dtanh  = Eq(Derivative(tanh(x),x), tanh(x).diff(x))
                self.dcoth  = Eq(Derivative(coth(x),x), coth(x).diff(x))
                self.dsech  = Eq(Derivative(sech(x),x), sech(x).diff(x))
                self.dcsch  = Eq(Derivative(csch(x),x), csch(x).diff(x))
        class cintegrals:
            """
            Sub formulary class for integrals.
            """
            def __init__(self):
                # Trigonometric Functions.
                self.isin   = Eq(Integral(sin(x),x), integrate(sin(x),x))
                self.icos   = Eq(Integral(cos(x),x), integrate(cos(x),x))
                self.itan   = Eq(Integral(tan(x),x), integrate(tan(x),x))
                self.icot   = Eq(Integral(cot(x),x), integrate(cot(x),x))
                self.isec   = Eq(Integral(sec(x),x), integrate(sec(x),x))
                self.icsc   = Eq(Integral(csc(x),x), integrate(sin(x),x))
                
                # Hyperbolic Functions.
                self.isinh  = Eq(Integral(sinh(x),x), integrate(sinh(x),x))
                self.icosh  = Eq(Integral(cosh(x),x), integrate(cosh(x),x))
                self.itanh  = Eq(Integral(tanh(x),x), integrate(tanh(x),x))
                self.icoth  = Eq(Integral(coth(x),x), integrate(coth(x),x))
                self.isech  = Eq(Integral(sech(x),x), integrate(sech(x),x))
                self.icsch  = Eq(Integral(csch(x),x), integrate(sinh(x),x))
        self.cderivatives = cderivatives()        
        self.cintegrals   = cintegrals()
        
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