# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
optics.py
Created on Fri Mar 11 12:53:36 2022

"""
from sympy.abc import*
from sympy import*
from libreflection import *
import mpmath as mp
import numpy as np
import scipy as sp
import libphyscon as pc

class optics(branch):
    """

    """
    _name = "optics"
    
    def define_symbols(self):
        """
        Common global symbols, functions.
        a: 
        F: 
        """
        global k,m
        global x,y,z,x0,y0,z0
        global x0min,x0max,y0min,y0max
        global l,t
        global _El, Eli, _I, Uapr
        
        k,m = symbols('k m', real=True, positive=True)
        x,y,z,x0,y0,z0 = symbols('x y z x_0 y_0 z_0', nonzero=True)
        x0min,x0max,y0min,y0max = symbols('x_0min x_0max y_0min y_0max', nonzero=True)
        l,t = symbols('lambda t', real=True, positive=True)
        
        _El   = Function('E')(x,y,z)
        _I    = Function('I')(x,y,z)
        Eli   = Function('E_i')(x0,y0)
        Uapr  = Function('U_apr')(x0,y0)
    
    def __init__(self, class_type='default'):
        """
        a: 
        F:
            
        integrate(exp(I*k*y0**2/(2*z))*exp(-I*k*y*y0/z), y0)
        
        Fraunhofer_Diff_Int:    Fraunhofer diffraction integral.
        Fresnel_Diff_IntEx:     Fraunhofer diffraction integral expanded version.
        Fresnel_Diff_Int:       Fresnel diffraction integral.
        Fresnel_Diff_intgd:     Integrand of Fresnel diffraction integral.
        """
        super().__init__()
        self.class_type = class_type
        self.define_symbols()
        
        # File settings
        self.input_dir  = "input/optics"
        self.output_dir = "output/optics"
        
        
        class subformulary:
            """
            Sub formulary class.
            
            Define global symbols in the outer class.
            """
            def __init__(self):
                self.rect = Lambda(x, Piecewise( (1, abs(x)<=0.5), (0, True) ))
        self.subformulary = subformulary() 


#### --- CLASSES ---


#### ABCD
        class ABCD(branch):
            """
            Ghatak2009 Chapter5. The Matrix Method in Paraxial Optics.
            
            O = Object.
            I = Image.
            
            D1 = u  Distance between object O and optical element.
            D2 = v  Distance between optical element and image I.
            
            Usage:
            ======    
                det(oopti.ABCD.T().rhs.doit()) = 1
            """
            global n1, n2, R1, R2                       # Material specific parameters.
            global x1, x2, alpha1, alpha2, D1, D2       # Optic configuration parameters.
            global lambda1, lambda2
            n1, n2, R1, R2 = symbols('n_1, n_2, R_1, R_2')
            x1, x2, alpha1, alpha2, D1, D2 = symbols('x_1, x_2, alpha_1, alpha_2, D_1 D_2')
            lambda1, lambda2 = symbols('lambda_1, lambda_2')
            
            def __init__(self):
                super().__init__()
                self.name = "ABCD Matrix Method"
                self.lambda1 = Eq(S('lambda_1'), n1*alpha1)
                self.lambda2 = Eq(S('lambda_2'), n2*alpha2)
                self.P       = Eq(S('P'), (n2-n1)/R1)             # (29)
                self.P1      = Eq(S('P_1'), (n-1)/R1)             # ~(50)
                self.P2      = Eq(S('P_2'), (1-n)/R2)             # ~(50)
                
                self.T = lambda D1=D1, n1=n1: Eq(S('T'), 
                            UnevaluatedExpr(Matrix(((1,    0),
                                                    (D1/n1, 1)))))  # (21)
                self.R = lambda n2=n2, n1=n1, R1=R1: Eq(S('R'), 
                            UnevaluatedExpr(Matrix(((1, -(n2-n1)/R1),
                                                    (0,  1)))))     # (31)
                self.SM = Eq(S('SM'), 
                            UnevaluatedExpr(Matrix((( b, -a),
                                                    (-d,  c)))))
#### --- SUBCLASSES ---


#### ----> Single Spherical Refracting Surface
                class single_spherical_refracting_surface(branch):
                    """
                    Single Spherical Refracting Surface.
                    Image = Transverse Matrix * Object
                    Transverse Matrix = Translation * Refraction * Translation
                    """
                    def __init__(self, parent):
                        super().__init__()
                        self.name = "Single Spherical Refracting Surface"
                        
                        self.transverse_matrix = \
                            Eq(Matrix(((lambda2),
                                       (x2))),
                               MatMul(parent.T(D2,n2).rhs.doit(), parent.R().rhs.doit(), parent.T(-D1,n1).rhs.doit(), 
                                      Matrix(((lambda1),
                                              (x1))))
                               ) # (37)

                        self.matrix = \
                            Eq(Matrix(((lambda2),
                                       (x2))),
                               MatMul(Matrix(((1 + parent.P.rhs*D1/n1, -parent.P.rhs),
                                       (0                   ,  1 - parent.P.rhs*D2/n2 ))),
                               Matrix(((lambda1),
                                       (x1))))                    
                               ) # (40)
                        self.M = Eq(S('m'), n1*D2/(n2*D1))

                        self.magnification = self.M
                    @staticmethod
                    def __doc__():
                        return "Single Spherical Refracting Surface."
                self.single_spherical_refracting_surface = single_spherical_refracting_surface(self)
                

#### ----> Coaxial Optical System
                class coaxial_optical_system(branch):
                    """
                    Coaxial Optical System.
                    Transverse Matrix = Translation * abcd matrix * Translation
                    """
                    def __init__(self, parent):
                        super().__init__()
                        self.name = "Coaxial Optical System"
                        
                        self.transverse_matrix = \
                            Eq(Matrix(((lambda2),
                                       (x2))),
                               MatMul(parent.T(D2,n2).rhs.doit(), parent.SM.rhs.doit(), parent.T(-D1,n1).rhs.doit(), 
                                      Matrix(((lambda1),
                                              (x1))))
                               ) # (41)
                        
                        self.M = Eq(S('M'), c - a*D2) # (45) Magnification
                        self.matrix = \
                            Eq(Matrix(((lambda2),
                                       (x2))),
                               MatMul(UnevaluatedExpr(Matrix(((1/self.M.rhs, -a),
                                                              (0,    self.M.rhs)))),
                                      Matrix(((lambda1),
                                              (x1))))
                               ) # (47)
                    @staticmethod
                    def __doc__():
                        return "Sub class for Coaxial Optical System."
                self.coaxial_optical_system = coaxial_optical_system(self)

                
#### ----> Thin Lens
                class thin_lens(branch):
                    """
                    Thin Lens formulas.
                    Transverse Matrix = Refraction2 * Translation * Refraction1
                    """
                    def __init__(self, parent):
                        super().__init__()
                        self.name = "Thin Lens"
                        
                        self.transverse_matrix = \
                            Eq(Matrix(((lambda2),
                                       (x2))),
                               MatMul(parent.R(1,n,R2).rhs.doit(), parent.T(t,n).rhs.doit(), parent.R(n,1,R1).rhs.doit(),
                               Matrix(((lambda1),
                                       (x1))))
                               ) # (41)
                        
                        self.system_matrixP1P2 = Eq(S('SM'),
                            UnevaluatedExpr(Matrix(((1, -parent.P1.rhs-parent.P2.rhs), (0, 1))))
                            )
                        
                        self.focal_length = Eq(f, 1/((n-1)*(1/R1 - 1/R2))) # (56)
                        
                    @staticmethod
                    def __doc__():
                        return "Sub class for Thin Lens."
                self.thin_lens = thin_lens(self)
                
                
                self.refraction_matrix  = self.R
                self.translation_matrix = self.T
                self.system_matrix      = self.SM
                
            @staticmethod
            def __doc__():
                return "Sub class with ABCD  matrix method in paraxial optics formulas from Ajoy Ghatak's Optics."
        self.MatrixMethod = self.ABCD = ABCD() 
                
            
            
#### Fiber Bragg Grating
        class Fiber_Bragg_Grating(branch):
            """
            Fiber Bragg Grating formulas from Ajoy Ghatak's Optics (Appendix C)
            """
            global Gamma, lambda_B, lambda_0, lambda_, n0, Delta_n
            Gamma, lambda_B, lambda_0, lambda_, n0, Delta_n = symbols('Gamma lambda_B lambda_0 Lambda n_0 Delta_n')           
            
            def __init__(self):
                super().__init__()
                self.name = "Fiber Bragg Grating"
                
                # Fundamental parameters
                self.lambda_B = self.Bragg_wavelength = Eq(lambda_B, 2*lambda_*n0)
                self.kappa = self.coupling_coefficient = Eq(kappa, (pi*Delta_n)/lambda_0)
                self.Gamma = Eq(Gamma, 4*pi*n0*(1/lambda_0 - 1/lambda_B))
                self.alpha = Eq(alpha, sqrt(kappa**2 - Gamma**2/4))
                # Ajoy Ghatak uses the word reflectivity instead of reflectance. Hecht uses reflectance.
                self.R = self.Reflectance = Eq( S('R'), kappa**2*sinh(alpha*L)**2 / (kappa**2*cosh(alpha*L)**2 - Gamma**2/4) ) # Ghatak2009 Appendix C Eq.3
                self.R_peak = Eq(S('R_peak'), tanh(pi*Delta_n*L/lambda_B)**2)
                    
            @staticmethod
            def __doc__():
                return "Sub class with Fiber Bragg Grating formulas from Ajoy Ghatak's Optics"
        self.Fiber_Bragg_Grating = self.FBG = Fiber_Bragg_Grating()        


#### Fraunhofer Diffraction
        class Fraunhofer(branch):
            """
            Sub class for Fraunhofer Diffraction.
            
            x0: x coordinate at the aperture at z=0.
            y0: y coordinate at the aperture at z=0.
            
            Reference: Degiorgio, Photonics a Short Course. (Eq.1.61)
            """
            def __init__(self):
                super().__init__()
                self.name = "Fraunhofer Diffraction"
                self.class_type = "Fraunhofer"
                self.integrand  = Uapr*Eli*exp(-I*2*pi/(l*z)*(x*x0+y*y0))
                self.integral   = Eq(_El, ( exp(I*k*z)/(I*l*z))*(exp(I*k/(2*z)*(x**2 + y**2)))*Integral(self.integrand, (y0, y0min, y0max), (x0,x0min,x0max)) )
                self.EField     = Eq(_El, self.integral.rhs)
                self.intensity  = Eq(_I, self.EField.rhs*conjugate(self.EField.rhs))
                # self.intensity= Eq(_I, abs(self.EField.rhs)**2)
        self.Fraunhofer = Fraunhofer()


#### Fresnel Diffraction
        class Fresnel(branch):
            """
            Sub class for Fresnel Diffraction.
            """
            def __init__(self):
                super().__init__()
                self.name = "Fresnel Diffraction"
                self.class_type = "Fresnel"
                self.integrand  = Uapr*Eli*exp( I*k/(2*z)*((x-x0)**2+(y-y0)**2) )
                self.integral   = Eq(_El, (exp(I*k*z)/(I*l*z))*Integral(self.integrand, (y0, y0min, y0max), (x0,x0min,x0max)) )
                # Fresnel Diffraction Expanded Integral
                self.Fresnel_Diff_IntEx  = Eq(_El, ( exp(I*k*z)/(I*l*z))*(exp(I*k/(2*z)*(x**2 + y**2)))*Integral(Uapr*Eli*exp( I*k/(2*z)*(x0**2+y0**2))*exp(-I*2*pi/(l*z)*(x*x0+y*y0)), (y0, y0min, y0max), (x0,x0min,x0max)) )
                self.EField     = Eq(_El, self.integral.rhs)
                self.intensity  = Eq(_I, exp(I*k*z)/(I*l*z)*conjugate(exp(I*k*z)/(I*l*z))* \
                                     (Integral(re(self.integrand), (y0, y0min, y0max), (x0,x0min,x0max))**2+ \
                                      Integral(im(self.integrand), (y0, y0min, y0max), (x0,x0min,x0max))**2 ))
        self.Fresnel = Fresnel()
             

#### Rayleigh-Sommerfeld Diffraction Integral
        class Rayleigh_Sommerfeld(branch):
            """
            Sub class for Rayleigh-Sommerfeld Diffraction Integral.
            """
            def __init__(self):
                super().__init__()
                self.class_type = "Rayleigh_Sommerfeld"
                self.name = "Rayleigh-Sommerfeld Diffraction Integral"
                self.integrand  = Uapr*Eli*z*exp(I*k*sqrt((x-x0)**2+(y-y0)**2+z**2))/((x-x0)**2+(y-y0)**2+z**2)
                self.integral   = Eq(_El, 1/(I*l)*Integral(self.integrand, (y0, y0min, y0max), (x0,x0min,x0max)) )
                self.EField     = Eq(_El, self.integral.rhs)
                self.intensity  = Eq(_I, 1/(I*l)*conjugate(1/(I*l))* \
                                    (Integral(re(self.integrand), (y0, y0min, y0max), (x0,x0min,x0max))**2+ \
                                     Integral(im(self.integrand), (y0, y0min, y0max), (x0,x0min,x0max))**2 ))
        self.Rayleigh_Sommerfeld = Rayleigh_Sommerfeld()

        
    @staticmethod
    def __doc__():
        return("Document of optics class.")
        
oopti = optics()