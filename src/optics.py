#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
optics.py
Created on Fri Mar 11 12:53:36 2022

"""
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


#### Fraunhofer Diffraction
        class Fraunhofer(branch):
            """
            Sub class for Fraunhofer Diffraction
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
            Sub class for Fresnel Diffraction
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
             
        
    @staticmethod
    def __doc__():
        return("Document of optics class.")
        
oopti = optics()