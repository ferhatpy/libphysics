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

exec(open("../src/libreflection.py").read())

class optics(branch):
    """

    """
    _name = "optics"
    class_type = "default"

    
    def define_symbols(self):
        """
        Global symbols, functions.
        a: 
        F: 
        """
        global k,m
        global x,y,z,x0,y0,z0
        global x0min,x0max,y0min,y0max
        global l,t
        global _El, Eli, _I, Uapr
        
        [k,m] = symbols('k m', real=True, positive=True)
        [x,y,z,x0,y0,z0] = symbols('x y z x_0 y_0 z_0', nonzero=True)
        [x0min,x0max,y0min,y0max] = symbols('x_0min x_0max y_0min y_0max', nonzero=True)
        [l,t] = symbols('lambda t', real=True, positive=True)
        
        _El   = Function('E')(x,y,z)
        _I    = Function('I')(x,y,z)
        Eli   = Function('E_i')(x0,y0)
        Uapr  = Function('U_apr')(x0,y0)
    
    def __init__(self):
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
        self.define_symbols()
        
        self.Rayleigh_Sommerfeld_Diff_Int = Eq(_El, 1/(I*l)*Integral(Uapr*Eli*z*exp(I*k*sqrt((x-x0)**2+(y-y0)**2+z**2))/((x-x0)**2+(y-y0)**2+z**2), (y0, y0min, y0max), (x0,x0min,x0max)) )
        self.Rayleigh_Sommerfeld_Diff_intgd = Uapr*Eli*z*exp(I*k*sqrt((x-x0)**2+(y-y0)**2+z**2))/((x-x0)**2+(y-y0)**2+z**2)
        self.Fraunhofer_Diff_Int = Eq(_El, ( exp(I*k*z)/(I*l*z))*(exp(I*k/(2*z)*(x**2 + y**2)))*Integral(Uapr*Eli*exp(-I*2*pi/(l*z)*(x*x0+y*y0)), (y0, y0min, y0max), (x0,x0min,x0max)) )
        self.Fresnel_Diff_IntEx  = Eq(_El, ( exp(I*k*z)/(I*l*z))*(exp(I*k/(2*z)*(x**2 + y**2)))*Integral(Uapr*Eli*exp( I*k/(2*z)*(x0**2+y0**2))*exp(-I*2*pi/(l*z)*(x*x0+y*y0)), (y0, y0min, y0max), (x0,x0min,x0max)) )
        self.Fresnel_Diff_Int    = Eq(_El, ( exp(I*k*z)/(I*l*z))*Integral(Uapr*Eli*exp( I*k/(2*z)*((x-x0)**2+(y-y0)**2)), (y0, y0min, y0max), (x0,x0min,x0max)) )
        self.Fresnel_Diff_intgd  = Uapr*Eli*exp( I*k/(2*z)*((x-x0)**2+(y-y0)**2))
        
        class subformulary:
            """
            Sub formulary class.
            
            Define global symbols in the outer class.
            """
            def __init__(self):
                # List of Moment of Inertia
                self.rect = Lambda(x, Piecewise( (1, abs(x)<=0.5), (0, True) ))
        self.subformulary = subformulary() 
        
        if self.class_type == "Rayleigh_Sommerfeld":
            self.El = Eq(_El, self.Rayleigh_Sommerfeld_Diff_Int.rhs)
            self.Int = Eq(_I, 1/(I*l)*conjugate(1/(I*l))* \
                          (Integral(re(self.Rayleigh_Sommerfeld_Diff_intgd), (y0, y0min, y0max), (x0,x0min,x0max))**2+ \
                           Integral(im(self.Rayleigh_Sommerfeld_Diff_intgd), (y0, y0min, y0max), (x0,x0min,x0max))**2 ))
        
        if self.class_type == "Fraunhofer":
            self.El = Eq(_El, self.Fraunhofer_Diff_Int.rhs)
            self.Int = Eq(_I, self.El.rhs*conjugate(self.El.rhs))
#            self.Int = Eq(_I, abs(self.El.rhs)**2)
            
        if self.class_type == "Fresnel":
            self.El = Eq(_El, self.Fresnel_Diff_Int.rhs)
            self.Int = Eq(_I, exp(I*k*z)/(I*l*z)*conjugate(exp(I*k*z)/(I*l*z))* \
                          (Integral(re(self.Fresnel_Diff_intgd), (y0, y0min, y0max), (x0,x0min,x0max))**2+ \
                           Integral(im(self.Fresnel_Diff_intgd), (y0, y0min, y0max), (x0,x0min,x0max))**2 ))
            
        
    @staticmethod
    def __doc__():
        return("Document of optics class.")
        
oopti = optics()
oopti.__init__()