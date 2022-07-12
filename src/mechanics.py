#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mechanics.py
Created on Fri Mar 11 12:48:37 2022

"""
from sympy import*
from sympy.diffgeom import *
from sympy.diffgeom.rn import *
from sympy.diffgeom.rn import R3_r, R3_s
from sympy.vector import CoordSys3D

from libreflection import *
import libphyscon as pc

exec(open("/media/veracrypt1/python/projects/libphysics/src/libreflection.py").read())
#exec(open("../src/libreflection.py").read())

class mechanics(branch):
    """

    """
    _name = "mechanics"
    class_type = {1:"scalar", 2:"vectorial"}[2]
        
    
    def define_symbols(self):
        """
        Global symbols, functions.
        a: 
        F: 
        """
        global a,F,r,t,v,w
        global k,m,M
        global r_t,v_t,a_t,p_t,L_t
        global x,y,z
        global theta, phi
        global C

        C = CoordSys3D('C') # Cartesian coordinate system.
        [a,F,r,t,v,w] = symbols('a F r t v w', real=True)
        [k,m,M]       = symbols('k m M', real=True, positive=True)
        [x,y,z]     = [Function('x')(t), Function('y')(t), Function('z')(t)]
        [theta,phi] = [Function('theta')(t), Function('phi')(t)]
        
        if self.class_type in ["vectorial"]:
            r_t = Function('r')(t)           # Position vector function.
            v_t = Function('v')(t)           # Velocity vector function.
            a_t = Function('a')(t)           # Acceleration vector function.
            p_t = Function('p')(t)           # Linear momentum vector function.
            L_t = Function('L')(t)           # Angular momentum vector function.
        
    def __init__(self):
        super().__init__()
        self.define_symbols()
        
        class subformulary:
            """
            Sub formulary class.
            
            Define global symbols in the outer class.
            """
            def __init__(self):
                # Transformations from Polar to Cartesian Coordinates.
                self.pol_to_cart_x = r_t*cos(theta)*C.i
                self.pol_to_cart_y = r_t*sin(theta)*C.j
                # Transformations from Cylindrical to Cartesian Coordinates.
                self.cyl_to_cart_x = r_t*sin(phi)*C.i
                self.cyl_to_cart_y = r_t*sin(phi)*C.j
                self.cyl_to_cart_z = z*C.k
                # Transformations from Spherical to Cartesian Coordinates.
                self.sph_to_cart_x = r_t*sin(theta)*cos(phi)*C.i
                self.sph_to_cart_y = r_t*sin(theta)*sin(phi)*C.j
                self.sph_to_cart_z = r_t*cos(theta)*C.k
                
                # List of Moment of Inertia
                self.Icm_sphere = S(2)/5*M*r**2
        self.subformulary = subformulary()
        
        if self.class_type == "scalar":
            self.a = a
            self.NewtonsLaw2 = Eq(F, m*a)
            self.HookesLaw   = Eq(F, -k*x)
            
        if self.class_type == "vectorial":
            """
            Example:
            omech.r -> Eq(r(t), x(t) + y(t) + z(t))
            omech.r.rhs -> x(t) + y(t) + z(t)
            omech.r.rhs.subs({x:_r*cos(theta)}) -> r(t)*cos(theta(t)) + y(t) + z(t)
            diff(omech.r.rhs.subs({x:_r*cos(theta)}),t) -> -r(t)*sin(theta(t))*Derivative(theta(t), t) + cos(theta(t))*Derivative(r(t), t) + Derivative(y(t), t) + Derivative(z(t), t)
            """
            [self.x, self.y, self.z] = [x, y, z]
            self.r = Eq( r_t, self.x + self.y + self.z)
            self.v = Eq( v_t, diff(self.r.rhs, t, evaluate=False) )
            self.a = Eq( a_t, diff(self.v.rhs, t, evaluate=False) )
#            self.p = Eq( p_t, m*self.v )
#            self.L = Eq( L_t, (self.r.rhs.args[0]*C.i+
#                               self.r.rhs.args[1]*C.j+
#                               self.r.rhs.args[2]*C.k).cross(p_t) )
            self.NewtonsLaw2 = Eq(F, m*self.a.rhs)
            self.HookesLaw   = Eq(F, -k*self.x)
        
        if self.class_type == "EulerLagrange":
            self.NewtonsLaw2 = Eq(F, m*a)
            self.HookesLaw   = Eq(F, -k*x)

        
    @staticmethod
    def __doc__():
        return("Document of mechanics class.")
        
omech = mechanics()
omech.__init__()