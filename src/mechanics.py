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
from sympy.physics.vector import *
from sympy.vector import CoordSys3D

from libreflection import *
import libphyscon as pc

#exec(open("/media/veracrypt1/python/projects/libphysics/src/libreflection.py").read())
exec(open("../src/libreflection.py").read())

class mechanics(branch):
    """

    """
    _name = "mechanics"
    class_type = {1:"scalar", 2:"vectorial"}[1]
        
    
    def define_symbols(self):
        """
        Global symbols, functions.
        a: 
        F: 
        """
        global a,r,t,v,w
        global F,Fx,Fy,Fz
        global k,m,M
        global r,v,a,p,L
        global _r,_v,_a,_F,_p,_L
        global x,y,z
        global theta, phi
        global C

        C           = CoordSys3D('C') # Cartesian coordinate system.
        t           = symbols('t', real=True)
        [k,m,M,w]   = symbols('k m M w', real=True, positive=True)
        [x,y,z]     = [Function('x')(t), Function('y')(t), Function('z')(t)]
        [theta,phi] = [Function('theta')(t), Function('phi')(t)]
       
        if self.class_type in ["scalar"]:
           [r,v,_a,F,p] = symbols('r v a F p', real=True) 
        
        if self.class_type in ["vectorial"]:
            [r,] = symbols('r,', real=True) 
            [rx,ry,rz] = symbols('r_x r_y r_z', real=True)
            [vx,vy,vz] = symbols('v_x v_y v_z', real=True)
            [ax,ay,az] = symbols('a_x a_y a_z', real=True)
            [Fx,Fy,Fz] = symbols('F_x F_y F_z', real=True)
            [px,py,pz] = symbols('p_x p_y p_z', real=True)
            [Lx,Ly,Lz] = symbols('L_x L_y L_z', real=True)
            _r = rx*C.i+ry*C.j+rz*C.k   # Function('r')(t)    # Position vector function.
            _v = vx*C.i+vy*C.j+vz*C.k   # Velocity vector function.
            _a = ax*C.i+ay*C.j+az*C.k   # Acceleration vector function.
            _F = Fx*C.i+Fy*C.j+Fz*C.k   # Force vector function.
            _p = px*C.i+py*C.j+pz*C.k   # Linear momentum vector function.
            _L = Lx*C.i+Ly*C.j+Lz*C.k   # Angular momentum vector function.
        
        
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
                self.pol_to_cart_x = r*cos(theta)
                self.pol_to_cart_y = r*sin(theta)
                # Transformations from Cylindrical to Cartesian Coordinates.
                self.cyl_to_cart_x = r*sin(phi)
                self.cyl_to_cart_y = r*sin(phi)
                self.cyl_to_cart_z = z
                # Transformations from Spherical to Cartesian Coordinates.
                self.sph_to_cart_x = r*sin(theta)*cos(phi)
                self.sph_to_cart_y = r*sin(theta)*sin(phi)
                self.sph_to_cart_z = r*cos(theta)
                
                # List of Moment of Inertia
                self.Icm_sphere = S(2)/5*M*r**2
        self.subformulary = subformulary()
        
        if self.class_type in ["scalar"]:
            self.a           = Eq(_a, diff(x, t, 2, evaluate=False))
            self.momentum    = Eq(p, m*v)
            self.NewtonsLaw2 = Eq(F, m*self.a.rhs)
            self.HookesLaw   = Eq(F, -k*x)
            
        if self.class_type in ["vectorial"]:
            """
            Example:
            omech.r -> Eq(r(t), x(t) + y(t) + z(t))
            omech.r.rhs -> x(t) + y(t) + z(t)
            omech.r.rhs.subs({x:_r*cos(theta)}) -> r(t)*cos(theta(t)) + y(t) + z(t)
            diff(omech.r.rhs.subs({x:_r*cos(theta)}),t) -> -r(t)*sin(theta(t))*Derivative(theta(t), t) + cos(theta(t))*Derivative(r(t), t) + Derivative(y(t), t) + Derivative(z(t), t)
            """
            [self.x, self.y, self.z] = [x, y, z]
            self.r = Eq(_r, self.x*C.i + self.y*C.j + self.z*C.k)
            self.v = Eq(_v, diff(self.r.rhs, t, evaluate=False))
            self.a = Eq(_a, diff(self.v.rhs, t, evaluate=False))
            self.p = Eq(_p, m*self.v.rhs)
            self.L = Eq(_L, self.r.rhs.cross(self.p.rhs))
            
            self.NewtonsLaw2 = Eq(_F, m*self.a.rhs)
            self.HookesLaw   = Eq(_F, -k*self.x)
        
        if self.class_type == "EulerLagrange":
            self.NewtonsLaw2 = Eq(F, m*a)
            self.HookesLaw   = Eq(F, -k*x)

        
    @staticmethod
    def __doc__():
        return("Document of mechanics class.")
        
omech = mechanics() # Create an omech object from mechanics class.
#omech.__init__()