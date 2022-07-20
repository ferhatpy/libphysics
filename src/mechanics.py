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

exec(open("/media/veracrypt1/python/projects/libphysics/src/libreflection.py").read())
#exec(open("../src/libreflection.py").read())

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
        global k,m,M
        global a,r,t,x,v,w,y,z
        global _a,_r,_x,_v
        global F,Fx,Fy,Fz,p,W
        global _F,_L,_p,_Tr,_W
        global dr
        global xi,xf,yi,yf,zi,zf
        global theta, phi
        global C

        C         = CoordSys3D('C') # Cartesian coordinate system.
        t         = symbols('t', real=True)
        k,m,M,w = symbols('k m M w', real=True, positive=True)
        xi,xf,yi,yf,zi,zf = symbols('x_i x_f y_i y_f z_i z_f', real=True)
        r         = Function('r')(t)
        x,y,z     = [Function('x')(t), Function('y')(t), Function('z')(t)]
        theta,phi = [Function('theta')(t), Function('phi')(t)]
       
        if self.class_type in ["scalar"]:
            _F,_W  = symbols('F W', real=True)
            _x  = Function('x')(t) 
            _v  = Function('v')(t)          # Velocity.
            _a  = Function('a')(t)          # Acceleration.
            _p  = Function('p')(t)
            
        if self.class_type in ["vectorial"]:
            _W, = symbols('W,', real=True)
#            rx   = Function('r_x', real=True)(t) # Possible time dependent definition.
            rx,ry,rz = symbols('r_x r_y r_z', real=True) # Components of position.
            vx,vy,vz = symbols('v_x v_y v_z', real=True) # Components of velocity.
            ax,ay,az = symbols('a_x a_y a_z', real=True) # Components of acceleration.
            Fx,Fy,Fz = symbols('F_x F_y F_z', real=True) # Components of force.
            Trx,Try,Trz = symbols('tau_x tau_y tau_z', real=True) # Components of torque.
            px,py,pz = symbols('p_x p_y p_z', real=True) # Components of linear momentum.
            Lx,Ly,Lz = symbols('L_x L_y L_z', real=True) # Components of angular momentum.
            dr = 1*C.i + 1*C.j + 1*C.k
            _r = rx*C.i+ry*C.j+rz*C.k       # Position vector.
            _v = vx*C.i+vy*C.j+vz*C.k       # Velocity vector.
            _a = ax*C.i+ay*C.j+az*C.k       # Acceleration vector.
            _F = Fx*C.i+Fy*C.j+Fz*C.k       # Force vector.
            _Tr = Trx*C.i+Try*C.j+Trz*C.k   # Torque vector.
            _p = px*C.i+py*C.j+pz*C.k       # Linear momentum vector.
            _L = Lx*C.i+Ly*C.j+Lz*C.k       # Angular momentum vector.

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
            self.x = x
            self.v = Eq(_v, diff(self.x, t,    evaluate=True))
            self.a = Eq(_a, diff(self.x, t, 2, evaluate=True))
            self.p = Eq(_p, m*self.v.rhs)
            self.F = self.NewtonsLaw2 = Eq(_F, m*self.a.rhs)
            self.HookesLaw = Eq(_F, -k*self.x)
            self.W = Eq(_W, Integral(self.F, (x,xi,xf)))
            
        if self.class_type in ["vectorial"]:
            """
            Example:
            omech.r -> Eq(r(t), x(t) + y(t) + z(t))
            omech.r.rhs -> x(t) + y(t) + z(t)
            omech.r.rhs.subs({x:_r*cos(theta)}) -> r(t)*cos(theta(t)) + y(t) + z(t)
            diff(omech.r.rhs.subs({x:_r*cos(theta)}),t) -> -r(t)*sin(theta(t))*Derivative(theta(t), t) + cos(theta(t))*Derivative(r(t), t) + Derivative(y(t), t) + Derivative(z(t), t)
            """
            self.x, self.y, self.z = [x, y, z]
            self.r = Eq(_r, self.x*C.i + self.y*C.j + self.z*C.k)
            self.v = Eq(_v, diff(self.r.rhs, t, evaluate=False))
            self.a = Eq(_a, diff(self.v.rhs, t, evaluate=False))
            self.p = Eq(_p, m*self.v.rhs)
            self.L = Eq(_L, self.r.rhs.cross(self.p.rhs))
            self.F = Eq(_F, m*self.a.rhs)
            self.NewtonsLaw2 = Eq(_F, m*self.a.rhs)
            self.Tr1 = Eq(_Tr, self.r.rhs.cross(diff(self.p.rhs, t, evaluate=False)))
            self.Tr2 = Eq(_Tr, self.r.rhs.cross(self.F.rhs))
            self.HookesLaw   = Eq(_F, -k*self.x)
            self.W   = Eq(_W, Integral(self.F.rhs.dot(dr), (z,zi,zf), (y,yi,yf), (x,xi,xf)))
        
        if self.class_type == "EulerLagrange":
            self.NewtonsLaw2 = Eq(F, m*a)
            self.HookesLaw   = Eq(F, -k*x)

        
    @staticmethod
    def __doc__():
        return("Document of mechanics class.")
        
omech = mechanics() # Create an omech object from mechanics class.
omech.__init__()