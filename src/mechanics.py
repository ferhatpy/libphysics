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
from sympy.plotting import plot_parametric
from sympy.vector import CoordSys3D
from sympy.physics.quantum.operator import *
from sympy.physics.quantum.qapply import *

from libreflection import *
import libphyscon as pc

exec(open("../src/libreflection.py").read())

class mechanics(branch):
    """

    """
    class_type = {1:"scalar", 2:"vectorial"}[1]
    
    def define_symbols(self):
        """
        Common global symbols, functions.
        a: 
        F: 
        """
        global C
        global C1,C2,C3
        global t,tau,p,s
        global alpha,beta,gamma
        global A,k,m,M,F0
        global G
        global a,r,x,v,w,w0,w1,w2,y,z
        global _a,_r,_x,_v
        global F,Fx,Fy,Fz,W
        global _F,_L,_p,_T,_Tr,_H,_W,_U
        global dr
        global xi,xf,yi,yf,zi,zf
        global x0,y0,z0,v0
        global theta, phi
        

        C         = CoordSys3D('C') # Cartesian coordinate system.
        C1,C2,C3  = symbols('C1 C2 C3')
        t,p,s     = symbols('t p s', real=True)
        t         = Symbol('t', real=True, positive=True)
        tau       = Symbol('tau', real=True, positive=True)
        alpha,beta,gamma = symbols('alpha beta gamma', real=True, positive=True)
        A,k,m,M,F0,w,w0,w1,w2 = symbols('A k m M F_0 w w_0 w_1 w_2', real=True, positive=True)
        xi,xf,yi,yf,zi,zf = symbols('x_i x_f y_i y_f z_i z_f', real=True)
        x0,y0,z0,v0 = symbols('x_0 y_0 z_0 v_0', real=True)
        G         = Function('G')(t)
        r         = Function('r')(t)
        x,y,z     = [Function('x')(t), Function('y')(t), Function('z')(t)]
        theta,phi = [Function('theta')(t), Function('phi')(t)]
       
        if self.class_type in ["scalar"]:
            _F,_W  = symbols('F W', real=True)
            _x = Function('x')(t)          # Position.
            _v = Function('v')(t)          # Velocity.
            _a = Function('a')(t)          # Acceleration.
            _p = Function('p')(t)          # Linear momentum.
            _H = Function('H')(t)          # Total energy.
            _T = Function('T')(t)          # Kinetic energy.
            _U = Function('U')(t)          # Potential energy.
            
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
            dr = 1*C.i + 1*C.j + 1*C.k      # Differential position element.
            _r = rx*C.i+ry*C.j+rz*C.k       # Position vector.
            _v = vx*C.i+vy*C.j+vz*C.k       # Velocity vector.
            _a = ax*C.i+ay*C.j+az*C.k       # Acceleration vector.
            _F = Fx*C.i+Fy*C.j+Fz*C.k       # Force vector.
            _H = Function('H')(t)           # Total energy.
            _T = Function('T')(t)           # Kinetic energy.
            _U = Function('U')(t)           # Potential energy.
            _Tr= Trx*C.i+Try*C.j+Trz*C.k    # Torque vector.
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
                # Harmonic oscillator substitutions.
                self.A = self.reduced_amplitude = F0/m
                self.beta = self.damping_parameter = gamma/(2*m)
                self.w0 = self.natural_frequency = sqrt(k/m)
                self.underdamping_criteria = {w0:sqrt(beta**2+w1**2)}
                self.critical_damping_criteria = {w0:beta}
                self.overdamping_criteria = {w0:sqrt(beta**2-w2**2)}
                
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
            self.G = self.GreensF = G # Green's function.
            self.x = self.position = x
            self.v = self.velocity = Eq(_v, diff(self.x, t, evaluate=True))
            self.a = self.acceleration = Eq(_a, diff(self.x, t, 2, evaluate=True))
            self.p = self.momentum = Eq(_p, m*self.v.rhs)
            self.F = self.NewtonsLaw2 = Eq(_F, m*self.a.rhs)
            self.HookesLaw = Eq(_F, -k*self.x)
            self.W = self.work = Eq(_W, Integral(self.F.rhs, (x,xi,xf)))
            self.T = self.kinetic_energy = Eq(_T, S(1)/2*m*self.v.rhs**2)
            self.U = self.potential_energy = Eq(_U, Integral(k*self.x, (x,0,self.x)))
            self.H = self.energy = Eq(_H, self.T.rhs + self.U.rhs)
            self.damped_harmonic_oscillator1 = Eq(m*diff(self.x, t, 2, evaluate=True)+gamma*diff(self.x, t, evaluate=True)+k*self.x, 0)
            self.damped_harmonic_oscillator2 = Eq(diff(self.x, t, 2, evaluate=True)+2*beta*diff(self.x, t, evaluate=True)+w0**2*self.x, 0)
            self.driven_oscillator1 = Eq(m*diff(self.x, t, 2, evaluate=True)+gamma*diff(self.x, t, evaluate=True)+k*self.x, F0*cos(w*t))
            self.driven_oscillator2 = Eq(diff(self.x, t, 2, evaluate=True)+2*beta*diff(self.x, t, evaluate=True)+w0**2*self.x, A*cos(w*t))
            self.driven_oscillator2_GreensF = Eq(diff(self.G, t, 2) + 2*beta*diff(self.G, t) + w0**2*self.G, 0)
            self.amplitude = None
            self.phase = None
            self.scaled_amplitude = None
            
        if self.class_type in ["vectorial"]:
            """
            Example:
            omech.r -> Eq(r(t), x(t) + y(t) + z(t))
            omech.r.rhs -> x(t) + y(t) + z(t)
            omech.r.rhs.subs({x:_r*cos(theta)}) -> r(t)*cos(theta(t)) + y(t) + z(t)
            diff(omech.r.rhs.subs({x:_r*cos(theta)}),t) -> -r(t)*sin(theta(t))*Derivative(theta(t), t) + cos(theta(t))*Derivative(r(t), t) + Derivative(y(t), t) + Derivative(z(t), t)
            """
            self.x, self.y, self.z = [x, y, z]
            self.r = self.position = Eq(_r, self.x*C.i + self.y*C.j + self.z*C.k)
            self.v = self.velocity = Eq(_v, diff(self.r.rhs, t, evaluate=False))
            self.a = self.acceleration = Eq(_a, diff(self.v.rhs, t, evaluate=False))
            self.p = self.momentum = Eq(_p, m*self.v.rhs)
            self.L = self.angular_momentum = Eq(_L, self.r.rhs.cross(self.p.rhs))
            self.F = self.NewtonsLaw2 = Eq(_F, m*self.a.rhs)
            self.Tr1 = self.Torque1 = Eq(_Tr, self.r.rhs.cross(diff(self.p.rhs, t, evaluate=False))) # Torque = r x dp/dt
            self.Tr2 = self.Torque1 = Eq(_Tr, self.r.rhs.cross(self.F.rhs)) # Torque = r x F
            self.HookesLaw   = Eq(_F, -k*self.x*C.i)
            self.W = self.work = Eq(_W, Integral(self.F.rhs.dot(dr), (z,zi,zf), (y,yi,yf), (x,xi,xf)))
            self.T = self.kinetic_energy  = Eq(_T, S(1)/2*m*self.v.rhs.dot(self.v.rhs))
            self.U = self.potential_energy= Eq(_U, Integral(k*self.r.rhs.dot(dr), (z,zi,zf), (y,yi,yf), (x,xi,xf)))
            self.H = self.energy = Eq(_H, self.T.rhs + self.U.rhs)
        
        if self.class_type == "EulerLagrange":
            """
            todo
            """
            self.x, self.y, self.z = [x, y, z]
            self.r = self.position = Eq(_r, self.x*C.i + self.y*C.j + self.z*C.k)
            self.v = self.velocity = Eq(_v, diff(self.r.rhs, t, evaluate=False))
            self.a = self.acceleration = Eq(_a, diff(self.v.rhs, t, evaluate=False))
            self.F = self.NewtonsLaw2 = Eq(_F, m*self.a.rhs)
            self.HookesLaw   = Eq(_F, -k*self.x)
        
    @staticmethod
    def __doc__():
        return("Document of mechanics class.")
        
omech = mechanics() # Create an omech object from mechanics class.
#omech.__init__()