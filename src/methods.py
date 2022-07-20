#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
template.py
Created on Fri Mar 11 12:53:36 2022

Find and replace template with desired class name.

SymPy 1.10.1 documentation » Reference Documentation » Matrices » Vector

"""
from sympy import*
from sympy.diffgeom import *
from sympy.diffgeom.rn import *
from sympy.diffgeom.rn import R3_r, R3_s
from sympy.vector import CoordSys3D

from libreflection import *
import libphyscon as pc

exec(open("../src/libreflection.py").read())

class methods(branch):
    """

    """
    _name = "methods"
    class_type = "default"

    
    def define_symbols(self):
        """
        Global symbols, functions.
        a: 
        F: 
        """
        global alpha,beta,gamma,phi,theta
        global a,b,c,d
        global k,m,t,w
        global vA,vB,vC,vD
        global Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz
        global x,y,z,r,t,w
        global C
        
        alpha,beta,gamma,phi,theta = symbols('alpha beta gamma phi theta', real=True)
        a,b,c,d   = symbols('a b c d', real=True)
        k,m,t,w   = symbols('k m t w', real=True, positive=True)
        vA,vB,vC,vD = symbols('A B C D', vector=True)
        Ax,Bx,Cx,Dx = symbols('A_x B_x C_x D_x', real=True)
        Ay,By,Cy,Dy = symbols('A_y B_y C_y D_y', real=True)
        Az,Bz,Cz,Dz = symbols('A_z B_z C_z D_z', real=True)
        x,y,z,r,t,w = symbols('x y z r t w', real=True)
        C = CoordSys3D('C')
    
    def __init__(self):
        super().__init__()
        self.define_symbols()
        
        self.vA = Ax*C.i + Ay*C.j + Az*C.k
        self.vB = Bx*C.i + By*C.j + Bz*C.k
        self.vC = Cx*C.i + Cy*C.j + Cz*C.k
        
        class subformulary:
            """
            Sub formulary class.
            
            Define global symbols in the outer class.
            """
            def __init__(self):
                pass
    
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
            self.DotProduct     = Eq(c, self.vA.dot(self.vB))
            self.DotProductTheta= Eq(c, self.vA.magnitude()*self.vB.magnitude()*cos(theta))
            self.CrossProduct   = Eq(self.vC, self.vA.cross(self.vB))
            self.CrossProductTheta = Eq(c, self.vA.magnitude()*self.vB.magnitude()*sin(theta))
            self.TripleProduct  = Eq(c, self.vA.dot(self.vB.cross(self.vC)))

        
    @staticmethod
    def __doc__():
        return("Document of methods class.")
        
ometh = methods()