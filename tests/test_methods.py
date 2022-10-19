#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_method.py

Find and replace template with desired class name.

Example: ostat
============= 
template.class_type = "micro_canonical_discrete_distinguihable"
template.__init__()
template.verbose = False
[mu,B] = symbols('mu B', real=True)
xreplaces = {g:1, engF:mu*B*(2*i-3), j:1, n:2}

# 4 lines
commands = ["xreplace", "ostat.Zsp", xreplaces]
template.process(commands)
Zsp = simplify(template.result.doit())
display(Zsp)

# Deep copy to create new instances of obranch class.
import copy
ostat2 = copy.deepcopy(ostat)
"""
import copy
import sys
import os
lstPaths = ["../src"]
for ipath in lstPaths:
    if ipath not in sys.path:
        sys.path.append(ipath)
from sympy import*
from sympy.vector import CoordSys3D
from libsympy import *
from methods import *


# print(sys.version)
# print(sys.path)

global R
R = CoordSys3D('R')


### Settings
class sets:
    """
    Setttings class.
        
    Instead of settings class, settings nametuble might be used.
    Settings = namedtuple("Settings", "type dropinf delta")
    sets = Settings(type="symbolic", dropinf=True, delta=0.1)
    """
    def __init__(self):
        pass
    
    input_dir  = "input/optics"
    output_dir = "output/optics"
    
    # Plotting settings
    plot_time_scale = {1:"xy", 2:"xz", 3:"yz"}[3]
    
    flow = [{100:"get_formulary", 150:"get_subformulary",
             200:"TripleVectorProduct", 300:"topic2",
             400:"topic3"}[i] 
            for i in [200]]

### Formulary
print("Test of the {0}.".format(sets.flow[0]))
if "get_formulary" in sets.flow:
#    omech = mechanics() # DO NOT create any instance.
    ometh.class_type = ""
    ometh.__init__()
    ometh.get_formulary()
    ometh.get_formulary(style="eq")
    
if "get_subformulary" in sets.flow:
    ometh.class_type = ""
    ometh.__init__()
    ometh.get_subformulary()    


### A Spin-1/2 Paramagnet
if "topic1" in sets.flow:
    print("A Spin-1/2 Paramagnet")
    
    ostat.class_type = "micro_canonical_discrete_distinguihable"
    ostat.__init__()
    ostat.verbose = False
    [mu,B] = symbols('mu B', real=True)
    xreplaces = {g:1, engF:mu*B*(2*i-3), j:1, n:2}
    display("Single particle partition function:", ostat.Zsp)
    """
    res = ostat.Zsp.xreplace(xreplaces)
    display(simplify(res.doit()))
    """
    ### 4 lines
    commands = ["xreplace", "ostat.Zsp", xreplaces]
    ostat.process(commands)
    Zsp1 = simplify(ostat.result.evalf()) # evalf() can give a result in a different format.
    Zsp2 = simplify(ostat.result.doit())
    Zsp3 = simplify(ostat.Zsp.evalf(subs=xreplaces).doit())
    display(Zsp1, Zsp2, Zsp3)
    
    commands = ["xreplace", "ostat.U", xreplaces]
    ostat.process(commands)
    U = simplify(ostat.result.doit())
    display(U)
    
    ### Get generated SymPy codes.
    print("Codes:\n", *ostat.get_codes())


elif "TripleVectorProduct" in sets.flow:
    """
    # (A x B) x C =! A x (B x C)
    [Ax,Bx,Cx] = symbols('A_x B_x C_x', real=True)
    [Ay,By,Cy] = symbols('A_y B_y C_y', real=True)
    [Az,Bz,Cz] = symbols('A_z B_z C_z', real=True)
    
    C = CoordSys3D('C')
    A = Ax*C.i + Ay*C.j + Az*C.k
    B = Bx*C.i + By*C.j + Bz*C.k 
    C = Cx*C.i + Cy*C.j + Cz*C.k
    lhs = A.cross(B).cross(C)
    rhs = A.cross(B.cross(C))
    print(simplify(lhs-rhs))
    print(lhs.equals(rhs))
    """
    ometh.__init__()
    ometh.verbose = True
    
    # 1. way
    vA,vB,vC = (ometh.vA, ometh.vB, ometh.vC)
    lhs = vA.cross(vB).cross(vC)
    rhs = vA.cross(vB.cross(vC))
    pprints("(A x B) x C=", lhs,
            "A x (B x C)=", rhs,
            "lhs.equals(rhs)", lhs.equals(rhs),
            "simplify(lhs-rhs)", simplify(lhs-rhs),
            output_style = {1:"display", 2:"pprint", 3:"print", 4:"latex"}[1],
            newline=True)
    
    # 2. Way
    vA = 3*C.i + 3*C.j + 3*C.k
    vB = 3*C.i + 3*C.j + 3*C.k
    xreplaces = {Ax:vA.components[C.i], Ay:vA.components[C.j], Az:vA.components[C.k],
                 Bx:vB.components[C.i], By:vB.components[C.j], Bz:vB.components[C.k],
                 Cx:vC.components[C.i], Cy:vC.components[C.j], Cz:vC.components[C.k]}
    commands = ["xreplace", "ometh.DotProduct", xreplaces]
    ometh.process(commands)
    
    # 3. Way
    ometh.DotProduct.evalf(subs=xreplaces)
    # 4. Way
    display(f"A= {vA}", f"B= {vB}",
            "Dot Product",
            ometh.DotProduct,
            ometh.DotProduct.xreplace(xreplaces),
            "Vector Product",
            ometh.CrossProduct.xreplace(xreplaces)
            )
    