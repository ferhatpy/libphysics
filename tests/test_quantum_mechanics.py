#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ## test_template.py

"""
test_quantum_mechanics.py

References:
    Leonardo Angelini - Solved Problems in Quantum Mechanics-Springer (2019)
    
sudo apt install pandoc # Converts a *.ipynb file to a *.pdf file via latex.


Example: ostat
============= 
template.class_type = "micro_canonical_discrete_distinguihable"
template.__init__()
template.verbose = False
[mu,B] = symbols('mu B', real=True)
substitutions = {g:1, engF:mu*B*(2*i-3), j:1, n:2}

# 4 lines
commands = ["xreplace", "ostat.Zsp", substitutions]
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
# The following is not compatible with jupyter-notebook.
# for ipath in lstPaths:
#    if os.path.join(os.path.dirname(__file__), ipath) not in sys.path:
#        sys.path.append(os.path.join(os.path.dirname(__file__), ipath))
from libsympy import *
from quantum_mechanics import *
# Execute jupyter-notebook related commands.
#exec(open('libnotebook.py').read())
# print(sys.version)
# print(sys.path)

# ### Settings

### Settings
class sets:
    """
    Setttings class.
        
    Instead of settings class, settings nametuble might be used.
    Settings = namedtuple("Settings", "type dropinf delta")
    sets = Settings(type="symbolic", dropinf=True, delta=0.1)
    """
    global dictflow, test_all
    
    def __init__(self):
        pass

    # File settings
    input_dir  = "input/quantum_mechanics"
    output_dir = "output/quantum_mechanics"
    
    # Plotting settings
    plot_time_scale = {1:"xy", 2:"xz", 3:"yz"}[3]
    
    # Execution settings.
    test_all = {0:False, 1:True}[0]
    dictflow = {100:"get_formulary", 150:"get_subformulary",
                200:"topic1", 300:"topic2",
                400:"topic3"}
    
    
    dictflow = dict(
        ch1 = {1:"p1.3",3:"p1.5",4:"p1.9",5:"p1.17"},
        ch2 = {4:"p2.4",7:"p2.7",9:"p2.9",50:"e2.5",11:"p2.11",12:"p2.12",
               232:"ch2.3.2",60:"e2.6",22:"p2.22",26:"ch2.6",233:"p2.33",41:"p2.41"},
        ch3 = {322:"p3.22", 330:"p3.30"},
        ch4 = {401:"p4.1",402:"e4.1",421:"ch4.2.1",411:"p4.11",4:"p4.12",404:"fig4.4",
               413:"p4.13",7:"p4.14",8:"p4.15",9:"ch4.3.1",10:"ch4.4.1",
               11:"e4.2",12:"p4.27",449:"p4.49",16:"p4.55 todo"},
        ch5 = {1:"p5.1 todo"},
        ch6 = {1:"p6.2"})
    flow = [dictflow["ch2"][i] for i in [9]]
    if test_all: flow = [dictflow[i] for i in dictflow.keys()]

print("Test of the {0}.".format(sets.flow))

# ### get_formulary

#### get_formulary
if "get_formulary" in sets.flow:
#    omech = mechanics() # DO NOT create any instance.
    omec.class_type = ""
#   omech.__init__('EulerLagrange')
    omec.__init__()
    omec.get_formulary()
    omec.get_formulary(style="eq")

# ### get_subformulary

#### get_subformulary    
if "get_subformulary" in sets.flow:
    omec.class_type = ""
    omec.__init__()
    omec.get_subformulary()    

if "p2.9" in sets.flow:
    # --- p2.9
    print("p2.9 <H>=?")
    oqmec.class_type = "position_space"
    oqmec.__init__()
    oqmec.verbose = True
    [A,a,m] = symbols('A a m', real=True, positive=True)
    psi = Wavefunction(A*x*(a-x), (x, 0, a))
    npsi= psi.normalize()
    substitutions = {oqmec.Psi:npsi.expr, xmin:0, xmax:a, oqmec.V:0}
    
    H = oqmec.H.xreplace(substitutions)
    exp_H = oqmec.exp_H.xreplace(substitutions)
    substitutions = {oqmec.Psi:npsi.expr, xmin:0, xmax:a, oqmec.V:0}
    exp_opH = oqmec.exp_opH.xreplace(substitutions)
    pprints("psi=", psi,
            "npsi=", npsi,
            "H=", H, H.doit(),
            "<npsi|H|npsi>=<npsi|p^2>/(2m)|npsi>=", exp_H, exp_H.doit(),
            "<npsi|H|npsi>=<npsi|p^2>/(2m)|npsi>= ERROR", exp_opH, exp_opH.doit(),
            )
    

if "p2.4" in sets.flow:
    # --- p2.4
    print("p.2.4")
    oqmec.class_type = "position_space"
    oqmec.__init__()
    oqmec.verbose = True
    
    psi = {1:Wavefunction(oqmec.subformulary.psi_infinite_qw.rhs, x).expr,
           2:oqmec.subformulary.psi_infinite_qw.rhs}[2]
    substitutions = {oqmec.Psi:psi, xmin:0, xmax:a}
    
    exp_x  = oqmec.exp_x.xreplace(substitutions)
    exp_x2 = oqmec.exp_x2.xreplace(substitutions)
    exp_px = oqmec.exp_px.xreplace(substitutions)
    exp_px2= oqmec.exp_px2.xreplace(substitutions)
    delta_x = oqmec.delta_x.xreplace(substitutions)
    delta_px = oqmec.delta_px.xreplace(substitutions)
    delta_XP = oqmec.delta_XP.xreplace(substitutions)
    min_deltaXP = simplify(delta_XP.doit()).subs({n:1})
    
    pprints("1. Way: SymPy derivative function used in operator definitions.",
            "psi=", psi,
            "psi*=", conjugate(psi),
            "<x>=",     exp_x,  exp_x.doit(),  simplify(exp_x.doit()),
            "<x^2>=",   exp_x2, exp_x2.doit(), simplify(exp_x2.doit()),
            "<p_x>=",   exp_px, exp_px.doit(), simplify(exp_px.doit()),
            "<p_x^2>=", exp_px2, exp_px2.doit(), simplify(exp_px2.doit()),
            "delta x=",  delta_x, delta_x.doit(), simplify(delta_x.doit()), 
            "delta p_x=",delta_px, delta_px.doit(), simplify(delta_px.doit()), 
            "deltaX*deltaP=", delta_XP, delta_XP.doit(), simplify(delta_XP.doit()),
            "At n=1, uncertaninty becomes minimum=", min_deltaXP,
            output_style="display")
    
    exp_x  = oqmec.exp_opX.xreplace(substitutions)
    exp_x2 = oqmec.exp_opX2.xreplace(substitutions)
    exp_px = oqmec.exp_opPx .xreplace(substitutions)
    exp_px2= oqmec.exp_opPx2.xreplace(substitutions)
    delta_x = oqmec.delta_opX.xreplace(substitutions)
    delta_px = oqmec.delta_opPx.xreplace(substitutions)
    delta_XP = oqmec.delta_opXopPx.xreplace(substitutions)
    min_deltaXP = simplify(delta_XP.doit()).subs({n:1})
    
    pprints("2. Way: DifferentialOperator used from sympy.physics.quantum.operator in operator definitions.",
            "psi=", psi,
            "psi*=", conjugate(psi),
            "<x>=",     exp_x,  exp_x.doit(),  simplify(exp_x.doit()),
            "<x^2>=",   exp_x2, exp_x2.doit(), simplify(exp_x2.doit()),
            "<p_x>=",   exp_px, exp_px.doit(), simplify(exp_px.doit()),
            "<p_x^2>=", exp_px2, exp_px2.doit(), simplify(exp_px2.doit()),
            "delta x=",  delta_x, delta_x.doit(), simplify(delta_x.doit()), 
            "delta p_x=",delta_px, delta_px.doit(), simplify(delta_px.doit()), 
            "deltaX*deltaP=", delta_XP, delta_XP.doit(), simplify(delta_XP.doit()),
            "At n=1, uncertaninty becomes minimum=", min_deltaXP,
            output_style="display")

### Chapter 1 The Wave Function
    
# ### p1.17

### p1.17
if "p1.17" in sets.flow:
    print("p.1.17")
    oqmec.class_type = "position_space"
    oqmec.__init__()
    oqmec.verbose = True
    
    [A,a,m] = symbols('A a m', real=True, positive=True)
    f = Piecewise((0, x < -a), (0, x > a), (A*(a**2-x**2), True))
    psi = Wavefunction(f, x)
    npsi = psi.normalize()
    solA = solve(psi.norm-1, A)[0]
    
    substitutions = {oqmec.Psi:npsi.expr, xmin:-a, xmax:a}
    expX_1 = oqmec.exp_x.evalf(subs=substitutions)
    expX_2 = oqmec.exp_x.evalf(subs=substitutions).doit()
    expX2_1 = oqmec.exp_x2.evalf(subs=substitutions)
    expX2_2 = oqmec.exp_x2.evalf(subs=substitutions).doit()
    
    exp_px_1 = oqmec.exp_px.evalf(subs=substitutions)
    exp_px_2 = oqmec.exp_px.evalf(subs=substitutions).doit()
    exp_px2_1 = oqmec.exp_px2.evalf(subs=substitutions)
    exp_px2_2 = oqmec.exp_px2.evalf(subs=substitutions).doit()
    exp_px2_3= oqmec.exp_px2.xreplace(substitutions).doit()
    
    commands = ["xreplace", "oqmec.exp_px2", substitutions]
    exp_px2_4 = oqmec.process(commands).doit()
    
    sigmaX = sqrt(expX2_2.rhs-expX_2.rhs**2)
    sigmaP = sqrt(exp_px2_2.rhs-exp_px_2.rhs**2)
    
    pprints("p1.17",
            "a)",
            "psi=", psi,
            "npsi=", npsi,
            "Normalization constant, A=", solA,
            "b)",
            "<x>=<psi|x|psi>=", oqmec.exp_x, expX_1, expX_2,
            "c)",
            "<p>=", oqmec.exp_px, exp_px_1, exp_px_2,
            "d)",
            "<x^2>=<psi|x^2|psi>=", oqmec.exp_x2, expX2_1, expX2_2,
            "e)",
            "<p^2>=", oqmec.exp_px2, exp_px2_1, exp_px2_3,
            "f)",
            "sigma_X=<x^2>-<x>^2=", sigmaX,
            "g)",
            "sigma_p_x=", sigmaP,
            "h)",
            "sigma_X*sigma_P=", sigmaX*sigmaP,
            output_style="display")

# ### p1.3

### p1.3
if "p1.3" in sets.flow:
    print(sets.flow)
    oqmec.class_type = "position_space"
    oqmec.__init__()
    oqmec.verbose = True
    
    [A,a,l] = symbols('A a lambda', real=True, positive=True)
    psi = Wavefunction(sqrt(A*exp(-l*(x-a)**2)), x)
    npsi = psi.normalize()
    normconst = psi.norm
    
    substitutions = {oqmec.Psi:npsi.expr, xmin:-Inf, xmax:Inf}
    expX_1 = oqmec.exp_x.evalf(subs=substitutions)
    expX_2 = oqmec.exp_x.evalf(subs=substitutions).doit()
    
    expX2_1 = oqmec.exp_x2.evalf(subs=substitutions)
    expX2_2 = oqmec.exp_x2.evalf(subs=substitutions).doit()
    
    deltaX2 = oqmec.delta_x2.evalf(subs=substitutions).doit()
    
    pprints("~p1.3:",
            "Probability distribution",
            "a)",
            "psi=", psi,
            "normalised psi=", npsi.expr,
            "normalization constant=", simplify(normconst),
            "normconst**2=", normconst**2,
            "|A|^2=", solve(normconst**2-1, A**2),
            
            "b)",
            "<x>",   expX_1,  expX_2,
            "<x^2>", expX2_1, expX2_2,
            
            "c)",
            "sigmaX2 = DeltaX2=", deltaX2,
            
            output_style="display")

#    [mu,B] = symbols('mu B', real=True)
#    
#    display("Single particle partition function:", ostat.Zsp)
#    
#    pprints("")
#
#    ### Magic 4 lines
#    commands = ["xreplace", "ostat.Zsp", substitutions]
#    ostat.result = ostat.process(commands).rhs
#    Zsp = simplify(ostat.result.doit())
#    display(Zsp)
#    
#    ### Alternative ways.
#    commands = ["xreplace", "ostat.Zsp", substitutions]
#    ostat.process(commands)
#    Zsp1 = simplify(ostat.result.evalf())   # 1. Way; evalf() can give a result in a different format.
#    Zsp2 = simplify(ostat.result.doit())    # 2. Way.
#    Zsp3 = simplify(ostat.Zsp.evalf(subs=substitutions))        # 3. Way, substituion without commands and execution.
#    Zsp4 = simplify(ostat.Zsp.evalf(subs=substitutions).doit()) # 4. Way, execution without commands.
#    Zsp5 = simplify(ostat.Zsp.xreplace(substitutions).doit())   # 5. Way.
#    display(Zsp1, Zsp2, Zsp3, Zsp4, Zsp5)
#    
#    commands = ["xreplace", "ostat.U", substitutions]
#    ostat.process(commands)
#    U = simplify(ostat.result.doit())
#    display(U)
    
    ### Get generated SymPy codes.
#    print("Codes:\n", *oqmec.get_codes())
