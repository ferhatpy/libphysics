#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_template.py

Find and replace template with desired class name.

Example: ostat
============= 
template.class_type = "micro_canonical_discrete_distinguihable"
template.__init__()
template.solver.verbose = False
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
lstPaths = ["../src", "../../libpython/src"]
for ipath in lstPaths:
    if ipath not in sys.path:
        sys.path.append(ipath)
from libsympy import *
from mechanics import *
from statistical_physics import *


# print(sys.version)
# print(sys.path)

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
             200:"topic1", 300:"topic2",
             400:"topic3"}[i] 
            for i in [200]]

### Formulary
print("Test of the {0}.".format(sets.flow[0]))
if "get_formulary" in sets.flow:
#    omech = mechanics() # DO NOT create any instance.
    omec.class_type = ""
    omec.__init__()
    omec.get_formulary()
    omec.get_formulary(style="eq")
    
if "get_subformulary" in sets.flow:
    omec.class_type = ""
    omec.__init__()
    omec.get_subformulary()    


### A Spin-1/2 Paramagnet
if "topic1" in sets.flow:
    print("A Spin-1/2 Paramagnet")
    
    ostat.class_type = "micro_canonical_discrete_distinguihable"
    ostat.__init__()
    ostat.solver.verbose = False
    
    [mu,B] = symbols('mu B', real=True)
    xreplaces = {g:1, engF:mu*B*(2*i-3), j:1, n:2}
    display("Single particle partition function:", ostat.Zsp)

    ### Magic 4 lines
    commands = ["xreplace", "ostat.Zsp", xreplaces]
    ostat.result = ostat.process(commands).rhs
    Zsp = simplify(ostat.result.doit())
    display(Zsp)
    
    ### Alternative ways.
    commands = ["xreplace", "ostat.Zsp", xreplaces]
    ostat.process(commands)
    Zsp1 = simplify(ostat.result.evalf())   # 1. Way; evalf() can give a result in a different format.
    Zsp2 = simplify(ostat.result.doit())    # 2. Way.
    Zsp3 = simplify(ostat.Zsp.evalf(subs=xreplaces))        # 3. Way, substituion without commands and execution.
    Zsp4 = simplify(ostat.Zsp.evalf(subs=xreplaces).doit()) # 4. Way, execution without commands.
    Zsp5 = simplify(ostat.Zsp.xreplace(xreplaces).doit())   # 5. Way.
    display(Zsp1, Zsp2, Zsp3, Zsp4, Zsp5)
    
    commands = ["xreplace", "ostat.U", xreplaces]
    ostat.process(commands)
    U = simplify(ostat.result.doit())
    display(U)
    
    ### Get generated SymPy codes.
    print("Codes:\n", *ostat.solver.get_codes())