#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ## test_nuclear_physics.py

"""
test_nuclear_physics.py connected to test_nuclear_physics.ipynb via "jupytext".
In ipynb notebook select File->Jupytext->Pair Notebook with Light Format.

References:
===========    
    Books:
    ======    
    Kenneth S. Krane - Introductory Nuclear Physics-John Wiley (1988)
    
    Problem Books:
    ==============    

    
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
import os
import sys
# Import path for library functions.
lstPaths = ["../src"]
for ipath in lstPaths:
    if ipath not in sys.path:
        sys.path.append(ipath)
# The following is not compatible with jupyter-notebook.
# for ipath in lstPaths:
#    if os.path.join(os.path.dirname(__file__), ipath) not in sys.path:
#        sys.path.append(os.path.join(os.path.dirname(__file__), ipath))
from libsympy import *
from sympy.abc import*
from quantum_mechanics import *
from nuclear_physics import *

import scipy.constants as pc
import scienceplots
plt.style.use(['science', 'notebook'])

# Execute jupyter-notebook related commands.
# exec(open('libnotebook.py').read())
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
    test_all = {0:False, 1:True}[1]
    dictflow = dict(
        ch1 = {1:"p1.1",3:"p1.5",4:"p1.9",5:"p1.17"},
        ch2 = {24:"p2.4",27:"p2.7",29:"p2.9",25:"e2.5",211:"p2.11",212:"p2.12",
               232:"ch2.3.2",26:"e2.6",222:"p2.22",26:"ch2.6",233:"p2.33",241:"p2.41"},
        ch3 = {322:"p3.22", 330:"p3.30"},
        ch4 = {401:"p4.1",402:"e4.1",421:"ch4.2.1",411:"p4.11",4:"p4.12",404:"fig4.4",
               413:"p4.13",7:"p4.14",8:"p4.15",9:"ch4.3.1",10:"ch4.4.1",
               11:"e4.2",12:"p4.27",449:"p4.49",16:"p4.55"},
        ch5 = {1:"p5.1"},
        ch6 = {61:"p6.1", 62:"p6.2", 611:"p6.11", 614:"p6.14", 615:"p6.15", 
               253:"c25.3"},
        ch9 = {91:"p9.1"})
    flow = [dictflow["ch9"][i] for i in [91]]
    if test_all: flow = flatten([list(dictflow[i].values()) for i in dictflow.keys()])

print("Test of the {0}.".format(sets.flow))

# ### get_formulary

#### get_formulary
if "get_formulary" in sets.flow:
    onucl.__init__()
    onucl.get_formulary()
    onucl.get_formulary(style="eq")

# ### get_subformulary

#### get_subformulary    
if "get_subformulary" in sets.flow:
    onucl.class_type = ""
    onucl.__init__()
    onucl.get_subformulary()

# ### p9.1

if "p9.1" in sets.flow:
    print("ENERGY RELEASE IN beta DECAY")
    onucl.__init__()
    onucl.verbose = True
    
    substitutions = {e:m_e, n:m_n, p:m_p, nu_e:m_nu_e, anu_e:m_nubar_e}
    react = onucl.Q(onucl.beta_negative_decay, substitutions)
    nums = {c:1, 
            m_e:pc.physical_constants["electron mass energy equivalent in MeV"][0],
            m_n:pc.physical_constants["neutron mass energy equivalent in MeV"][0],
            m_p:pc.physical_constants["proton mass energy equivalent in MeV"][0],
            m_nubar_e:0.013}
    display(react.subs(nums))


