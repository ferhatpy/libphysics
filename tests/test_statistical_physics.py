#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ## test_statistical_physics.py

"""
test_statistical_physics.py

omec.class_type = ""
omec.__init__()
omec.verbose = True
commands = ["solve", "NewtonsLaw2", a]
print(omec.process(commands))
"""
import copy
import sys
import os
lstPaths = ["../src"]
for ipath in lstPaths:
    if ipath not in sys.path:
        sys.path.append(ipath)
from libsympy import *
from statistical_physics import *
# print(sys.version)
# print(sys.path)
# Execute jupyter-notebook related commands.
#exec(open('libnotebook.py').read())


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
    input_dir  = "input/statistical_physics"
    output_dir = "output/statistical_physics"
    
    # Plotting settings
    plot_time_scale = {1:"xy", 2:"xz", 3:"yz"}[3]
    
    # Execution settings.
    test_all = {0:False, 1:True}[0]
    dictflow = {100:"get_formulary", 
                310:"1D_1/2_paramagnet_way1", 311:"1D_1/2_paramagnet_way2", 
                331:"1D_simple_harmonic_oscillator", 332:"",
                430:"monoatomic_ideal_gas",
                710:"ideal_gas_canonical"}
    flow = [dictflow[i] for i in [310]]
    if test_all: flow = [dictflow[i] for i in dictflow.keys()]

print("Test of the {0}.".format(sets.flow))

# ### get_formulary

### get_formulary
if "get_formulary" in sets.flow:
    ostat.class_type = "micro_canonical_discrete_distinguihable"
    ostat.__init__()
    ostat.get_formulary()
    
    ostat.class_type = "micro_canonical_discrete_indistinguihable"
    ostat.__init__()
    ostat.get_formulary()
    
    ostat.class_type = "micro_canonical_continuous_indistinguihable"
    ostat.__init__()
    ostat.get_formulary()

# ## 3 Paramagnets and Oscillators

# ### A Spin-1/2 Paramagnet Way1

#----A Spin-1/2 Paramagnet Way1
if "1D_1/2_paramagnet_way1" in sets.flow:
    print("A Spin-1/2 Paramagnet Way1")
    
    ostat.class_type = "micro_canonical_discrete_distinguihable"
    ostat.__init__()
    ostat.verbose = True
    [mu,B] = symbols('mu B', real=True)
    xreplaces = {g:1, engF:mu*B*(2*i-3), j:1, n:2}
    display("Single particle partition function:", ostat.Zsp)
    """
    res = ostat.Zsp.xreplace(xreplaces)
    display(simplify(res.doit()))
    """
    commands = ["xreplace", "ostat.Zsp", xreplaces]
    ostat.process(commands)
    Zsp = simplify(ostat.result.doit())
    display(Zsp)
    
    commands = ["xreplace", "ostat.U", xreplaces]
    ostat.process(commands)
    U = simplify(ostat.result.doit())
    display(U)
    
    commands = ["xreplace", "ostat.Cv", xreplaces]
    ostat.process(commands)
    Cv = simplify(ostat.result.doit())
    display(Cv)
    
    commands = ["xreplace", "ostat.S", xreplaces]
    ostat.process(commands)
    S = simplify(ostat.result.doit())
    display(S)
    
    commands = ["xreplace", "ostat.F", xreplaces]
    ostat.process(commands)
    F = simplify(ostat.result.doit())
    display(F)
    
    commands = ["xreplace", "ostat.M", xreplaces]
    ostat.process(commands)
    M = simplify(ostat.result.doit())
    display(M)
    
    simplify(ostat.M.evalf(subs=xreplaces).doit())

# ### A Spin-1/2 Paramagnet Way2

### A Spin-1/2 Paramagnet Way2
if "1D_1/2_paramagnet_way2" in sets.flow:
    print("A Spin-1/2 Paramagnet Way2")
    
    ostat.class_type = "micro_canonical_discrete_distinguihable"
    ostat.__init__()
    ostat.verbose = True
    [mu,B] = symbols('mu B', real=True)
    xreplaces = {g:1, engF:mu*B*(2*i-3), j:1, n:2}
    display("Single particle partition function:", ostat.Zsp)
    
    Zsp = simplify(ostat.Zsp.evalf(subs=xreplaces).doit())
    U   = simplify(  ostat.U.evalf(subs=xreplaces).doit())
    Cv  = simplify( ostat.Cv.evalf(subs=xreplaces).doit())
    S   = simplify(  ostat.S.evalf(subs=xreplaces).doit())
    F   = simplify(  ostat.F.evalf(subs=xreplaces).doit())
    M   = simplify(  ostat.M.evalf(subs=xreplaces).doit())

    # list(map(display, [Zsp,U,Cv,S,F,M]))
    display(Zsp,U,Cv,S,F,M)

# ### An Array of 1-D Simple Harmonic Oscillators    

### An Array of 1-D Simple Harmonic Oscillators
if "1D_simple_harmonic_oscillator" in sets.flow:
    print("1D_simple_harmonic_oscillator")
    print("An Array of 1-D Simple Harmonic Oscillators")
    
    ostat.class_type = "micro_canonical_discrete_distinguihable"
    ostat.__init__()
    ostat.verbose = False
    [h,nu,theta] = symbols('h nu theta', real=True, positive=True)
    xreplaces = {g:1, engF:(i+S(1)/2)*h*nu, j:0, n:inf, (h*nu)/kB:theta}
    display("Single particle partition function:", ostat.Zsp)

    commands = ["xreplace", "ostat.Zsp", xreplaces]
    ostat.result = ostat.process(commands).rhs
    Zsp = simplify(ostat.result.doit())
    display(Zsp)
    
    commands = ["xreplace", "ostat.U", xreplaces]
    ostat.result = ostat.process(commands).rhs
    U = simplify(ostat.result.doit())
    display(U)
    
    commands = ["xreplace", "ostat.Cv", xreplaces]
    ostat.result = ostat.process(commands).rhs
    Cv = simplify(ostat.result.doit())
    display(Cv)

# ## An Array of 3-D Simple Harmonic Oscillators todo

# ## 4 Indistinguishable Particles and Monatomic Ideal Gases  

# ### Monoatomic Ideal Gas

### Monoatomic Ideal Gas
if "monoatomic_ideal_gas" in sets.flow:    
    print("Monoatomic Ideal Gas")
    
    ostat.class_type = "micro_canonical_continuous_indistinguihable"
    ostat.__init__()
    ostat.verbose = False
    [h,nu,theta] = symbols('h nu theta', real=True, positive=True)
    xreplaces = {i:eng, g:4*m*pi*V*(2*m*eng)**(S(1)/2)/(h**3), engF:eng}
    display("Single particle partition function:", ostat.Zsp)

    commands = ["xreplace", "ostat.Zsp", xreplaces]
    ostat.result = ostat.process(commands).rhs
    Zsp = simplify(ostat.result.doit())
    display(Zsp)
    
    commands = ["xreplace", "ostat.U", xreplaces]
    ostat.result = ostat.process(commands).rhs
    U = simplify(ostat.result.doit())
    display(U)
    
    commands = ["xreplace", "ostat.S", xreplaces]
    ostat.process(commands)
    S = simplify(ostat.result.doit())
    display(S)

# ## 7 Electrons in Metals 

# ### The Ideal Gas in the Canonical Ensemble

#---The Ideal Gas in the Canonical Ensemble
if "ideal_gas_canonical" in sets.flow:
    print("The Ideal Gas in the Canonical Ensemble")
    
    ostat.class_type = "canonical"
    ostat.__init__()
    ostat.verbose = True
    
    ostat.ZN  = Eq( ostat.ZN.lhs, ostat.subformulary.Z_Ideal_Gas)
    display(ostat.F)
