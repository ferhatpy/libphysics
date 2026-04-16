# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# ## test_template.py

"""
test_template.py

Installation:
=============
sudo apt install pandoc             # Converts a *.ipynb file to a *.pdf file via latex.
sudo pip3 install nbextensions      # Jupyter-notebook extension.

Converting from .py to .ipynb format
====================================
In a notebook, do the following in a cell:
%pip install jupytext
!jupytext --to notebook <name_of_script_file>.py

References:
===========    


Find and replace template with desired class name.

Example: otemp
============= 
template.class_type = "micro_canonical_discrete_distinguihable"
template.__init__()
template.solver.verbose = False
[mu,B] = symbols('mu B', real=True)
xreplaces = {g:1, engF:mu*B*(2*i-3), j:1, n:2}

# 4 lines
commands = ["xreplace", "otemp.Zsp", xreplaces]
template.process(commands)
Zsp = simplify(template.result.doit())
display(Zsp)

# Deep copy to create new instances of obranch class.
import copy
otemp2 = copy.deepcopy(otemp)


assert euler_equations(D(x, t)**2/2, {x}) == [Eq(-D(x, t, t), 0)]
"""
import copy
from libphysics.libsympy import *
from libphysics.template import *
# Execute jupyter-notebook related commands.
#exec(open('libnotebook.py').read())
# print(sys.version)
# print(sys.path)

# ### Settings

### Settings
class sets:
    """
    Setttings class.
        
    Instead of settings class, settings nametuple might be used.
    Settings = namedtuple("Settings", "type dropinf delta")
    sets = Settings(type="symbolic", dropinf=True, delta=0.1)
    print(set.type)
    """
    global dictflow, test_all
    
    def __init__(self):
        pass

    # File settings
    input_dir  = "input/template"
    output_dir = "output/template"
    
    # Plotting settings
    plot_time_scale = {1:"xy", 2:"xz", 3:"yz"}[3]
    
    # Execution settings.
    test_all = {0:False, 1:True}[0]
    dictflow = {100:"get_formulary", 150:"get_subformulary",
                200:"topic1", 300:"topic2",
                400:"topic3"}
    flow = [dictflow[i] for i in [200]]
    if test_all: flow = [dictflow[i] for i in dictflow.keys()]

print("Test of the {0}.".format(sets.flow))

# ### get_formulary

#### get_formulary
if "get_formulary" in sets.flow:
#    omech = mechanics() # DO NOT create any instance.
    otemp.class_type = ""
    otemp.__init__()
    otemp.get_formulary()
    otemp.get_formulary(style="eq")

# ### get_subformulary

#### get_subformulary    
if "get_subformulary" in sets.flow:
    otemp.class_type = ""
    otemp.__init__()
    otemp.get_subformulary()    

# ### A Spin-1/2 Paramagnet

### A Spin-1/2 Paramagnet
if "topic1" in sets.flow:
    print("A Spin-1/2 Paramagnet")
    otemp.__init__(class_type = "micro_canonical_discrete_distinguihable")
    otemp.verbose = True
    
    [mu,B] = symbols('mu B', real=True)
    xreplaces = {g:1, engF:mu*B*(2*i-3), j:1, n:2}
    display("Single particle partition function:", otemp.Zsp)

    ### Magic 4 lines
    commands = ["xreplace", "otemp.Zsp", xreplaces]
    otemp.result = otemp.process(commands).rhs
    Zsp = simplify(otemp.result.doit())
    display(Zsp)
    
    ### Alternative ways.
    commands = ["xreplace", "otemp.Zsp", xreplaces]
    otemp.process(commands)
    Zsp1 = simplify(otemp.result.evalf())   # 1. Way; evalf() can give a result in a different format.
    Zsp2 = simplify(otemp.result.doit())    # 2. Way.
    Zsp3 = simplify(otemp.Zsp.evalf(subs=xreplaces))        # 3. Way, substituion without commands and execution.
    Zsp4 = simplify(otemp.Zsp.evalf(subs=xreplaces).doit()) # 4. Way, execution without commands.
    Zsp5 = simplify(otemp.Zsp.xreplace(xreplaces).doit())   # 5. Way. THe most beautiful way.
    display(Zsp1, Zsp2, Zsp3, Zsp4, Zsp5)
    
    commands = ["xreplace", "otemp.U", xreplaces]
    otemp.process(commands)
    U = simplify(otemp.result.doit())
    display(U)
    
    #### Get generated SymPy codes.
    print("Codes:\n", *otemp.solver.get_codes())
    
    #### Print calculations, results, etc.
    pprints("Z_{sp}1=", Zsp1, 
            "U=", U,
            output_style = otemp.output_style )
    
    list(map(display, [Zsp1, U]))
    
    print(multiline_latex(Zsp1.lhs, Zsp1.rhs))
