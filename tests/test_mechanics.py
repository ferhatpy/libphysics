#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_mechanics.py

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
    
    input_dir  = "input/mechanics"
    output_dir = "output/mechanics"
    
    # Plotting settings
    plot_time_scale = {1:"xy", 2:"xz", 3:"yz"}[3]
    
    flow = [{100:"get_formulary", 150:"get_subformulary",
             200:"simple_harmonic_oscillator", 
             2321:"coordinate_systems", 2322:"moving_particle"}[i] 
            for i in [200]]

### Formulary
print("Test of the {0}.".format(sets.flow))
if "get_formulary" in sets.flow:
    omech.class_type = "scalar"
    omech.__init__()
    omech.get_formulary()
    omech.get_formulary(style="eq")
    
    omech.class_type = "vectorial"
    omech.__init__()
    omech.get_formulary()    
    
    omech.class_type = "EulerLagrange"
    omech.__init__()
    omech.get_formulary()    

if "get_subformulary" in sets.flow:
    omech.__init__()
    omech.get_subformulary()

if "simple_harmonic_oscillator" in sets.flow:
    """       
    Example: Solve a from F = ma
    """
#    omech = mechanics() # DO NOT create any instance.
    omech.class_type = "scalar"
    omech.__init__()
    omech.solver.verbose = True
    commands = ["solve", "NewtonsLaw2", a]
    omech.process(commands)


    """
    Example: Solve position of a spring.
    F = ma, F = -kx
    -kx = ma
    -kx = m d^2 x/dt^2
    w = sqrt(k/m)
    x(t) = C1*sin(wt) + C2*sin(wt)
    """
    omech.__init__()
    omech.solver.verbose = True
    display("Newton's 2nd Law", omech.NewtonsLaw2, 
            "Hooke's Law", omech.HookesLaw)
    commands = ["Eq", "NewtonsLaw2", "HookesLaw"]
    omech.process(commands)
    commands = ["subs", "omech.result", [(a, diff(x, t, 2, evaluate=False))]]
    res = omech.process(commands)
    simp = simplify(res.lhs/m)
    omech.result = Eq(simp, 0)
    commands = ["subs", "omech.result", [(k/m, w**2)]]
    omech.process(commands)
    commands = ["dsolve", "omech.result", x]
    omech.process(commands)
    
    print("Codes:\n", *omech.solver.get_codes())
    
if "coordinate_systems" in sets.flow:
    print("Coordinate Systems")
    
    print("Polar Coordinates")
    omech.class_type = "vectorial"
    omech.__init__()
    omech.solver.verbose = False
    

    xreplaces = {x:omech.subformulary.pol_to_cart_x,
                 y:omech.subformulary.pol_to_cart_y,
                 z:0} # C.k
    
    # xreplaces = {x:r_t*cos(theta)*C.i,
    #              y:r_t*sin(theta)*C.j,
    #              z:0}
    display(omech.r, omech.v, omech.a)
    display(xreplaces)

    commands = ["xreplace", "omech.r", xreplaces]
    x = omech.x.evalf(subs=xreplaces).doit()
    y = omech.y.evalf(subs=xreplaces).doit()
    z = omech.z.evalf(subs=xreplaces).doit()
    r = omech.process(commands)
    v = omech.v.evalf(subs=xreplaces).doit()
    a = omech.a.evalf(subs=xreplaces).doit()
    display(x,y,z,r,v,a)
    
    print("Components of r")
    [display(r.rhs.args[i]) for i in range(2)]
    print("Components of v")
    [display(v.rhs.args[i]) for i in range(2)]
    print("Components of a")
    [display(a.rhs.args[i]) for i in range(2)]
            
"""
def main():    
if __name__ == '__main__':
    main()
"""