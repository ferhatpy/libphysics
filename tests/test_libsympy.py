#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_libsympy.py

"""
# Import path for library functions.
import sys
lstPaths = ["../src"]
for ipath in lstPaths:
    if ipath not in sys.path:
        sys.path.append(ipath)
from libsympy import *

func_names = [{100:"plot_sympfunc", 200:"plot_list", 350:"read_latex_file",
              300:"pprints", 500:"substitute"}[i] 
                for i in [350,]]
input_path = "input/libsympy/"
print(func_names, "\n")


#----Algebra
if "substitute" in func_names:
    print("substitute([x,y**2,z], {x: 1, y: 2, z: 3})-> [x,y,z]=", 
          substitute([x,y**2,z], {x: 1, y: 2, z: 3}))

#----Converters
if "read_latex_file" in func_names:
    file_path = input_path+'read_latex_file.tex'
    formulas, tex = read_latex_file(file_path)
    
    for iformula in formulas:
        print(iformula)    

#----Differential Equations


#----Linear Algebra


#----Plotting
if "plot_list" in func_names:
    xs1 = np.linspace(0, 2*np.pi, 90)
    ys1 = [sin(ix) for ix in xs1]
    ys2 = [cos(ix) for ix in xs1]
    list2Ds = [[xs1, ys1], [xs1, ys2]]
    plot_list(list2Ds, plabels=["sin(x)","cos(x)"],
              xlabel="$x$", ylabel="$y$",
              pxscale="linear", pyscale="linear",
              pgrid=True, paxis=True)

if "plot_sympfunc" in func_names:
    m = Symbol('m', real = True)
    ms=[1,2,3]
    psi = m*(exp(x) - exp(-x))
    psis = [psi.subs({m:im})*conjugate(psi.subs({m:im})) for im in ms]
    plot_sympfunc(psis, (-5.5,5.5,101), plabels = ms, xlabel="$x$", ylabel="$|\psi|^2$", paxis = True)
    plt.show()

    f = x**2
    plot_sympfunc([f,],(-4,4,101))
    plt.show()


#----Printing
if "pprints" in func_names:
    output_styles=["display","pprint","print","latex"]
    alist = [1,2,3,]
    astring = "A string"
    asympy = cosh(sqrt(x))
    for i in output_styles:
        pprints(f"Output Style: {i}",
                "alist=", alist, 
                "astring=", astring,
                "asympy=", asympy,
                output_style = i)
        print("\n")