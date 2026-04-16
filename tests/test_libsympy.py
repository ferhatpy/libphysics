# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
test_libsympy.py

"""
import scipy.constants as pc
from libphysics.libsympy import *
from libphysics.quantum_mechanics import *

func_names = [{100:"plot_sympfunc", 200:"plot_list", 350:"read_latex_file",
              300:"pprints", 500:"substitute", 600:"table_function",
              700:"plot_energy_levels"}[i] 
                for i in [700,]]
input_path = "input/libsympy/"
print(func_names, "\n")


#----Algebra
if "substitute" in func_names:
    print("substitute([x,y**2,z], {x: 1, y: 2, z: 3})-> [x,y,z]=", 
          substitute([x,y**2,z], {x: 1, y: 2, z: 3}))

#----Converters
if "read_latex_file" in func_names:
    file_path = input_path + 'read_latex_file.tex'
    formulas, tex = read_latex_file(file_path)
    
    for iformula in formulas:
        print(iformula)    

#----Differential Equations


#----Linear Algebra


#----Plotting
#----> plot_energy_levels
if "plot_energy_levels" in func_names:
    # --- Example Usage ---
    freq = 1e10
    data_list = {
        1e8:[[0.1, 0, 0.000206783384846193],
        [0.1, 1, 0.00229529557179274],
        [0.1, 2, 0.0229757018902605],
        [0.99, 1, 0.000620371041749169],
        [0.99, 2, 0.00103402136239374]],
        
        1e10:[[0.1, 0, 0.0206783384846193],
        [0.1, 1, 0.229529557179274],
        [0.1, 2, 2.29757018902605],
        [0.99, 0, 0.0206783384846193],
        [0.99, 1, 0.0620371041749169],
        [0.99, 2, 0.103402136239374]],
        
        1e11:[[0.1, 0, 0.206783384846193],
        [0.1, 1, 2.29529557179274],
        [0.1, 2, 22.9757018902605],
        [0.99, 0, 0.206783384846193],
        [0.99, 1, 0.620371041749169],
        [0.99, 2, 1.03402136239374]],
        
        1e12:[[0.1, 0, 2.06783384846193],
        [0.1, 1, 22.9529557179274],
        [0.1, 2, 229.757018902605],
        [0.99, 0, 2.06783384846193],
        [0.99, 1, 6.20371041749169],
        [0.99, 2, 10.3402136239374]],
        
        1e15:[[0.1, 0, 2067.83384846193],
        [0.1, 1, 22952.9557179274],
        [0.1, 2, 229757.018902605],
        [0.99, 0, 2067.83384846193],
        [0.99, 1, 6203.71041749169],
        [0.99, 2, 10340.2136239374]]
        }[freq]

    constYs = np.array([10,20])*1000*pc.Boltzmann/pc.electron_volt
    plot_energy_levels(data_list, constYs, title=freq)
    
#----> plot_list
if "plot_list" in func_names:
    xs1 = np.linspace(0, 2*np.pi, 90)
    ys1 = [sin(ix) for ix in xs1]
    ys2 = [cos(ix) for ix in xs1]
    list2Ds = [[xs1, ys1], [xs1, ys2]]
    plot_list(list2Ds, plabels=["sin(x)","cos(x)"],
              xlabel="$x$", ylabel="$y$",
              pxscale="linear", pyscale="linear",
              pgrid=True, paxis=True)

#----> plot_sympfunc
if "plot_sympfunc" in func_names:
    m = Symbol('m', real = True)
    ms=[1,2,3]
    psi = m*(exp(x) - exp(-x))
    psis = [psi.subs({m:im})*conjugate(psi.subs({m:im})) for im in ms]
    plot_sympfunc(psis, (-5.5,5.5,101), plabels = ms, xlabel="$x$", ylabel=r"$|\psi|^2$", paxis = True)
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
        
if "table_function" in func_names:
    print("table_function")
    print(table_function(oqmec.qho.En(n), n=[0,1,2]))
    display(*table_function(oqmec.qho.En(n), n=[0,1,2]))
    print_latex(table_function(oqmec.qho.En(n), n=[0,1,2]))