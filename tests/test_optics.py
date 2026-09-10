# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# ## test_optics

"""
test_optics.py connected to test_optics.ipynb via "jupytext".
In ipynb notebook select File->Jupytext->Pair Notebook with Light Format.

omec.__init__()
omec.verbose = True
commands = ["solve", "NewtonsLaw2", a]
print(omec.process(commands))

References:
===========    
Abedin, Islam, Haider, 2007, Computer simulation of Fresnel diffraction from rectangular apertures and obstacles using the Fresnel integral
Gerrard, A., Burch, J.M., 1975, Introduction to matrix methods in optics
"""

import copy
import sys
import os

def find_local():
    current = os.getcwd()
    while True:
        candidate = os.path.join(current, "libphysics", "src")
        if os.path.exists(candidate):
            return candidate
        parent = os.path.dirname(current)
        if parent == current:
            return None
        current = parent
local_path = find_local()
sys.path.insert(0, local_path)
print("-> Using local libphysics")

"""lstPaths = ["../../src"]
libphysics\src
for ipath in lstPaths:
    if ipath not in sys.path:
        sys.path.append(ipath)"""
import matplotlib
matplotlib.use('Qt5Agg')   # or 'Qt5Agg', 'Qt6Agg', 'GTK3Agg'
import matplotlib.pyplot as plt
from libsympy import *
from optics import *
import torchsympy
import importlib
importlib.reload(torchsympy)
# from numba import jit
import torch
from torchquad import GaussLegendre
# Execute jupyter-notebook related commands.
#exec(open('libnotebook.py').read())
print(sys.version); print(sys.path)


# ### Settings

#### Settings
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
    
    input_dir  = "input/optics"
    output_dir = "output/optics"
    
    # Plotting settings
    plot_time_scale = {1:"xy", 2:"xz", 3:"yz"}[3]
    
    # Execution settings.
    test_all = {0:False, 1:True}[0]
    usecupy = {0:False, 1:True}[0]
    usetorchsympy = {0:False, 1:True}[1]
    dictflow = dict(
        ch1 = {100:"get_formulary", 150:"get_subformulary",
               200:"ABCD_2_thin_lens", 202:"ABCD_microscope",
               250:"diffraction_rectangular", 252:"Fraunhofer_3D",
               350:"FBG_Reflection"})
    flow = [dictflow["ch1"][i] for i in [200]]
    if test_all: flow = flatten([list(dictflow[i].values()) for i in dictflow.keys()])

print("Test of the {0}.".format(sets.flow))
print("Using cupy: {0}".format(sets.usecupy))

# ### get_formulary

#### get_formulary
if "get_formulary" in sets.flow:
    oopti.__init__()
    oopti.get_formulary()

# ### get_subformulary

#### get_subformulary
if "get_subformulary" in sets.flow:
    oopti.__init__()
    oopti.get_subformulary() 
    
    
#### ABCD
#----> ABCD_2_lens
if "ABCD_2_thin_lens" in sets.flow:
    """
    Gerrard1975, 2.7 Illustrative Problems, Problem 5 p48.
    
    """
    # 
    print("ABCD_2_thin_lens")
    print("")
    oopti.__init__()
    oopti.verbose = False
    
    abcd = oopti.ABCD

    L1 = abcd.thin_lens.abcd(f1).rhs.doit()
    L2 = abcd.thin_lens.abcd(f2).rhs.doit()
    D = abcd.T(d).rhs.doit()
    G = abcd.T(g).rhs.doit()
    B = abcd.T(b).rhs.doit()
    
    # Analitic Calculation
    TM_ = MatMul(B, L2, D, L1, G)
    TM_ = Eq(M, UnevaluatedExpr(TM_))
    TM = TM_.rhs.doit()
    solb = solve(TM[0,1], b)[0]
    d_expr = f1 + E + f2
    solb_E = simplify(solb.subs(d, d_expr)) # Find the distance between image and the 2. lens.
    M = 1/TM[1,1]
    list(map(display, [L1, L2, D, G, B, TM_, TM, solb, d_expr]))
    
    # Numerical Calculation
    num_values = {f1:0.08, f2:-0.12, d:0.06, g:0.24}
    TM_num = TM.evalf(subs=num_values)
    b_num = solb.evalf(subs=num_values) # Distance between image and the 2. lens.
    M_num = 1/TM_num[1,1]               # Magnification
    list(map(display, [num_values, TM_num, b_num, M_num]))

#----> ABCD_microscope    
if "ABCD_microscope" in sets.flow:
    print("ABCD_microscope")
    print("")
    oopti.__init__()
    oopti.verbose = False
    
    pprints(
        "System matrix=", oopti.ABCD.microscope.system_matrix,
        "Trasnfer matrix=", oopti.ABCD.microscope.transfer_matrix,
        "The distance between image and the ocular, b=", oopti.ABCD.microscope.imaging_condition,
        "Magnification", oopti.ABCD.microscope.magnification)


# ### DIFFRACTION RECTANGULAR

#### DIFFRACTION RECTANGULAR
if "diffraction_rectangular" in sets.flow:
    """
    from scipy.integrate import quad
    from scipy.special import jn
    integrand = lambda x, r: np.exp(-x** 2 ) * np.exp(np.pi* 1j *(-x)) * jn( 0 , r*x) * x
    intensity = lambda r: np.sqrt(quad( lambda x: np.real(integrand(x, r)),0,5 )[0]** 2 + quad(lambda x: np.imag(integrand(x, r)), 0 , 5)[0]**2)
    
    import matplotlib.pyplot as plt
    t = np.linspace(0, 3)
    plt.plot(t, np.vectorize(intensity)(t))
    """
    
    print("Diffraction from a Rectangular Aperture")
    class_type = {1:"Rayleigh_Sommerfeld", 2:"Fraunhofer", 3:"Fresnel",
                  4:"FresnelJit"}[3]
    oopti.__init__(class_type)
    oopti.verbose = False
    
    [Lx,Ly] = symbols('L_x L_y', real=True, positive=True)
    rectX = oopti.subformulary.rect(x0/Lx)
    rectY = oopti.subformulary.rect(y0/Ly)
    xreplaces = {Uapr:1, Eli:1, 
                 x0min:-S(1)/2*Lx, x0max:S(1)/2*Lx,
                 y0min:-S(1)/2*Ly, y0max:S(1)/2*Ly}
#    xreplaces = {Uapr:1, Eli:1, 
#                 x0min:-Lx, x0max:Lx,
#                 y0min:-Ly, y0max:Ly}
    
    # Numerical Calculations
    # config in mm.
    # Screen mesh size.  The integrand is evaluated for every quadrature point
    # and every screen point at once, i.e. on an array of
    # (quadrature points) x n x n numbers, so the memory needed grows with n**2:
    # n=200 needs ~1.5 GB and n=1000 would need ~38 GB in a single array.
    # torchsympy splits such a mesh into chunks automatically; chunk_size_params
    # below sets the chunk explicitly.  See LARGE_GRIDS.md and
    # examples/large_screen_diffraction.py for the details and for a separable
    # formulation that is ~400x faster (and more accurate) for this aperture.
    n = 1000
    # Number of screen points evaluated simultaneously; keeps the memory bounded
    # independently of n.  Lower it if you still run out of memory.
    chunk_size_params = 4096
    screen_factor = 2.5;
    config = {0:0, 1:"LakshminarayananFig11_3",
              61:"Abedin2005Fig6a", 62:"Abedin2005Fig6b",
              63:"Abedin2005Fig6c", 64:"Abedin2005Fig6d", 
              72:"Abedin2005Fig7b"}[61]
    
    # Aperture size (mm). Square or rectangle.
    [nLx, nLy] = {0:[1, 1],
                "LakshminarayananFig11_3":[0.11, 0.11],
                "Abedin2005Fig6a":[2,2], "Abedin2005Fig6b":[2,2],
                "Abedin2005Fig6c":[2,2], "Abedin2005Fig6d":[2,2], 
                "Abedin2005Fig7b":[2,2]}[config]
    
    # Wavelength in (mm).
    nl = {0:1, 
          "LakshminarayananFig11_3":560e-6,
          "Abedin2005Fig6a":632e-6, "Abedin2005Fig6b":632e-6,
          "Abedin2005Fig6c":632e-6, "Abedin2005Fig6d":632e-6,
          "Abedin2005Fig7b":1264e-6}[config]
    
    # Screen distance in (mm).
    nz = {0:0.5, 
          "LakshminarayananFig11_3":3, 
          "Abedin2005Fig6a":400, "Abedin2005Fig6b":800, # Abedin2005Fig6a 400
          "Abedin2005Fig6c":1700,"Abedin2005Fig6d":8000, 
          "Abedin2005Fig7b":400}[config]
    brightness = {"Rayleigh_Sommerfeld":1, "Fraunhofer":0.1, "Fresnel":1}[oopti.class_type];
    lsX = np.linspace(-nLx*screen_factor, nLx*screen_factor, n)
    lsY = np.linspace(-nLy*screen_factor, nLy*screen_factor, n)
    X,Y = np.meshgrid(lsX, lsY)
    #    X,Y = np.mgrid[-nLx:nLx:0.0015, -nLy:nLy:0.0015]
    subs = {z:nz, l:nl, Lx:nLx, Ly:nLy, k:2*pi/l}
    print(config, " ", oopti.class_type)
    
    # Symbolic Calculations
    commands = ["xreplace", "oopti."+class_type+".integral", xreplaces]
    oopti.process(commands)
    diffr_form = oopti.result.doit()
    display(diffr_form)
    
    commands = ["xreplace", "oopti."+class_type+".EField", xreplaces]
    oopti.process(commands)
    El = oopti.result.doit()
    display(El)
    
    commands = ["xreplace", "oopti."+class_type+".intensity", xreplaces]
    oopti.process(commands)
    Int = oopti.result
    if oopti.class_type == "Fraunhofer":
        Int = Int.doit()
        Int = simplify(Int.rewrite(sin))
        oopti.result = Int
    display(Int)
    
    commands = ["subs", "oopti.result", subs]
    Int = oopti.process(commands)
    print(multiline_latex(Int.lhs, Int.rhs))
    lt = torchsympy.TorchSymPy()
    
    
    #----> Fraunhofer
    if oopti.class_type == "Fraunhofer":
        # Method 1
#        fInt = lambdify([x,y], Int.rhs)
#        Z = fInt(X,Y)
        
        # Method 2 - Convert symbolic expression to numerical expression. 
#        fInt = lambda ix,iy: Int.rhs.xreplace({x:ix, y:iy}).doit().evalf()
        # fInt = lambdify([x,y], Int.rhs.xreplace({x:x, y:y}).doit().evalf(), "scipy")
        if not sets.usetorchsympy:
            fInt = lambdify([x,y], Int.rhs.evalf(quad='osc'), "scipy")
            Z = np.vectorize(fInt)(X,Y)
        else:
            print("using torchsympy...")
            """Z = lt.eval_numeric(Int, params_values=[X, Y],
                solver="batched", N=1001, method="gauss-legendre",
                chunk_size_points=4096, dtype=torch.complex128)"""
            
            Z = lt.eval_numeric(Int, params_values=[X, Y], method=GaussLegendre(),
                solver="vectorized", N=5001, #high node number is needed
                chunk_size_params=chunk_size_params)
    
    
    #----> Rayleigh_Sommerfeld, Fresnel
    if oopti.class_type in ["Rayleigh_Sommerfeld", "Fresnel"]:
        if not sets.usecupy:
            if not sets.usetorchsympy:
                # fInt = lambda ix,iy: Int.rhs.xreplace({x:ix, y:ix}).doit().evalf() # it takes longer time.
                fInt = lambdify([x,y], Int.rhs.xreplace({x:x, y:y}).doit().evalf(quad='osc'), "scipy")
                Z = np.vectorize(fInt)(X,Y)
            else:
                print("using torchsympy...")
                # Z = lt.eval_numeric(Int, params_values=[X, Y],
                #     solver="batched", N=1001, method="gauss-legendre",
                #     chunk_size_points=4096, dtype=torch.complex128)
                
                # Z = lt.eval_numeric(Int, params_values=[X, Y], method=GaussLegendre(),
                #     solver="vectorized", N=501) #high node number is needed 
                
                # N is shared between both integration axes here, so N=5001
                # really means ~70 nodes per axis; method=GaussLegendre() (as in
                # the Fraunhofer branch above) is markedly more accurate than
                # the default Simpson rule at the same cost.
                Z = lt.eval_numeric(Int, params_values=[X, Y],
                    solver="vectorized", N=5001, #high node number is needed
                    chunk_size_params=chunk_size_params)
        else:
            import cupy as cp
            
            """
            # Assuming Int.rhs.xreplace({x:x, y:y}).doit().evalf(quad='osc') is the symbolic expression:
            fInt = lambdify([x, y], Int.rhs.xreplace({x: x, y: y}).doit().evalf(quad='osc'), "scipy")
            # Define a GPU-based vectorized function
            def fInt_cupy(x, y):
                x_np, y_np = cp.asnumpy(x), cp.asnumpy(y)  # Convert CuPy arrays to NumPy
                result = fInt(x_np, y_np)  # Apply the original SciPy-based function
                return cp.asarray(result)  # Convert the result back to CuPy
            # Perform element-wise computationss
            Z = fInt_cupy(X, Y).get()
            """
            
            fInt = lambdify([x,y], Int.rhs.xreplace({x:x, y:y}).doit().evalf(quad='osc'), "scipy")
            Z = cp.asarray(fInt(cp.asnumpy(X), cp.asnumpy(Y))).get()


    #----> Plotting 2D Diffraction Intensity
    if torch.is_tensor(Z):
        Z = Z.detach().cpu().numpy()
    fig = plt.figure(figsize=(5, 5))
    ax1 = fig.add_subplot(111)
    ax1.imshow(Z, cmap=plt.cm.gray, interpolation ='bilinear', origin='lower', vmin=np.min(Z), vmax=brightness*np.max(Z))
    ax1.set_xticks(np.linspace(0, n, 5))
    ax1.set_xticklabels([-nLx*screen_factor, -nLx*screen_factor*0.5, 0,
                          nLx*screen_factor*0.5, nLx*screen_factor])
    ax1.set_yticks(np.linspace(0, n, 5))
    ax1.set_yticklabels([-nLy*screen_factor, -nLy*screen_factor*0.5, 0,
                          nLy*screen_factor*0.5, nLy*screen_factor])
    os.makedirs(sets.output_dir, exist_ok=True)
    plt.savefig("{0}/{1}_{2}_{3}.{4}".format(\
                sets.output_dir, sets.flow[0], oopti.class_type, config, "png"), format="png", dpi=600, bbox_inches='tight')
    plt.show()

    #----> Plotting 3D Fraunhofer Diffraction Intensity
    if oopti.class_type == "Fraunhofer":
        if "Fraunhofer_3D" in sets.flow:
            plot3d(Int.rhs, (x,-1,1), (y,-1,1))
        

#### FIBER BRAGG GRATING
if "FBG_Reflection" in sets.flow:
    print("Fiber Bragg Grating")
    print("R versus lambda_0, Ghatak2009 Appendix C Eq.3")
    oopti.__init__()
    oopti.verbose = False
    
    sym_replaces = {
                lambda_B: oopti.FBG.lambda_B.rhs,   # Most fundamental: λ_B = 2Λn₀
                kappa: oopti.FBG.kappa.rhs,         # Depends on λ_B and Δn
                Gamma: oopti.FBG.Gamma.rhs,         # Depends on λ_B and n₀
                alpha: oopti.FBG.alpha.rhs          # Depends on κ and Γ
                }

    # Perform recursive substitution
    # todo PUT INTO NOTES
    R = oopti.FBG.R.rhs
    for _ in range(3):  # 3 iterations sufficient for this dependency chain
        R = R.subs(sym_replaces)
    R1 = Eq(oopti.FBG.R.lhs, R)
    
    # Same as above but walrus operator is used.
    R = oopti.FBG.R.rhs
    R = Eq(oopti.FBG.R.lhs, [R := R.xreplace(sym_replaces) for _ in range(3)][-1])
    
    lambda_c = 1300e-9              # central wavelength
    period = 10
    nGaAs, nAlAs = [3.4059, 2.9086] # refractive index
    tGaAs, tAlAs = [lambda_c/(4*nGaAs), lambda_c/(4*nAlAs)] # thickness of the layers
    n0_ = (nGaAs + nAlAs)/2         # effective refractive index
    Delta_n_ = (nGaAs - nAlAs)/sqrt(2)
    num_Lambda = tGaAs + tAlAs     # the period of the z-dependent variation
    tot_length = num_Lambda * period   # Length of the fiber of # period periodic FBG structure.
    num_replaces = {lambda_ : num_Lambda,
                    n0:n0_,
                    Delta_n: Delta_n_,
                    L: tot_length,
                    pi: np.pi
                    }
    
    # 1. Way takes long time for computing.
    """
    R_num = lambda ilambda0: R.xreplace(num_replaces).rhs.evalf(subs={lambda_0:ilambda0})
    llist = np.linspace(0.75, 2, 400) # 9e-9, 40e-9, 40
    Rlist = [R_num(il) for il in llist*1e-6]
    """
        
    # 2. Way takes short time for completing.
    # Create a vectorized function using lambdify with complex-safe sqrt
    Rnum = R.xreplace(num_replaces).rhs
    # R_func = lambdify(lambda_0, Rnum, modules="numpy") # gives <lambdifygenerated-4>:2: RuntimeWarning: invalid value encountered in sqrt
    R_func = lambdify(lambda_0, Rnum, [{"sqrt": np.lib.scimath.sqrt}, "numpy"])
    
    # Compute wavelengths (convert to meters for calculation, nanometers for plotting)
    llist_um = np.linspace(0.75, 2, 1000)  # Wavelengths in micrometers
    llist_m = llist_um * 1e-6              # Convert to meters
    
    # Evaluate R in a vectorized manner and take real components
    Rlist = R_func(llist_m)
    # Rlist = np.real(R_func(llist_m))
    
    # Plotting the results
    file_dir = f"output"+"/"+oopti.classname+"/"+oopti.FBG.classname
    os.makedirs(file_dir, exist_ok=True)
    file_name = f"FBG_R_l_PR_T={period}"
    file_path = file_dir+"/"+file_name
    
    plt.figure(figsize=(8, 5))
    plt.plot(llist_um, Rlist, linestyle='-', color='b', label=r'$R$')
    plt.xlabel(r'$\lambda (\mu m)$', fontsize=14)
    plt.ylabel(r'$R$ (%)', fontsize=14)
    plt.title(r'Plot of $R$ vs $\lambda$', fontsize=16)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend(fontsize=12)
    plt.tight_layout()
    
    # Save the plot as a high-resolution image
    plt.savefig(f"{file_path}.png", format="png", dpi=300)
    plt.show()
    
    # Save results to a file.
    # Save data to file (llist_um in micrometers, Rlist in %)
    np.savetxt(f"{file_path}.txt", 
               np.column_stack( (llist_um, np.real(Rlist)) ),
               header='Wavelength(um) Reflectance(%)')
    