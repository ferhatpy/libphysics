#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
libsympy.py

Usage
=====

In jupyter-notebook:
python2.#
execfile('/media/hdd/python/projects/libpython/src/libsympy.py')

In a python module, jupyter-notebook.
# Import path for library functions.
import sys
lstPaths = ["/media/hdd/python/projects/libpython/src"]
for ipath in lstPaths:
    if ipath not in sys.path:
        sys.path.append(ipath)
from libsympy import *
"""
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
import mpmath as mp
import numpy as np
from pylab import *
from sympy import *
from sympy.abc import x,y,z,t
from sympy.integrals.manualintegrate import manualintegrate
from sympy.plotting import *
from sympy.solvers.ode import *
from sympy.vector import *

# Initiate rendering Latex, HTML, Math etc.
from IPython import get_ipython
from IPython.display import display, HTML, Latex, Math
from sympy.interactive import printing
printing.init_printing()


#----Global definitions
# Integer symbols
[i,j,k,n,s] = symbols('i j k n s', integer=True, real=True)

# Real symbols.
# A = Symbol('A', real=True)
lst_symbols_reals = ['x','y','z','r','t','T','q',
                     'a','b','c',
                     'alpha','beta','gamma','theta','phi',
                     'C1','C2','C3']
[x,y,z,r,t,T,q,
 a,b,c,
 alpha,beta,gamma,theta,phi,
 C1,C2,C3] = [Symbol(isr, real=True) for isr in lst_symbols_reals]  # SymPy symbol definitions.

# Functions
lst_functions = ['f','g','h']
[f,g,h] = [Function(ifun) for ifun in lst_functions] # SymPy function definitions.
# V=symbols('V', cls=Function)


#----Algebra
def eliminate(system, symbols):
    """
    todo: Find reference from https://stackoverflow.com/
    
    Eliminates given variables from a system of equations.
    
    Args:
        system: List of SymPy equations.
        symbols: List of SymPy symbols to eliminate. Must be shorter than list of equations.
            
    Returns:
        A list of SymPy equations with N_symbols fewer equations than the input. The new
        equations have had the given symbols eliminated.
        
    Example:
        eliminate(flatten([Eqs_BCs12, Eqs_BCs23]), [A2p,B2p,A2m,B2m])
    """
    new_system = system.copy()
    for symbol in symbols:
        solvedfor = None
        elimidx = None
        for idx, eq in enumerate(new_system):
            if symbol in eq.free_symbols:
                if solvedfor is None:
                    solvedfor = solve(eq, symbol, dict=True)[0]
                    elimidx = idx
                else:
                    new_system[idx] = eq.subs(solvedfor)
        del new_system[elimidx]
    return new_system

def condense(system, symbols, elim):
    """
    todo: Find reference from https://stackoverflow.com/
    
    Solves a system for a set of variables while eliminating others.
    Depends on function eliminate.
    
    Args:
        system: List of SymPy equations. Must have len(system) == len(symbols) + len(elim).
        symbols: List of SymPy symbols to solve for.
        elim: List of SymPy symbols to eliminate.
            
    Returns:
        A list of len(symbols) SymPy equations with the given variables eliminated.
        Equations have symbols to be solved for isolated on the left hand side.
    
    Example:
        condense(flatten([Eqs_BCs12, Eqs_BCs23]) , [A1p,B1p,A1m,B1m,A2p,B2p,A2m,B2m], [A2p,B2p,A2m,B2m])
    
    """
    return [Eq(k, v) for (k, v) in solve(eliminate(system, elim), symbols).items()]

def substitute(pexpressions, psubstitutions):
    """
    substitute([x,y,z], {x: 1, y: 2, z: 3})
    
    pexpressions = [x, y, z]
    psubstitutions = {x: 1, y: 2, z: 3} OR
    psubstitutions = [(x, 1), (y, 2), (z, 3)]
    substitute(pexpressions, psubstitutions)
    """
    res = []
    for ieq in pexpressions:
        ieq = ieq.subs(psubstitutions)
        res.append(ieq)
    return res


#----Differential Equations
def solve_odes(equations, func=y, output_style="display"):
    """
    x,t,k,a = symbols("x, t, k, a")
    [y,p,q,r] = [Function('y')(x), Function('p'), Function('q'), Function('r')(t)]

    diff_equations = {"diffeq":Eq(y.diff(x,2)+3*y.diff(x)+2*y,1/(1+exp(x)))}
    libsympy.solve_odes(diff_equations)
    """
    for key,value in equations.items():
        (label,eq) = (key,value)
        sol = dsolve(eq, func, check=True)
        pprints(label,  eq,
                "Solution=", sol,
                "Simplified solution=", simplify(sol),
                "Checking=", checkodesol(eq,sol),
                output_style = output_style)
        

#----Functions
def get_iterated_functions(f, fixed_vals={C1:0, C2:0}, prm=alpha, 
                           param_vals=np.arange(1,2,0.2)):
    """
    fixed_vals = {A:1, w0:1}
    fixed_func = function.subs(fixed_vals)
    ifunc = lambda i:fixed_func.subs(beta,i) # Lambda function
    funcs = list(map(ifunc, np.arange(0.1,1.2,0.1))) # All functions.
    
    Returns a list of functions.
    
    Usage:
        fixed_vals = {A:1, w0:1}
        param_vals = np.arange(0.1,1.2,0.1)
        display(get_iterated_functions(omech.scaled_amplitude, fixed_vals, beta, np.arange(0.1,1.2,0.1)))
    """
    fixed_func = f.subs(fixed_vals)          # Substitute fixed numerical values to symbols.
    ifunc = lambda i:fixed_func.subs(prm, i) # Construct a lambda function with an independent parameter prm.
    funcs = list(map(ifunc, param_vals))
    return(funcs)
        
def get_piecewise():
    """
    todo 
    Returns a piecewise function of the function.
    
    def getfV(self, fV, xs):
        # Plots potential function
        xintervals = []
        for i in range(len(xs)):
            if i == 0:
                xintervals.append((fV[i], x < xs[i]))
            else:
                xintervals.append((fV[i], (xs[i-1] < x) & (x < xs[i])))
                
            if i == len(xs)-1:
                xintervals.append((0, True)) # fV[i]
                
        res = Piecewise(*xintervals)
        return(res)
    """
    xintervals = []
    [fV, xs] = [self.V, self.xs]
    for i in range(len(xs)):
        if i == 0:
            xintervals.append((fV[i], x < xs[i]))
        else:
            xintervals.append((fV[i], (xs[i-1] < x) & (x < xs[i])))
            
        if i == len(xs)-1:
            xintervals.append((0, True)) # fV[i]
            
    res = Piecewise(*xintervals)
    return(res)
    
#def plot_V(sub_num):
#    # todo
#    plot(get_V().subs(sub_num))
    
def rect(x):
    res = Lambda((x), Heaviside(x + 1/2) - Heaviside(x - 1/2))
    return(res)


#----Linear Algebra
def eq_convert_to_matrix(eqs, x):
    # 
    [A,b] = linear_eq_to_matrix(eqs, x)
    return([A, Matrix(x), b])

def eq_solve(eqs, x, ptype="LUsolve"):
    if ptype == "LUsolve":
        [A,x,b] = eq_convert_to_matrix(eqs, x)
        x = A.LUsolve(b)
    elif ptype == "nonlinsolve":
        # todo not working
        nonlinsolve(eqs, x)
    return(x)
    
def solve_2x2_matrix(lhsM, coeffsM, rhsM):
    """
    |A1-| = |A  B||A1+|
    |B1-|   |C  D||B1+|

    Solve unknown A,B,C,D for known A1-,B1-,A1+,B1+
    
    eqs = []
    eqs.append( Eq(solA1m, A*solA1p + B*solB1p) )
    eqs.append( Eq(solB1m, C*solA1p + D*solB1p) )
    
    Parameters
    ----------
    lhsM : 2 x 1 matrix
        DESCRIPTION.
    coeffsM : 2 x 2 matrix
        DESCRIPTION.
    rhsM : 2 x 1 matrix
        DESCRIPTION.

    Example
    -------
    [A,B,C,D] = symbols('A,B,C,D')
    lhsM    = Matrix([[A1m], [B1m]])
    coeffsM = Matrix([[A,B], [C,D]])
    rhsM    = Matrix([[A1p], [B1p]])
    solve_2x2_matrix(lhsM, coeffsM, rhsM)
    
    # nonlinsolve(flatten([Eqs_BCs12, Eqs_BCs23]), [A2p,B2p,A2m,B2m]) gives extra solution???.
    solAB12 = solve(flatten([Eqs_BCs12, Eqs_BCs23]), [A1p,B1p,A1m,B1m], dict=True)
    """
    eqs = []
    rhs1 = coeffsM*rhsM
#        rhs2 = coeffsM.inv()*lhsM
    for i in range(2):
        ieq1 = Eq(lhsM[i], rhs1[i])
#            ieq2 = Eq(rhsM[i], rhs2[i])
        eqs.append(ieq1)
#            eqs.append(ieq2)
    
    res = solve(flatten(eqs), flatten(coeffsM.tolist()), dict=True)
    return(res)    


#----Plotting
def plot_list(plist2Ds, plabels=[1,2,3], 
              xlabel="$x$", ylabel="$y$",
              pxscale="linear", pyscale="linear", 
              pgrid=False, paxis=False):
    """
    Plots functions within a specified region.
    
    Parameters
    ==========
    
    plist2Ds: A lists of xs & ys. [[xs1, ys1], [xs1, ys2]]

    Examples
    ========

    xs1 = np.linspace(0, 2*np.pi, 90)
    ys1 = [sin(ix) for ix in xs1]
    ys2 = [cos(ix) for ix in xs1]
    list2Ds = [[xs1, ys1], [xs1, ys2]]
    plot_func(list2Ds, plabels=["sin(x)","cos(x)"],
              xlabel="$x$", ylabel="$y$",
              pxscale="linear", pyscale="linear",
              pgrid=False, paxis=False)
    """
    fig, ax = plt.subplots()
    linestyles = ['-',':','-.','--']
    
    for i, ilist in enumerate(plist2Ds):
        [xs, ys] = [ilist[0], ilist[1]]
        ax.plot(xs, ys, color='black', 
                linestyle=linestyles[i % 4], 
                linewidth=2, 
                label=plabels[i].__str__())
    
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 18,
            }
    
    ax.set_xlabel(xlabel, fontdict=font)
    ax.set_ylabel(ylabel, fontdict=font)
    ax.set_xscale(pxscale)
    ax.set_yscale(pyscale)
    ax.legend(loc='best', fontsize=14)
    ax.grid(pgrid)
    
    # Plot axis of the origin.
    if paxis:
        plt.axvline(0, color='k')
        plt.axhline(0, color='k')
    
    fig.tight_layout()


def plot_sympfunc(pfuncs, prange=(-1,1,500), plabels=[1,2,3], xlabel="$x$", ylabel="$y$", 
                  pxscale="linear", pyscale="linear", pgrid=False, paxis=True):
    """
    Plots sympy functions within a specified region.
    
    Parameters
    ==========
    
    pfuncs: A list of sympy functions.
    prange: (xmin,xmax,# of points)

    Examples
    ========

    f=x**2
    plot_sympfunc([f,], (-4,4,101))
    plot_sympfunc([[R_nl(i, j, b*a, Z=1/a).subs({a**(3/2):1}).evalf().subs({b:x}) for j in range(i)] for i in range(1,4)], (0, 18, 100), xlabel="$r/a$", ylabel="$R_{nl}(r)$")
    """
    fig, ax = plt.subplots()
    [xmin, xmax, points] = prange
    xs = np.linspace(xmin, xmax, points)
    linestyles = ['-',':','-.','--']
    
    for i, ifunc in enumerate(pfuncs):
        # Lambdify sympy function and obtain f(x) for processing with numpy.
        # lambdify(x, pfunc)(ix)->f(x)
        f = lambdify(x, ifunc, 'numpy') # 'mpmath'
        ys=[f(ix) for ix in xs]
#        ys=[lambdify(x,ifunc,'numpy')(ix) for ix in xs]
        ax.plot(xs, ys, color='black', linestyle=linestyles[i % 4], linewidth=2, label=plabels[i].__str__())
    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)
    ax.set_xscale(pxscale)
    ax.set_yscale(pyscale)
    ax.legend(loc='best')
    ax.grid(pgrid)
    # Plot axis of the origin.
    if paxis:
        plt.axvline(0, color='k')
        plt.axhline(0, color='k')
    plt.rcParams["text.usetex"]
    fig.tight_layout()
    
def plot_save(pfilepath="output", ppad_inches=0.05, pformats=("png","pdf","svg")):
    """
    Saves plot into an image file.
    ppad_inches=0.05 gives nice padding.
    """
    if "eps" in pformats:
        plt.savefig(pfilepath+".eps",format='eps',dpi=1200,bbox_inches='tight',pad_inches=ppad_inches)
    elif "png" in pformats:
        plt.savefig(pfilepath+".png",format='png',dpi=200, bbox_inches='tight',pad_inches=ppad_inches)
    elif "pdf" in pformats:
        plt.savefig(pfilepath+".pdf",format='pdf',dpi=1200,bbox_inches='tight',pad_inches=ppad_inches)
    elif "svg" in pformats:
        plt.savefig(pfilepath+".svg",format="svg")


#----Printing
def pprints(func, *funcs, **kwargs):
    """
    Reference: https://butterflyofdream.wordpress.com/
    
    Parameters
    ==========
    **kwargs
    output_style : display, pprint, print, latex
    
    Examples
    ========
    pprints("f(x)", f,
            "collect(f,x)=", collect(f,x),
            output_style = {1:"display", 2:"pprint", 3:"print", 4:"latex"}[1],
            newline=True)
    """
    # Get **kwargs optional parameters.
    output_style = kwargs.get("output_style", "display")
    newline = kwargs.get("newline", False)
    
    if output_style == "display":
        display(func)
        if funcs is None: return
        for i, f in enumerate(funcs):
            display(f)
            if newline==True and (i % 2)==True: print()
    
    elif output_style == "pprint":
        pprint(func)
        if funcs is None: return
        for i, f in enumerate(funcs):
            pprint(f)
            if newline==True and (i % 2)==True: print()
    
    elif output_style == "print":
        print(func)
        if funcs is None: return
        for i, f in enumerate(funcs):
            print(f)
            if newline==True and (i % 2)==True: print()
    
    elif output_style == "latex":
        print(func)
        if funcs is None: return
        for i, f in enumerate(funcs):
            # todo select by type 
            if type(f) == type(""):
                print(r"$\rm {0}$".format(f))
            # elif type(f) == type([]):
            #     print(*f, sep = "\n")
            else:
                print(r"\[{0}\]".format(print_latex(f)))
                # print(r"$\displaystyle {0}$".format(print_latex(f)))
            if newline==True and (i % 2)==True: print(r"\\")
            # if newline==True : print(r"\\")

def print_matrix_elements(mat, **kwargs):
    """
    M = Matrix([[1,2],[3,4]])
    todo apply latex template.
    """
    output_style = kwargs.get("output_style", "latex")
    
    M = Matrix([[1,2],[3,4]])
    for i in range(mat.rows):
        for j in range(mat.cols):
            if output_style == "display":
                display(f"m{i+1}{j+1} = ", mat[i,j])
            elif output_style == "pprint":
                pprint(f"m{i+1}{j+1} = {mat[i,j]}")
            elif output_style == "print":
                print(f"m{i+1}{j+1} = {mat[i,j]}")
            elif output_style == "latex":
                print("\\begin{multline}")
                print_latex(f"m{i+1}{j+1} = ")
                print_latex(mat[i,j])
                print("\\end{multline}")
    

#----Statistics
def meanT(pfunc, pvar=t, pint=[t,t+T]):
    res = 1/T*integrate(pfunc, (pvar, pint[0], pint[1]))
    return res            


print("libsympy is loaded.")

"""
if sys.version_info[0]<3:
    execfile('./../../projects/libpython/nbcommon.py')
elif sys.version_info[0]>=3:
    exec(compile(open(filename, "rb").read(), filename, 'exec'))
"""