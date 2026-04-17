# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
libsympy.py

Requirements
============
pip3 install latex2sympy2 # This also installs antlr4-python3-runtime==4.11 
pip3 install antlr4-python3-runtime==4.11 # Only for using parse_latex method.
              

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
import mpmath as mp
import numpy as np
import scipy.constants as pc

try:
    import matplotlib.pyplot as plt
    from matplotlib.colors import LightSource
    from pylab import arange, linspace, figure, subplot, show, xlabel, ylabel, title, legend, grid, savefig
except ImportError:
    plt = None
    LightSource = None
from sympy import (
    symbols, Symbol, Function, Eq, solve, simplify, checkodesol, dsolve, Piecewise, Lambda, Heaviside, linear_eq_to_matrix,
    Matrix, flatten, nonlinsolve, Mul, lambdify, pprint, integrate, collect
)
from sympy.abc import x, y, z, t
from sympy.integrals.manualintegrate import manualintegrate
from sympy.plotting import plot
from sympy.physics.quantum import TensorProduct
from sympy.vector import gradient, divergence, curl, vector_integrate
from sympy.vector.parametricregion import ParametricRegion
from sympy import sympify

# Initiate rendering Latex, HTML, Math etc.
from IPython import get_ipython
from IPython.display import display, HTML, Latex, Math
from sympy.parsing.latex import parse_latex
from sympy.parsing.sympy_parser import parse_expr
from sympy.interactive import printing
printing.init_printing()


# Sets global defaults for all plots
plt.rcParams.update({
    'font.family': "serif",
    'font.size': 14,          # Global font size
    'axes.labelsize': 15,     # X and Y label size
    'xtick.labelsize': 14,    # X tick size
    'ytick.labelsize': 14,    # Y tick size
    'legend.fontsize': 12,    # Legend size
})

# Import External Libraries
# from latex2sympy2 import latex2sympy

#----Global definitions
# Integer symbols
i, j, k, n, s = symbols('i j k n s', integer=True, real=True)

# Real symbols.
lst_symbols_reals = ['x','y','z','r','t','T', 'q',
                     'a','b','c', 'g', 'L', 'W',
                     'alpha','beta','gamma','theta','phi', 'tau', 'omega',
                     'C1','C2','C3', 'A', 'M' ,'N', 'B', 'V', 'E', 'H', 'Q', 'Z',
                     'V0', 'x0', 'y0', 'z0', 'v0', 'k0', 'w0',
                     'Ax','Ay','Az', 'Bx','By','Bz', 'Cx','Cy','Cz', 'Dx','Dy','Dz']
x, y, z, r, t, T, q, a, b, c, g, L, W, alpha, beta, gamma, theta, phi, tau, omega, C1, C2, C3, A, M, N, B, V, E, H, Q, Z, V0, x0, y0, z0, v0, k0, w0, Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, Dx, Dy, Dz = [Symbol(isr, real=True) for isr in lst_symbols_reals]  # SymPy symbol definitions.

__all__ = [
    'i', 'j', 'k', 'n', 's',
    'x', 'y', 'z', 'r', 't', 'T', 'q', 'a', 'b', 'c', 
    'alpha', 'beta', 'gamma', 'theta', 'phi', 'tau', 'omega',
    'C1', 'C2', 'C3', 'A', 'M', 'N', 'B', 'V', 'E', 'H', 'Q', 'Z',
    'V0', 'x0', 'y0', 'z0', 'v0', 'k0', 'w0',
    'Ax','Ay','Az', 'Bx','By','Bz', 'Cx','Cy','Cz', 'Dx','Dy','Dz',
    'f', 'g', 'h',
    'eliminate', 'condense', 'substitute', 'solve_odes', 'read_latex_file',
    'table_function', 'get_iterated_functions', 'get_piecewise', 'rect',
    'eq_convert_to_matrix', 'eq_solve', 'solve_2x2_matrix', 'tensor_to_product',
    'plot_energy_levels', 'plot_list', 'plot_sympfunc', 'plot_save',
    'pprints', 'print_matrix_elements', 'meanT',
    'display', 'HTML', 'Latex', 'Math'
]


# Functions
lst_functions = ['f','g','h']
f, g, h = [Function(ifun) for ifun in lst_functions] # SymPy function definitions.
# V=symbols('V', cls=Function)


#----Algebra
from typing import Any, List, Dict

def eliminate(system: List[Any], symbols: List[Any]) -> List[Any]:
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

def condense(system: List[Any], symbols: List[Any], elim: List[Any]) -> List[Any]:
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

def substitute(pexpressions: List[Any], psubstitutions: Any) -> List[Any]:
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
def solve_odes(equations: Dict[Any, Any], func: Any = y, output_style: str = "display") -> None:
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


#----Converters
"""
expr = parse_latex(r"\frac{1}{2}")
print(expr)
"""
def read_latex_file(file_path: str) -> Any:
    import re
    def convert_latex_expression(latex_expression):
        # Replace any occurrence of \hat{} with 'hat'
        hat_pattern = re.compile(r'\\hat\{(.+?)\}')
        latex_expression = re.sub(hat_pattern, r'\\hat', latex_expression)
        return latex_expression
    
    try:
        with open(file_path, 'r') as file:
            latex_content = file.read()

            # Extracting equations using regex without dollar signs
            equation_pattern = re.compile(r'(?<=\$).*?(?=\$)')
            tex = re.findall(equation_pattern, latex_content)

            # Perform additional conversions on each equation
            converted_equations = [convert_latex_expression(eq) for eq in tex]

            # Converting LaTeX expressions to Sympy objects using parse_latex
            res = []
            for ieq in converted_equations:
                try:
                    res.append(parse_latex(ieq))
                except:
                    continue
            return(res, tex)

    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
        return None

#----Functions
from itertools import product
def table_function(expr: Any, eval_: bool = False, verbose: bool = False, **kwargs: Any) -> List[Any]:
    """
    Generates a table for a sympy expression with any number of variables.
    The order of arguments in the function call determines the nesting order.
    
    # 1. Define SymPy symbols
    q, n, m = symbols('q n m')
    
    # 2. Create a complex expression
    # Formula: (q * n) + m
    my_expr = q * n + m
    
    # 3. Call the function with 3 variables
    # Loop order: q (outer), then n (middle), then m (inner)
    table_results = table_function(
        my_expr, 
        q=[0.1, 0.2], 
        n=[1, 2], 
        m=[10, 20]
    )
    
    # 4. Print the table formatted nicely
    print(f"{'q':>5} | {'n':>5} | {'m':>5} | {'Result':>10}")
    print("-" * 35)
    for row in table_results[1:]:  # Skip the header for custom formatting
        print(f"{row[0]:>5} | {row[1]:>5} | {row[2]:>5} | {row[3]:>10.2f}")
    """
    # 1. Find the actual symbol objects used inside the expression
    # This ensures we match 'n' even if it has special assumptions (like integer=True)
    expr_symbols = {s.name: s for s in expr.free_symbols}
    
    # 2. Map the names from kwargs to the actual symbols found in the expression
    ordered_names = list(kwargs.keys())
    value_lists = list(kwargs.values())
    
    # Identify which actual symbol objects to substitute
    # If the name isn't in the expression, we create a dummy symbol
    subs_symbols = [expr_symbols.get(name, Symbol(name)) for name in ordered_names]
    
    # 3. Build the table
    header = ordered_names + ["Result"]
    table = [header]
    
    for values in product(*value_lists):
        # Create the substitution dictionary: {Actual_Symbol_Object: Value}
        subs_dict = dict(zip(subs_symbols, values))
        
        # .subs() replaces the symbol
        # .doit() evaluates the Hermite polynomial and Factorial
        # .evalf() converts the final expression to a decimal
        if eval_:
            # result = expr.subs(subs_dict).doit().evalf()
            result = expr.xreplace(subs_dict).doit().evalf()
        else:
            # result = expr.subs(subs_dict)
            result = expr.xreplace(subs_dict)
        
        if verbose:
            display(result)
        
        table.append(list(values) + [result])
        
    return table


def get_iterated_functions(
    f: Any,
    fixed_vals: Dict[Any, Any] = {C1:0, C2:0},
    prm: Any = alpha,
    param_vals: Any = np.arange(1,2,0.2)
) -> List[Any]:
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
        
def get_piecewise(fV: Any, xs: Any) -> Any:
    """
    Returns a piecewise function (potential function) for given fV and xs.
    Usage: get_piecewise(fV, xs)
    """
    xintervals = []
    for i in range(len(xs)):
        if i == 0:
            xintervals.append((fV[i], x < xs[i]))
        else:
            xintervals.append((fV[i], (xs[i-1] < x) & (x < xs[i])))
        if i == len(xs) - 1:
            xintervals.append((0, True))
    res = Piecewise(*xintervals)
    return res

    
#def plot_V(sub_num):
#    # todo
#    plot(get_V().subs(sub_num))
    
def rect(x: Any) -> Any:
    res = Lambda((x), Heaviside(x + 1/2) - Heaviside(x - 1/2))
    return(res)


#----Linear Algebra
def eq_convert_to_matrix(eqs: Any, x: Any) -> List[Any]:
    # 
    [A,b] = linear_eq_to_matrix(eqs, x)
    return([A, Matrix(x), b])

def eq_solve(eqs: Any, x: Any, ptype: str = "LUsolve") -> Any:
    if ptype == "LUsolve":
        [A,x,b] = eq_convert_to_matrix(eqs, x)
        x = A.LUsolve(b)
    elif ptype == "nonlinsolve":
        # todo not working
        nonlinsolve(eqs, x)
    return(x)
    
def solve_2x2_matrix(lhsM: Any, coeffsM: Any, rhsM: Any) -> Any:
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

def tensor_to_product(expr: Any) -> Any:
    """
    from sympy import TensorProduct, Mul
    
    Replace TensorProduct(a, b, ...) with ordinary multiplication a*b*...
    expr = TensorProduct(a, b) + TensorProduct(a, b, c)
    converted = tensor_to_product(expr)
    """
    return expr.replace(
        lambda x: isinstance(x, TensorProduct),
        lambda x: Mul(*x.args)
    )


#----Plotting
from typing import Optional, Sequence

def plot_energy_levels(
    data: Sequence[Any],
    constYs: Optional[Sequence[Any]] = None,
    title: str = "",
    line_width: float = 0.025
) -> None:
    """
    Plots horizontal lines at y-points grouped by x-coordinates.
    data: List of lists [[x, n, y], ...]
    line_width: The horizontal span of the 'little lines'
    
    
    """
    plt.figure(figsize=(9, 7))
    
    # Extract unique x values to help with axis formatting
    unique_x = sorted(list(set(row[0] for row in data)))
    
    for x, n, y in data:
        # Draw a horizontal line centered at x
        # xmin and xmax define the start and end of the 'little line'
        plt.hlines(y=y, xmin=x - line_width, xmax=x + line_width, 
                   color='royalblue', linewidth=2.0)
        
        # Optional: Label each line with its 'n' value
        plt.text(x + line_width + 0.01, y, f"n={n}", 
                 verticalalignment='center', fontsize=10, color='black')
    
    if type(constYs) is not type(None):
        for i, iy in enumerate(constYs):
            plt.hlines(y=iy, xmin=0, xmax=1, 
                       color=['blue','red'][i], 
                       linewidth=1.0,
                        label=[r'$T_c$=10 K', r'$T_h$=20 K'][i])
            plt.text(0.9, iy+iy*0.2, [rf"$T_c$=10 K", rf"$T_h$=20 K"][i] , 
                     verticalalignment='center', fontsize=10, color='black')

    # Formatting the plot
    plt.xticks(unique_x)
    plt.xlabel("Deformation Parameter q")
    plt.ylabel("Energy (meV)")
    # plt.title(rf"{title:g}")
    
    # If your y-values vary greatly (e.g., 0.2 to 229), 
    # a log scale often makes the plot much more readable:
    plt.yscale({1:'log', 2:'linear'}[2]) 
    
    plt.legend(title=rf"$\nu = {title:g}$ Hz")
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.savefig(f"output/libsympy/E_levels_nu={title:g}Hz.pdf", format='pdf', dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.show()
    
    
def plot_list(
    plist2Ds: Sequence[Any],
    plabels: Sequence[Any] = [1, 2, 3],
    xlabel: str = "$x$",
    ylabel: str = "$y$",
    pxscale: str = "linear",
    pyscale: str = "linear",
    pgrid: bool = False,
    paxis: bool = False
) -> None:
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
    plot_list(list2Ds, plabels=["sin(x)","cos(x)"],
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


def plot_sympfunc(
    pfuncs: Sequence[Any],
    prange: Any = (-1, 1, 500),
    plabels: Sequence[Any] = [1, 2, 3],
    xlabel: str = "$x$",
    ylabel: str = "$y$",
    pxscale: str = "linear",
    pyscale: str = "linear",
    pgrid: bool = False,
    paxis: bool = True,
    ptitle: str = ""
) -> None:
    """
    Plots sympy functions within a specified region.
    
    Parameters
    ==========
    
    pfuncs: A list of sympy functions.
    prange: (xmin,xmax,# of points)

    Examples
    ========

    f=x**3
    plot_sympfunc([f,], (-4,4,101))
    plot_sympfunc([[R_nl(i, j, b*a, Z=1/a).subs({a**(3/2):1}).evalf().subs({b:x}) for j in range(i)] for i in range(1,4)], (0, 18, 100), xlabel="$r/a$", ylabel="$R_{nl}(r)$")
    """
    fig, ax = plt.subplots()
    [xmin, xmax, points] = prange
    xs = np.linspace(xmin, xmax, points)
    linestyles = ['-',':','-.','--']
    colors = ["black", "red", "blue", "green"]    
    
    for i, ifunc in enumerate(pfuncs):
        # Lambdify sympy function and obtain f(x) for processing with numpy.
        # lambdify(x, pfunc)(ix)->f(x)
        f = lambdify(x, ifunc, 'numpy') # 'mpmath'
        ys=[f(ix) for ix in xs]
#        ys=[lambdify(x,ifunc,'numpy')(ix) for ix in xs]
        ax.plot(xs, ys, color=colors[i % 4], linestyle=linestyles[i % 4], linewidth=2, label=plabels[i].__str__())
    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)
    ax.set_xscale(pxscale)
    ax.set_yscale(pyscale)
    ax.set_title(ptitle)
    ax.legend(loc='best', frameon=True, edgecolor='black')
    ax.grid(pgrid)
    
    # Plot axis of the origin.
    if paxis:
        plt.axvline(0, color='k')
        plt.axhline(0, color='k')
    plt.rcParams["text.usetex"]
    fig.tight_layout()
    plt.show()
    
    
def plot_save(
    pfilepath: str = "output",
    ppad_inches: float = 0.05,
    pformats: Sequence[str] = ("png", "pdf", "svg")
) -> None:
    """
    Saves plot into an image file.
    ppad_inches=0.05 gives nice padding.
    
    plot_save(pfilepath="output/libdiscretemath_DFT_"+seq_name, pformats="pdf")
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
def pprints(func: Any, *funcs: Any, **kwargs: Any) -> None:
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
            # If type is a string 
            if type(f) == type(""):
                print(r"$\rm {0}$".format(f))
            # elif type(f) == type([]):
            #     print(*f, sep = "\n")
            else:
                print(r"\[{0}\]".format(print_latex(f)))
                # print(r"$\displaystyle {0}$".format(print_latex(f)))
            if newline==True and (i % 2)==True: print(r"\\")
            # if newline==True : print(r"\\")

def print_matrix_elements(mat: Any, **kwargs: Any) -> None:
    """
    Usage
    =====
    M = Matrix([[1,2],[3,4]])
    print_matrix_elements(M)
    print_matrix_elements(M, output_style="display")
    
    todo apply latex template.
    """
    output_style = kwargs.get("output_style", "latex")

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
def meanT(pfunc: Any, pvar: Any = t, pint: Any = [t, t+T]) -> Any:
    res = 1/T*integrate(pfunc, (pvar, pint[0], pint[1]))
    return res            


print("libsympy is loaded.")

"""
if sys.version_info[0]<3:
    execfile('./../../projects/libpython/nbcommon.py')
elif sys.version_info[0]>=3:
    exec(compile(open(filename, "rb").read(), filename, 'exec'))
"""
