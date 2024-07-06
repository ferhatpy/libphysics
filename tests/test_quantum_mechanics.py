#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ## test_quantum_mechanics.py

"""
test_quantum_mechanics.py connected to test_quantum_mechanics.ipynb 
via "jupytext" light pairing.

References:
===========
    Leonardo Angelini, Solved Problems in Quantum Mechanics, Springer, 2019
    David J. Griffiths, Introduction to Quantum Mechanics, Pearson, 2005
    Zwiebach, MIT 8.05 Quantum Physics II, Fall 2013
    
sudo apt install pandoc # Converts a *.ipynb file to a *.pdf file via latex.


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
# import scienceplots
# plt.style.use(['science', 'notebook'])

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
    test_all = {0:False, 1:True}[0]
    dictflow = dict(
        ch1 = {1:"p1.3",3:"p1.5",4:"p1.9",5:"p1.17"},
        ch2 = {24:"p2.4",27:"p2.7",29:"p2.9",25:"e2.5",211:"p2.11",212:"p2.12",
               232:"ch2.3.2",26:"e2.6",222:"p2.22",26:"ch2.6",233:"p2.33",241:"p2.41"},
        ch3 = {322:"p3.22", 330:"p3.30"},
        ch4 = {401:"p4.1",402:"e4.1",421:"ch4.2.1",411:"p4.11",4:"p4.12",404:"fig4.4",
               413:"p4.13",7:"p4.14",8:"p4.15",9:"ch4.3.1",10:"ch4.4.1",
               11:"e4.2",12:"p4.27",449:"p4.49",16:"p4.55"},
        ch5 = {1:"p5.1"},
        ch6 = {61:"p6.1", 62:"p6.2", 611:"p6.11", 614:"p6.14", 615:"p6.15", 
               253:"c25.3"},
        ch7 = {701:"e7.1"})
    flow = [dictflow["ch6"][i] for i in [253]]
    if test_all: flow = flatten([list(dictflow[i].values()) for i in dictflow.keys()])

print("Test of the {0}.".format(sets.flow))

# ### get_formulary

#### get_formulary
if "get_formulary" in sets.flow:
    omec.__init__()
    omec.get_formulary()
    omec.get_formulary(style="eq")

# ### get_subformulary

#### get_subformulary    
if "get_subformulary" in sets.flow:
    omec.class_type = ""
    omec.__init__()
    omec.get_subformulary()

# ## Chapter 6 Time-Independent Perturbation Theory

# ### ----> 25.3 The Anharmonic Oscillator (Zwiebach)

#----> 25.3 The Anharmonic Oscillator (Zwiebach)
if "c25.3" in sets.flow:
    print("c25.3 The Anharmonic Oscillator (Zwiebach)")
    oqmec.__init__("position_space")
    oqmec.verbose = False
    
    # Apply nondegenerate perturbation theory upto 3rd order.
    # Calculate perturbation contributions to energy and wavefunctions.
    l = var('lambda_')
    psi0c = oqmec.qho.nb
    psi0  = oqmec.qho.nk
    En0   = oqmec.qho.En
    
    # x^4 : Anharmonic oscillator
    oqmec.Hp = Hp = S(1)/4*hbar*w*(oqmec.qho.a + oqmec.qho.ad)**4 # Zwiebach works without mass m.
    j_ = 4 # Power of the perturbation term x^j, j={2,3,4}.
    
    """
    j=4 case
    En1_s1 = oqmec.En_ND_PT(1, psi0c, psi0, Hp, En0, k2min=n, k2max=n)
    En2_s1 = oqmec.En_ND_PT(2, psi0c, psi0, Hp, En0, k2min=n-4, k2max=n+4)
    En3_s1 = oqmec.En_ND_PT(3, psi0c, psi0, Hp, En0, k2min=n-4, k2max=n+4)
    psi1n_s1 = oqmec.psin_ND_PT(1, psi0c, psi0, Hp, En0, k1min=n-4, k1max=n+4)
    psi2n_s1 = oqmec.psin_ND_PT(2, psi0c, psi0, Hp, En0, k1min=n-4, k1max=n+4, k2min=n-4, k2max=n+4)
    # psi3n_s1 = oqmec.psin_ND_PT(3, psi0c, psi0, Hp, En0, k1min=n-4, k1max=n+4, k2min=n-4, k2max=n+4, k3min=n-4, k3max=n+4)
    """    
    En1_s1 = oqmec.En_ND_PT(1, psi0c, psi0, Hp, En0, k2min=n, k2max=n)
    En2_s1 = oqmec.En_ND_PT(2, psi0c, psi0, Hp, En0, k2min=n-j_, k2max=n+j_)
    En3_s1 = oqmec.En_ND_PT(3, psi0c, psi0, Hp, En0, k2min=n-j_, k2max=n+j_)
    psi1n_s1 = oqmec.psin_ND_PT(1, psi0c, psi0, Hp, En0, k1min=n-j_, k1max=n+j_)
    psi2n_s1 = oqmec.psin_ND_PT(2, psi0c, psi0, Hp, En0, k1min=n-j_, k1max=n+j_, k2min=n-j_, k2max=n+j_)
    # psi3n_s1 = oqmec.psin_ND_PT(3, psi0c, psi0, Hp, En0, k1min=n-j_, k1max=n+j_, k2min=n-j_, k2max=n+j_, k3min=n-j_, k3max=n+j_)
    
    # Simplify results.
    En1,En2,En3 = simplify(En1_s1), simplify(En2_s1), simplify(En3_s1)
    En = En0().rhs + l*En1.rhs + l**2*En2.rhs + l**3*En3.rhs
    psi1n = simplify(psi1n_s1)
    psi2n = simplify(psi2n_s1)
    # psi3n = simplify(psi3n_s1)
    psin = psi0() + l*psi1n.rhs + l**2*psi2n.rhs #+ l**3*psi3n.rhs
    
    pprints("oqmec.Hp=", oqmec.Hp,
             "En1=", En1_s1, En1,
             "En2=", En2_s1, En2,
             "En3=", En3_s1, En3,
             "psi1n=", psi1n_s1, psi1n,
             "psi2n=", psi2n_s1, psi2n,
             # "psi3n=", psi3n_s1, psi3n
             )
    
    # Get energy versus n functions.
    substitutions = {hbar:1, w:1, g:1, m:1}
    fEn0 = lambdify(n, En0().rhs.xreplace(substitutions))
    fEn1 = lambdify(n, En1.rhs.xreplace(substitutions))
    fEn2 = lambdify(n, En2.rhs.xreplace(substitutions))
    # fEn3 = lambdify(n, En3.rhs.xreplace(substitutions))
    
    # Total energy with perturbation contributions upto 3rd order
    fEn = lambda n,l: fEn0(n) + l*fEn1(n) + l**2*fEn2(n) #+ l**3*fEn3(n)
    pprints("E(n)=", fEn(n,l))
    
    # Plot energy versus n.
    n_ = np.arange(0,10,1) # Quantum level index.
    l_ = 0.01 # Perturbation constant.
    plt.scatter(n_, fEn0(n_), label=r'$E_n^0$', color='red')
    # plt.scatter(n, l*fEn1(n))
    # plt.scatter(n, l**2*fEn2(n))
    # plt.scatter(n, l**3*fEn3(n))
    plt.scatter(n_, fEn(n_,l_), label=r'$E_n$', color='blue')
    plt.xlabel(r"n")
    plt.ylabel(r"E(n)")
    plt.legend(loc='upper left', frameon=True, edgecolor='black')
    plt.xlim(-0.5,11)
    plt.show()

    # Convert ket states to wavefunctions.
    # Unperturbed wavefunction of harmonic oscillator.
    wfpsi0 = lambda n=n: oqmec.qho.psix(n).rhs.xreplace(substitutions)
    wfpsi1 = lambda n: oqmec.ket_to_wavefunction(n, j:=4, psi0, wfpsi0, psi1n.rhs.xreplace(substitutions))
    wfpsi2 = lambda n: oqmec.ket_to_wavefunction(n, j:=4, psi0, wfpsi0, psi2n.rhs.xreplace(substitutions))
    # wfpsi3 = lambda n=n: ket_to_wavefunction(n, psi0, wfpsi0, l**3*psi3n)
    # Total perturbed wavefunction.
    wfpsi = lambda n,l: wfpsi0(n) + l*wfpsi1(n) + l**2*wfpsi2(n) #+ l**3*wfpsi3(n)
    
    # Plot wavefunctions
    # plot(psi0n(6), (x,-5*pi,5*pi))
    # libsympy.plot_sympfunc([wfpsi0(3),], xlabel=r"$x$", ylabel=r"$\psi_0(x)$", prange=(-5,5,500))
    n_=8 # Excited level index.
    libsympy.plot_sympfunc([wfpsi0(n_), l_*wfpsi1(n_), l_**2*wfpsi2(n_)],
                           plabels=[rf"$\psi_{n_}^0$", rf"$\lambda\psi_{n_}^1$", rf"$\lambda^2\psi_{n_}^2$", rf"$\lambda^3\psi_{n_}^3$"],
                           xlabel=r"$x$", ylabel=r"$\psi(x)$", prange=(-5,5,500))
    
    libsympy.plot_sympfunc([wfpsi0(n_), wfpsi(n_,l_)],
                           plabels=[rf"$\psi_{n_}^0$", rf"$\psi_{n_}$"],
                           xlabel=r"$x$", ylabel=r"$\psi(x)$", prange=(-5,5,500))

# ### ----> p6.2

#----> p6.2
if "p6.2" in sets.flow:
    print("p6.2 Perturbation of quantum harmonic oscillator (1+e)*k")
    oqmec.__init__("position_space")
    oqmec.verbose = True
    
    oqmec.Hp = S(1)/2*epsilon*k*oqmec.qho.x2op.rhs
    # oqmec.Hp = oqmec.qho.x2op*oqmec.qho.x2op
    substitutions = {Hp:oqmec.Hp, nb:oqmec.qho.nb(), nk:oqmec.qho.nk()}
    En1_s1 = oqmec.En1_ND_PT_bk.rhs.xreplace(substitutions)
    En1_s2 = qapply(En1_s1)
    En1 = En1_s2.xreplace({k:m*w**2}).collect(epsilon/2*hbar*w)
    
    pprints("H'=1/2*e*k*x^2=", oqmec.Hp,
            "substitutions=", substitutions,
            "En1=", En1_s1,
            "En1=", En1_s2,
            "En1=", En1,
            output_style = "display")


# ### ----> p6.11

#----> 6.11 Harmonic Oscillator: Cubic Correction 
if "p6.11" in sets.flow:
    print("6.11 Harmonic Oscillator: Cubic Correction (Angelini2019)")
    oqmec.__init__("position_space")
    oqmec.verbose = True
    A,m,w = symbols('A m w', real=True)
    psi0c = oqmec.qho.nb
    psi0  = oqmec.qho.nk
    En0   = oqmec.qho.En
    Hp = A*oqmec.qho.xop.rhs**3
    En1 = oqmec.En_ND_PT(1, psi0c, psi0, Hp, En0)
    En2 = oqmec.En_ND_PT(2, psi0c, psi0, Hp, En0, k2min=n-3, k2max=n+3)

    pprints("Hp", Hp,
            "En1=", En1,
            "En2=", En2,
            output_style = "display")


# ### ----> p6.14

#----> 6.14 Charged Harmonic Oscillator in an Electric Field
if "p6.14" in sets.flow:
    print("6.14 Charged Harmonic Oscillator in an Electric Field (Angelini2019)")
    oqmec.__init__("position_space")
    oqmec.verbose = True
    
    # En^1 1st order correction to energy.
    print("En^1=")
    epsilon,m,q,w = symbols('epsilon m q w', real=True)
    oqmec.Hp = -q*epsilon*oqmec.qho.xop.rhs
    substitutions = {Hp:oqmec.Hp, nb:oqmec.qho.nb(), nk:oqmec.qho.nk()}
    En1_s1 = oqmec.En1_ND_PT_bk.rhs.xreplace(substitutions)
    En1_s2 = qapply(En1_s1)
    En1 = oqmec.En_ND_PT(1, oqmec.qho.nb, oqmec.qho.nk, oqmec.Hp, oqmec.qho.En)
    En2 = oqmec.En_ND_PT(2, oqmec.qho.nb, oqmec.qho.nk, oqmec.Hp, oqmec.qho.En, k2min=n-2, k2max=n+2)
    psi1n = oqmec.psin_ND_PT(1, oqmec.qho.nb, oqmec.qho.nk, oqmec.Hp, oqmec.qho.En, k1min=n-2, k1max=n+2)
    
    pprints("En1_s1=", En1_s1,
            "En1_s2", En1_s2,
            "En1=", En1, simplify(En1),
            "En2=", En2, simplify(En2),
            "psi1n", psi1n,
            "psi1n", simplify(psi1n))


# ### ----> p6.15

#----> 6.15 Harmonic Oscillator: Second Harmonic Potential I
if "p6.15" in sets.flow:
    print("p6.15 Harmonic Oscillator: Second Harmonic Potential I (Angelini2019)")
    oqmec.__init__("position_space")
    oqmec.verbose = True
    
    psi0c = oqmec.qho.nb
    psi0  = oqmec.qho.nk
    En0   = oqmec.qho.En
    Hp = S(1)/2*m*alpha**2*oqmec.qho.x2op.rhs
    En1 = oqmec.En_ND_PT(1, psi0c, psi0, Hp, En0, k2min=n, k2max=n)
    En2 = oqmec.En_ND_PT(2, psi0c, psi0, Hp, En0, k2min=n-4, k2max=n+4)
    En = Eq(S("E_n"), En0(n).rhs + En1.rhs + En2.rhs)
    
    pprints("oqmec.Hp=", oqmec.Hp,
             "En1=", En1, simplify(En1),
             "En2=", En2, simplify(En2),
             "En=", En, simplify(En),
             output_style="display")

# ## Chapter 7 The Variational Principle

# ### ----> 7.1

#----> 7.1
if "e7.1" in sets.flow:
    print("Griffiths2005 e7.1")
    oqmec.__init__("position_space")
    oqmec.verbose = True
    varfx = Wavefunction(A*exp(-b*x**2), x)
    nvarfx = varfx.normalize().simplify()
    Vx = S(1)/2*m*w**2*x**2
    xreplaces = {xmin:-oo, xmax:oo, Psi:nvarfx.expr, V:Vx}
    expH = oqmec.exp_H.xreplace(xreplaces)
    expHs = expH.doit()
    solb = solve(diff(expHs.rhs, b), b)[1]
    expHmin = expHs.subs({b:solb})
    
    pprints(
        "V(x)=", Vx,
        "<H>=", expH,
        "<H>=", expH.doit(),
        "b=", solb,
        "<H>min=", expHmin,
        output_style="display")


