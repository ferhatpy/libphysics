# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
mechanics.py
Created on Fri Mar 11 12:48:37 2022

Important Points:
    Use diff not Derivative
    self.v = self.velocity = Eq(S('v'), diff(self.x, t, evaluate=False))
"""
from __future__ import annotations
from itertools import combinations_with_replacement
from typing import List, Tuple, Any, Dict, Union, Optional, Callable
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d
from sympy import (Derivative, Symbol, Matrix, ratsimp, Derivative as D, Eq, Function, symbols, sqrt, sin, cos, exp, sinh, cosh, tanh, diff, integrate, solve, dsolve, latex, Rational, S, pi, I, Integral, Sum, DiracDelta, Heaviside, sign)
from sympy.vector import CoordSys3D
from .libsympy import *
from libphysics.libreflection import branch
import libphysics.libphyscon as pc

# Global symbols
i, j, n = symbols('i j n', integer=True)
C1, C2, C3 = symbols('C1 C2 C3')
alpha, beta, gamma, tau = symbols('alpha beta gamma tau', real=True, positive=True)
g, t = symbols('g t', real=True, positive=True)
p = symbols('p')
k, m, w, w0, w1, w2, A, M = symbols('k m w w0 w1 w2 A M', real=True, positive=True)
F0, m1, m2, m3, m4, m5 = symbols('F0 m1 m2 m3 m4 m5', real=True, positive=True)
xi, xf, yi, yf, zi, zf = symbols('x_i x_f y_i y_f z_i z_f', real=True)
x0, y0, z0, v0 = symbols('x0 y0 z0 v0', real=True)

theta, phi = [Function('theta')(t), Function('phi')(t)]
q, u, v = [Function('q')(t), Function('u')(t), Function('v')(t)]
r, x, y, z = [Function('r')(t), Function('x')(t), Function('y')(t), Function('z')(t)]
px, py, pz = [Function('p_x')(t), Function('p_y')(t), Function('p_z')(t)]

C = CoordSys3D('C')  # Cartesian coordinate system.
dr = 1*C.i + 1*C.j + 1*C.k # Default, will be updated in __init__ if needed

_G = Function('G')(t, tau)
_Gw = Function('Gtilde')(w)

# EulerLagrange-specific symbols (initialized here for IDE/type-checker support;
# define_symbols() overwrites them with the same values when class_type='EulerLagrange')
_lst_funcs_el = ['x1','x2','x3','x4','x5','y1','y2','y3','y4','y5','z1','z2','z3','z4','z5']
x1, x2, x3, x4, x5 = [Function(_n)(t) for _n in _lst_funcs_el[:5]]
y1, y2, y3, y4, y5 = [Function(_n)(t) for _n in _lst_funcs_el[5:10]]
z1, z2, z3, z4, z5 = [Function(_n)(t) for _n in _lst_funcs_el[10:]]
xdot, ydot, zdot = Function('xdot')(t), Function('ydot')(t), Function('zdot')(t)
pxdot, pydot, pzdot = Function('pdot_x')(t), Function('pdot_y')(t), Function('pdot_z')(t)
q_i, p_i = Function('q_i')(t), Function('p_i')(t)
q_idot, p_idot = Function('qdot_i')(t), Function('pdot_i')(t)
H = Function('H')(q_i, p_i, t)
L = Function('L')(q_i, q_idot, t)
T = Function('T')(q_i, q_idot, t)
V = Function('V')(q_i)
Gamma, lambda_B, lambda_0, lambda_, n0, Delta_n = symbols('Gamma lambda_B lambda_0 Lambda n_0 Delta_n')

__all__ = [
    'mechanics', 'omech',
    'i', 'j', 'n', 'C1', 'C2', 'C3', 'alpha', 'beta', 'gamma', 'tau', 'g', 't', 'p',
    'k', 'm', 'w', 'w0', 'w1', 'w2', 'A', 'M', 'F0', 'm1', 'm2', 'm3', 'm4', 'm5',
    'xi', 'xf', 'yi', 'yf', 'zi', 'zf', 'x0', 'y0', 'z0', 'v0',
    'theta', 'phi', 'q', 'u', 'v', 'r', 'x', 'y', 'z', 'px', 'py', 'pz', 'C', 'dr', '_G', '_Gw',
    'x1', 'x2', 'x3', 'x4', 'x5', 'y1', 'y2', 'y3', 'y4', 'y5', 'z1', 'z2', 'z3', 'z4', 'z5',
    'xdot', 'ydot', 'zdot', 'pxdot', 'pydot', 'pzdot',
    'q_i', 'p_i', 'q_idot', 'p_idot', 'H', 'L', 'T', 'V',
    'Gamma', 'lambda_B', 'lambda_0', 'lambda_', 'n0', 'Delta_n',
]

class mechanics(branch):
    """
    
    """
    _name = "mechanics"
        
    def define_symbols(self) -> None:
        """
        Common global symbols, functions, objects.
        """
        # Symbols are now defined at the module level for better IDE support.
        pass
        
        # Common definitions.
        if self.class_type in ["scalar", "vectorial"]:
            global _G, _Gw
            # G         = Lambda((t,tau), Function('G')(t,tau)) # Callable Green's function.
            # Gw        = Lambda((w), Function('Gtilde')(w))    # Callable Green's tilde function.
            _G  = Function('G')(t,tau)    # Incallable Green's function.
            _Gw = Function('Gtilde')(w)   # Incallable Green's tilde function.
        
        if self.class_type in ["EulerLagrange"]:
            global x1,x2,x3,x4,x5, y1,y2,y3,y4,y5, z1,z2,z3,z4,z5
            lst_funcs = ['x1','x2','x3','x4','x5',
                         'y1','y2','y3','y4','y5',
                         'z1','z2','z3','z4','z5']
            [x1,x2,x3,x4,x5,
             y1,y2,y3,y4,y5,
             z1,z2,z3,z4,z5] = [Function(ifun)(t) for ifun in lst_funcs]
            global xdot, ydot, zdot
            global pxdot, pydot, pzdot
            global q_i,p_i,q_idot,p_idot
            xdot, ydot, zdot = [Function('xdot')(t), Function('ydot')(t), Function('zdot')(t)]
            pxdot, pydot, pzdot = [Function('pdot_x')(t), Function('pdot_y')(t), Function('pdot_z')(t)]
            q_i,p_i       = [Function('q_i')(t), Function('p_i')(t)]  # Generalized coordinates
            q_idot,p_idot = [Function('qdot_i')(t), Function('pdot_i')(t)] 
            global f,H,L,T,V
            f       = Function('f')(u,t)
            H,L,T,V = [Function('H')(q_i,p_i,t),    # Hamiltonian
                       Function('L')(q_i,q_idot,t), # Lagrangian
                       Function('T')(q_i,q_idot,t), # Kinetic energy
                       Function('V')(q_i)]          # Potential enerhy
        

    def __init__(self, class_type: str = 'scalar'):
        super().__init__()
        self.class_type = class_type
        self.define_symbols()
        
        # File settings
        self.input_dir: str = "input/mechanics"
        self.output_dir: str = "output/mechanics"
        
        class subformulary:
            """
            Sub formulary class.
            
            Define global symbols in the outer class.
            """
            def __init__(self) -> None:
                # Transformations from Polar to Cartesian Coordinates.
                self.pol_to_cart_x = r*cos(theta)
                self.pol_to_cart_y = r*sin(theta)
                
                # Transformations from Cylindrical to Cartesian Coordinates.
                self.cyl_to_cart_x = r*sin(phi)
                self.cyl_to_cart_y = r*sin(phi)
                self.cyl_to_cart_z = z
                
                # Transformations from Spherical to Cartesian Coordinates.
                self.sph_to_cart_x = r*sin(theta)*cos(phi)
                self.sph_to_cart_y = r*sin(theta)*sin(phi)
                self.sph_to_cart_z = r*cos(theta)
                
        self.subformulary = subformulary()
        
        class moment_of_inertia(branch):
            """
            List of Moment of Inertia
            """
            def __init__(self) -> None:
                super().__init__()
                self.name: str = "moment_of_inertia"
                
                # Define common symbols
                M, r, a, b, L, R1, R2 = symbols('M r a b L R1 R2')
                
                # Solid objects (about center of mass unless specified)
                self.I_sphere_cm: Any = S(2)/5 * M * r**2          # Solid sphere
                self.I_cylinder_cm: Any = S(1)/2 * M * r**2        # Solid cylinder (central axis)
                self.I_rod_cm: Any = S(1)/12 * M * L**2            # Thin rod (center)
                self.I_rod_end: Any = S(1)/3 * M * L**2            # Thin rod (one end)
                self.I_rect_plate_cm: Any = S(1)/12 * M * (a**2 + b**2)  # Rectangular plate (center)
                self.I_rect_plate_edge: Any = S(1)/3 * M * a**2    # Rectangular plate (edge)
                self.I_cube_cm: Any = S(1)/6 * M * a**2            # Cube (center)
                self.I_cube_edge: Any = S(2)/3 * M * a**2          # Cube (edge)
                self.I_cone_axis_cm: Any = S(3)/10 * M * r**2      # Solid cone (central axis)
                self.I_cone_perpendicular_cm: Any = S(3)/20 * M * (r**2 + 4*L**2)  # Solid cone (perpendicular axis)
                
                # Hollow objects (thin shells or composite)
                self.I_hollow_sphere_cm: Any = S(2)/3 * M * r**2  # Thin spherical shell (center)
                self.I_hollow_cylinder_cm: Any = M * r**2          # Thin cylindrical shell (central axis)
                self.I_hollow_ring_cm: Any = M * r**2              # Thin ring/hoop (central axis)
                self.I_hollow_rect_frame_cm: Any = S(1)/12 * M * (a**2 + b**2)  # Rectangular frame (center)
                self.I_hollow_rect_frame_edge: Any = S(1)/3 * M * a**2  # Rectangular frame (edge)
                self.I_hollow_cylinder_annulus_cm: Any = M * (R1**2 + R2**2)/2  # Hollow cylinder (R1=inner, R2=outer)
                self.I_hollow_sphere_annulus_cm: Any = S(2)/5 * M * (R2**5 - R1**5)/(R2**3 - R1**3)  # Hollow sphere (R1=inner, R2=outer)
                
            def __doc__(self) -> None:
                print("List of Moment of Inertia")
        self.moment_of_inertia = moment_of_inertia()

    ####    1) scalar        
        if self.class_type in ["scalar"]:
            #----> Green's Function Methods
            self.G: Any = _G     # Green's function.
            self.Greens_function: Any = _G
            self.Gw: Any = _Gw    # Green's tilde function.
            self.Greensw_function: Any = _Gw
            self.IFT_Gw: Any = Eq(self.G, 1/sqrt(2*pi)*Integral(self.Gw*exp(I*w*(t-tau)), (w)))
            self.IFT_Dirac_delta: Any = Eq(DiracDelta(t-tau), 1/(2*pi)*Integral(exp(I*w*(t-tau)), (w)))
            self.inverse_Fourier_transform_Gw: Any = self.IFT_Gw
            self.inverse_Fourier_transform_Dirac_delta: Any = self.IFT_Dirac_delta
             
            #----> Mechanics Methods
            self.x: Any = x
            self.position: Any = x
            self.v: Any = Eq(S('v'), diff(self.x, t, evaluate=False))
            self.velocity: Any = self.v
            self.a: Any = Eq(S('a'), diff(self.x, t, 2, evaluate=False))
            self.acceleration: Any = self.a
            self.p: Any = Eq(S('p'), m*self.v.rhs)
            self.momentum: Any = self.p
            self.F: Any = Eq(S('F'), m*self.a.rhs)
            self.force: Any = self.F
            self.NewtonsLaw2: Any = self.F
            self.W: Any = Eq(S('W'), Integral(self.F.rhs, (x,xi,xf)))
            self.work: Any = self.W
            self.T: Any = Eq(S('T'), S(1)/2*m*self.v.rhs**2)
            self.kinetic_energy: Any = self.T
            self.U: Any = Eq(S('U'), Integral(k*self.x, (x,0,self.x)))
            self.potential_energy: Any = self.U
            self.H: Any = Eq(S('H'), self.T.rhs + self.U.rhs)
            self.energy: Any = self.H
            self.HookesLaw: Any = Eq(S('F'), -k*self.x)
             
            self.amplitude: Any = None
            self.phase: Any = None
            self.scaled_amplitude: Any = None

####    2) vectorial            
        if self.class_type in ["vectorial"]:
            """
            Example:
            omech.r -> Eq(r(t), x(t) + y(t) + z(t))
            omech.r.rhs -> x(t) + y(t) + z(t)
            omech.r.rhs.subs({x:r*cos(theta)}) -> r(t)*cos(theta(t)) + y(t) + z(t)
            diff(omech.r.rhs.subs({x:r*cos(theta)}),t) -> -r(t)*sin(theta(t))*diff(theta(t), t, evaluate=False) + cos(theta(t))*diff(r(t), t, evaluate=False) + diff(y(t), t, evaluate=False) + diff(z(t), t, evaluate=False)
            """
            self.x, self.y, self.z = [x, y, z]
            self.r = self.position = Eq(S('r_x')*C.i+S('r_y')*C.j+S('r_z')*C.k, 
                                        self.x*C.i + self.y*C.j + self.z*C.k)
            self.v = self.velocity = Eq(S('v_x')*C.i+S('v_y')*C.j+S('v_z')*C.k, 
                                        diff(self.r.rhs, t, evaluate=False))
            self.a = self.acceleration = Eq(S('a_x')*C.i+S('a_y')*C.j+S('a_z')*C.k,
                                            diff(self.v.rhs, t, evaluate=False))
            self.p = self.momentum = Eq(S('p_x')*C.i+S('p_y')*C.j+S('p_z')*C.k , m*self.v.rhs)
            self.L = self.angular_momentum = Eq(S('L_x')*C.i+S('L_y')*C.j+S('L_z')*C.k , 
                                                self.r.rhs.cross(self.p.rhs))
            self.F = self.force = self.NewtonsLaw2 = Eq(S('F_x')*C.i+S('F_y')*C.j+S('F_z')*C.k, 
                                                        m*self.a.rhs)
            self.Tr1 = self.Torque1 = Eq(S('tau_x')*C.i+S('tau_y')*C.j+S('tau_z')*C.k, 
                                         self.r.rhs.cross(diff(self.p.rhs, t, evaluate=False))) # Torque = r x dp/dt
            self.Tr2 = self.Torque1 = Eq(S('tau_x')*C.i+S('tau_y')*C.j+S('tau_z')*C.k, 
                                         self.r.rhs.cross(self.F.rhs))                          # Torque = r x F
            self.W = self.work = Eq(S('W'), Integral(self.F.rhs.dot(dr), (z,zi,zf), (y,yi,yf), (x,xi,xf)))
            self.T = self.kinetic_energy  = Eq(S('T'), S(1)/2*m*self.v.rhs.dot(self.v.rhs))
            self.U = self.potential_energy= Eq(S('U'), Integral(k*self.r.rhs.dot(dr), (z,zi,zf), (y,yi,yf), (x,xi,xf)))
            self.H = self.energy = Eq(S('H'), self.T.rhs + self.U.rhs)
            self.HookesLaw = Eq(S('F_x')*C.i+S('F_y')*C.j+S('F_z')*C.k, -k*self.x)
        
####    3) EulerLagrange
        if self.class_type in ["EulerLagrange"]:
            """
            todo Write explanation
            """
            self.x, self.y, self.z = [x, y, z]
            self.f, self.q, self.u, self.v = [f,q,u,v]
            self.T = self.kinetic_energy   = Eq(S('T'), T)
            self.V = self.potential_energy = Eq(S('V'), V)
            
            #----> Lagrangian Mechanics
            self.L = L
            self.Lag = self.Lagrangian = Eq(L, self.T.rhs-self.V.rhs)
            self.Eulers_equation = Eq( Sum((-1)**n*D(D(f,u), (t,n), ), (n,1,m)), 0 )
            self.Lagrange_equations_I = Eq( D( D(L, q_idot), t) - D(L, q_i), 0 )
#            self.Lagrange_equations_I = UnevaluatedExpr(Eq( D( D(S('L'), q_idot), t) - D(S('L'), q_i), 0 ))
#            self.u = dynamicsymbols('u')   # gives time dependent function
#            u = IndexedBase('u')           # Generates an error.
#            f = Function('f')(u[n],x)
#            self.Eulers_equations = Eq( Sum((-1)**n*diff(diff(f,u[n]), (x,n)), (n,0,oo)), 0)
            
            #self.Lagrange_equations_II with constraints todo
            
            #----> Hamiltonian Mechanics
            self.H   = Eq(S('H'), Sum(p_i*q_idot, (i,1,n)) - L)
            self.F_i = Eq(S('F_i'), D(L, q_i))
            self.q_i = q_i
            self.p_i = Eq(p_i, D(L, q_idot))
#            self.q_idot = self.Hamiltons_equations_I  = Eq(var('qdot_i'),  D(H, p_i))
            self.q_idot = self.Hamiltons_equations_I  = Eq(q_idot,  D(H, p_i))
            # self.p_idot = self.Hamiltons_equations_II = Eq(var('pdot_i'), -D(H, q_i))
            self.p_idot = self.Hamiltons_equations_II = Eq(p_idot, -D(H, q_i))
            
        # Common text definitions.
        self.Hamiltonian = self.H
        
        
        
#### --- CLASSES ---



#### Harmonic Oscillator
        class oscillator(branch):
            """
            Sub Class for Harmonic, Damped, Driven Oscillator.
            """
            global Gamma, lambda_B, lambda_0, lambda_, n0, Delta_n
            Gamma, lambda_B, lambda_0, lambda_, n0, Delta_n = symbols('Gamma lambda_B lambda_0 Lambda n_0 Delta_n')           
            
            def __init__(self):
                super().__init__()
                self.name = "oscillator"
                
                # Harmonic oscillator substitutions.
                self.A = self.reduced_amplitude = F0/m
                self.beta = self.damping_parameter = gamma/(2*m)
                self.w0 = self.natural_frequency = sqrt(k/m)
                self.underdamping_criteria = {w0:sqrt(beta**2+w1**2)}
                self.critical_damping_criteria = {w0:beta}
                self.overdamping_criteria = {w0:sqrt(beta**2-w2**2)}

                # Damped Oscillator
                self.damped_harmonic_oscillator1 = Eq(m*diff(x, t, 2, evaluate=False) + gamma*diff(x, t, evaluate=False) + k*x, 0)
                self.damped_harmonic_oscillator2 = Eq(diff(x, t, 2, evaluate=False) + 2*beta*diff(x, t, evaluate=False) + w0**2*x, 0)
                
                # Driven Oscillator
                self.driven_oscillator1 = Eq(m*diff(x, t, 2, evaluate=False) + gamma*diff(x, t, evaluate=False) + k*x, F0*cos(w*t))
                self.driven_oscillator2 = Eq(diff(x, t, 2, evaluate=False) + 2*beta*diff(x, t, evaluate=False) + w0**2*x, A*cos(w*t))
                self.driven_oscillator3 = Eq(m*D(x, t, 2) + 2*gamma*m*D(x, t,1) + m*w0**2*x, F0*cos(w*t))
                # self.driven_oscillator_GreenF = Eq( m*(1/sqrt(2*pi))*diff(ometh.IFT_Gw,t,2), DiracDelta(t-tau)) #todo
                    
                # Green Function for Driven Oscillator
                self.G_driven_oscillator_weak_damping = 1/(m*sqrt(w0**2-gamma**2))*exp(-gamma*(t-tau))*sin(sqrt(w0**2-gamma**2)*(t-tau))
                self.G_driven_oscillator_strong_damping = 1/(m*sqrt(gamma**2-w0**2))*exp(-gamma*(t-tau))*sinh(sqrt(gamma**2-w0**2)*(t-tau))
                self.G_driven_oscillator_critical_damping = (t-tau)/m*exp(-gamma*(t-tau))
            
            @staticmethod
            def __doc__():
                return "Sub Class for Harmonic, Damped, Driven Oscillator."
        self.oscillator = self.osc = oscillator()
        
    
#### Global Methods
#----> Eulers_equation_sympy
    @staticmethod
    def Eulers_equation_sympy(L: Any, funcs: Any = (), vars: Any = ()) -> Tuple[List[Any], List[Any]]:
        r"""
        Modified verbose version of sympy.calculus.euler.euler_equations.
        
        Find the Euler-Lagrange equations [1] for a given Lagrangian.
    
               ∞                                  
            _____                                
            ╲                                    
             ╲            n + 1                 
     dF       ╲      n   d                      
    ───── =   ╱  (-1) ⋅─────────(f(u(x), x)) = 0
     du      ╱           n   (a)                   
            ╱          dx  du_(n)(x)                
            ‾‾‾‾‾
            n = 0
                   
    
        Parameters
        ==========
    
        L : Expr
            The Lagrangian that should be a function of the functions listed
            in the second argument and their derivatives.
    
            For example, in the case of two functions `f(x,y)`, `g(x,y)` and
            two independent variables `x`, `y` the Lagrangian would have the form:
    
                .. math:: L\left(f(x,y),g(x,y),\frac{\partial f(x,y)}{\partial x},
                          \frac{\partial f(x,y)}{\partial y},
                          \frac{\partial g(x,y)}{\partial x},
                          \frac{\partial g(x,y)}{\partial y},x,y\right)
    
            In many cases it is not necessary to provide anything, except the
            Lagrangian, it will be auto-detected (and an error raised if this
            couldn't be done).
    
        funcs : Function or an iterable of Functions
            The functions that the Lagrangian depends on. The Euler equations
            are differential equations for each of these functions.
    
        vars : Symbol or an iterable of Symbols
            The Symbols that are the independent variables of the functions.
    
        Returns
        =======
    
        eqns : list of Eq
            The list of differential equations, one for each function.
    
        Examples
        ========
    
        >>> from sympy import Symbol, Function
        >>> from sympy.calculus.euler import euler_equations
        >>> x = Function('x')
        >>> t = Symbol('t')
        >>> L = (x(t).diff(t))**2/2 - x(t)**2/2
        >>> euler_equations(L, x(t), t)
        [Eq(-x(t) - Derivative(x(t), (t, 2)), 0)]
        >>> u = Function('u')
        >>> x = Symbol('x')
        >>> L = (u(t, x).diff(t))**2/2 - (u(t, x).diff(x))**2/2
        >>> euler_equations(L, u(t, x), [t, x])
        [Eq(-Derivative(u(t, x), (t, 2)) + Derivative(u(t, x), (x, 2)), 0)]
        
        x = Symbol('x')
        u = Function('u')(x)
        f = Function('f')(u,x)
        omech.Eulers_equation_sympy(f,[u],x)
    
        References
        ==========
    
        [1] http://en.wikipedia.org/wiki/Euler%E2%80%93Lagrange_equation
    
        """
    
        funcs = tuple(funcs) if iterable(funcs) else (funcs,)
    
        if not funcs:
            funcs = tuple(L.atoms(Function))
        else:
            for f in funcs:
                if not isinstance(f, Function):
                    raise TypeError('Function expected, got: %s' % f)
    
        vars = tuple(vars) if iterable(vars) else (vars,)
    
        if not vars:
            vars = funcs[0].args
        else:
            vars = tuple(sympify(var) for var in vars)
    
        if not all(isinstance(v, Symbol) for v in vars):
            raise TypeError('Variables are not symbols, got %s' % vars)
    
        for f in funcs:
            if not vars == f.args:
                raise ValueError("Variables %s don't match args: %s" % (vars, f))
    
        order = max(len(d.variables) for d in L.atoms(Derivative)
                    if d.expr in funcs)
        
        eqns, steps = [],[]
        for f in funcs:
            eq = diff(L, f)
            # Verbose calculation steps.
            sL = Function('L')(*funcs)
            steps.append(Eq(Derivative(sL,f), diff(L, f)))
            for i in range(1, order + 1):
                for p in combinations_with_replacement(vars, i):
                    eq = eq + S.NegativeOne**i*diff(L, diff(f, *p), *p)
            eqns.append(Eq(eq, 0))
    
        return(eqns, steps)

#----> Eulers_equation_1D
    @staticmethod
    def Eulers_equation_1D(f: Any, uns: List[Any], ivar: Any) -> Tuple[Any, List[Any]]:
        """
        Eulers_equation_1D(f,uns,ivar)
        f: f(x, u(x), u_x(x), ...) is known as the density of the functional F.
           u_x = du/dx
        uns: Dependent functions and their derivatives, [u(x), u_x(x)].
        ivars: Independent variable x or t.
        
        u = Function('u')(x)
        f = Function('f')(u,x)
        omech.Eulers_equation_1D(f,[u],x)
        omech.Eulers_equation_1D(f,[u,u.diff(x)],x)
        omech.Eulers_equation_1D(f,[u,u.diff(x,2)],x)  
        """
        _sum = 0
        steps = []
        for n in range(0,len(uns)):
            df_du = diff(f, uns[n])
            # Verbose calculation steps.
            sL = Function('L')(*uns)
            steps.append(Eq(D(sL,uns[n]), df_du))
            _sum = _sum + S.NegativeOne**n*diff(df_du, (ivar, n), evaluate=False)
        
        eqn = Eq(_sum, 0)
        return(eqn, steps)

#----> Hamiltons_equations
    def Hamiltons_equations(
        self,
        pL: Any,
        lst_qi: List[Any],
        lst_qidot: List[Any],
        lst_pi: List[Any],
        lst_pidot: List[Any]
    ) -> Tuple[List[Any], List[Any]]:
        """
        Usage:
        ======    
        lst_qi    = [x,y,z]
        lst_qidot = [xdot, ydot, zdot]
        lst_pi    = [px,py,pz]
        lst_pidot = [pxdot, pydot, pzdot]
        pL = omech.L
        
        [lst_qidot, lst_pidot] = omech.Hamiltons_equations(pL, lst_qi, lst_qidot, lst_pi, lst_pidot)
        display(lst_qidot, lst_pidot)
        
        1. Calculate generalize momenta by taking derivative of Lagrangian with respect to q_idot.
        2. Solve q_idots from generalize momenta equations.
        3. Replace q_idots in Lagrangian with corresponding generalize momenta.
        4. Replace pi*qidot in Hamiltonian with expressions written in terms of generalize momenta.
        5. Calculate qidot, p_idot, p_idot by Hamilton's equations. 
        """

        #  1. Calculate generalize momenta by taking derivative of Lagrangian with respect to q_idot.
        dim = len(lst_qidot)
        eq_pis, sol_qidots, res_qidots, res_pidots = [],[],[],[]
        sub_qidots = dict()
        for i in range(dim):
            """
            eq_px = omech.p_i.xreplace({L:omech.L.rhs, q_idot:xdot, p_i:px}).doit() -> p_x = m*xdot
            eq_py = omech.p_i.xreplace({L:omech.L.rhs, q_idot:ydot, p_i:py}).doit()
            eq_pz = omech.p_i.xreplace({L:omech.L.rhs, q_idot:zdot, p_i:pz}).doit()
            """
            ieq_pi = self.p_i.xreplace({L:self.L.rhs, q_idot:lst_qidot[i], p_i:lst_pi[i]}).doit()
            eq_pis.append(ieq_pi)
   
        # 2. Solve q_idots from generalize momenta equations.
        for i in range(dim):
            """
            sol_x_dot = solve(eq_px, xdot)[0]
            sol_y_dot = solve(eq_py, ydot)[0]
            sol_z_dot = solve(eq_pz, zdot)[0]
            """
            isol = solve(eq_pis[i], lst_qidot[i])[0]
            sol_qidots.append(isol)
        
        # 3. Replace q_idots in Lagrangian with corresponding generalize momenta.
        for i in range(dim):
            # sub_qidots = {xdot:sol_x_dot, ydot:sol_y_dot, zdot:sol_z_dot}
            sub_qidots[lst_qidot[i]] =  sol_qidots[i]
        self.L = self.L.subs(sub_qidots)
        
        # 4. Replace pi*qidot in Hamiltonian with expressions written in terms of generalize momenta.
        # piqidot = Matrix([[px,py,pz]]).dot(Matrix([[sol_x_dot,sol_y_dot,sol_z_dot]]))
        piqidot = Matrix([lst_pi]).dot(Matrix([sol_qidots]))
        sub_H = {n:1, L:self.L.rhs, p_i*q_idot:piqidot}
        self.H = Eq(self.H.lhs, ratsimp(self.H.rhs.xreplace(sub_H).doit()))
        
        # 5. Calculate qidot, p_idot by Hamilton's equations.
        for i in range(dim):
            """
            xdot = omech.q_idot.xreplace({H:omech.H.rhs, p_i:px, q_idot:xdot})
            ydot = omech.q_idot.xreplace({H:omech.H.rhs, p_i:py, q_idot:ydot})
            zdot = omech.q_idot.xreplace({H:omech.H.rhs, p_i:pz, q_idot:zdot})
            # zdot = omech.Hamiltons_equations_I.xreplace({H:omech.H.rhs, p_i:pz, q_idot:zdot})
            pxdot = omech.p_idot.xreplace({H:omech.H.rhs, q_i:x, p_idot:pxdot})
            pydot = omech.p_idot.xreplace({H:omech.H.rhs, q_i:y, p_idot:pydot})
            pzdot = omech.p_idot.xreplace({H:omech.H.rhs, q_i:z, p_idot:pzdot})
            """
            iqidot = self.q_idot.xreplace({H:self.H.rhs, p_i:lst_pi[i], q_idot:lst_qidot[i]}).doit()
            ipidot = self.p_idot.xreplace({H:self.H.rhs, q_i:lst_qi[i], p_idot:lst_pidot[i]}).doit()
            res_qidots.append(iqidot)
            res_pidots.append(ipidot)
            
        if self.verbose:
            libsympy.pprints(
                pL, lst_qi, lst_qidot, lst_pi, lst_pidot, omech.p_i,
                eq_pis, sub_qidots, self.L,
                sub_H, self.H,
                self.Hamiltons_equations_I,
                self.Hamiltons_equations_II,
                *res_qidots,
                *res_pidots,
                output_style=self.output_style)

        return(res_qidots, res_pidots)
        
    @staticmethod
    def __doc__():
        return("Document of mechanics class.")
        
omech = mechanics("scalar") # Create an omech object from mechanics class.
# omech = mechanics("scalar")
