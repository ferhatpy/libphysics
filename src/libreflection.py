#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
libreflection.py
Created on Fri Mar 11 12:53:02 2022
todo: place into book explain step by step

https://github.com/ferhatpy/libphysics

omech.__init__()
omech.verbose = True
commands = ["solve", "NewtonsLaw2", a]
print(omech.process(commands))
"""
import sys
lstPaths = ["../../libpython/src"]
for ipath in lstPaths:
    if ipath not in sys.path:
        sys.path.append(ipath)
#from collections import namedtuple
from sympy import *
import libsympy

"""
Getting Class Name:
omech.__class__.__name__ --> 'mechanics'
self.__class__.__name__ --> 'mechanics'
type(self).__name__     --> 'mechanics'

getattr(globals()['mechanics'](), 'NewtonsLaw2')
getattr(globals()['mechanics'](), 'NewtonsLaw2')('sample arg')
globals()["solve"](omech.NewtonsLaw2, a)
"""

#---- branch
class branch:
    """
    Example:
    ========
    omech = mechanics()
    omech.x              # Symbolic x position as a function of time.
    omech.x_t            # Numeric x position as a function of time.
    omech.NewtonsLaw2
    """
    def __init__(self):
        self.classname = type(self).__name__
        # self.classname = self.__class__.__name__ # same as above
        self.codes = []
        self.newline = True     # enable newline break in the output
        self.verbose = False    # enable verbose output of intermediate steps.
        self.output_style = {1:"display", 2:"pprint", 3:"print", 4:"latex"}[1]
        # self.solver = solver() # Assign a solver to my branch.
        
#    def process(self, commands):
#        # omech.process(commands)
#        self.result = self.process(commands, self.classname) 
#        return(self.result)
    
    def get_codes(self): # todo
        return(self.codes)
    
    def get_formulary(self, style="name-eq", verbose=True):
        """
        Display all instance variables of a <branch> object as a formulary list.
        
        [(i,j) for (i,j) in vars(ocar).items()] -> [('brand', 'A'), ('year', 2022)]
        vars(ocar).items() -> dict_items([('brand', 'A'), ('year', 2022)])
        [getattr(ocar,ikey) for (ikey, ival) in vars(ocar).items()] -> ['A', 2022]
        
        attribute_values = [getattr(ocar,ikey) for (ikey, ival) in vars(ocar).items()]
        display(*omech.get_formulary())
        
        todo: 
            1) export to a latex file.
            2) 
        
        
        Parameters
        ==========
    
        expr : .
    
        vars : .
    
        kwargs : ``style``, ``newline`` 
        
        Examples
        ========
        display(*omech.get_subformulary())
        """
        if style == "eq":
            res = [getattr(self,ikey) for (ikey, ival) in vars(self).items()]
            if verbose:
                libsympy.pprints(*res,
                           output_style=self.output_style,
                           newline=self.newline)
        if style == "name-eq":
            """
            name = ikey; ival = eq
            """
            res = [(ikey,ival) for (ikey,ival) in vars(self).items()]
            if verbose:
                for ikey,ival in res:
                    libsympy.pprints(ikey, ival,
                               output_style=self.output_style,
                               newline=self.newline)
        if style == "mathematica":
            res = []
            for (ikey, ival) in vars(self).items():
                try:
                    icode = mathematica_code(ival)
                    if verbose:
                        libsympy.pprints(ikey, icode,
                                   output_style=self.output_style,
                                   newline=self.newline)
                except:
                    pass 
                res.append((ikey, icode))
        
        return(res)
        
    def get_subformulary(self, style="name-eq", verbose=True):
        """
        Display all instance variables of subformulary class of a <branch> object as a formulary list.
        vars(omech.subformulary).items()
        
        Examples
        ========
        display(*omech.get_subformulary())
        """
        
        if style == "eq":
            res = [getattr(self.subformulary,ikey) for (ikey, ival) in vars(self.subformulary).items()]
            if verbose:
                libsympy.pprints(*res,
                           output_style=self.output_style,
                           newline=self.newline)
        if style == "name-eq":
            """
            name = ikey; ival = eq
            """
            res = [(ikey,ival) for (ikey,ival) in vars(self.subformulary).items()]
            if verbose:
                for ikey,ival in res:
                    libsympy.pprints(ikey, ival,
                               output_style=self.output_style,
                               newline=self.newline)
        return(res)
        
    def get_symbols(self):
        """
        todo 
        List all symbols used in the branch and defined under define_symbols() function.
        """
        pass

    def process(self, commands):
        """
        Structure
        ---------
        Sentence with 3 words: verb - subject - object
        Sentence with 4 words: verb - subject - object - args
        verb = dsolve, solve, Eq, Subs, subs, etc.
        
        getattr(expr, verb)
        getattr(globals()[classname], method)
        getattr(globals()['omech'], 'NewtonsLaw2') --> F = m*d^2/dt^2
        getattr(globals()['mechanics'](), "NewtonsLaw2") --> F = m*d^2/dt^2
        OR
        vars(omech)['NewtonsLaw2'] --> vars(classname)[method] --> F = m*d^2/dt^2
        
        globals()["solve"](omech.NewtonsLaw2, omech.a.rhs)
        
        Equate Newton's 2nd Law to Hooke's Law
        --------------------------------------
        commands = [verb, subject, object] --> ["Eq", "NewtonsLaw2", "HookesLaw"]
        commands = [verb, subject, object] --> ["solve", "NewtonsLaw2", a]
        
        Example: Solve a from F = ma
        ----------------------------
        [a,F] = symbols('a F', real=True)
        [k,m,t,w] = symbols('k m t w', real=True, positive=True)
        x = Function('x')(t)
        omech.__init__()
        omech.verbose = True
        commands = ["solve", "NewtonsLaw2", omech.a.rhs]
        print(omech.process(commands))
        
        """
        
#        commands = "sentence1. sentence2. sentence3."
#        commands = commands.split('.')
#        commands.__delitem__(-1)
        if len(commands)==3:
            # verb - subject - object
            [verb, subject, obj] = [commands[0], commands[1], commands[2]]
        elif len(commands)==4:
            # verb - subject - object - args
            [verb, subject, obj, args] = [commands[0], commands[1], commands[2], commands[3]]
        
        # kaldik dallanma yap verbose_type = python_code, text_explanation     
        if self.verbose:
            output = [#"{0}({1}, {2})".format(cmd.__name__, expr, params),
                      ' '.join(map(str, commands))]
            libsympy.pprints(*output,
                       output_style=self.output_style,
                       newline=self.newline)
        
        if verb in ["dsolve", "solve"]:
            """
            dsolve
            ======
            dsolve(eq, f(x), hint) --> Solve ordinary differential equation eq for function f(x), using method hint.
            commands = ["dsolve", "omech.result", x]
            dsolve(Derivative(f(x), x, x) + 9*f(x), f(x))
            
            solve
            =====
            solve(f, *symbols, **flags)
            commands = ["solve", "NewtonsLaw2", a]
            
            solve(F=ma, a) --> command[0](command[1], command[2])
            solve(Eq(F, a*m),a)
            
            cmd(expr, params, args)
            dsolve_system(eqs, ics={f(0): 1, g(0): 0})
            
            """
            # Check omech.result (object.method) pattern in a command.
            if subject.find('.') != -1:
                (classname, method) = subject.split('.')
#                expr = getattr(globals()[classname], method)
                expr = vars(self)[method]
            else:
#                getattr(globals()['mechanics'](), 'NewtonsLaw2')
#                expr = getattr(globals()[self](), subject)
                expr = vars(self)[subject]
            cmd  = globals()[verb] # dsolve, etc.
            params = obj
            
            if len(commands)==3:
                res = cmd(expr, params)
                # kaldik strcode yerine print_python denenmelidir.
                strcode = "{0}({1}, {2})".format(cmd.__name__, expr, params)
            elif len(commands)==4:
                res = cmd(expr, params, args)
                strcode = "{0}({1}, {2}, {3})".format(cmd.__name__, expr, params, args)
            
        elif verb in ["Eq",]:
            """
            Eq
            ==
            Eq(y, x + x**2)
            commands = ["verb", "subject", "object"]
            commands = ["Eq", "NewtonsLaw2", "HookesLaw"]
            
            The following lines generates an error vector calculations.
            The solutions is to use global reflection methods instead of class
            based reflection methods.
            cmd = globals()[verb] # Eq
            lhs = vars(self)[subject].rhs
            rhs = vars(self)[obj].rhs
            """
            cmd = globals()[verb] # Eq
            lhs = getattr(globals()[self.classname](), subject).rhs
            rhs = getattr(globals()[self.classname](), obj).rhs
            
            expr = lhs-rhs
            params = 0
            res = cmd(expr, 0)
            
            strcode = "{0}({1}, {2})".format(cmd.__name__, expr, 0)
            
        elif verb in ["laplace_transform",]:
            """
            laplace_transform
            ==
            lap_trans = Eq(laplace_transform(omech.driven_oscillator2.lhs, t, p),
                           laplace_transform(omech.driven_oscillator2.rhs, t, p, noconds=True))
            
            commands = ["laplace_transform", "driven_oscillator2", (t,p)]
            omech.process(commands)
            """
            cmd = globals()[verb] # laplace_transform
#            lhs = getattr(globals()[classname](), subject).lhs
            lhs = vars(self)[subject].lhs
#            rhs = getattr(globals()[classname](), subject).rhs
            rhs = vars(self)[subject].rhs
            param1, param2 =  (obj[0], obj[1])
            res = Eq(cmd(lhs, param1, param2), cmd(rhs, param1, param2, noconds=True))
            
            strcode  = "Laplace transform of the {0} equation.\n".format(subject)
            strcode += "Eq({0}({1}, {3}, {4}), {0}({2}, {3}, {4}, noconds=True))".format(cmd.__name__, lhs, rhs, param1, param2)
            
        elif verb in ["Subs",]:
            """
            Subs
            ====
            Subs(expr, x, x0) represents the expression resulting from substituting x with x0 in expr.
            
            """
            cmd = globals()[verb] # Sub
            params = commands[3]
            (expr, x, x0) = (subject, obj, params)
            res = cmd(expr, x, x0)
            
            strcode = "{0}({1}, {2}, {3})".format(cmd.__name__, expr, x, x0)
            
        elif verb in ["subs","xreplace"]:
            """
            subs
            ====
            expr = x**3 + 4*x*y - z
            expr.subs([(x, 2), (y, 4), (z, 0)])
            expr.subs([(x, a)])
            commands = ["subs", "omech.result", [(a, diff(x, t, 2, evaluate=False))]]
            ostat.process(commands)
            
            xreplace
            ========
            Replace occurrences of objects within the expression.
            
            xreplaces = {g:1, engF:mu*B*(2*i-3), j:1, n:2}
            commands = ["xreplace", "ostat.Zsp", xreplaces]
            ostat.process(commands)
            
            """
            (classname, method) = subject.split('.') # 'obj.result' -> ['obj', 'result']
#            expr = getattr(globals()[classname], method)
            expr = vars(self)[method]
            cmd  = getattr(expr, verb) # expr.subs
            params = obj
            res = cmd(params)
            
            strcode = "{0}({1}, {2})".format(expr, cmd.__name__, params)
        
         
        self.codes.append(strcode+'\n')            
        if self.verbose:print(strcode)
        libsympy.pprints(res, output_style=self.output_style)
        self.result = res
        return(res)
        
"""    
#---- solver
class solver:
    def __init__(self):
        self.codes = []
        self.newline = True     # enable newline break in the output
        self.verbose = False    # enable verbose output of intermediate steps.
        self.output_style = {1:"display", 2:"pprint", 3:"print", 4:"latex"}[1]
    
    def get_codes(self): # todo
        return(self.codes)
    
    def process(self, commands, classname):
        pass
"""        