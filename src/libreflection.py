#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
libreflection.py
Created on Fri Mar 11 12:53:02 2022
todo: place into book explain step by step

omec.__init__()
omec.solver.verbose = True
commands = ["solve", "NewtonsLaw2", a]
print(omec.process(commands))
"""
import sys
lstPaths = ["../../libpython/src"]
for ipath in lstPaths:
    if ipath not in sys.path:
        sys.path.append(ipath)
#from collections import namedtuple
from sympy import print_python
import libsympy as ls


class branch:
    """
    Example:
    ========
    omec = mechanics()
    omec.NewtonsLaw2
    
    getattr(globals()['mechanics'](), 'NewtonsLaw2')
    getattr(globals()['mechanics'](), 'NewtonsLaw2')('sample arg')
    globals()["solve"](om.NewtonsLaw2,a)
    """
    
    def __init__(self):
        self.classname = type(self).__name__
#        self.classname = self.__class__.__name__ # same as above
        self.solver = solver()
        
    def process(self, commands):
        # omec.process(commands)
        self.result = self.solver.process(commands, self.classname) 
        return(self.result)
    
    
    def get_formulary(self, style="name-eq", verbose=True):
        """
        Display all instance variables of a <branch> object as a formulary list.
        
        [(i,j) for (i,j) in vars(ocar).items()] -> [('brand', 'A'), ('year', 2022)]
        vars(ocar).items() -> dict_items([('brand', 'A'), ('year', 2022)])
        [getattr(ocar,ikey) for (ikey, ival) in vars(ocar).items()] -> ['A', 2022]
        
        attribute_values = [getattr(ocar,ikey) for (ikey, ival) in vars(ocar).items()]
        display(*omec.get_formulary())
        
        todo: export to a latex file.
        
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
                ls.pprints(*res,
                           output_style=self.solver.output_style,
                           newline=self.solver.newline)
        if style == "name-eq":
            """
            name = ikey; ival = eq
            """
            res = [(ikey,ival) for (ikey,ival) in vars(self).items()]
            if verbose:
                for ikey,ival in res:
                    ls.pprints(ikey, ival,
                               output_style=self.solver.output_style,
                               newline=self.solver.newline)
        return(res)
        
    def get_subformulary(self, style="name-eq", verbose=True):
        """
        Display all instance variables of subformulary class of a <branch> object as a formulary list.
        vars(omec.subformulary).items()
        
        Examples
        ========
        display(*omech.get_subformulary())
        """
        
        if style == "eq":
            res = [getattr(self.subformulary,ikey) for (ikey, ival) in vars(self.subformulary).items()]
            if verbose:
                ls.pprints(*res,
                           output_style=self.solver.output_style,
                           newline=self.solver.newline)
        if style == "name-eq":
            """
            name = ikey; ival = eq
            """
            res = [(ikey,ival) for (ikey,ival) in vars(self.subformulary).items()]
            if verbose:
                for ikey,ival in res:
                    ls.pprints(ikey, ival,
                               output_style=self.solver.output_style,
                               newline=self.solver.newline)
        return(res)
        
    def get_symbols(self):
        """
        todo 
        List all symbols used in the branch and defined under define_symbols() function.
        """
        pass


class solver:
    def __init__(self):
        self.codes = []
        self.newline = True     # enable newline break in the output
        self.verbose = False    # enable verbose output of intermediate steps.
        self.output_style = {1:"display", 2:"pprint", 3:"print", 4:"latex"}[1]
    
    def get_codes(self): # todo
        return(self.codes)
    
    def process(self, commands, classname):
        """
        Structure
        ---------
        verb - subject - object
        Equate Newton's 2nd Law to Hooke
        commands = ["Eq", "NewtonsLaw2", "HookesLaw"]
        
        verb = dsolve, solve, Eq, Subs, subs, etc.
        
        commands = [verb, subject, object] -> ["solve", "NewtonsLaw2", a]
        getattr(globals()['mechanics'](), "solve") -> getattr(globals()[classname], method)
        globals()["solve"](om.NewtonsLaw2, a) -> getattr(globals()[classname](), commands[1])
        
        
        Example: Solve a from F = ma
        ----------------------------
        [a,F] = symbols('a F', real=True)
        [k,m,t,w] = symbols('k m t w', real=True, positive=True)
        x = Function('x')(t)
        omec.__init__()
        omec.solver.verbose = True
        commands = ["solve", "NewtonsLaw2", a]
        print(omec.process(commands))
        
        """
#        commands = "sentence1. sentence2. sentence3."
#        commands = commands.split('.')
#        commands.__delitem__(-1)
        [verb, subject, obj] = [commands[0], commands[1], commands[2]]
        
        # kaldik dallanma yap verbose_type = python_code, text_explanation     
        if self.verbose:
            output = [#"{0}({1}, {2})".format(cmd.__name__, expr, params),
                      ' '.join(map(str, commands))]
            ls.pprints(*output,
                       output_style=self.output_style,
                       newline=self.newline)
        
        if verb in ["dsolve", "solve"]:
            """
            dsolve
            ======
            dsolve(eq, f(x), hint) -> Solve ordinary differential equation eq for function f(x), using method hint.
            commands = ["dsolve", "om.result", x]
            
            solve
            =====
            solve(f, *symbols, **flags)
            commands = ["solve", "NewtonsLaw2", a]
            
            solve(F=ma, a) ~> command[0](command[1], command[2])
            solve(Eq(F, a*m),a)
            
            """
            # Check om.result (object.method) pattern in a command.
            if subject.find('.') != -1:
                (classname, method) = subject.split('.')
                expr = getattr(globals()[classname], method)
            else:
                # getattr(globals()['mechanics'](), 'NewtonsLaw2')
                expr = getattr(globals()[classname](), subject) # f
                
            cmd  = globals()[verb] # dsolve
            params = obj
            res = cmd(expr, params)
            
            # kaldik strcode yerine print_python denenmelidir.
            strcode = "{0}({1}, {2})".format(cmd.__name__, expr, params)
            self.codes.append(strcode+'\n')
            
            if self.verbose:
                print(strcode)
            
        elif verb in ["Eq",]:
            """
            Eq
            ==
            Eq(y, x + x**2)
            commands = ["Eq", "NewtonsLaw2", "HookesLaw"]
            
            """
            cmd = globals()[verb] # Eq
            lhs = getattr(globals()[classname](), subject).rhs
            rhs = getattr(globals()[classname](), obj).rhs
            expr = lhs-rhs
            params = 0
            res = cmd(expr, 0)
            
            strcode = "{0}({1}, {2})".format(cmd.__name__, expr, 0)
            self.codes.append(strcode+'\n')
            
            if self.verbose:
                ls.pprints(subject,
                           getattr(globals()[classname](), subject),
                           obj,
                           getattr(globals()[classname](), obj),
                           output_style=self.output_style)
                print(strcode)
        
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
            self.codes.append(strcode+'\n')
            
            if self.verbose:
                print(strcode)
        
        elif verb in ["subs","xreplace"]:
            """
            subs
            ====
            expr = x**3 + 4*x*y - z
            expr.subs([(x, 2), (y, 4), (z, 0)])
            expr.subs([(x, a)])
            commands = ["subs", "om.result", [(a, diff(x, t, 2, evaluate=False))]]
            osm.process(commands)
            
            xreplace
            ========
            Replace occurrences of objects within the expression.
            
            xreplaces = {g:1, engF:mu*B*(2*i-3), j:1, n:2}
            commands = ["xreplace", "osm.Zsp", xreplaces]
            osm.process(commands)
            
            """
            (classname, method) = subject.split('.') # 'obj.result' -> ['obj', 'result']
            expr = getattr(globals()[classname], method)
            cmd  = getattr(expr, verb) # expr.subs
            params = obj
            res = cmd(params)
            
            strcode = "{0}({1}, {2})".format(expr, cmd.__name__, params)
            self.codes.append(strcode+'\n')
            
            if self.verbose:
                print(strcode)
        
        ls.pprints(res, output_style=self.output_style)
        return(res)