from sympy import *
from sympy.physics.quantum import *
from sympy.physics.quantum.sho1d import *

n = symbols('n', integer=True)
k = symbols('k', integer=True)
psi = SHOKet(n)
expr = psi / 2
print(expr)
print(qapply(expr))
