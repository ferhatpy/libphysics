from sympy import *
from sympy.physics.quantum import *
from sympy.physics.quantum.sho1d import *

n = symbols('n', integer=True)
a = LoweringOp('a')
ad = RaisingOp('a')
expr = (a + ad)**4
psi = SHOKet(n)
res = qapply(expr * psi)
print(res)
