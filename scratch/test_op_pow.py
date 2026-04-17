from sympy import *
from sympy.physics.quantum import *
from sympy.physics.quantum.sho1d import *

a = LoweringOp('a')
ad = RaisingOp('a')
expr = (a + ad)**4
print(expr)
print(type(expr))

psi = SHOKet(0)
bra = SHOBra(0)
res = qapply(bra * expr * psi)
print(res)
