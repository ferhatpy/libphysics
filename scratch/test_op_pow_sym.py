from sympy import *
from sympy.physics.quantum import *
from sympy.physics.quantum.sho1d import *

n = symbols('n', integer=True)
a = LoweringOp('a')
ad = RaisingOp('a')
expr = (a + ad)**4
psi = SHOKet(n)
bra = SHOBra(n)
res = qapply(bra * expr * psi)
print(res)
