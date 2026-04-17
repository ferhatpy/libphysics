from sympy import *
from sympy.physics.quantum import *
from sympy.physics.quantum.sho1d import *

# Emulate libphysics environment
n = symbols('n', positive=True, integer=True)
hbar, w = symbols('hbar w', real=True, positive=True)
m = symbols('m', real=True, positive=True)

class QHO:
    def __init__(self):
        self.a = LoweringOp('a')
        self.ad = RaisingOp('a')
        self.nk = lambda n=n: SHOKet(n)
        self.nb = lambda n=n: SHOBra(n)

qho = QHO()
Hp = S(1)/4*hbar*w*(qho.a + qho.ad)**4

def V_ij(i, j):
    return qapply(qho.nb(i) * qapply(Hp * qho.nk(j)))

res = V_ij(n, n)
print(res)
print(qapply(res))
