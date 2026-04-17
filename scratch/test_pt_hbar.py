from sympy import *
from sympy.physics.quantum import *
from sympy.physics.quantum.sho1d import *
from sympy.physics.quantum.constants import hbar

# Emulate libphysics environment
n = symbols('n', positive=True, integer=True)
w = symbols('w', real=True, positive=True)
a = LoweringOp('a')
ad = RaisingOp('a')
Hp = S(1)/4*hbar*w*(a + ad)**4

def psi0c(i): return SHOBra(i)
def psi0(j): return SHOKet(j)

total = qapply(psi0c(n)*qapply(Hp*psi0(n)))
print(total)
