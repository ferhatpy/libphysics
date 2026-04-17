from libphysics.quantum_mechanics import *
from sympy import *

oqmec = quantum_mechanics()
n = symbols('n', positive=True, integer=True)
hbar, w = symbols('hbar w', real=True, positive=True)
psi0c = oqmec.qho.nb
psi0  = oqmec.qho.nk
En0   = oqmec.qho.En

oqmec.Hp = Hp = S(1)/4*hbar*w*(oqmec.qho.a + oqmec.qho.ad)**4

try:
    En1_s1 = oqmec.En_ND_PT(1, psi0c, psi0, Hp, En0, k2min=n, k2max=n)
    print("Success:", En1_s1)
except TypeError as e:
    print("Caught TypeError:", e)
    import traceback
    traceback.print_exc()
except Exception as e:
    print("Caught Exception:", type(e), e)
    import traceback
    traceback.print_exc()
