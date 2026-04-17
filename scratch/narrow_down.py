from libphysics.quantum_mechanics import *
from sympy import *
import traceback

oqmec = quantum_mechanics()
n = symbols('n', positive=True, integer=True)
hbar, w = symbols('hbar w', real=True, positive=True)

Hp = S(1)/4*hbar*w*(oqmec.qho.a + oqmec.qho.ad)**4
psi = oqmec.qho.nk(n)

print("Hp type:", type(Hp))
print("psi type:", type(psi))

try:
    print("Trying qapply(Hp * psi)...")
    res = qapply(Hp * psi)
    print("Success!")
except TypeError:
    print("Caught expected TypeError!")
    traceback.print_exc()
except Exception as e:
    print(f"Caught other exception: {type(e).__name__}: {e}")
    traceback.print_exc()
